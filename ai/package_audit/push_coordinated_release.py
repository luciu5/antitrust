#!/usr/bin/env python3
"""Push the seven frozen release candidates without changing main/master.

Default: validate locally and print the plan (no network).
--check-remote: also check GitHub, without pushing.
--push: preflight all seven remotes, then create the exact candidate branches.

Only manifest commits are uploaded. Uncommitted files are left untouched.
Existing remote candidate branches must already equal the selected SHA, unless
--previous-manifest selects their reviewed previous SHAs. Updates must be fast
forwards from those SHAs. Explicit Git leases reject concurrent changes rather
than overwriting them. Cross-repository pushes cannot be atomic: after a failure,
rerun the same command to resume, preserving already matching refs.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import subprocess
import sys
from pathlib import Path

from check_release_manifest import ManifestError, PACKAGE_NAMES, load_manifest


WORKSPACE = Path(__file__).resolve().parents[2]
ORDER = PACKAGE_NAMES


class PushError(Exception):
    """A failed precondition or Git operation; stop before further pushes."""


def git_env() -> dict[str, str]:
    env = os.environ.copy()
    for key in (
        "GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_COMMON_DIR",
        "GIT_NAMESPACE", "GIT_OBJECT_DIRECTORY", "GIT_ALTERNATE_OBJECT_DIRECTORIES",
    ):
        env.pop(key, None)
    env.update(GIT_OPTIONAL_LOCKS="0", GIT_NO_REPLACE_OBJECTS="1", GIT_TERMINAL_PROMPT="0")
    return env


def git(repo: Path, *args: str, timeout: int = 120) -> str:
    try:
        result = subprocess.run(
            ["git", "-C", str(repo), *args], env=git_env(),
            capture_output=True, text=True, timeout=timeout, check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise PushError(f"{repo.name}: Git operation failed ({type(error).__name__}).") from error
    if result.returncode:
        # Git credential helpers/servers can write sensitive diagnostics. Do
        # not put their raw output in the terminal transcript or saved report.
        operation = args[2] if args[0] == "-c" else args[0]
        raise PushError(f"{repo.name}: git {operation} failed (exit {result.returncode}); check access and authentication.")
    return result.stdout.strip()


def remote_url(repository: str, transport: str) -> str:
    return f"https://github.com/{repository}.git" if transport == "https" else f"git@github.com:{repository}.git"


def local_plan(args: argparse.Namespace) -> list[dict]:
    manifest = load_manifest(args.manifest)
    previous = load_manifest(args.previous_manifest) if args.previous_manifest else None
    handoff = json.loads(args.handoff.read_text())
    validation = json.loads(args.validation.read_text())
    release_id = manifest["release_id"]
    if previous and previous["release_id"] != release_id:
        raise PushError("Previous and selected manifests must identify the same release branches.")
    if not re.fullmatch(r"coordinated-\d{4}-\d{2}-\d{2}", release_id):
        raise PushError("Unexpected coordinated release ID.")
    if handoff.get("release_id") != release_id or validation.get("release_id") != release_id:
        raise PushError("Manifest, handoff and validation release IDs differ.")
    proposals = handoff.get("proposals", [])
    checks = validation.get("verification", [])
    if len(proposals) != len(ORDER) or {p["package"] for p in proposals} != set(ORDER):
        raise PushError("Handoff must contain exactly one proposal for each of the seven packages.")
    if len(checks) != len(ORDER) or {r["package"] for r in checks} != set(ORDER):
        raise PushError("Validation must cover exactly the seven packages.")
    proposals = {p["package"]: p for p in proposals}
    checks = {r["package"]: r for r in checks}
    if handoff.get("publication_order") != list(ORDER):
        raise PushError("Handoff publication order differs from the supported dependency order.")
    # Run the existing source/version/dependency checker in a clean Git env.
    checker = Path(__file__).with_name("check_release_manifest.py")
    result = subprocess.run(
        [sys.executable, str(checker), str(args.manifest), "--repositories-root", str(args.repositories_root)],
        env=git_env(), capture_output=True, text=True, check=False,
    )
    if result.returncode:
        raise PushError("Local manifest/source validation failed. Run check_release_manifest.py for details.")
    plan = []
    for name in ORDER:
        candidate, proposal, check = manifest["packages"][name], proposals[name], checks[name]
        sha = candidate["sha"]
        branch = f"release/{release_id}"
        ref = f"refs/heads/{branch}"
        expected_target = "master" if name in ("antitrust", "trade", "competitiontoolbox") else "main"
        if (proposal.get("repository") != candidate["repository"] or
            proposal.get("candidate_sha") != sha or proposal.get("version") != candidate["version"] or
            proposal.get("local_candidate_ref") != ref or
            proposal.get("proposed_remote_candidate_branch") != branch or
            proposal.get("proposed_target_branch") != expected_target):
            raise PushError(f"{name}: handoff target/ref/version differs from the frozen manifest.")
        target_sha = proposal.get("observed_target_sha", "")
        if not re.fullmatch(r"[0-9a-f]{40}", target_sha):
            raise PushError(f"{name}: missing exact reviewed mainline SHA.")
        if (check.get("candidate_sha") != sha or check.get("version") != candidate["version"] or
            check.get("package_check") != "OK" or check.get("failures") != 0):
            raise PushError(f"{name}: no successful package check for this exact candidate.")
        repo = args.repositories_root / name
        if git(repo, "rev-parse", "--verify", ref) != sha:
            raise PushError(f"{name}: local candidate branch moved; refusing to push.")
        git(repo, "merge-base", "--is-ancestor", target_sha, sha)
        previous_sha = None
        if previous:
            old = previous["packages"][name]
            if old["repository"] != candidate["repository"] or old["version"] != candidate["version"]:
                raise PushError(f"{name}: previous manifest repository/version differs.")
            previous_sha = old["sha"]
            git(repo, "merge-base", "--is-ancestor", previous_sha, sha)
        plan.append(dict(package=name, repository=candidate["repository"], version=candidate["version"],
                         sha=sha, ref=ref, target_ref=f"refs/heads/{expected_target}", target_sha=target_sha,
                         previous_sha=previous_sha,
                         url=remote_url(candidate["repository"], args.transport), status="pending"))
    return plan


def remote_refs(repo: Path, item: dict) -> dict[str, str]:
    output = git(repo, "ls-remote", "--refs", item["url"], item["target_ref"], item["ref"])
    refs = {}
    for line in output.splitlines():
        fields = line.split()
        if len(fields) != 2 or not re.fullmatch(r"[0-9a-f]{40}", fields[0]):
            raise PushError(f"{item['package']}: malformed remote ref response.")
        sha, ref = fields
        if ref not in (item["target_ref"], item["ref"]) or ref in refs:
            raise PushError(f"{item['package']}: unexpected/duplicate remote ref.")
        refs[ref] = sha
    return refs


def preflight_remote(repo: Path, item: dict) -> None:
    refs = remote_refs(repo, item)
    if refs.get(item["target_ref"]) != item["target_sha"]:
        raise PushError(f"{item['package']}: remote mainline moved or is missing; review and refresh the handoff first.")
    existing = refs.get(item["ref"])
    if existing != item["sha"]:
        expected = item.get("previous_sha")
        if existing != expected:
            raise PushError(f"{item['package']}: remote release branch differs from the reviewed previous state; it will not be overwritten.")
    item["remote_sha"] = existing
    item["status"] = "already-present" if existing == item["sha"] else "ready"


def push_one(repo: Path, item: dict) -> None:
    # Only compare-and-create, or compare-and-fast-forward from a locally
    # verified ancestor. A different remote SHA fails the exact lease.
    item["status"] = "push-outcome-unknown"
    git(repo, "-c", "push.followTags=false", "push", "--porcelain", "--no-follow-tags",
        "--recurse-submodules=no", f"--force-with-lease={item['ref']}:{item.get('remote_sha') or ''}",
        item["url"], f"{item['sha']}:{item['ref']}", timeout=600)
    item["status"] = "pushed-unverified"
    if remote_refs(repo, item).get(item["ref"]) != item["sha"]:
        raise PushError(f"{item['package']}: post-push SHA verification failed; inspect the remote before resuming.")
    item["status"] = "pushed"


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest", type=Path, default=WORKSPACE / "artifacts/coordinated_release_manifest.json")
    p.add_argument("--handoff", type=Path, default=WORKSPACE / "artifacts/coordinated_release_handoff.json")
    p.add_argument("--validation", type=Path, default=WORKSPACE / "artifacts/coordinated_release_validation.json")
    p.add_argument("--previous-manifest", type=Path, help="allow only reviewed fast-forward updates from these exact previous SHAs")
    p.add_argument("--repositories-root", type=Path, default=WORKSPACE.parent)
    p.add_argument("--transport", choices=("https", "ssh"), default="https")
    modes = p.add_mutually_exclusive_group()
    modes.add_argument("--check-remote", action="store_true", help="read-only preflight against GitHub")
    modes.add_argument("--push", action="store_true", help="push the exact seven remote release branches")
    p.add_argument("--report", type=Path, help="push result JSON (default: timestamped file in release evidence)")
    return p


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    plan = []
    error = None
    report = None
    report_created = False
    completed = False
    try:
        if args.push:
            stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
            report = args.report or WORKSPACE / f"artifacts/coordinated_release_evidence/push-{stamp}.json"
            report.parent.mkdir(parents=True, exist_ok=True)
            # Fail before any remote mutation if reporting is not writable.
            with report.open("x") as output:
                json.dump({"status": "preflight-started"}, output)
            report_created = True
        plan = local_plan(args)
        print("Only frozen manifest commits are selected; working-tree edits are untouched.", flush=True)
        for item in plan:
            print(f"{item['package']} {item['version']}: {item['sha']} -> {item['repository']} {item['ref']}", flush=True)
        if not (args.check_remote or args.push):
            print("Local plan verified. Use --check-remote for read-only GitHub checks, or --push to upload.")
            return 0
        # No uploads until all seven remote destinations pass preflight.
        for item in plan:
            preflight_remote(args.repositories_root / item["package"], item)
            print(f"Remote preflight: {item['package']} {item['status']}", flush=True)
        if args.push:
            for item in plan:
                if item["status"] == "already-present":
                    continue
                # Recheck immediately before this repo's upload. This does
                # not lock mainline; promotion still needs a later head check.
                preflight_remote(args.repositories_root / item["package"], item)
                if item["status"] == "already-present":
                    continue
                push_one(args.repositories_root / item["package"], item)
                print(f"Verified push: {item['package']} {item['sha']}", flush=True)
            print("All seven remote release branches match the manifest. Main/master were not updated.")
        completed = True
        return 0
    except (PushError, ManifestError, OSError, ValueError, KeyError, TypeError) as caught:
        error = str(caught)
        print(f"Release push stopped: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        error = "Interrupted. Rerun the same command to verify remote state and resume."
        print(f"Release push stopped: {error}", file=sys.stderr)
        return 130
    finally:
        if report_created:
            report.write_text(json.dumps(dict(
                completed_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                status="complete" if completed else "stopped", error=error, packages=plan,
            ), indent=2) + "\n")
            print(f"Push report: {report}", flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
