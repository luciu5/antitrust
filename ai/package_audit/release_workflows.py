#!/usr/bin/env python3
"""Plan, push and certify the frozen coordinated release in a GitHub terminal.

Default is an offline plan. --push-and-certify advances the seven reviewed
release branches and waits for their automatic full release-tier checks.
Requires Python 3.10+, Git and an authenticated GitHub CLI (`gh`). No mainline
branch, tag, release or package registry is changed. Matching certification
runs are reused on restart. Failed runs require --retry-failed to retry.
Push events work before promotion; manual workflow dispatch would require
workflow definitions on protected default branches and is not used here.
Explicit, SHA-bound user waivers in the handoff skip only the named workflow;
reports retain waived results separately from successful certifications.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path
import subprocess
import sys
import time
from urllib.parse import quote, urlencode

import push_coordinated_release as push


WORKFLOWS = {
    "antitrust": ("qa.yml", (), {"test_tier": "nightly"}),
    "coordination": ("R-CMD-check.yaml", ("antitrust",), {}),
    "trade": ("qa.yml", ("antitrust", "coordination"), {"test_tier": "nightly"}),
    "vertical": ("R-CMD-check.yaml", ("antitrust", "coordination", "trade"), {}),
    "competitiontoolbox": ("dependency-stack-qa.yml", ("antitrust", "coordination", "trade", "vertical"), {"test_tier": "NIGHTLY"}),
    "financial": ("qa.yml", ("antitrust",), {"extended": True}),
    "antitrustBayes": ("qa.yml", ("antitrust", "coordination", "trade"), {"run_stan_tests": False}),
}


def text_value(value):
    return str(value).lower() if isinstance(value, bool) else str(value)


def certification_inputs(name, packages):
    _, deps, extra = WORKFLOWS[name]
    values = {"release_certification": True} if deps else {}
    values.update({f"{dep}_release_sha": packages[dep]["sha"] for dep in deps})
    values.update(extra)
    return values


def run_name(values):
    return "certification-" + "-".join(text_value(value) for value in values.values())


def run_name_expression(name):
    """Source expression mirrored by the helper's exact run identity check."""
    _, deps, extra = WORKFLOWS[name]
    keys = (["release_certification"] if deps else []) + [f"{dep}_release_sha" for dep in deps] + list(extra)
    template = "certification-" + "-".join("{" + str(i) + "}" for i in range(len(keys)))
    fields = ", ".join("inputs." + key for key in keys)
    return "${{ github.event_name == 'push' && startsWith(github.ref, 'refs/heads/release/') && format('release-certification-{0}', github.sha) || github.event_name == 'workflow_dispatch' && format('" + template + "', " + fields + ") || github.event_name }}"


def gh_api(path, payload=None):
    command = ["gh", "api", "--hostname", "github.com", "--method", "POST" if payload is not None else "GET", path]
    if payload is not None:
        command.extend(["--input", "-"])
    try:
        result = subprocess.run(command, input=json.dumps(payload) if payload is not None else None,
                                capture_output=True, text=True, timeout=120, check=False)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise push.PushError(f"GitHub CLI operation failed ({type(error).__name__}); check gh access/authentication.") from error
    if result.returncode:
        # Credential helper and remote stderr are deliberately not recorded.
        raise push.PushError(f"GitHub CLI request failed (exit {result.returncode}); check gh authentication and Actions permissions.")
    return json.loads(result.stdout) if result.stdout.strip() else {}


def endpoint(item, suffix):
    return f"repos/{item['repository']}/{suffix}"


def verify_branch(item):
    branch = item["ref"].removeprefix("refs/heads/")
    response = gh_api(endpoint(item, "git/ref/heads/" + quote(branch, safe="/")))
    if response.get("object", {}).get("sha") != item["sha"]:
        raise push.PushError(f"{item['package']}: GitHub release branch differs from the frozen candidate.")


def matches(item, run):
    return (run.get("head_sha") == item["sha"] and run.get("event") == "push"
            and run.get("display_title") == item["run_name"]
            and isinstance(run.get("workflow_id"), int)
            and (not item.get("workflow_id") or run["workflow_id"] == item["workflow_id"])
            and run.get("head_branch") == item["ref"].removeprefix("refs/heads/")
            and run.get("path", "").split("@", 1)[0] == ".github/workflows/" + item["workflow"])


def matching_run(item):
    query = urlencode({"event": "push", "branch": item["ref"].removeprefix("refs/heads/"),
                       "head_sha": item["sha"], "per_page": 100})
    data = gh_api(endpoint(item, f"actions/runs?{query}"))
    runs = [r for r in data.get("workflow_runs", []) if matches(item, r)]
    return max(runs, key=lambda r: r["id"]) if runs else None


def remember_run(item, run):
    item.update(run_id=run["id"], workflow_id=run["workflow_id"], run_attempt=run.get("run_attempt", 1), url=run["html_url"], status=run["status"], conclusion=run.get("conclusion"))


def apply_waivers(plan, handoff):
    waivers = handoff.get("certification_waivers", {})
    if not isinstance(waivers, dict) or set(waivers) - {"competitiontoolbox"}:
        raise push.PushError("Only the explicitly waived competitiontoolbox workflow may be skipped.")
    for item in plan:
        waiver = waivers.get(item["package"])
        if waiver is None:
            continue
        if (not isinstance(waiver, dict) or waiver.get("candidate_sha") != item["sha"]
            or waiver.get("workflow") != ".github/workflows/" + item["workflow"]
            or waiver.get("authorized_by") != "user"
            or not isinstance(waiver.get("reason"), str) or not waiver["reason"].strip()):
            raise push.PushError("Toolbox waiver must name this exact candidate/workflow and the user's reason.")
        item["waiver"] = dict(waiver)


def certify(plan, retry_failed=False, poll_seconds=30, timeout_seconds=10800):
    # Every branch is checked before optional retries. Initial checks are
    # triggered by the release push, even if no workflow exists on main/master.
    for item in plan:
        verify_branch(item)
    started = time.monotonic()
    for item in plan:
        if item.get("waiver"):
            item.update(status="waived", conclusion=None)
            print(f"Certification: {item['package']} WAIVED: {item['waiver']['reason']}", flush=True)
            continue
        verify_branch(item)
        existing = matching_run(item)
        if existing:
            remember_run(item, existing)
            if existing["status"] == "completed" and existing.get("conclusion") != "success" and retry_failed:
                item["status"] = "retry-outcome-unknown"
                item["retry_after_attempt"] = existing.get("run_attempt", 1)
                gh_api(endpoint(item, f"actions/runs/{existing['id']}/rerun-failed-jobs"), {})
                item["status"] = "retry-requested"
        else:
            item["status"] = "waiting-for-push-run"
        verify_branch(item)
        print(f"Certification: {item['package']} {item['status']}", flush=True)
    # Runs are independent because every dependency source has already been
    # pushed. A failed package is reported while other evidence completes.
    while True:
        pending = False
        for item in plan:
            if item.get("status") in ("completed", "waived"):
                continue
            run = matching_run(item)
            if run:
                # A queued retry can briefly expose the prior failed attempt.
                if item.get("status") == "retry-requested" and run.get("run_attempt", 1) <= item["retry_after_attempt"]:
                    pending = True
                    continue
                remember_run(item, run)
            if item.get("status") != "completed":
                pending = True
            print(f"{item['package']}: {item.get('conclusion') or item['status']} {item.get('url', '')}", flush=True)
        if not pending:
            return all(item.get("status") == "waived" or item.get("conclusion") == "success" for item in plan)
        if time.monotonic() - started >= timeout_seconds:
            raise push.PushError("Certification wait timed out. Rerun to reuse matching runs and continue monitoring.")
        time.sleep(poll_seconds)


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ("manifest", "handoff", "validation"):
        p.add_argument("--" + name, type=Path, default=push.WORKSPACE / f"artifacts/coordinated_release_{name}.json")
    p.add_argument("--previous-manifest", type=Path, default=push.WORKSPACE / "artifacts/coordinated_release_previous_manifest.json")
    p.add_argument("--repositories-root", type=Path, default=push.WORKSPACE.parent)
    p.add_argument("--transport", choices=("https", "ssh"), default="https")
    modes = p.add_mutually_exclusive_group()
    modes.add_argument("--push-and-certify", action="store_true", help="push reviewed updates and monitor automatic exact-pin full-tier checks")
    modes.add_argument("--certify", action="store_true", help="monitor automatic checks for candidates already pushed")
    p.add_argument("--retry-failed", action="store_true", help="retry failed jobs of matching certification runs")
    p.add_argument("--report", type=Path, help="new JSON report path; existing reports are preserved")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    push_args = ["--manifest", str(args.manifest), "--handoff", str(args.handoff),
                 "--validation", str(args.validation), "--repositories-root", str(args.repositories_root),
                 "--transport", args.transport, "--previous-manifest", str(args.previous_manifest)]
    plan, error, completed, report = [], None, False, None
    try:
        plan = push.local_plan(push.parser().parse_args(push_args))
        packages = push.load_manifest(args.manifest)["packages"]
        for item in plan:
            name = item["package"]
            item.update(workflow=WORKFLOWS[name][0], inputs=certification_inputs(name, packages))
            item["run_name"] = "release-certification-" + item["sha"]
        apply_waivers(plan, json.loads(args.handoff.read_text()))
        for item in plan:
            name = item["package"]
            if item.get("waiver"):
                print(f"{name}: WAIVED for {item['sha']}: {item['waiver']['reason']}", flush=True)
                continue
            print(f"{name}: {item['sha']} {item['workflow']} {json.dumps(item['inputs'])}", flush=True)
        if not (args.certify or args.push_and_certify):
            print("Offline plan verified. Run --push-and-certify in a terminal with GitHub/gh access.")
            return 0
        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        destination = args.report or push.WORKSPACE / f"artifacts/coordinated_release_evidence/certify-{stamp}.json"
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("x") as output:
            output.write("{}\n")
        report = destination
        # Repository run listing works even when the new workflow exists only
        # on the release branch; workflow_dispatch/default-branch lookup does not.
        for item in plan:
            if item.get("waiver"):
                continue
            gh_api(endpoint(item, "actions/runs?per_page=1"))
        if args.push_and_certify and push.main(push_args + ["--push"]):
            raise push.PushError("Candidate push stopped; rerun to verify/resume the partial release push and its automatic checks.")
        completed = certify(plan, retry_failed=args.retry_failed)
        if not completed:
            error = "One or more certification runs failed; see their recorded GitHub URLs."
            print(error, file=sys.stderr)
        return 0 if completed else 1
    except (push.PushError, push.ManifestError, OSError, ValueError, KeyError, TypeError) as caught:
        error = str(caught)
        print(f"Release certification stopped: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        error = "Interrupted; rerun to reuse matching certification runs."
        print(error, file=sys.stderr)
        return 130
    finally:
        if report:
            waived = [item["package"] for item in plan if item.get("waiver")]
            status = ("complete_with_waivers" if waived else "complete") if completed else "stopped"
            report.write_text(json.dumps(dict(status=status, waived_packages=waived, error=error, packages=plan), indent=2) + "\n")
            print(f"Certification report: {report}", flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
