"""Local Git integration gates for release pushes; never contact GitHub.

Run: python3 -m unittest discover -s ai/package_audit -p 'test_push_*.py'
Each test names a material release-integrity or data-preservation failure.
"""

import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import push_coordinated_release as release


class ReleasePushTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="release-push-test-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.local = self.root / "repos"
        self.remotes = self.root / "remotes"
        self.local.mkdir()
        self.remotes.mkdir()
        self.env_patch = patch.dict(os.environ, {
            "GIT_CONFIG_NOSYSTEM": "1", "GIT_CONFIG_GLOBAL": os.devnull,
            "GIT_CONFIG_COUNT": "0", "GIT_AUTHOR_NAME": "Release test",
            "GIT_AUTHOR_EMAIL": "release-test@example.invalid",
            "GIT_COMMITTER_NAME": "Release test",
            "GIT_COMMITTER_EMAIL": "release-test@example.invalid",
        })
        self.env_patch.start()
        self.addCleanup(self.env_patch.stop)
        self.manifest = dict(schema_version=1, release_id="coordinated-2026-11-01", packages={})
        self.handoff = dict(release_id=self.manifest["release_id"], publication_order=list(release.ORDER), proposals=[])
        self.validation = dict(release_id=self.manifest["release_id"], verification=[])
        self.ref = "refs/heads/release/coordinated-2026-11-01"
        self.original = {}
        for name in release.ORDER:
            repo, remote = self.local / name, self.remotes / name
            repo.mkdir()
            self.git(repo, "init", "--quiet")
            self.git(repo, "config", "core.autocrlf", "false")
            target = "master" if name in ("antitrust", "trade", "competitiontoolbox") else "main"
            self.git(repo, "symbolic-ref", "HEAD", f"refs/heads/{target}")
            (repo / "DESCRIPTION").write_text(f"Package: {name}\nVersion: 1.0.0\n")
            (repo / "payload").write_text("baseline\n")
            self.git(repo, "add", ".")
            self.git(repo, "commit", "--quiet", "-m", "baseline")
            base = self.git(repo, "rev-parse", "HEAD")
            self.git(self.remotes, "init", "--quiet", "--bare", str(remote))
            self.git(repo, "push", "--quiet", str(remote), f"{base}:refs/heads/{target}")
            (repo / "payload").write_text("reviewed release\n")
            self.git(repo, "add", "payload")
            self.git(repo, "commit", "--quiet", "-m", "reviewed candidate")
            sha = self.git(repo, "rev-parse", "HEAD")
            self.git(repo, "update-ref", self.ref, sha)
            self.git(repo, "tag", "must-not-push")
            # Active branches and both staged/unstaged work must survive.
            (repo / "payload").write_text("unrelated staged work\n")
            self.git(repo, "add", "payload")
            (repo / "payload").write_text("unrelated unstaged work\n")
            self.original[name] = (self.git(repo, "rev-parse", "HEAD"),
                                   (repo / ".git/index").read_bytes(), (repo / "payload").read_bytes())
            self.manifest["packages"][name] = dict(repository=f"luciu5/{name}", version="1.0.0", sha=sha, minimums={}, remotes=[])
            self.handoff["proposals"].append(dict(package=name, repository=f"luciu5/{name}", version="1.0.0",
                candidate_sha=sha, local_candidate_ref=self.ref,
                proposed_remote_candidate_branch="release/coordinated-2026-11-01",
                proposed_target_branch=target, observed_target_sha=base))
            self.validation["verification"].append(dict(package=name, candidate_sha=sha, version="1.0.0", package_check="OK", failures=0))
        self.paths = {}
        for name, data in [("manifest", self.manifest), ("handoff", self.handoff), ("validation", self.validation)]:
            self.paths[name] = self.root / f"{name}.json"
            self.paths[name].write_text(json.dumps(data))
        self.url_patch = patch.object(release, "remote_url", side_effect=lambda repo, transport: str(self.remotes / repo.split("/")[-1]))
        self.url_patch.start()
        self.addCleanup(self.url_patch.stop)

    def git(self, repo, *args):
        result = subprocess.run(["git", "-C", str(repo), *args], env=release.git_env(),
                                capture_output=True, text=True, check=True)
        return result.stdout.strip()

    def invoke(self, *flags, report="result.json"):
        argv = ["--repositories-root", str(self.local), "--report", str(self.root / report)]
        for name, path in self.paths.items():
            argv.extend([f"--{name}", str(path)])
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return release.main(argv + list(flags))

    def remote_branches(self, name):
        return self.git(self.remotes / name, "for-each-ref", "--format=%(refname) %(objectname)")

    def assert_preserved(self):
        for name, (head, index, payload) in self.original.items():
            repo = self.local / name
            self.assertEqual(self.git(repo, "rev-parse", "HEAD"), head)
            self.assertEqual((repo / ".git/index").read_bytes(), index)
            self.assertEqual((repo / "payload").read_bytes(), payload)

    def test_default_plan_is_offline_and_remote_check_does_not_upload(self):
        """Prevent accidental network access/uploads from plan and read-only modes."""
        before = {name: self.remote_branches(name) for name in release.ORDER}
        with patch.object(release, "remote_refs", side_effect=AssertionError("unexpected network")):
            self.assertEqual(self.invoke(), 0)
        self.assertEqual(self.invoke("--check-remote"), 0)
        self.assertEqual(before, {name: self.remote_branches(name) for name in release.ORDER})
        self.assert_preserved()

    def test_push_preserves_dirty_work_and_targets_and_is_resumable(self):
        """Prevent wrong SHAs, tag leakage, local data loss, and duplicate-run failure."""
        self.assertEqual(self.invoke("--push"), 0)
        for proposal in self.handoff["proposals"]:
            name = proposal["package"]
            remote = self.remotes / name
            self.assertEqual(self.git(remote, "rev-parse", self.ref), proposal["candidate_sha"])
            self.assertEqual(self.git(remote, "rev-parse", proposal["proposed_target_branch"]), proposal["observed_target_sha"])
            self.assertEqual(self.git(remote, "tag", "--list"), "")
        with patch.object(release, "push_one", side_effect=AssertionError("duplicate upload")):
            self.assertEqual(self.invoke("--push", report="resume.json"), 0)
        report = json.loads((self.root / "resume.json").read_text())
        self.assertTrue(all(p["status"] == "already-present" for p in report["packages"]))
        self.assert_preserved()

    def test_last_repo_conflict_prevents_any_earlier_upload(self):
        """Prevent a late preflight conflict from creating a partial release stack."""
        name = release.ORDER[-1]
        base = self.handoff["proposals"][-1]["observed_target_sha"]
        self.git(self.local / name, "push", "--quiet", str(self.remotes / name), f"{base}:{self.ref}")
        before = {n: self.remote_branches(n) for n in release.ORDER}
        self.assertEqual(self.invoke("--push"), 1)
        self.assertEqual(before, {n: self.remote_branches(n) for n in release.ORDER})
        self.assert_preserved()

    def test_moved_local_candidate_and_remote_mainline_are_rejected(self):
        """Prevent publishing an unvalidated candidate or silently ignoring target drift."""
        name = release.ORDER[0]
        base = self.handoff["proposals"][0]["observed_target_sha"]
        sha = self.manifest["packages"][name]["sha"]
        self.git(self.local / name, "update-ref", self.ref, base)
        with patch.object(release, "remote_refs", side_effect=AssertionError("local preflight must fail first")):
            self.assertEqual(self.invoke("--push"), 1)
        self.git(self.local / name, "update-ref", self.ref, sha)
        self.git(self.local / name, "push", "--quiet", str(self.remotes / name), f"{sha}:refs/heads/master")
        before = {n: self.remote_branches(n) for n in release.ORDER}
        self.assertEqual(self.invoke("--push", report="drift.json"), 1)
        self.assertEqual(before, {n: self.remote_branches(n) for n in release.ORDER})

    def test_concurrent_branch_creation_cannot_be_overwritten(self):
        """Prevent the preflight/push race from replacing another actor's branch."""
        real_push = release.push_one
        name = release.ORDER[0]
        base = self.handoff["proposals"][0]["observed_target_sha"]
        def race(repo, item):
            self.git(repo, "push", "--quiet", str(self.remotes / name), f"{base}:{self.ref}")
            return real_push(repo, item)
        with patch.object(release, "push_one", side_effect=race):
            self.assertEqual(self.invoke("--push"), 1)
        self.assertEqual(self.git(self.remotes / name, "rev-parse", self.ref), base)
        for other in release.ORDER[1:]:
            self.assertNotIn(self.ref, self.remote_branches(other))
        self.assert_preserved()

    def test_existing_report_is_not_overwritten_and_no_upload_occurs(self):
        """Prevent reruns from destroying previous audit evidence."""
        path = self.root / "result.json"
        path.write_text("previous audit evidence\n")
        before = {n: self.remote_branches(n) for n in release.ORDER}
        self.assertEqual(self.invoke("--push"), 1)
        self.assertEqual(path.read_text(), "previous audit evidence\n")
        self.assertEqual(before, {n: self.remote_branches(n) for n in release.ORDER})

    def test_partial_upload_failure_stops_dependents_and_resumes(self):
        """Prevent downstream uploads after failure or loss of partial-run evidence."""
        real_push = release.push_one
        failed_name = release.ORDER[2]
        def fail_third(repo, item):
            if item["package"] == failed_name:
                raise release.PushError("simulated third-repository connection failure")
            real_push(repo, item)
        with patch.object(release, "push_one", side_effect=fail_third):
            self.assertEqual(self.invoke("--push"), 1)
        report = json.loads((self.root / "result.json").read_text())
        self.assertEqual(report["status"], "stopped")
        self.assertEqual([p["status"] for p in report["packages"][:2]], ["pushed", "pushed"])
        for name in release.ORDER[2:]:
            self.assertNotIn(self.ref, self.remote_branches(name))
        self.assertEqual(self.invoke("--push", report="resume.json"), 0)
        report = json.loads((self.root / "resume.json").read_text())
        self.assertEqual([p["status"] for p in report["packages"][:2]], ["already-present", "already-present"])
        self.assertEqual([p["status"] for p in report["packages"][2:]], ["pushed"] * 5)
        self.assert_preserved()

    def test_lost_push_acknowledgement_is_reported_as_unknown_and_resumes(self):
        """Prevent misreporting an accepted push as untouched after connection loss."""
        real_git = release.git
        def lose_ack(repo, *args, **kwargs):
            output = real_git(repo, *args, **kwargs)
            if repo.name == release.ORDER[0] and args[0] == "-c":
                raise release.PushError("simulated lost acknowledgement after server accepted push")
            return output
        with patch.object(release, "git", side_effect=lose_ack):
            self.assertEqual(self.invoke("--push"), 1)
        report = json.loads((self.root / "result.json").read_text())
        self.assertEqual(report["status"], "stopped")
        self.assertEqual(report["packages"][0]["status"], "push-outcome-unknown")
        self.assertEqual(self.git(self.remotes / release.ORDER[0], "rev-parse", self.ref),
                         self.manifest["packages"][release.ORDER[0]]["sha"])
        self.assertEqual(self.invoke("--push", report="resume.json"), 0)
        report = json.loads((self.root / "resume.json").read_text())
        self.assertEqual(report["packages"][0]["status"], "already-present")
        self.assert_preserved()

    def prepare_revision(self):
        """Build a descendant candidate without altering the dirty worktrees."""
        previous = self.root / "previous-manifest.json"
        previous.write_text(json.dumps(self.manifest))
        for row in self.handoff["proposals"]:
            name = row["package"]
            repo = self.local / name
            old = row["candidate_sha"]
            self.git(repo, "push", "--quiet", str(self.remotes / name), f"{old}:{self.ref}")
            tree = self.git(repo, "rev-parse", f"{old}^{{tree}}")
            new = self.git(repo, "commit-tree", tree, "-p", old, "-m", "CI-only descendant candidate")
            self.git(repo, "update-ref", self.ref, new)
            self.manifest["packages"][name]["sha"] = new
            row["candidate_sha"] = new
            next(r for r in self.validation["verification"] if r["package"] == name)["candidate_sha"] = new
        for name, data in [("manifest", self.manifest), ("handoff", self.handoff), ("validation", self.validation)]:
            self.paths[name].write_text(json.dumps(data))
        return previous

    def test_revision_requires_exact_previous_manifest_and_fast_forwards(self):
        """Prevent accidental replacement while allowing reviewed CI repair commits."""
        previous = self.prepare_revision()
        self.assertEqual(self.invoke("--push"), 1)
        self.assertEqual(self.invoke("--push", "--previous-manifest", str(previous), report="revision.json"), 0)
        for name in release.ORDER:
            self.assertEqual(self.git(self.remotes / name, "rev-parse", self.ref), self.manifest["packages"][name]["sha"])
        self.assert_preserved()

    def test_previous_manifest_cannot_authorize_rewriting_history(self):
        """Prevent an unrelated selected commit from exploiting an explicit lease."""
        previous = self.prepare_revision()
        name = release.ORDER[0]
        repo = self.local / name
        tree = self.git(repo, "rev-parse", self.ref + "^{tree}")
        unrelated = self.git(repo, "commit-tree", tree, "-m", "unrelated root")
        self.git(repo, "update-ref", self.ref, unrelated)
        self.manifest["packages"][name]["sha"] = unrelated
        self.handoff["proposals"][0]["candidate_sha"] = unrelated
        self.validation["verification"][0]["candidate_sha"] = unrelated
        for key, data in [("manifest", self.manifest), ("handoff", self.handoff), ("validation", self.validation)]:
            self.paths[key].write_text(json.dumps(data))
        with patch.object(release, "remote_refs", side_effect=AssertionError("ancestry must fail locally")):
            self.assertEqual(self.invoke("--push", "--previous-manifest", str(previous)), 1)

    def test_revision_lease_rejects_a_concurrent_remote_advance(self):
        """Prevent a concurrent reviewed-branch update from being overwritten."""
        previous = self.prepare_revision()
        old_manifest = json.loads(previous.read_text())
        real_push = release.push_one
        other = []
        # The concurrent commit is written into the local object database, so
        # transfer it without modifying a branch before updating the remote.
        def race_with_transfer(repo, item):
            old = old_manifest["packages"][item["package"]]["sha"]
            tree = self.git(repo, "rev-parse", f"{old}^{{tree}}")
            actor = self.git(repo, "commit-tree", tree, "-p", old, "-m", "another actor's update")
            self.git(repo, "push", "--quiet", str(self.remotes / item["package"]), f"{actor}:{item['ref']}")
            other.append(actor)
            real_push(repo, item)
        with patch.object(release, "push_one", side_effect=race_with_transfer):
            self.assertEqual(self.invoke("--push", "--previous-manifest", str(previous)), 1)
        self.assertEqual(self.git(self.remotes / release.ORDER[0], "rev-parse", self.ref), other[0])
        self.assert_preserved()


if __name__ == "__main__":
    unittest.main()
