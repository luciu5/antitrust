"""Prevent uncertified sources, default-branch mutations and false CI success."""

import contextlib
import io
import json
import unittest
from unittest.mock import patch

import push_coordinated_release as push
import release_workflows as workflows
import test_push_coordinated_release as fixtures


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.f = fixtures.ReleasePushTests(methodName="runTest")
        self.f.setUp()
        self.addCleanup(self.f.doCleanups)
        argv = ["--repositories-root", str(self.f.local)]
        for name, path in self.f.paths.items():
            argv.extend(["--" + name, str(path)])
        self.argv = argv + ["--previous-manifest", str(self.f.paths["manifest"])]
        self.plan = push.local_plan(push.parser().parse_args(self.argv))
        self.runs = {name: [] for name in push.ORDER}
        self.posts = []
        self.bad_branch = None
        for item in self.plan:
            name = item["package"]
            item.update(workflow=workflows.WORKFLOWS[name][0], inputs=workflows.certification_inputs(name, self.f.manifest["packages"]))
            item["run_name"] = "release-certification-" + item["sha"]
            self.runs[name].append(dict(id=push.ORDER.index(name) + 1, head_sha=item["sha"],
                event="push", display_title=item["run_name"], workflow_id=push.ORDER.index(name) + 100,
                head_branch=item["ref"].removeprefix("refs/heads/"),
                path=".github/workflows/" + item["workflow"] + "@" + item["ref"].removeprefix("refs/heads/"),
                html_url=f"https://github.com/luciu5/{name}/actions/runs/1", status="completed", conclusion="success", run_attempt=1))

    def api(self, path, payload=None):
        name = path.split("/")[2]
        item = next(p for p in self.plan if p["package"] == name)
        if payload is not None:
            self.posts.append((path, payload))
            if path.endswith("/rerun-failed-jobs"):
                run = self.runs[name][-1]
                run.update(run_attempt=run["run_attempt"] + 1, conclusion="success")
            else:
                raise AssertionError("Unexpected mutation " + path)
            return {}
        if "/git/ref/heads/" in path:
            return {"object": {"sha": "f" * 40 if name == self.bad_branch else item["sha"]}}
        if "/runs?" in path:
            return {"workflow_runs": self.runs[name]}
        if "/actions/workflows/" in path:
            raise AssertionError("No default-branch workflow lookup is permitted")
        raise AssertionError("Unexpected read " + path)

    def certify(self, **kwargs):
        with patch.object(workflows, "gh_api", side_effect=self.api), contextlib.redirect_stdout(io.StringIO()):
            return workflows.certify(self.plan, poll_seconds=0, timeout_seconds=2, **kwargs)

    def test_offline_plan_never_calls_github_or_pushes(self):
        """Prevent default invocation from mutating GitHub or touching the network."""
        with patch.object(workflows, "gh_api", side_effect=AssertionError("network")), \
             patch.object(push, "main", side_effect=AssertionError("upload")), \
             contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(workflows.main(self.argv), 0)

    def test_exact_push_runs_need_no_dispatch_or_default_branch_workflow(self):
        """Certify release-only workflows without changing protected default branches."""
        self.assertTrue(self.certify())
        self.assertEqual(self.posts, [])
        self.assertEqual(self.plan[0]["inputs"], {"test_tier": "nightly"})
        self.assertEqual(self.plan[4]["inputs"]["test_tier"], "NIGHTLY")
        self.assertIs(self.plan[5]["inputs"]["extended"], True)
        self.assertIs(self.plan[6]["inputs"]["run_stan_tests"], False)
        self.assertTrue(self.certify())
        self.assertEqual(self.posts, [])
        self.f.assert_preserved()

    def test_last_branch_mismatch_prevents_all_retries(self):
        """Prevent remote mutations when any release source differs from review."""
        self.bad_branch = push.ORDER[-1]
        with self.assertRaises(push.PushError):
            self.certify()
        self.assertEqual(self.posts, [])

    def test_bootstrap_or_wrong_source_run_is_not_mistaken_for_certification(self):
        """Prevent false green from a successful run with wrong source or tier marker."""
        run = self.runs[push.ORDER[0]][0]
        for field, value in [("display_title", "bootstrap development refs"), ("event", "workflow_dispatch"), ("head_sha", "f" * 40)]:
            original = run[field];run[field] = value
            with patch.object(workflows, "gh_api", side_effect=self.api), contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaises(push.PushError):
                    workflows.certify(self.plan, poll_seconds=0, timeout_seconds=0)
            run[field] = original
        self.assertEqual(self.posts, [])

    def test_failed_run_is_not_reported_green_and_retry_tracks_new_attempt(self):
        """Prevent false success or an endless wait after a very fast failed-job retry."""
        self.assertTrue(self.certify())
        self.runs[push.ORDER[2]][0]["conclusion"] = "failure"
        self.posts.clear()
        self.assertFalse(self.certify())
        self.assertEqual(self.posts, [])
        self.assertTrue(self.certify(retry_failed=True))
        self.assertEqual(len(self.posts), 1)
        self.assertTrue(self.posts[0][0].endswith("/rerun-failed-jobs"))
        self.assertEqual(self.plan[2]["run_attempt"], 2)

    def test_lost_retry_acknowledgement_reuses_completed_attempt_on_restart(self):
        """Prevent duplicate paid CI after GitHub accepted a failed-job retry."""
        self.runs[push.ORDER[0]][0]["conclusion"] = "failure"
        def lose_ack(path, payload=None):
            value = self.api(path, payload)
            if payload is not None:
                raise push.PushError("simulated lost retry acknowledgement")
            return value
        with patch.object(workflows, "gh_api", side_effect=lose_ack), contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(push.PushError):
                workflows.certify(self.plan, retry_failed=True, poll_seconds=0, timeout_seconds=2)
        self.assertEqual(self.plan[0]["status"], "retry-outcome-unknown")
        self.assertTrue(self.certify())
        self.assertEqual(len(self.posts), 1)

    def toolbox_waiver(self):
        item = next(p for p in self.plan if p["package"] == "competitiontoolbox")
        return {"certification_waivers": {"competitiontoolbox": {
            "candidate_sha": item["sha"], "workflow": ".github/workflows/" + item["workflow"],
            "authorized_by": "user", "reason": "Skip the private vertical dependency certification."}}}

    def test_waived_private_workflow_is_not_queried_retried_or_counted_passed(self):
        """Prevent a waived private-dependency job from requiring secrets or claiming success."""
        workflows.apply_waivers(self.plan, self.toolbox_waiver())
        self.runs["competitiontoolbox"][0]["conclusion"] = "failure"
        def without_private_actions(path, payload=None):
            if "/competitiontoolbox/actions/" in path:
                raise AssertionError("Waived workflow must not be read or retried")
            return self.api(path, payload)
        with patch.object(workflows, "gh_api", side_effect=without_private_actions), contextlib.redirect_stdout(io.StringIO()):
            self.assertTrue(workflows.certify(self.plan, retry_failed=True, poll_seconds=0, timeout_seconds=2))
        item = next(p for p in self.plan if p["package"] == "competitiontoolbox")
        self.assertEqual(item["status"], "waived")
        self.assertIsNone(item["conclusion"])
        self.assertEqual(self.posts, [])

    def test_waiver_cannot_hide_other_failures_or_cover_another_commit(self):
        """Prevent a one-workflow exception from hiding a Bayes failure or future source change."""
        data = self.toolbox_waiver()
        workflows.apply_waivers(self.plan, data)
        self.runs["antitrustBayes"][0]["conclusion"] = "failure"
        self.assertFalse(self.certify())
        data["certification_waivers"]["competitiontoolbox"]["candidate_sha"] = "f" * 40
        with self.assertRaises(push.PushError):
            workflows.apply_waivers(self.plan, data)
        with self.assertRaises(push.PushError):
            workflows.apply_waivers(self.plan, {"certification_waivers": {"antitrustBayes": {}}})

    def test_completed_report_distinguishes_waived_from_successful_certification(self):
        """Prevent a zero exit status from being reported as seven successful hosted checks."""
        self.f.handoff.update(self.toolbox_waiver())
        self.f.paths["handoff"].write_text(json.dumps(self.f.handoff))
        report = self.f.root / "waived-report.json"
        with patch.object(workflows, "gh_api", side_effect=self.api), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(workflows.main(self.argv + ["--certify", "--report", str(report)]), 0)
        data = json.loads(report.read_text())
        self.assertEqual(data["status"], "complete_with_waivers")
        self.assertEqual(data["waived_packages"], ["competitiontoolbox"])
        item = next(p for p in data["packages"] if p["package"] == "competitiontoolbox")
        self.assertEqual(item["status"], "waived")
        self.assertIsNone(item["conclusion"])


if __name__ == "__main__":
    unittest.main()
