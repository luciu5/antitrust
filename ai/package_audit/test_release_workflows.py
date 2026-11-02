"""Prevent uncertified sources, duplicate dispatches and false CI success."""

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
            item["run_name"] = workflows.run_name(item["inputs"])

    def api(self, path, payload=None):
        name = path.split("/")[2]
        item = next(p for p in self.plan if p["package"] == name)
        if payload is not None:
            self.posts.append((path, payload))
            if path.endswith("/dispatches"):
                self.assertEqual(payload["inputs"], item["inputs"])
                self.assertEqual(payload["ref"], item["ref"].removeprefix("refs/heads/"))
                self.runs[name].append(dict(id=len(self.runs[name]) + 1, head_sha=item["sha"],
                    event="workflow_dispatch", display_title=item["run_name"],
                    workflow_id=push.ORDER.index(name) + 1, head_branch=item["ref"].removeprefix("refs/heads/"),
                    path=".github/workflows/" + item["workflow"] + "@" + item["ref"].removeprefix("refs/heads/"),
                    html_url=f"https://github.com/luciu5/{name}/actions/runs/1", status="completed", conclusion="success", run_attempt=1))
                return {"workflow_run_id": self.runs[name][-1]["id"], "html_url": self.runs[name][-1]["html_url"]}
            elif path.endswith("/rerun-failed-jobs"):
                run = self.runs[name][-1]
                run.update(run_attempt=run["run_attempt"] + 1, conclusion="success")
            else:
                raise AssertionError("Unexpected mutation " + path)
            return {}
        if "/git/ref/heads/" in path:
            return {"object": {"sha": "f" * 40 if name == self.bad_branch else item["sha"]}}
        if "/runs?" in path:
            return {"workflow_runs": self.runs[name]}
        if "/actions/runs/" in path:
            run_id = int(path.rsplit("/", 1)[1])
            return next(run for run in self.runs[name] if run["id"] == run_id)
        if "/actions/workflows/" in path:
            return {"id": push.ORDER.index(name) + 1, "state": "active", "path": ".github/workflows/" + item["workflow"]}
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

    def test_exact_full_tiers_dispatch_once_and_reuse_successful_runs(self):
        """Prevent accidental FAST/bootstrap certification and duplicate run costs."""
        self.assertTrue(self.certify())
        self.assertEqual(len(self.posts), 7)
        self.assertEqual(self.plan[0]["inputs"], {"test_tier": "nightly"})
        self.assertEqual(self.plan[4]["inputs"]["test_tier"], "NIGHTLY")
        self.assertIs(self.plan[5]["inputs"]["extended"], True)
        self.assertIs(self.plan[6]["inputs"]["run_stan_tests"], False)
        self.assertTrue(self.certify())
        self.assertEqual(len(self.posts), 7)
        self.f.assert_preserved()

    def test_last_branch_mismatch_prevents_all_dispatches(self):
        """Prevent partial certification of a stack with an unreviewed remote source."""
        self.bad_branch = push.ORDER[-1]
        with self.assertRaises(push.PushError):
            self.certify()
        self.assertEqual(self.posts, [])

    def test_bootstrap_run_is_not_mistaken_for_exact_certification(self):
        """Prevent a successful workflow with wrong dependency inputs being reused."""
        self.assertTrue(self.certify())
        for name in push.ORDER:
            self.runs[name][0]["display_title"] = "bootstrap development refs"
        for item in self.plan:
            item.pop("dispatch_receipt", None)
        self.posts.clear()
        self.assertTrue(self.certify())
        self.assertEqual(len(self.posts), 7)

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

    def test_lost_dispatch_acknowledgement_reuses_accepted_run_on_restart(self):
        """Prevent duplicate paid CI after GitHub accepted a dispatch before timeout."""
        def lose_ack(path, payload=None):
            value = self.api(path, payload)
            if payload is not None:
                raise push.PushError("simulated lost dispatch acknowledgement")
            return value
        with patch.object(workflows, "gh_api", side_effect=lose_ack), contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(push.PushError):
                workflows.certify(self.plan, poll_seconds=0, timeout_seconds=2)
        self.assertEqual(self.plan[0]["status"], "dispatch-outcome-unknown")
        self.assertTrue(self.certify())
        self.assertEqual(len(self.posts), 7)


if __name__ == "__main__":
    unittest.main()
