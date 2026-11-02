#!/usr/bin/env python3
"""Statically validate and smoke-test the seven-package release workflows.

This guards material CI failures: duplicate YAML keys, invalid job-level GitHub
contexts, malformed Bash, Rscript before R setup, and certification paths that
accept mutable or missing release SHA inputs. Release selectors are exercised
with local event/environment stubs; no remote refs are fetched.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError as error:  # pragma: no cover - environment diagnostic
    raise SystemExit("PyYAML is required to validate GitHub workflow files") from error


PACKAGES = (
    "antitrust",
    "antitrustBayes",
    "competitiontoolbox",
    "coordination",
    "financial",
    "trade",
    "vertical",
)
JOB_ENV_CONTEXTS = {"github", "needs", "strategy", "matrix", "vars", "secrets", "inputs"}
WORKFLOW_ENV_CONTEXTS = {"github", "vars", "secrets", "inputs"}
EXPRESSION_RE = re.compile(r"\$\{\{(.*?)\}\}", re.DOTALL)
CONTEXT_RE = re.compile(r"\b([A-Za-z_][A-Za-z0-9_]*)\s*(?:\.|\[)")
SHA_RE = re.compile(r"^[0-9a-fA-F]{40}$")
SHA_GUARD_RE = re.compile(r"\{40\}")


class DuplicateKeyLoader(yaml.SafeLoader):
    """Safe YAML loader that rejects duplicate keys at every mapping depth."""


def construct_unique_mapping(loader: DuplicateKeyLoader, node: Any, deep: bool = False) -> dict:
    loader.flatten_mapping(node)
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        try:
            duplicate = key in result
        except TypeError as error:
            raise yaml.constructor.ConstructorError(
                "while constructing a mapping", node.start_mark,
                "found an unhashable mapping key", key_node.start_mark,
            ) from error
        if duplicate:
            raise yaml.constructor.ConstructorError(
                "while constructing a mapping", node.start_mark,
                f"found duplicate key {key!r}", key_node.start_mark,
            )
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


DuplicateKeyLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, construct_unique_mapping
)


def fail(message: str) -> None:
    raise ValueError(message)


def load_workflow(path: Path) -> dict[str, Any]:
    try:
        document = yaml.load(path.read_text(encoding="utf-8"), Loader=DuplicateKeyLoader)
    except (OSError, yaml.YAMLError) as error:
        fail(f"{path}: invalid or duplicate-key YAML: {error}")
    if not isinstance(document, dict) or not isinstance(document.get("jobs"), dict):
        fail(f"{path}: workflow must define a jobs mapping")
    return document


def event_config(document: dict[str, Any]) -> dict[str, Any]:
    # PyYAML's YAML 1.1 resolver treats the GitHub key "on" as boolean True.
    value = document.get("on", document.get(True, {}))
    return value if isinstance(value, dict) else {}


def strings(value: Any):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for nested in value.values():
            yield from strings(nested)
    elif isinstance(value, (list, tuple)):
        for nested in value:
            yield from strings(nested)


def check_context_mapping(mapping: Any, allowed: set[str], label: str) -> None:
    if not isinstance(mapping, dict):
        return
    for key, value in mapping.items():
        for text in strings(value):
            for expression in EXPRESSION_RE.findall(text):
                for match in CONTEXT_RE.finditer(expression):
                    context = match.group(1)
                    if context not in allowed:
                        fail(f"{label}.{key} uses unavailable GitHub context {context!r}")


def step_list(job: dict[str, Any]) -> list[dict[str, Any]]:
    steps = job.get("steps", [])
    return steps if isinstance(steps, list) else []


def validate_workflow(path: Path, document: dict[str, Any]) -> None:
    check_context_mapping(
        document.get("env"), WORKFLOW_ENV_CONTEXTS, f"{path}: workflow env"
    )
    for job_name, job in document["jobs"].items():
        if not isinstance(job, dict):
            fail(f"{path}: job {job_name!r} must be a mapping")
        check_context_mapping(
            job.get("env"), JOB_ENV_CONTEXTS, f"{path}: job {job_name} env"
        )
        steps = step_list(job)
        r_setup = [
            index for index, step in enumerate(steps)
            if str(step.get("uses", "")).startswith("r-lib/actions/setup-r@")
        ]
        for index, step in enumerate(steps):
            script = step.get("run")
            if not isinstance(script, str):
                continue
            shell = str(step.get("shell", "")).strip().lower()
            runs_on = " ".join(str(item) for item in (job.get("runs-on", "") if isinstance(job.get("runs-on", ""), list) else [job.get("runs-on", "")])).lower()
            bash_script = shell == "bash" or (not shell and ("ubuntu" in runs_on or "linux" in runs_on))
            if bash_script:
                result = subprocess.run(
                    ["bash", "-n"], input=script, text=True, capture_output=True, check=False
                )
                if result.returncode:
                    fail(f"{path}: job {job_name} step {index + 1} fails bash -n: {result.stderr.strip()}")
            calls_r = shell.startswith("rscript") or re.search(r"(^|[\s;&|])Rscript(?:\s|$)", script) is not None
            if calls_r and (not r_setup or r_setup[0] >= index):
                fail(f"{path}: job {job_name} step {index + 1} calls Rscript before setup-r")

    dispatch = event_config(document).get("workflow_dispatch", {})
    inputs = dispatch.get("inputs", {}) if isinstance(dispatch, dict) else {}
    release_inputs = [name for name in inputs if name.endswith("_release_sha")]
    if not release_inputs:
        return
    if "release_certification" not in inputs:
        fail(f"{path}: release SHA inputs exist without a release_certification dispatch input")

    bindings: list[tuple[str, str]] = []
    scripts: list[str] = []
    for job in document["jobs"].values():
        if not isinstance(job, dict):
            continue
        scopes = [document.get("env", {}), job.get("env", {})]
        steps = step_list(job)
        scopes.extend(step.get("env", {}) for step in steps if isinstance(step, dict))
        for scope in scopes:
            if isinstance(scope, dict):
                bindings.extend((str(key), str(value)) for key, value in scope.items())
        scripts.extend(
            step["run"] for step in steps
            if isinstance(step, dict) and isinstance(step.get("run"), str)
        )

    for field in release_inputs:
        package = field.removesuffix("_release_sha").upper().replace("-", "_")
        aliases = {
            name for name, value in bindings
            if name in {f"{package}_RELEASE_SHA", f"{package}_RELEASE_SHA_INPUT"}
            and f"inputs.{field}" in value
        }
        if not aliases:
            fail(f"{path}: {field} is not bound to a raw release SHA environment variable")
        guarded = any(
            "RELEASE_CERTIFICATION" in script
            and SHA_GUARD_RE.search(script)
            and any(re.search(rf"\b{re.escape(alias)}\b", script) for alias in aliases)
            for script in scripts
        )
        if not guarded:
            fail(f"{path}: {field} lacks a full-SHA certification guard on its raw input")


def parse_dcf(path: Path) -> dict[str, str]:
    fields: dict[str, str] = {}
    current: str | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        if line[0].isspace() and current:
            fields[current] += " " + line.strip()
        elif ":" in line:
            current, value = line.split(":", 1)
            fields[current.strip()] = value.strip()
    return fields


def description_remotes(package_root: Path) -> dict[str, str]:
    entries = parse_dcf(package_root / "DESCRIPTION").get("Remotes", "").split(",")
    result = {}
    for entry in filter(None, (item.strip() for item in entries)):
        match = re.fullmatch(r"luciu5/([^/@]+)@([0-9a-fA-F]{40})", entry)
        if not match or match.group(1) in result:
            fail(f"{package_root.name}: malformed or duplicate full-SHA Remote {entry!r}")
        result[match.group(1)] = match.group(2).lower()
    return result


def manifest_refs(package_root: Path) -> dict[str, str]:
    path = package_root / "ci/release-dependencies.json"
    if not path.is_file():
        return {}
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        fail(f"{path}: release dependencies must be a JSON object")
    return {str(key): str(sha).lower() for key, sha in value.items()}


def selector_expected(package: str, package_root: Path) -> dict[str, str] | None:
    remotes = description_remotes(package_root)
    pins = manifest_refs(package_root)
    if package == "antitrustBayes":
        return pins
    if package == "competitiontoolbox":
        return remotes
    if package == "coordination":
        return {"ref": remotes["antitrust"]}
    if package == "financial":
        return {"antitrust_ref": pins["antitrust"]} if all(
            SHA_RE.fullmatch(value) for value in pins.values()
        ) else None
    if package == "trade":
        return {key: remotes[key] for key in ("antitrust", "coordination")}
    if package == "vertical":
        combined = {"antitrust": remotes["antitrust"], "trade": remotes["trade"], **pins}
        return {
            "antitrust_ref": combined["antitrust"],
            "coordination_ref": combined["coordination"],
            "trade_ref": combined["trade"],
        } if all(SHA_RE.fullmatch(value) for value in combined.values()) else None
    return None


SELECTORS = {
    "antitrustBayes": ("qa.yml", "fast", "dependency-refs", "bash"),
    "competitiontoolbox": (
        "dependency-stack-qa.yml", "dependency-stack", "Initialize isolated library and provenance", "bash"
    ),
    "coordination": ("R-CMD-check.yaml", "check", "antitrust-ref", "bash"),
    "financial": ("qa.yml", "check", "antitrust-ref", "bash"),
    "trade": ("qa.yml", "check", "release-pins", "rscript"),
    "vertical": ("R-CMD-check.yaml", "check", "dependency-refs", "bash"),
}


def key_values(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key] = value
    return result


def selector_environment(package: str, temp: Path) -> dict[str, str]:
    env = os.environ.copy()
    env.update({
        "RUNNER_TEMP": str(temp),
        "GITHUB_OUTPUT": str(temp / "output"),
        "GITHUB_ENV": str(temp / "env"),
        "GITHUB_STEP_SUMMARY": str(temp / "summary"),
        "GITHUB_EVENT_NAME": "push",
        "GITHUB_REF": "refs/heads/release/validator-smoke",
        "GITHUB_HEAD_REF": "",
        "GITHUB_BASE_REF": "",
        "RELEASE_CERTIFICATION": "false",
        "RELEASE_PIN_BUILD": "true",
        "AUTOMATIC_RELEASE_PINS": "true",
        "ANTITRUST_RELEASE_SHA": "",
        "ANTITRUST_RELEASE_SHA_INPUT": "",
        "COORDINATION_RELEASE_SHA": "",
        "COORDINATION_RELEASE_SHA_INPUT": "",
        "TRADE_RELEASE_SHA": "",
        "VERTICAL_RELEASE_SHA": "",
        "ANTITRUST_REQUESTED_REF": "master",
        "COORDINATION_REQUESTED_REF": "main",
        "TRADE_REQUESTED_REF": "master",
        "ANTITRUST_REF": "master",
        "COORDINATION_REF": "main",
        "TRADE_REF": "master",
        "VERTICAL_REF": "main",
        "ANTITRUST_DEVELOPMENT_REF": "master",
        "TEST_TIER": "FAST",
        "RELEASE_BRANCH_BUILD": "true",
    })
    if package == "coordination":
        env.update({
            "EVENT_NAME": "push",
            "EVENT_REF": env["GITHUB_REF"],
            "EVENT_HEAD_REF": "",
            "EVENT_BASE_REF": "",
        })
    if package == "competitiontoolbox":
        bin_dir = temp / "bin"
        bin_dir.mkdir()
        stub = bin_dir / "git"
        stub.write_text("#!/bin/sh\nprintf '" + ("0" * 40) + "\\n'\n", encoding="utf-8")
        stub.chmod(0o755)
        env["PATH"] = str(bin_dir) + os.pathsep + env.get("PATH", "")
    return env


def smoke_selector(sources_root: Path, package: str) -> str:
    package_root = sources_root / package
    workflow_name, job_name, selector, shell_kind = SELECTORS[package]
    document = load_workflow(package_root / ".github/workflows" / workflow_name)
    steps = step_list(document["jobs"][job_name])
    step = next(
        (item for item in steps if item.get("id") == selector or item.get("name") == selector),
        None,
    )
    if not step or not isinstance(step.get("run"), str):
        fail(f"{package}: automatic release selector {selector!r} is missing")
    expected = selector_expected(package, package_root)
    if shell_kind == "rscript" and not shutil.which("Rscript"):
        return f"{package}: skipped selector execution (Rscript is unavailable)"

    with tempfile.TemporaryDirectory(prefix="release-selector-smoke-") as temp_name:
        temp = Path(temp_name)
        env = selector_environment(package, temp)
        script = step["run"]
        script = script.replace("${{ github.run_id }}", "validator-run")
        script = script.replace("${{ github.run_attempt }}", "1")
        command = ["Rscript", "-"] if shell_kind == "rscript" else ["bash", "-e"]
        result = subprocess.run(
            command, input=script, cwd=package_root, env=env,
            text=True, capture_output=True, timeout=30, check=False,
        )
        transcript = result.stdout + result.stderr
        outputs = key_values(Path(env["GITHUB_OUTPUT"]))
        exported = key_values(Path(env["GITHUB_ENV"]))
        if expected is None:
            if result.returncode == 0:
                fail(f"{package}: release selector accepted a mutable or pending pin")
            return f"{package}: pending release manifest fails closed as expected"
        if result.returncode:
            fail(f"{package}: release selector smoke failed: {transcript.strip()}")
        actual = exported if package == "competitiontoolbox" else outputs
        actual_keys = {
            "antitrustBayes": {"antitrust_ref": "antitrust", "coordination_ref": "coordination", "trade_ref": "trade"},
            "competitiontoolbox": {f"{name.upper()}_REF": name for name in expected},
            "coordination": {"ref": "ref"},
            "financial": {"antitrust_ref": "antitrust_ref"},
            "trade": {"antitrust": "antitrust", "coordination": "coordination"},
            "vertical": {"antitrust_ref": "antitrust_ref", "coordination_ref": "coordination_ref", "trade_ref": "trade_ref"},
        }[package]
        for output_key, expected_key in actual_keys.items():
            if actual.get(output_key, "").lower() != expected[expected_key].lower():
                fail(
                    f"{package}: release selector emitted {output_key}={actual.get(output_key)!r}; "
                    f"expected the exact local pin {expected[expected_key]}"
                )
    return f"{package}: release selector chose exact staged pins"


def validate_all(sources_root: Path) -> None:
    workflows = []
    for package in PACKAGES:
        package_root = sources_root / package
        if not package_root.is_dir():
            fail(f"missing staged source root: {package_root}")
        paths = sorted(
            path for path in (package_root / ".github/workflows").glob("*")
            if path.suffix in {".yml", ".yaml"}
        )
        if not paths:
            fail(f"{package}: no GitHub workflow YAML files found")
        for path in paths:
            document = load_workflow(path)
            validate_workflow(path, document)
            workflows.append(path)
    print(f"Static workflow checks passed for {len(workflows)} files across {len(PACKAGES)} roots")
    smoke_errors = []
    for package in SELECTORS:
        try:
            print(smoke_selector(sources_root, package))
        except (OSError, ValueError, subprocess.SubprocessError, json.JSONDecodeError) as error:
            smoke_errors.append(f"{package}: {error}")
            print(f"{package}: selector smoke failed: {error}", file=sys.stderr)
    if smoke_errors:
        fail("automatic selector smoke failures: " + "; ".join(smoke_errors))


def main() -> int:
    default_root = Path(__file__).resolve().parents[3]
    sources_root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else default_root
    try:
        validate_all(sources_root)
    except (OSError, ValueError, subprocess.SubprocessError, json.JSONDecodeError) as error:
        print(f"release workflow validation failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
