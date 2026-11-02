#!/usr/bin/env python3
"""Check an exact, local coordinated-release candidate manifest.

The checker reads candidate commit objects with ``git show`` and never checks
out a branch, modifies a repository, or contacts a remote. The manifest itself
must remain outside the antitrust candidate commit to avoid a self-referential
core SHA.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any


PACKAGE_NAMES = (
    "antitrust",
    "coordination",
    "trade",
    "vertical",
    "competitiontoolbox",
    "financial",
    "antitrustBayes",
)
SHA_RE = re.compile(r"^[0-9a-fA-F]{40}$")
VERSION_RE = re.compile(r"^[0-9]+(?:\.[0-9]+)*$")
DESCRIPTION_DEPENDENCY_FIELDS = ("Depends", "Imports", "LinkingTo", "Suggests")


class ManifestError(Exception):
    """Raised for an invalid manifest or candidate tree."""


def git(repo: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(repo), *args],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode:
        detail = result.stderr.strip() or result.stdout.strip()
        raise ManifestError(f"git {' '.join(args)} failed in {repo}: {detail}")
    return result.stdout.rstrip("\n")


def parse_dcf(raw: str, source: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    current: str | None = None
    for line_number, line in enumerate(raw.splitlines(), start=1):
        if not line.strip():
            continue
        if line[0].isspace():
            if current is None:
                raise ManifestError(f"{source}:{line_number}: continuation without a field")
            fields[current] += " " + line.strip()
            continue
        if ":" not in line:
            raise ManifestError(f"{source}:{line_number}: malformed DESCRIPTION line")
        current, value = line.split(":", 1)
        current = current.strip()
        if not current or current in fields:
            raise ManifestError(f"{source}:{line_number}: empty or duplicate field {current!r}")
        fields[current] = value.strip()
    return fields


def version_tuple(value: str, label: str) -> tuple[int, ...]:
    if not VERSION_RE.fullmatch(value):
        raise ManifestError(f"{label} must be a numeric R package version, got {value!r}")
    return tuple(int(part) for part in value.split("."))


def parse_package_minima(description: dict[str, str], package: str) -> dict[str, str]:
    minima: dict[str, str] = {}
    for field in DESCRIPTION_DEPENDENCY_FIELDS:
        value = description.get(field, "")
        if not value.strip():
            continue
        for clause in value.split(","):
            clause = clause.strip()
            match = re.match(
                r"^([A-Za-z][A-Za-z0-9.]*)\s*(?:\(\s*([^\s)]+)\s+([^\s)]+)\s*\))?$",
                clause,
            )
            if not match:
                raise ManifestError(
                    f"{package} DESCRIPTION has an unrecognized {field} entry: {clause!r}"
                )
            dependency, operator, version = match.groups()
            if dependency in PACKAGE_NAMES and operator is not None:
                if operator != ">=":
                    raise ManifestError(
                        f"{package} DESCRIPTION uses {operator} for internal dependency "
                        f"{dependency}; release minima must use >="
                    )
                if dependency in minima:
                    raise ManifestError(
                        f"{package} DESCRIPTION declares internal dependency "
                        f"{dependency} more than once"
                    )
                minima[dependency] = version
    return minima


def parse_remotes(description: dict[str, str], package: str) -> dict[str, str]:
    value = description.get("Remotes", "")
    remotes: dict[str, str] = {}
    for entry in filter(None, (part.strip() for part in value.split(","))):
        match = re.fullmatch(r"luciu5/([A-Za-z][A-Za-z0-9.]*)@([^\s,]+)", entry)
        if not match:
            raise ManifestError(
                f"{package} DESCRIPTION Remote must use luciu5/<package>@<full SHA>: {entry!r}"
            )
        dependency, ref = match.groups()
        if dependency in remotes:
            raise ManifestError(f"{package} DESCRIPTION repeats Remote {dependency}")
        if dependency not in PACKAGE_NAMES:
            raise ManifestError(f"{package} DESCRIPTION Remote {dependency} is outside the release stack")
        if not SHA_RE.fullmatch(ref):
            raise ManifestError(
                f"{package} DESCRIPTION Remote {dependency} must be a full 40-character SHA, "
                f"got {ref!r}"
            )
        remotes[dependency] = ref.lower()
    return remotes


def load_manifest(path: Path) -> dict[str, Any]:
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ManifestError(f"cannot read manifest {path}: {error}") from error
    if not isinstance(manifest, dict):
        raise ManifestError("manifest root must be a JSON object")
    if manifest.get("schema_version") != 1:
        raise ManifestError("schema_version must equal 1")
    if not isinstance(manifest.get("release_id"), str) or not manifest["release_id"].strip():
        raise ManifestError("release_id must be a nonempty string")
    packages = manifest.get("packages")
    if not isinstance(packages, dict) or set(packages) != set(PACKAGE_NAMES):
        missing = sorted(set(PACKAGE_NAMES) - set(packages or {}))
        extra = sorted(set(packages or {}) - set(PACKAGE_NAMES))
        raise ManifestError(f"packages must contain the seven release packages; missing={missing}, extra={extra}")

    for name in PACKAGE_NAMES:
        candidate = packages[name]
        if not isinstance(candidate, dict):
            raise ManifestError(f"packages.{name} must be an object")
        if candidate.get("repository") != f"luciu5/{name}":
            raise ManifestError(f"packages.{name}.repository must equal 'luciu5/{name}'")
        if not isinstance(candidate.get("version"), str):
            raise ManifestError(f"packages.{name}.version must be a string")
        version_tuple(candidate["version"], f"packages.{name}.version")
        sha = candidate.get("sha")
        if not isinstance(sha, str) or not SHA_RE.fullmatch(sha):
            raise ManifestError(
                f"packages.{name}.sha must be an explicit full 40-character commit SHA; "
                f"mutable refs such as main/refactor and pending placeholders are rejected"
            )
        candidate["sha"] = sha.lower()
        minimums = candidate.get("minimums", {})
        remotes = candidate.get("remotes", [])
        if not isinstance(minimums, dict) or not all(
            isinstance(dep, str) and dep in PACKAGE_NAMES and isinstance(version, str)
            for dep, version in minimums.items()
        ):
            raise ManifestError(f"packages.{name}.minimums must map internal package names to versions")
        for dependency, minimum in minimums.items():
            version_tuple(minimum, f"packages.{name}.minimums.{dependency}")
        if not isinstance(remotes, list) or not all(
            isinstance(dep, str) and dep in PACKAGE_NAMES for dep in remotes
        ):
            raise ManifestError(f"packages.{name}.remotes must list internal package names")
        if len(set(remotes)) != len(remotes):
            raise ManifestError(f"packages.{name}.remotes contains duplicates")
    return manifest


def validate_candidate(
    package: str,
    candidate: dict[str, Any],
    all_candidates: dict[str, Any],
    repositories_root: Path,
) -> None:
    repository = repositories_root / package
    if not repository.is_dir():
        raise ManifestError(f"local repository directory is missing for {package}: {repository}")
    expected_sha = candidate["sha"]
    resolved_sha = git(repository, "rev-parse", "--verify", f"{expected_sha}^{{commit}}").lower()
    if resolved_sha != expected_sha:
        raise ManifestError(f"{package} SHA resolved to {resolved_sha}, expected {expected_sha}")

    description_raw = git(repository, "show", f"{expected_sha}:DESCRIPTION")
    description = parse_dcf(description_raw, f"{package}@{expected_sha}:DESCRIPTION")
    actual_version = description.get("Version")
    if actual_version != candidate["version"]:
        raise ManifestError(
            f"{package}@{expected_sha} DESCRIPTION version is {actual_version!r}; "
            f"manifest expects {candidate['version']!r}"
        )

    actual_minimums = parse_package_minima(description, package)
    expected_minimums = candidate.get("minimums", {})
    if actual_minimums != expected_minimums:
        raise ManifestError(
            f"{package}@{expected_sha} dependency minima differ: "
            f"DESCRIPTION={actual_minimums}, manifest={expected_minimums}"
        )
    for dependency, minimum in actual_minimums.items():
        candidate_version = all_candidates[dependency]["version"]
        if version_tuple(candidate_version, f"packages.{dependency}.version") < version_tuple(
            minimum, f"{package} minimum for {dependency}"
        ):
            raise ManifestError(
                f"{package} requires {dependency} >= {minimum}, but the selected candidate is "
                f"{candidate_version}"
            )

    actual_remotes = parse_remotes(description, package)
    expected_remotes = candidate.get("remotes", [])
    if set(actual_remotes) != set(expected_remotes):
        raise ManifestError(
            f"{package}@{expected_sha} Remotes entries differ: "
            f"DESCRIPTION={sorted(actual_remotes)}, manifest={sorted(expected_remotes)}"
        )
    for dependency, ref in actual_remotes.items():
        expected_dependency_sha = all_candidates[dependency]["sha"]
        if ref != expected_dependency_sha:
            raise ManifestError(
                f"{package}@{expected_sha} Remote {dependency} pins {ref}; "
                f"manifest candidate is {expected_dependency_sha}"
            )

    ci_path = "ci/release-dependencies.json"
    tracked_paths = set(git(repository, "ls-tree", "-r", "--name-only", expected_sha).splitlines())
    if ci_path in tracked_paths:
        try:
            ci_pins = json.loads(git(repository, "show", f"{expected_sha}:{ci_path}"))
        except json.JSONDecodeError as error:
            raise ManifestError(f"{package} has malformed {ci_path}: {error}") from error
        if not isinstance(ci_pins, dict) or not ci_pins:
            raise ManifestError(f"{package} {ci_path} must contain dependency SHA mappings")
        for dependency, ref in ci_pins.items():
            if dependency not in PACKAGE_NAMES or dependency == package:
                raise ManifestError(f"{package} {ci_path} has an invalid dependency: {dependency}")
            if ref != all_candidates[dependency]["sha"]:
                raise ManifestError(f"{package} {ci_path} does not match the frozen {dependency} SHA")


def validate_manifest_not_in_core_commit(
    manifest_path: Path,
    antitrust_repository: Path,
    antitrust_sha: str,
) -> None:
    resolved_manifest = manifest_path.resolve()
    resolved_repository = antitrust_repository.resolve()
    try:
        relative_path = resolved_manifest.relative_to(resolved_repository).as_posix()
    except ValueError:
        return
    tracked_paths = set(git(antitrust_repository, "ls-tree", "-r", "--name-only", antitrust_sha).splitlines())
    if relative_path in tracked_paths:
        raise ManifestError(
            f"release manifest {relative_path} is committed in the antitrust candidate; "
            "keep the exact-SHA manifest outside the core commit to avoid self-reference"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="exact coordinated-release manifest JSON")
    parser.add_argument(
        "--repositories-root",
        type=Path,
        default=Path(__file__).resolve().parents[3],
        help="directory containing the seven local Git repositories (default: package workspace parent)",
    )
    args = parser.parse_args()
    try:
        manifest = load_manifest(args.manifest)
        repositories_root = args.repositories_root.resolve()
        packages = manifest["packages"]
        for name in PACKAGE_NAMES:
            validate_candidate(name, packages[name], packages, repositories_root)
        validate_manifest_not_in_core_commit(
            args.manifest,
            repositories_root / "antitrust",
            packages["antitrust"]["sha"],
        )
    except ManifestError as error:
        print(f"release manifest check failed: {error}", file=sys.stderr)
        return 1
    print(f"release manifest {manifest['release_id']} is valid")
    for name in PACKAGE_NAMES:
        candidate = manifest["packages"][name]
        print(f"{name} {candidate['version']} {candidate['sha']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
