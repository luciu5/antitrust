#!/usr/bin/env bash
set -euo pipefail

# Commit only the reviewed test-audit edits. Existing vertical production edits
# and the untracked trade artifact are deliberately excluded.
for spec in trade:refactor vertical:main coordination:main antitrustBayes:main financial:main antitrust:refactor; do
  repo=${spec%%:*}
  branch=${spec#*:}
  dir="/mnt/d/Projects/$repo"
  test "$(git -C "$dir" branch --show-current)" = "$branch"
  git -C "$dir" diff --cached --quiet
  git -C "$dir" diff --check
done

commit_audit() {
  local repo=$1
  local message=$2
  shift 2
  local dir="/mnt/d/Projects/$repo"
  git -C "$dir" add -A -- "$@"
  if [[ "$repo" == antitrust ]]; then
    # A raw .patch records blank context lines as a leading space. Git's
    # whitespace check flags those patch-text lines even though the source
    # diffs themselves pass. Check every other staged audit file here.
    git -C "$dir" diff --cached --check -- . \
      ':(exclude)audit/testthat-2026-10-01/patches/*.patch'
  else
    git -C "$dir" diff --cached --check
  fi
  git -C "$dir" commit -m "$message"
}

commit_audit trade \
  "test: tier trade validation and consolidate contract checks" \
  .github/workflows/qa.yml tests/testthat

commit_audit vertical \
  "test: remove obsolete vertical parity checks" \
  tests/testthat/test-parity-fixtures.R \
  tests/testthat/test-subset-support.R

commit_audit coordination \
  "test: remove redundant coordination smokes" \
  tests/testthat/test-price-leadership-blp.R \
  tests/testthat/test-price-leadership-invariants.R \
  tests/testthat/test-s4-dispatch.R

commit_audit antitrustBayes \
  "test: tier Bayesian fixtures and repair Stan smoke" \
  .github/workflows/qa.yml \
  doc/market_intercept_enterprise/README.md \
  doc/market_intercept_enterprise/fixture_checks.R \
  doc/test_audit_inventory.csv \
  doc/test_suite_audit.md \
  tests/testthat

commit_audit financial \
  "test: move expensive financial parity checks to extended validation" \
  tests/testthat/test-bank-parity.R \
  tests/testthat/test-debt-parity.R \
  tests/testthat/test-skeleton.R

commit_audit antitrust \
  "test: rationalize ecosystem suites and fix auction owner handling" \
  AGENTS.md R/Retention.R TEST_SUITE_AUDIT.md audit tests/testthat
