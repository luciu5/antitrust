# Patches for read-only package worktrees

These patches were generated after a per-block inventory, substantive review,
and deletion red-team pass. `git apply --check` passed in each source repo on
2026-10-01. They do not touch production files or existing `vertical` edits.

The user applied all four patches to the source worktrees on 2026-10-01. The
following is the reproducible application command for a clean checkout:

```bash
for entry in trade:refactor vertical:main coordination:main financial:main; do
  repo=${entry%%:*}; branch=${entry#*:}
  dir="/mnt/d/Projects/$repo"
  test "$(git -C "$dir" branch --show-current)" = "$branch" || exit 1
  git -C "$dir" apply --check "/mnt/d/Projects/antitrust/audit/testthat-2026-10-01/patches/$repo.patch" || exit 1
  git -C "$dir" apply "/mnt/d/Projects/antitrust/audit/testthat-2026-10-01/patches/$repo.patch" || exit 1
done
```

The temporary-copy runs passed: trade fast/extended/nightly 18.5/59.6/106.7 s,
vertical 9.5 s, coordination 23.1 s, all zero failures/errors. After application,
the actual trade fast suite, vertical `test_dir`, and coordination suite passed.
Financial's package cannot load here without `cubature`; its patch passed
syntax checks and standalone fixture calculations. See the top-level report
for exact verification commands and unresolved limits.
