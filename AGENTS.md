Work only on refactor; never modify/merge/rebase master.
Preserve uncommitted local changes.
Legacy S4 economics are the behavioral source of truth.
calibrate, update, respecify, counterfactual, and simulate have distinct semantics.
Registries list complete implemented models, not arbitrary demand × conduct combinations.
respecify() must not secretly recalibrate.
Synthetic-market code must respect multi-product ownership and test FOCs.
Test economic behavior, not just object construction.
trade depends on/reuses antitrust; synthetic-market design and realization are antitrust-owned.
Report tests and substantive changes; never claim a test was run when it wasn't.
Every permanent test must name a plausible material economic, numerical, API, or
cross-package failure that other tests would miss. Prefer independent FOCs and
identities to snapshots or construction checks. Label migration parity A
(permanent), B (temporary with a sunset), or C (legacy trap). Deep-test distinct
economic logic; use representative smoke checks for syntactic variants. Route
expensive statistical validation to an explicit extended tier.
