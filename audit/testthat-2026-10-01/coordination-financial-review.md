# Coordination and financial test audit (2026-10-01)

## Scope and evidence

Read `coordination` and `financial` current `main` branches. Inventory of every `test_that()` block, classification, observed runtime (coordination), proposed tier, and deletion red-team rationale: `inventory_before.csv`. Read implementation for Stackelberg Logit/CES, core-fringe, PLE/BLP, Grim Trigger, revenue retention, signed rates, financial debt/bank equilibrium, liquidity risk, posterior mapping, and shared Beta moments. Both source repositories remained untouched because they are outside the configured writable roots.

| Package | Before | Proposed patch | Routine active blocks | Baseline evidence |
|---|---:|---:|---:|---|
| coordination | 12 files, 78 blocks | 12 files, 75 blocks | 75 | `devtools::test()` against loaded antitrust refactor: 0 failures, 24 warnings, 26.1 s wall / 21.9 s summed block time |
| financial | 7 files, 24 blocks | 6 files, 23 blocks | 21; 2 gated extended | `devtools::test()` could not load: missing `cubature`. CRAN install failed both sandboxed and approved due unreachable package index. Runtime unmeasured. Standalone Beta integral and 275-row tail fixture passed. |

The coordination patch was applied only to a temporary copy and retested: 75 blocks, 0 failures/errors, 19 warnings, 23.1 s wall. Ordinary runtime improvement is small and noisy, so no stable speedup is claimed. The financial patch parsed and `git apply --check` passed; its suite was not run. Patch files are `patches/coordination.patch` and `patches/financial.patch` (apply from each corresponding repository root with `git apply`).

## Package conclusions

### Coordination

The suite is relatively large for the package (78 blocks), but much of that mass protects genuinely distinct economics. The strongest tests independently compute Stackelberg follower reactions and leader reduced-profit gradients, core-fringe endpoint formulas, retention-weighted portfolio profit FOCs, and a BLP fringe firm's finite-difference profit derivative. Signed and zero input rates catch plausible price-domain mistakes; real-rate tests are unusually valuable. The antitrust endpoint comparisons are A-grade permanent integration benchmarks, though they belong in full CI. The hottest files are Stackelberg Logit (8.0 s), PLE BLP (3.5 s), Stackelberg CES (2.9 s), and retention (2.8 s) by summed per-block time. Runtime is reasonable; no broad migration to extended validation is justified.

Weak areas: repeated finite-output constructor checks, integration-node slot assertions that duplicate each other, S4 metadata introspection, and Grim Trigger tests that recompute the same fit repeatedly while failing to assert the firm-level repeated-game inequality independently. The suite often checks residuals returned by the implementation; the best tests compare against direct profit gradients. Retain that independent pattern. A future Grim Trigger test should calculate the firm-level threshold from reported Coord/Defect/Punish payoffs and discount factors and compare to `IC` for both pre and post ownership, especially multiproduct firms. Revenue-retention exit should check the retained portfolio FOC, not merely finite prices and an ownership matrix.

Tier 1: representative Logit/CES Stackelberg FOCs, exact core-fringe formulas and role FOCs, retention-weighted profits, Grim Trigger IC, PLE fringe FOC, signed-rate edge, public constructors and dispatch. Tier 2: broader gamma/conduct combinations, antitrust endpoint equivalence, ownership/cost/exit combinations, high-tariff limit, BLP 2D integration state. No current Tier 3 requirement.

### Financial

This package is under-tested economically despite 24 blocks. It has a large legacy oracle grid but few independent bank/debt end-to-end equilibrium identities. The highest-value existing tests are the observed 3-bank chain-connected ownership block bug, the asset-weighted psi bug, closed-form liquidity against an independent Monte Carlo comparator, diversification ordering, the Beta moment numerical integral, and the DebtFit default probability direction. Retain those aggressively.

The six-case debt golden master is B-grade temporary migration assurance, especially since `zero_uniform` expects a legacy BBsolve non-convergence warning. That warning is a C-grade legacy behavior trap and should not become a permanent correctness rule. The exact same-seed 200,000-draw Monte Carlo test pins RNG implementation more than economics; move it to extended validation. Keep the stored independent MC-versus-closed-form comparison in routine tests. The `FinancialFit` virtual-class skeleton adds little once DebtFit/BankFit public lifecycle tests run.

Major gaps: bank merger outcomes are checked mainly for finite deltas and object shape, without a public-equilibrium multi-bank FOC, balance-sheet accounting identity, or independently computed liquidity penalty in a small market. The debt test's closed form is compared to a captured fixture, while `calcShares(..., outSideOnly=TRUE)` is compared to a *different* capture; this misses a direct economic threshold relationship. Debt post-merger FOCs, probability-mass accounting, zero-debt limit under an explicit normalization, and welfare effects need independent tests. A reported `posterior_prob_positive` is checked only to be in [0,1]; a two/three-draw hand-aggregated benchmark would catch incorrect posterior mapping.

Tier 1: bank ownership FOC/chain regression, asset weighting, Beta moments/integration, normal liquidity closed form versus independent MC, debt default direction and lifecycle. Tier 2: bank hybrid exact-vs-normal and diversification, merger scenario accounting, posterior multi-draw aggregation. Tier 3: 200k-draw same-seed comparator and six-case legacy debt grid until replaced by independent economic benchmarks. The patch gates the latter two via `FINANCIAL_EXTENDED_TESTS=1`; R syntax was checked, but full execution requires `cubature`.

## Proposed changes and red-team disposition

| Current cluster | Disposition | Why and adversarial check |
|---|---|---|
| coordination `test-price-leadership-invariants.R`, repeat-fit package-location assertion | DELETE AS REDUNDANT (patch) | Creates two identical fits from the same code path; there is no actual package-location contrast. Plausible state leakage would be caught by BLP seeded integration and repeated public behavior tests. |
| coordination `test-price-leadership-blp.R`, final finite-output smoke | DELETE AS REDUNDANT (patch) | The preceding 2D BLP block checks shares, finite outcomes and an independent fringe profit FOC; parity file also constructs BLP. Red-team: no unique economic assertion in the removed block. |
| coordination `test-s4-dispatch.R`, `getClass()` exists-once smoke | DELETE AS REDUNDANT (patch) | `getClass()` does not establish uniqueness. Constructor and dispatch tests catch broken usability. Retained S4 inheritance and antitrust method ownership tests protect the real boundary. |
| coordination Grim Trigger structure/high-low/payoff tests | REWRITE/CONSOLIDATE (proposal only) | Repeated fixture fitting and loose payoff orderings can miss incorrect aggregation of product payoffs to a multiproduct firm. Replace with independent firm-level IC threshold, both ownership states. Preserve post-merger branch; do not delete it. |
| coordination BLP integration metadata tests | KEEP BUT CONSOLIDATE (proposal only) | Several assertions check the same draw count/weights and finite price. Red-team restored 2D node-state and FOC checks because integration errors can alter economics silently. |
| coordination retention exit smoke | REWRITE (proposal only) | Ownership matrix and finite prices can pass with wrong retention-weighted FOCs. Add direct profit-gradient check after exit. |
| coordination core-fringe endpoint parity | KEEP, FULL CI | A permanent A-grade integration benchmark: antitrust simultaneous/MonCom outcomes independently constrain endpoint logic. |
| financial `test-skeleton.R` | DELETE AS LOW VALUE (patch) | Class metadata is already exercised by DebtFit lifecycle. Red-team: BankFit lifecycle also checks public fit/simulation, so no independent failure remains. |
| financial 200k-draw exact Monte Carlo repeat | MOVE TO EXTENDED (patch) | Exact same RNG stream is migration assurance (B). Red-team: the current MC implementation itself could regress; keep in extended, retain routine stored independent MC-vs-closed-form check. |
| financial six-case debt legacy grid | MOVE TO EXTENDED (patch) | Broad B-grade migration assurance; no independent FOC and one C-grade legacy warning. Red-team: density/debt interactions are otherwise thinly covered, so retain in extended until independent tests exist. |
| financial closed-form liquidity vs captured MC | KEEP | A-grade independent algorithm; 5% band is wide but justified by stored MC noise. |
| financial single/multi-market bank comparison | KEEP, FULL CI | Diversification direction and hybrid-vs-normal difference are distinct economically. |
| financial 3-bank chain block-solve | KEEP | Actual serious wrong-FOC bug; catches overlapping-block algorithm. The private helper is justified because the bug is in ownership FOC algebra and no public FOC test currently isolates it. |
| financial Beta h-focal vs numeric integration | KEEP | Independent mathematical benchmark shared by debt and banking. |

### Parity/oracle classification

- **A permanent:** coordination hard endpoint formulas and antitrust endpoint equivalence; PLE hard-coded alpha (review its calibration convention periodically); financial closed-form liquidity vs independent MC, zero-psi/large-phi identities, diversification ordering, Beta numerical integral.
- **B migration:** coordination repeat-fit/package-location and most constructor/extraction checks; financial exact legacy closed-form value replication, Beta fixture grid, same-seed MC, debt six-case grid. The latter grid remains extended because debt economic coverage is thin.
- **C legacy trap:** financial `zero_uniform` cell's required BBsolve non-convergence warning. A valid future solver fix would break that assertion. Replace it with a meaningful zero-debt limiting FOC before retiring the grid.

## Future test rule

Require each new test to name the plausible material failure it detects. Prefer hand calculations, direct profit derivatives, accounting identities, and cross-package behavior over slots, snapshots, implementation branches, and repeated finite-output checks. Label migration tests with a removal condition. Deep-test distinct economic mechanisms; use a small dispatch matrix for syntactic variants. Put expensive statistical or Monte Carlo validation behind an explicit extended tier while keeping a fast independent benchmark in routine CI.
