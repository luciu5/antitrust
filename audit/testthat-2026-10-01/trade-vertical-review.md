# Trade and vertical test audit

Sources: `trade` refactor and `vertical` main working trees on 2026-10-01. Production implementation read: trade registries, architecture, promotion, fitted tariff accounts, synthetic policy adapter; vertical registry, lifecycle, ownership classification, counterfactuals, output, margins and welfare methods. Every `test_that()` block is recorded in `inventory_before.csv` with measured baseline seconds. The active vertical tree has pre-existing changes; no source-file edits were made there or in trade.

## Measured baseline and staged result

| Package | Before | Baseline runtime | Staged blocks | Staged runtime | Verification |
|---|---:|---:|---:|---:|---|
| trade | 19 files, 143 blocks | 113.689 s all-tests, 0 failures, 1 skip, 113 warnings | 18 files, 137 blocks plus tier helper | fast 18.482 s; extended 59.571 s; nightly 106.651 s | all tiers 0 failures/errors; nightly 0 skips |
| vertical | 5 files, 19 blocks | 10.542 s after antitrust fix, 0 failures, 46 warnings | 5 files, 17 blocks | 9.462 s | 0 failures/errors |

Before antitrust's vector ownership correction, vertical had 3 second-score failures. They demonstrated the value of cross-package tests. `antitrust/R/Retention.R` now normalizes `Auction2ndLogit` vector ownership through `ownerToMatrix`; all 19 vertical baseline blocks passed after that fix.

Trade runtime is concentrated: BLP provided integration calibration took 35.478 s and GH/Monte Carlo recovery took 47.104 s. The full fitted 14-game tariff matrix took 6.448 s; exact promotion route matrix 2.742 s. The prior `trade/tests/testthat/TESTING.md` claim that trade had no slow duplication was empirically false.

## Economic and API value

Trade's strongest permanent tests are the independent Cournot tariff analytical solution, log-linear plant FOCs, capacity KKT and zero-output corner; mixed-tariff physical-profit FOCs for Bertrand/Cournot; physical cost plus producer/government accounting identities; hard-coded native tariff/quota numerical oracles; and the antitrust/coordination-to-trade promotion and output-game matrices. Registry tests correctly list complete implemented combinations rather than constructing arbitrary demand × conduct cells. `calibrate`, `specify`, `simulate`, `update`, `respecify`, and policy promotion have separate test anchors.

Vertical's strongest tests are hybrid ownership classification by relationships, per-tier HHI, sequential cost/ownership promotion, second-score diversion residual, and inactive-product boundary conventions. Its 19 blocks contain many hundreds of assertions but are only about 10 s to run. The subset matrix has three chain levels × direct and fitted lifecycles; it is useful, but repeated output-shape assertions add less value than the first deep case. The nested and second-score cases reflect distinct equations and should stay.

## Parity/oracle classification

- **A permanent:** trade's three hard-coded tariff/quota/supplied-parameter oracles; zero-tariff antitrust/trade integration; independently checked tariff-accounting/FOC matrices. Vertical's selected second-score and hybrid price, welfare and FOC anchors can be A until an independent benchmark replaces them.
- **B migration assurance:** trade's six supplied `sim()` path comparisons and direct constructor-versus-lifecycle parity; vertical's full pre-extraction fixture comparisons across flat, nested and chain-level cases. They currently protect the refactor but should be reviewed after canonicalization. The trade supplied-parameter matrix is staged for nightly.
- **C legacy trap:** the exact frozen fixture SHA and exact UPP/CMCR error-message snapshots in vertical. The SHA assertion is staged for deletion; exact error strings remain in source and should be rewritten as public error-class/capability checks when vertical becomes writable.

## Staged deletions and adversarial check

| Existing cluster | Disposition | Why a material bug remains detectable |
|---|---|---|
| trade `test-price-start.R` | delete | AIDS hard-coded post-price oracle and custom-start result in `test-trade-architecture.R` catch solver-start changes that alter economics. Private `@priceStart` value alone promises little. |
| trade `test-structural-fit-contract.R`: arbitrary slot constructor | delete | Real `calibrate()`/`specify()` and shared StructuralFit inheritance/dispatch exercise actual fits. |
| same file: duplicate named `simulate()` | delete | Stronger `antitrust::simulate(fit)` versus `trade::simulate(fit)` comparison remains. |
| same file: dummy sibling fallback | delete | Antitrust owns/tests the shared fallback; a trade wrapper cannot silently change it without failing antitrust or cross-package dispatch tests. |
| trade `test-trade-model-registry.R`: private accessor delegation | delete | Public `respecify()` translation outcomes and transition registry behavior remain. A private accessor can change safely. |
| trade `test-tariff-game-fit.R`: fractional ownership | delete | The block was always skipped because `antitrust::specify()` rejected the fixture first. An active test in `test-fitted-tariff-equivalence.R:248` mutates a valid source owner matrix and reaches `as_trade_fit()`, confirming the public `trade_tariff_unsupported_ownership` contract. |
| vertical fixture SHA | delete | No economic assertion depends on this literal source commit; numerical fixtures remain. |
| vertical subset nested-second-score rejection | delete | The same public rejection appears in `test-parity-fixtures.R`; the registered model set also excludes this combination. |

The staged trade RDS test uses a real calibrated fit, serializes it, then compares tariff prices/costs after restoration. This protects saved user models better than serializing an arbitrary slot list.

## Tiering and remaining gaps

The staged trade patch adds `TRADE_TEST_TIER`: fast on pushes/local runs, extended on PRs, nightly on schedule. Extended keeps the exhaustive promotion/output matrices and BLP supplied-integration calibration. Nightly keeps both BLP recovery rules and the temporary supplied-parameter migration matrix. The suite remains cumulative by tier. All three staged tiers were verified, including the nightly all-block run.

Gaps remain. Trade's BLP recovery fixture generates shares with antitrust internals, so an independent external or hand-coded small benchmark would strengthen it. The fitted 14-game matrix compares to public antitrust/coordination solvers that share components with trade; it tests the adapter and cost transformation, not independent equilibrium truth. Retain separate numerical FOCs and physical accounting identities. Trade has 113 warnings in the baseline, which obscure new convergence warnings. Vertical's `subset_foc()` uses package `calcMargins()` within its residual, so the dense assertion grid is partly self-referential. Add a hand-derived vertical bargaining and welfare benchmark with independently computed demand derivatives. Vertical's `update()` test asserts only object class; it should show changed observed margins recalibrate economic primitives while `simulate()` changes only the scenario. A small cross-package test of vertical tariff/quota policy should be added once the current uncommitted policy implementation is complete.

## Patch instructions

Both patches were generated from temporary clones and passed `git apply --check` against the untouched source worktrees:

- `patches/trade.patch`: from `trade` refactor, `git apply patches/trade.patch`
- `patches/vertical.patch`: from `vertical` main, `git apply patches/vertical.patch`

Neither patch touches production code. The vertical patch touches only two test files and does not disturb pre-existing uncommitted vertical production changes. Commit each package separately after applying and verifying in writable repository worktrees.
