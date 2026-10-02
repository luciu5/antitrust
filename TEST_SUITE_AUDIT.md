# ioverse `testthat` audit — 2026-10-01

## Scope, method, and status

This audit covers `antitrust` and `trade` on `refactor`, and `vertical`, `coordination`, `antitrustBayes`, and `financial` on `main`. The four latter repositories have no local or remote `refactor` branch. `antitrust`, `trade`, `coordination`, and `financial` were already current when pulled; `antitrustBayes` fast-forwarded one documentation commit. `vertical` was fetched but not pulled because it contains pre-existing tracked edits. No scoped `DESCRIPTION` imports `iopolicy`.

Source tips inspected were `antitrust` `33a374f`, `trade` `4f13747`, `vertical` `8fe0be5`, `coordination` `24eea76`, `antitrustBayes` `1c2cfa9`, and `financial` `70f78fe` before audit edits. The exact `vertical` worktree also contains the pre-existing uncommitted policy changes noted below.

I inventoried all 833 original `test_that()` blocks in 131 files, measured per-block runtime where dependencies permitted, read the implementation behind the economic, lifecycle, numerical, and cross-package clusters, classified parity checks, and recorded a deletion red-team pass before edits. The machine-readable [before inventory](audit/testthat-2026-10-01/inventory_before.csv) contains package, branch, location, description, runtime, behavior, duplicate cluster, type, value, disposition, parity class, and deletion risk. The [actual after inventory](audit/testthat-2026-10-01/inventory_after_actual.csv) records checked-out source at report generation time. Detailed package notes and red-team records are in [the audit folder](audit/testthat-2026-10-01/).

The `antitrust` and `antitrustBayes` edits were made directly. The other four package worktrees are mounted read-only for the agent, including escalated commands; the user applied their reviewed [patches](audit/testthat-2026-10-01/patches/) from a writable shell. **The after column is the actual checked-out state at audit time.** The agent could not create package commits because each `.git` directory is mounted read-only for its tools. The existing `vertical` production edits and `trade` untracked artifact were preserved.

## A. Executive summary

| Package / branch | Before files / blocks | After files / blocks | Baseline wall time | Verified revised wall time |
|---|---:|---:|---:|---:|
| `antitrust` / `refactor` | 39 / 283 | 35 / 263 | fast 36.8 s | fast 22.8 s; extended 91.9 s; nightly 211.0 s |
| `trade` / `refactor` | 19 / 143 | 18 / 137 | all 113.7 s | actual fast 24.4 s and extended pass; temporary-copy extended 59.6 s, nightly 106.7 s |
| `vertical` / `main` | 5 / 19 | 5 / 17 | 10.5 s after antitrust fix | actual `test_dir` 11.7 s; temporary-copy 9.5 s |
| `coordination` / `main` | 12 / 78 | 12 / 75 | 26.1 s | actual 25.1 s; temporary-copy 23.1 s |
| `antitrustBayes` / `main` | 49 / 286 | 48 / 272 | routine 33.5 s | routine 32.2 s; extended fixtures 6.9 s; focused Stan 6.8 s warm / 54.4 s cold |
| `financial` / `main` | 7 / 24 | 6 / 23 | unmeasured | unmeasured: `cubature` unavailable |
| **Total** | **131 / 833** | **124 / 787** | approximately 221 s across five measurable suites | approximately 116 s for five routine suites; different test commands and runs, not a controlled aggregate benchmark |

The final accounting is 782 original blocks retained, 20 old antitrust S4/structural blocks consolidated into three, one trade saved-fit block rewritten, six Bayesian benchmark-fixture checks moved to a separate extended script, 24 old blocks deleted without replacement, and one new cross-package regression added. Eleven retained blocks change routine tier: two in `antitrust`, seven in `trade`, and two in `financial`. These categories overlap only where an old block was replaced, so the reliable size metric is **833 to 787 ordinary blocks**, a 46-block net reduction. The largest runtime gain comes from tiering two expensive trade BLP tests, not from deleting useful numerical validation.

The suite **was bloated selectively**: `antitrust` had duplicated lifecycle and infrastructure assertions; `trade` ran two BLP calibration checks costing 82.6 seconds together on every routine pass; `antitrustBayes` included a 104-cell synthetic benchmark-harness check and basic class/structure snapshots. A wholesale 30–50% deletion would be unjustified. Many dense suites protect distinct economic mechanisms, especially `coordination` and the antitrust BLP, synthetic-market, and sequential-counterfactual clusters.

## B. Package assessments and final tiers

### `antitrust`

The strongest tests use independent Logit/CES/Linear calculations, Bertrand and Cournot master-derived numerical anchors, mixed-retention profit FOCs, BLP derivative/heterogeneity limits, multi-product synthetic-market FOCs, and sequential state/cost identities. `supportedModels()` enumerates complete implemented combinations rather than a Cartesian demand × conduct table. The `calibrate`/`specify`/`update`/`respecify` distinction is protected; in particular, `respecify()` must preserve source structural primitives without silently recalibrating.

The weakest area was repeated S4 metadata and lifecycle smoke across four files. I consolidated it, removed repository CSV and test-helper self-tests, and added a vector-owner auction regression after the `vertical` suite exposed a real cross-package failure. The retained `test-legacy-model-parity.R` is mixed: some cells encode rare model economics, while many `sim()` route comparisons are **B: temporary migration assurance**. The latter remain in nightly validation until refactor canonicalization. Fast tests include independent FOCs/API/state; extended includes broader solver, ALM and stochastic checks; nightly contains rare-family and migration matrices. Gap: AIDS/PCAIDS and some rare Stackelberg/auction results still rely heavily on legacy parity instead of independent FOCs or welfare benchmarks.

### `trade`

The strongest tests independently check tariff-incidence and physical-cost accounting, Bertrand/Cournot tariff FOCs, capacity KKT and zero-output behavior, hard-coded tariff/quota outcomes, and real antitrust/coordination promotion. Its registered model tests correctly distinguish implemented combinations. The applied changes delete vacuous getters, dummy S4 construction, private delegation, and an always-skipped ownership case; an active public fractional-ownership rejection test remains. They replace a dummy RDS test with a saved fitted model whose simulated tariff prices and costs are compared after restoration.

Fast tests keep representative economic and API paths. Full CI keeps exhaustive output-game/promotion matrices and provided-point BLP calibration. Nightly keeps GH-versus-Monte-Carlo BLP recovery and B-grade `sim()` migration parity. Gap: BLP calibration fixtures derive their target shares through antitrust internals, so recovery is not an independent external numerical oracle; 113 baseline warnings also obscure new solver warnings.

### `vertical`

The suite checks hybrid ownership classification, three chain levels, subset/inactive-product behavior, second-score diversion, and selected prices/welfare against pre-extraction fixtures. The cross-package suite caught the antitrust vector-ownership regression; its value is unusually high relative to its 19 blocks. The applied changes remove a frozen source SHA and a duplicated unsupported-model error, while preserving numerical fixtures through the migration.

The subset tests contain hundreds of assertions, but some FOC residuals use the same `calcMargins()` implementation as the solver, making them partly circular. Add an independently calculated small-market vertical bargaining and welfare benchmark. Its current uncommitted tariff-policy implementation is not covered by a completed trade/vertical integration test; preserve those edits and add one when the feature settles. Tier 1 keeps one deep case per chain-level/ownership mechanism; Tier 2 keeps selected migration fixtures and cross-package scenarios; Tier 3 is unnecessary once independent welfare benchmarks replace old fixtures.

### `coordination`

This suite is dense but mostly valuable. Leader/follower Stackelberg FOCs, core-fringe endpoint formulas, retention-weighted portfolio FOCs, BLP fringe finite-difference gradients, and signed/zero input-rate tests all guard economically distinct logic. The applied changes remove three pure duplicate or metadata smokes; the actual suite passes. Keep representative Logit/CES FOCs, IC logic, real-rate boundaries, and public dispatch in Tier 1; keep broader 2D BLP integration and antitrust endpoint equivalence in Tier 2. No existing Tier 3 is essential.

The main gap is Grim Trigger: several tests inspect payoff ordering or package residuals, but no simple independent multiproduct, firm-level repeated-game IC threshold is computed from Coord/Defect/Punish payoffs. The exit/retention smoke similarly needs an independent post-exit portfolio gradient. These are proposed rewrites, not safe deletions.

### `antitrustBayes`

The largest test count is not simply waste. High-value checks include analytical five-conduct Logit/CES FOCs, dense Gaussian posterior-quality/intercept conditioning, exact-`k` model evidence and role probabilities, Student-t likelihood identities, and actual `simulate()` integration into antitrust and coordination. Eight low-value class/structure/benchmark-harness blocks were deleted; six manifest and DGP checks moved to an explicit extended fixture script. The moved checks include the 104-cell synthetic-market smoke: the independent deletion review found that the archived recovery runner could no longer serve as its substitute. Default suite retains 272 deterministic/public-contract blocks. Keep those in Tier 1, selected CmdStan compile and short posterior sanity tests in Tier 2, and recovery/model comparison in Tier 3.

A serious gap remains: 12 Stan sampling blocks are opt-in and skipped by default, so routine success is not statistical recovery evidence. The published 124-task enterprise recovery runner is stale after CmdStanR migration: it requires `rstan`, looks for package `stanmodels` that no longer exists, and reads `fit@stanfit`. The six moved manifest/DGP checks remain executable, but the HMC/recovery study itself cannot currently run. The audit also exposed a skipped test that passed an RStan-era `algorithm="Fixed_param"` option into a CmdStanR backend; that test was changed to short HMC and reverified. Until a current recovery harness exists, posterior recovery and Bayes-factor model selection need a new independently established validation experiment.

The workflow currently runs the Stan sampling job on a schedule or explicit dispatch, not on pull requests. Thus the proposed Tier 2 sampling checks are an extended gate in the present CI configuration; deterministic Bayes checks still run on pull requests.

### `financial`

This is the clearest case of **under-testing despite a legacy oracle grid**. The best checks are the previously observed chain-connected three-bank ownership bug, asset-weighted `psiRule`, Beta partial-moment numerical integration, liquidity closed form against independent Monte Carlo, zero-penalty limits, and public DebtFit/BankFit lifecycle. The applied changes delete a package-load skeleton and move exact same-seed 200,000-draw Monte Carlo and the six-case debt legacy grid to Tier 3. The grid's required BBsolve non-convergence warning is **C: a legacy-behavior trap**; a valid solver improvement would break it.

Tier 1 should center on independent bank FOCs, liquidity and balance-sheet identities, Beta moments, and debt default/zero-debt behavior. Tier 2 should test bank scenario and posterior multi-draw aggregation. Tier 3 retains migration grid and large Monte Carlo only until stronger independent checks replace them. Missing now: public multi-bank equilibrium FOC, small-market balance-sheet identity, direct debt threshold/FOC and welfare benchmark, and a hand-aggregated two/three-draw posterior probability check. The full suite could not load here because `cubature` is absent and the package index was unreachable; only R syntax, patch applicability, standalone Beta integral points, and a legacy tail fixture were verified.

`financial` has no scheduled extended CI job at present. Until one is added, maintainers must explicitly run `FINANCIAL_EXTENDED_TESTS=1` before releases; otherwise the two gated migration checks would be effectively dormant.

## C. Changes, disposition, replacement, and removal risk

| Current location / cluster | Disposition | Rationale and replacement | Removal risk after red-team |
|---|---|---|---|
| `antitrust/test-api-coverage.R` | DELETE AS LOW VALUE | CSV/export sync does not prove an API works; behavior and registry tests remain. | Low: new untested exports are a review issue, not caught meaningfully by this CSV. |
| `antitrust/test-capability-audit.R` | DELETE AS REDUNDANT | Private flags and finite-output cases duplicate public quality, entry, and translation tests. | Low: those public routes fail if capability rejection changes. |
| `antitrust/test-tier-routing.R` | DELETE AS LOW VALUE | Tests a simple test helper, not package behavior. | Low: CI tier runs reveal skipped blocks. |
| `antitrust/test-s4-lifecycle-dispatch.R` plus `test-structural-fit.R` | KEEP BUT CONSOLIDATE | Twenty blocks become three tests of shared slots, real dispatch, unsupported sibling failure, and `stats::simulate()`. Model API and sequential tests retain economics. | Moderate before review; independent reviewer found no lost material contract. |
| `antitrust/test-legacy-model-parity.R`: BLP message parity | DELETE AS LOW VALUE, C | Exact prose against a private old path has no economic contract. | Low; result and failure paths remain. |
| Same file: BLP alias | MOVE TO FULL CI, B | Backward-compatible conduct mapping remains protected off fast tier. | Low for local runs; CI catches it. |
| `antitrust/test-stochastic-reproducibility.R`: generated draws | MOVE TO FULL CI | Fixed-draw RNG protection stays fast; seeded generation stays in extended. | Low for local runs; CI catches it. |
| `antitrust/R/Retention.R` and `test-revenue-retention.R` | REWRITE/ADD ECONOMIC REGRESSION | Normalize valid vector ownership through `ownerToMatrix`; test two multiproduct firms, uniform and mixed pre/post retention with active subset. | Fixes observed vertical second-score failure. |
| `trade/test-price-start.R` | DELETE AS LOW VALUE | Private start-vector snapshot; hard-coded tariff oracle and custom-start economic result remain. | Low. |
| `trade/test-structural-fit-contract.R`: dummy slots/dispatch/fallback/RDS | KEEP BUT CONSOLIDATE / REWRITE | Replace dummy RDS with real fitted tariff economics after serialization; retain extension and actual dispatch checks. | Low after stronger fit test. |
| `trade/test-trade-model-registry.R`: private lookup delegation | DELETE AS REDUNDANT | Public transition and registry tests remain. | Low. |
| `trade/test-tariff-game-fit.R`: fractional ownership skip | DELETE AS OBSOLETE | It always skipped before reaching trade; active `test-fitted-tariff-equivalence.R` checks public rejection. | Low. |
| `trade` BLP and exhaustive tariff/promotion matrices | MOVE TO FULL CI / EXTENDED | Seven gated blocks preserve expensive economics with cumulative tier routing; no deletion. | No permanent coverage loss; delayed local feedback. |
| `vertical/test-parity-fixtures.R`: literal SHA | DELETE AS LOW VALUE, C | Preserve numerical fixtures; delete frozen historical source hash. | None for economics. |
| `vertical/test-subset-support.R`: duplicate nested-second-score error | DELETE AS REDUNDANT | Same public rejection in parity fixtures and registry. | Low. |
| `coordination/test-price-leadership-invariants.R`: repeat-fit identity | DELETE AS REDUNDANT, B | Two runs of same code path are not an independent extraction oracle. | Low. |
| `coordination/test-price-leadership-blp.R`: finite smoke | DELETE AS REDUNDANT | Preceding 2D FOC and other construction tests are stronger. | Low. |
| `coordination/test-s4-dispatch.R`: class-exists smoke | DELETE AS LOW VALUE | Actual inheritance/dispatch remains. | Low. |
| `antitrustBayes/test-classes.R`, `test-helpers.R` snapshots | DELETE AS REDUNDANT | Real fit/simulation and numerical Stan-data contracts remain. | Low after posterior and public API tests. |
| `antitrustBayes/test-market-intercept-enterprise-harness.R` | MOVE six manifest/DGP checks TO EXTENDED; DELETE one worker-path smoke | Fixture checks become executable `doc/market_intercept_enterprise/fixture_checks.R`; the 104-cell synthetic smoke remains because the archived runner cannot replace it. The path-only check adds no separate economic protection. | Moderate: fixture logic retained; runner itself is stale and must be rebuilt. |
| `antitrustBayes/test-posterior-predictive.R`: saved-fit lazy-cache assertion | DELETE AS LOW VALUE | Always skipped without an externally supplied RDS; public posterior-prediction tests remain. | Low; this did not provide routine protection. |
| `antitrustBayes/test-rho-studentt.R` formula-only limit and duplicate `test-stan-core-policy.R` option check | DELETE AS LOW VALUE / REDUNDANT | Neither added independent implementation coverage; density and core-policy tests remain. | Low. |
| `financial/test-skeleton.R` | DELETE AS LOW VALUE | Real DebtFit/BankFit lifecycle checks supersede virtual-class smoke. | Low. |
| `financial/test-bank-parity.R` exact RNG repeat and `test-debt-parity.R` legacy grid | MOVE TO EXTENDED VALIDATION, B | Keep costly migration checks behind `FINANCIAL_EXTENDED_TESTS=1` until independent economic benchmarks exist. | Moderate if Tier 3 is not run; explicit release/nightly scheduling needed. |

All other blocks are classified in the inventory. The broader `antitrust` legacy parity file and `vertical` pre-extraction numerical fixtures were **not** deleted: rare-family economics still lacks independent replacements. Likewise, `coordination` Grim Trigger and `financial` debt tests are candidates for stronger rewrites, not deletion today.

## D. Important uncovered failures

1. **Bayesian recovery is not operationally validated.** The enterprise runner uses obsolete RStan internals; repair it around current CmdStanR fits and verify a small known-truth recovery/model-evidence study before treating Tier 3 as real.
2. **Financial equilibria lack independent public FOCs/accounting.** A wrong solver can return finite, plausible merger deltas while legacy fixture comparisons remain green.
3. **Vertical welfare and subset residuals are partly circular.** Derive one small-market bargaining/welfare calculation independently, including multi-product ownership and inactive products.
4. **Coordination Grim Trigger lacks a firm-level IC oracle.** Calculate the discount threshold from independent per-firm profits before and after merger.
5. **Rare antitrust conduct families still lean on master parity.** Replace B-grade migration cells with analytical FOCs as refactor becomes canonical. Keep independent Bertrand/Cournot master anchors.
6. **Trade BLP recovery shares come from its dependency's own kernel.** Add an independent small integration/calibration benchmark; retain the existing full-CI cost/FOC tests.
7. **Cross-package version drift remains dangerous.** The antitrust vector-owner regression broke three vertical tests while antitrust fast passed. Test package stacks against the same refactor SHA before release.

## E. Permanent test policy

A new test must state the plausible material failure it would catch and why existing tests would miss it. Prefer a hand calculation, independent numerical benchmark, economic identity, direct profit FOC, accounting identity, or public cross-package behavior. Test each distinct economic mechanism deeply; use a small registry/dispatch matrix for syntactic variants. Do not retain private snapshots, trivial slot construction, repeated finite-output checks, or a broad grid merely for coverage. Label migration parity **A permanent**, **B remove after canonicalization/replacement**, or **C legacy trap**, and give B tests a sunset condition. Keep deterministic algebra in fast tests, broader integration in full CI, and costly sampling/recovery/Monte Carlo in explicit extended validation. A skipped or stale validation runner is a coverage gap, not a passing test.

The concise version is also added to `antitrust/AGENTS.md` so future agents see it while changing the kernel.

## Verification and remaining action

- `antitrust` refactor: fast, extended, and nightly `testthat::test_local()` passed with zero failures/errors after the main edit batch; the added post-state retention assertion passed in a focused rerun. Its `vertical` main integration suite changed from three second-score failures to 19/19 passing blocks before pruning.
- `trade` refactor: actual worktree fast tier passed, 137 blocks and zero failures/errors in 24.4 s. Its actual extended tier also exited successfully with no test failures; temporary-copy extended and nightly tiers passed in 59.6 and 106.7 s.
- `vertical` main: actual worktree `pkgload::load_all()` plus `testthat::test_dir()` passed 17 blocks with zero failures/errors in 11.7 s. `testthat::test_local()` produced two summary-dispatch failures; that loader-specific discrepancy is unresolved and is not counted as a passing `test_local()` run.
- An installed-package vertical check is blocked independently: its pre-existing uncommitted `DESCRIPTION` Collate entry names `R/VerticalHorizontalFits.R`, but that file is absent. The audit did not alter this in-progress production change.
- `coordination` main: actual worktree passed 75 blocks, zero failures/errors, 25.1 s.
- `antitrustBayes` main: final routine `test_local()` passed 272 blocks in 32.2 s with zero failures/errors, 12 opt-in skips and 17 expected warnings. The separate six-block manifest/DGP fixture suite passed in 6.9 s. The corrected one-market CmdStanR test passed four blocks with zero failures/errors in 54.4 s on a cold compile and 6.8 s cached. These are short-chain checks, not posterior recovery validation.
- `financial` main: full suite unrun because `cubature` could not be installed in this environment. The patch parses and applies cleanly; a standalone Beta integral check passed.

The four patches have been applied. Package commits require a writable user shell because the agent cannot write Git metadata. The reviewed [commit script](audit/testthat-2026-10-01/commit_audit.sh) checks branches and staged-index cleanliness, then stages only audit paths and commits package by package. The pre-existing `vertical` production changes were outside this audit's edits and are not part of those commits.
