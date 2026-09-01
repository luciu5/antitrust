# `refactor` migration report

## Branch and design

The work is on branch `refactor`, with `master` treated as the behavioral
oracle.  The new fitted-state representation is the `AntitrustFit` S4 wrapper
around an existing calibrated S4 object.  Simulation results remain the
existing S4 classes.

## Registered generalized model combinations

| Demand | Conduct | Existing class | `calibrate()` | `specify()` / `simulate()` | Legacy paths tested |
|---|---|---|---:|---:|---|
| Linear | Bertrand | `Linear` | yes | yes | `linear()`, `sim()` |
| AIDS | Bertrand | `AIDS` | yes | yes | `aids()`, `sim()` |
| LogLin | Bertrand | `LogLin` | yes | yes | `loglinear()`, `sim()` |
| Linear | Cournot | `Cournot` | yes | simulate only | `cournot()` |
| LogLin | Cournot | `Cournot` | yes | simulate only | `cournot()` |
| Linear | Stackelberg | `Stackelberg` | yes | simulate only | `stackelberg()` |
| LogLin | Stackelberg | `Stackelberg` | yes | simulate only | `stackelberg()` |
| Logit | Bertrand | `Logit` | yes | yes | `logit()`, `sim()` |
| Logit | Cournot | `LogitCournot` | yes | yes | `logit.cournot()`, `sim()` |
| Logit | second-score auction | `Auction2ndLogit` | yes | yes | `auction2nd.logit()`, `sim()` |
| Logit | bargaining | `BargainingLogit` | yes | yes | `bargaining.logit()`, `sim()` |
| Logit | second-score bargaining | `Bargaining2ndLogit` | yes | yes | `bargaining2nd.logit()`, `sim()` |
| CES | Bertrand | `CES` | yes | yes | `ces()`, `sim()` |
| CES | Cournot | `CESCournot` | yes | yes | `ces.cournot()`, `sim()` |
| CES | second-score auction | `Auction2ndCES` | yes | yes | `auction2nd.ces()`, `sim()` |
| CES | bargaining | `BargainingCES` | yes | yes | `bargaining.ces()`, `sim()` |
| CES | second-score bargaining | `Bargaining2ndCES` | yes | yes | `bargaining2nd.ces()`, `sim()` |
| Nested Logit | Bertrand | `LogitNests` | yes | yes | `logit.nests()`, `sim()` |
| Nested CES | Bertrand | `CESNests` | yes | yes | `ces.nests()`, `sim()` |
| LogitCap | Bertrand | `LogitCap` | yes | yes | `logit.cap()`, `sim()` |
| PCAIDS | Bertrand | `PCAIDS` | yes | simulate only | `pcaids()` |
| Nested PCAIDS | Bertrand | `PCAIDSNests` | yes | simulate only | `pcaids.nests()` |
| Capacity-constrained second-score auction | specialized auction | `Auction2ndCap` | yes | simulate only | `auction2nd.cap()` |
| Logit | vertical bargaining | `VertBargBertLogit` | yes | simulate only | `vertical.barg()` |
| BLP | Bertrand | `LogitBLP` | no observed-data calibrator | yes | `sim()` |
| BLP | Cournot | `CournotBLP` | no observed-data calibrator | yes | `sim()` |
| Logit | Bertrand | `LogitALM` | yes (`variant = "alm"`) | simulate only | `logit.alm()` |
| Logit | Cournot | `LogitCournotALM` | yes (`variant = "alm"`) | simulate only | `logit.cournot.alm()` |
| CES | Bertrand | `CESALM` | yes (`variant = "alm"`) | simulate only | `ces.alm()` |
| CES | Cournot | `CESCournotALM` | yes (`variant = "alm"`) | simulate only | `ces.cournot.alm()` |
| Nested Logit | Bertrand | `LogitNestsALM` | yes (`variant = "alm"`) | simulate only | `logit.nests.alm()` |
| LogitCap | Bertrand | `LogitCapALM` | yes (`variant = "alm"`) | simulate only | `logit.cap.alm()` |
| Logit | second-score auction | `Auction2ndLogitALM` | yes (`variant = "alm"`) | simulate only | `auction2nd.logit.alm()` |
| CES | second-score auction | `Auction2ndCESALM` | yes (`variant = "alm"`) | simulate only | `auction2nd.ces.alm()` |
| Logit | bargaining | `BargainingLogitALM` | yes (`variant = "alm"`) | simulate only | `bargaining.logit.alm()` |
| CES | bargaining | `BargainingCESALM` | yes (`variant = "alm"`) | simulate only | `bargaining.ces.alm()` |

PCAIDS parameter specification is intentionally not generalized in this
slice: its existing constructor identifies slopes from known elasticities,
and nested PCAIDS additionally estimates nesting parameters from margins.
Both variants use the shared `simulate()` boundary after calibration.

The capacity-constrained second-score auction is also migrated through a
specialized registry entry. It is not represented as a generic demand/supply
composition: `calibrate()` delegates seller-cost and buyer-value recovery to
`auction2nd.cap()`, while `simulate()` updates ownership, capacity changes,
reserves, and expected prices through the existing auction methods. Supplied
parameter construction is not exposed because `parmsStart` is an optimizer
starting value, not a complete structural-parameter loading interface.

General linear and log-linear Cournot are migrated as separate complete
quantity-game entries. `simulate()` recalculates plant quantities through the
existing `calcQuantities()` method after changing ownership and plant-level
cost deltas; it does not combine a Bertrand demand implementation with a
generic supply module. Direct `specify()` is not exposed because the legacy
Cournot constructor has no separate supplied-demand-parameter loading path.

Linear and log-linear Stackelberg are migrated as separate complete
quantity-game entries. Calibration retains the existing leader-aware demand
and cost identification routine. Simulation accepts the post-merger leader
matrix and other existing post-merger plant-structure fields through `...`,
then invokes the legacy Stackelberg quantity solver. Direct `specify()` is not
exposed because the legacy constructor identifies plant cost parameters from
margins.

Standard Logit vertical bargaining is migrated as a complete two-sided model
entry. `calibrate()` delegates downstream demand and upstream/downstream
bargaining-power identification to `vertical.barg()` using pre-merger
ownership on both sides. `simulate()` accepts upstream and downstream post
ownership in a named list, updates integration-dependent bargaining power and
both cost-delta vectors, and invokes the existing vertical price system.
Nested and second-score vertical variants remain legacy-only because their
scenario and equilibrium methods require a separate migration slice.

BLP is deliberately marked calibration-ineligible: the repository has a
parameterized BLP construction path but no observed-data BLP estimator that
belongs behind `calibrate()`.

## Compatibility and parity coverage

The test foundation includes hard-coded master-oracle values for representative
Logit-Bertrand and Logit-Cournot cases, including demand parameters, recovered
costs, pre/post prices, shares, margins, elasticities, diversion, UPP, CMCR,
and CV.  PCAIDS and nested PCAIDS parity tests additionally compare their
known-elasticity/nesting parameters and inherited AIDS outputs. The
architecture tests then compare old and new paths for all
registered families, including supplied-parameter construction and repeated
counterfactuals where applicable.

The BLP regression coverage independently checks the contraction fixed point
for both outside-good and no-outside-good markets.  During this audit,
`master` was found to force an outside option during contraction even when
the fitted object used a no-outside normalization.  The refactor branch now
uses the same normalization in `sim()`, `calcSlopes()`, and `calcShares()`;
the ordinary outside-good path is unchanged.

The existing full `devtools::test()` suite passes after each committed phase.
The final suite retains only the repository's prior warnings: missing-data ALM
boundary warnings, CES optimizer diagnostics, and the existing bargaining CES
optimizer diagnostic.

Direct branch comparison was also run from clean archives of `master` and
`refactor`.  Both full test suites passed, and both built packages returned
`Status: OK` from `R CMD check --no-manual --no-vignettes`.  The shared
`ai/examples` scripts had matching exit status: the three BLP examples that
complete produced matching numeric results, while the historical
price-leadership examples failed on both branches because `ple` is not
available in the package namespace.  The BLP informational messages emitted
by legacy `sim()` are preserved by the compatibility wrapper; performance
timings remain inherently variable.

## Intentionally preserved behavior

* Existing S4 classes, slots, methods, solvers, defaults, normalizations, and
  fallback behavior remain the source of truth.
* `calibrate()` captures the legacy same-ownership warning used by its
  pre-merger placeholder and stores it in `AntitrustFit@diagnostics` rather
  than printing it as a spurious user-facing merger warning.
* `simulate()` keeps the existing model-specific `mcDelta` and outside-good
  conventions.  It does not impose one supply formula on all demand systems.
* BLP keeps its existing nonlinear price solver; AG is not silently substituted
  for it.
* `sim()` remains exported with its historical argument signature and falls
  back to the original implementation for unsupported/non-registered paths.

## Known discrepancies and suspected issues not changed

No unexplained economic-output discrepancy remains in the migrated generalized
combinations covered by the tests.  The following historical constraints or
oddities were observed and intentionally left unchanged:

* The legacy parameterized LogLin path can reject a supplied slope matrix when
  the diversion matrix mechanically derived from it has rows above zero.  This
  was observed while probing parameters recovered by a legacy `loglinear`
  call; the validation and error behavior were not broadened.
* BLP has no observed-data calibration function in the current repository, so
  `calibrate(demand = "blp", ...)` remains unsupported rather than inventing a
  new estimator.
* The `master` BLP contraction path did not honor its no-outside-good
  normalization when recovering `meanval` from shares summing to one.  This
  was reproduced independently and fixed on `refactor`; the correction is
  covered by the BLP contraction tests.
* The regression example `ai/examples/test_sim_regressions.R` currently stops
  at the historical LogitCap requirement for a finite `meanval`; this is
  unchanged and is not treated as a refactor fix.
* Specialized ALM and nested/second-score vertical constructors are not
  generalized `sim()` demand/conduct entries.  Their legacy functions were
  not removed or mechanically wrapped.  General linear and log-linear
  Cournot/Stackelberg and standard Logit vertical bargaining are registered
  for the new fit pipeline, while the legacy `sim()` route remains unchanged.

These are migration notes, not economic corrections.
