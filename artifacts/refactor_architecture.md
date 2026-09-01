# `antitrust` calibration/simulation architecture

This document describes the architecture on the `refactor` branch.  The
existing model-specific S4 classes remain the economic implementation and the
existing simulation result objects remain the result API.

## Model identity and registry

`model_spec(demand, conduct, variant = "standard")` returns a small S3 list
with normalized `demand`, `conduct`, `variant`, and `id` fields.  Names such
as `Logit`, `CESNests`, `auction_2nd`, and `LogitBLP` are normalized to the
vocabulary used by the registry.  Existing ALM calibration variants are
selected with `variant = "alm"` (or an ALM demand alias).

`supportedModels()` is generated from one internal registry in
`R/ModelRegistry.R`.  Each entry records the normalized demand/conduct pair,
variant, existing S4 class, legacy calibrator mapping, and whether the entry
supports `calibrate()`, `specify()`, and `simulate()`.

The registry describes the combinations already exposed by generalized
`sim()`, plus the existing ALM constructors: Linear, AIDS, PCAIDS (including
nested PCAIDS), LogLin, Logit, CES, nested Logit/CES, LogitCap, BLP, Cournot,
second-score auction, bargaining, and second-score bargaining where the legacy
package supports those pairings. PCAIDS is registered as its own Bertrand
demand family even though its result classes inherit from AIDS; its
known-elasticity and nested-parameter calibration remains in the existing
PCAIDS-specific methods.
ALM entries point to their complete legacy model-specific implementations;
they are not assembled from a generic supply module.  Specialized capacity
auctions, vertical models, and general quantity-game constructors remain
outside this generalized registry and retain their legacy APIs.

## Fitted state

The branch uses approach A from the design proposal: `AntitrustFit` is a small
S4 wrapper containing:

* the normalized model specification;
* the existing calibrated S4 model object;
* recovered or supplied structural parameters;
* observed inputs; and
* captured calibration warnings/messages and dispatch metadata.

This wrapper cleanly marks the state before a particular counterfactual.  It
does not replace or alter the existing S4 inheritance hierarchy.  `simulate()`
returns the existing S4 model class (`Logit`, `CES`, `BargainingLogit`, and so
on), so downstream methods such as `summary()`, `elast()`, `diversion()`,
`upp()`, `cmcr()`, `CV()`, and `HypoMonTest()` continue to operate on their
established result types.

## Calibration and parameter specification

`calibrate()` dispatches through the registry to the existing model-specific
constructor, including the ALM entries selected by `variant = "alm"`.  For
PCAIDS, `knownElast` and `mktElast` are explicit calibration inputs and nested
PCAIDS also receives its `nests` and observed margins. The constructor remains
responsible for its own calibration equations, parameter
identification, normalization, cost recovery, solver, and warnings.  The
wrapper supplies pre-merger ownership as the required legacy post-ownership
placeholder; that placeholder is replaced before a counterfactual is
simulated.

`specify()` dispatches supplied structural parameters through the legacy
parameterized construction path.  It validates and loads the parameters
without recalibrating demand.  Linear and LogLin use an explicit `quantities`
argument because their historical constructors distinguish quantities from
shares.  BLP supports this route for both Bertrand and Cournot, but there is no
observed-data BLP calibrator in the current package, so its registry entries
mark `calibrate = FALSE`.

## Simulation

`simulate(fit, ownerPost, mcDelta, subset, ...)` applies the scenario to the
fitted model, recovers post-counterfactual marginal costs with the model's
existing `calcMC()` method, and calls that model's existing equilibrium
method.  The dispatch preserves important differences:

* Logit/CES Bertrand use their existing nonlinear or AG price solvers.
* Logit/CES Cournot use their existing quantity-game price methods.
* Bargaining retains its bargaining-specific equilibrium equations and solver
  choice, including AG support where the legacy constructor provides it.
* Second-score auction and second-score bargaining retain their direct or
  model-specific equilibrium methods.
* AIDS and PCAIDS recompute their ownership-dependent price-delta equation at
  simulation time; nested PCAIDS retains its model-specific nesting
  calibration.

There is no generic supply module that mechanically combines arbitrary demand
and conduct modules.  The registry selects a complete existing model
implementation, and scenario arguments are passed only to the methods that
understand them.  For example, the legacy second-score cost conventions are
not reinterpreted as a universal marginal-cost formula.

The same fit can be simulated repeatedly:

```r
fit <- calibrate(
  demand = "logit", conduct = "cournot",
  prices = prices, shares = shares, margins = margins,
  ownerPre = ownerPre
)

r1 <- simulate(fit, ownerPost = mergerAB)
r2 <- simulate(fit, ownerPost = mergerAC)
r3 <- simulate(fit, ownerPost = mergerBC)
```

Legacy constructors remain available.  Migrated `logit()`, `logit.cournot()`,
and generalized `sim()` calls are compatibility wrappers around the new
boundary where registry coverage and parity tests are complete; non-registered
legacy constructors continue through their original code paths.
