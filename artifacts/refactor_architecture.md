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
package supports those pairings. Linear and log-linear general Cournot and
Stackelberg are also registered as complete quantity-game entries. Standard
Logit vertical bargaining is registered as complete two-sided model entries for
both downstream Bertrand and downstream second-score conduct.
The capacity-constrained second-score auction
is registered as a complete specialized model entry, rather than as a demand
module mechanically combined with auction conduct. PCAIDS is registered as its own Bertrand
demand family even though its result classes inherit from AIDS; its
known-elasticity and nested-parameter calibration remains in the existing
PCAIDS-specific methods.
ALM entries point to their complete legacy model-specific implementations;
they are not assembled from a generic supply module. Nested vertical
second-score variants remain outside this generalized registry and retain
their legacy APIs.

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

Demand-transition validation remains model-specific.  In ordinary output
markets, CES uses `gamma > 1` and nested CES uses `sigma_g > gamma > 1`; these
are distinct restrictions from nested-Logit `lambda_g` in `(0, 1]`.  The
translation layer carries the source model's `output` orientation into the
target supplied-parameter construction and reports both orientations in its
diagnostics.  Input-market transitions are accepted only where a verified
target construction and input-market restrictions exist; unsupported cases
fail explicitly rather than applying output-market bounds.

## Updating and respecifying

`update(fit, ...)` is a genuine recalibration operation.  The fit stores a
canonical baseline `calibrate()` call, and `update()` replaces named observed
inputs or specification arguments before dispatching to `calibrate()` again.
This means a conduct change uses the target conduct's calibration equations;
it does not coerce the source S4 object.  A fit created only by `specify()` has
no observed-data identifying call and therefore rejects `update()`.

`respecify(fit, ...)` is narrower.  An explicit transition registry allows
standard Logit and CES Bertrand/Cournot transitions, and tested nested
Logit/CES translations where the target primitives and domain are valid.
Same-demand conduct transitions retain the portable demand primitives, while
flat and nested Logit/CES transitions use a separate deterministic demand
translation: baseline prices, ownership, active products, and accounting are
retained; target mean values are solved analytically to match baseline shares;
and target curvature is either retained or supplied explicitly. The target
supply state and marginal costs are then reconstructed through `specify()`
without using source margins to recalibrate demand.

This is a local translation, not a claim that price-level Logit and
log-price CES are globally equivalent.  Their counterfactual predictions may
diverge away from the preserved baseline.  Arbitrary demand-family,
variant, or specialized-model conversions remain rejected until an
economically valid transition rule and target construction path are
established.  Transition metadata records retained, recomputed, and
invalidated quantities, while local-translation diagnostics also report share
and quantity discrepancies, deterministic parameter mappings, output/input
orientation, and local matrix discrepancies.

`respecify()` is not a calibration operation.  Its result keeps the source
calibration call only as provenance in `source_calibration_args` and clears
`calibration_args`; consequently `update(respecify_fit)` fails rather than
silently recalibrating the target model from source margins.  Use `update()`
when the intended question is a fresh calibration under the target
specification.

### Transition taxonomy

The transition registry distinguishes structural restrictions, algebraic
translations, conditional translations, first-order Linear representations,
and first-order LogLin representations. Structural restrictions can discard
nesting while retaining portable primitives. Conditional translations require
explicitly supplied nests and nesting parameters. Algebraic Logit/CES
translations require the target curvature when it is not retained. First-order
Linear and LogLin translations derive slopes and intercepts analytically from
the fitted baseline elasticity matrix. No transition silently estimates a
missing primitive from margins; if a target model requires calibration, use
`update()` instead.

## Simulation

## Counterfactual boundary

`counterfactual()` is a lightweight container for simultaneous post-fit
changes. It stores only supplied environment fields such as `ownership`,
`costs`, `exit`, `capacity`, `bargaining`, and Stackelberg leadership/product
state. Model identity fields (`demand`, `conduct`, and `variant`) are rejected
there and remain the responsibility of `update()` or `respecify()`.

`simulate(fit, cf)` validates the requested fields against registry-derived
capabilities, translates them to the existing model-specific slots, and
performs one equilibrium solve. Legacy argument calls remain supported and
are translated internally. Results retain lightweight `counterfactual`
metadata as an attribute; the original fit is never mutated. The same
counterfactual can therefore be reused with multiple compatible fits.

`simulate(fit, ownerPost, mcDelta, subset, ...)` applies the scenario to the
fitted model, recovers post-counterfactual marginal costs with the model's
existing `calcMC()` method, and calls that model's existing equilibrium
method.  The dispatch preserves important differences:

* Logit/CES Bertrand use their existing nonlinear or AG price solvers.
* Logit/CES Cournot use their existing quantity-game price methods.
* Linear and log-linear general Cournot and Stackelberg recalculate plant
  quantities with their existing quantity solvers before calculating prices.
  Their plant-level cost and (for Stackelberg) leader-status conventions are
  retained, and they are not exposed through `specify()`.
* Bargaining retains its bargaining-specific equilibrium equations and solver
  choice, including AG support where the legacy constructor provides it.
* Second-score auction and second-score bargaining retain their direct or
  model-specific equilibrium methods. The capacity-constrained auction keeps
  its seller-cost distribution, buyer-valuation, reserve, and capacity
  calculations together.
* AIDS and PCAIDS recompute their ownership-dependent price-delta equation at
  simulation time; nested PCAIDS retains its model-specific nesting
  calibration.
* Logit vertical bargaining updates upstream/downstream ownership,
  integration-dependent bargaining power, and two-sided cost deltas before
  invoking its existing vertical price system. The downstream second-score
  entry retains the separate `VertBarg2ndLogit` price and share methods.
  Nested vertical second-score remains legacy-only.

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
