# Structural calibration/simulation refactor plan

This branch starts from the current `master` S4 implementation.  The package's
legacy constructors both recover structural parameters/costs and solve the
post-counterfactual equilibrium; the new API will wrap those existing model
methods without combining demand and conduct mathematics.

## Repository-grounded sequence

1. **Oracle coverage.**  Record stable outputs from representative Logit
   Bertrand and Logit Cournot markets, including parameters, cost recovery,
   equilibria, shares, margins, elasticities, diversion, UPP, CMCR, and CV.
   Run the existing full test suite before production changes.
2. **Registry.**  Add normalized model specifications and one internal table
   for the currently supported demand/conduct combinations.  Existing
   constructors and `sim()` remain unchanged at this point.
3. **Prototype.**  Add an `AntitrustFit` S4 wrapper plus `calibrate()` and
   `simulate()` for Logit-Bertrand and Logit-Cournot.  Calibration delegates to
   the existing model-specific constructors; simulation clones the calibrated
   S4 object, applies the scenario, recovers post-merger costs, and calls the
   existing price solver.  This preserves the distinct Bertrand and Cournot
   methods and keeps result objects in their existing S4 classes.
4. **Compatibility.**  Add `specify()` for supplied structural parameters and
   migrate legacy Logit constructors and `sim()` only after direct parity tests
   pass.  Preserve legacy argument names, defaults, warnings, classes, and
   result slots.
5. **Incremental migration.**  Extend the same dispatch boundary one coherent
   registered family at a time.  Each family gets focused parity tests, a full
   suite run, and its own commit.  Structurally distinct supply models are not
   mechanically combined with demand models.
6. **Documentation/cleanup.**  Once parity is established, make the registry
   authoritative for dispatch and supported-model reporting, document both
   workflows and repeated simulations, and record intentionally preserved
   historical behavior or suspected bugs in a migration report.

## Initial design decision

The fitted state will use a small `AntitrustFit` S4 wrapper around an existing
S4 model object.  This retains all current downstream methods (`summary`,
`elast`, `diversion`, `upp`, `cmcr`, `CV`, and so on) on simulation results,
while making the pre-counterfactual state explicit.  The wrapper will carry a
simple normalized model specification and diagnostics; it will not introduce a
new economic supply abstraction.
