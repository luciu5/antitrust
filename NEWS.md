# antitrust 0.99.43.9000

- Preserve plain Logit auction value indices, incremental cost shocks, and active products when advancing sequential counterfactuals. This maintains the solved post-state as the next baseline without recalibration.
- Vectorize nested `elast()` for BLP over draws, and return zero instead of NaN for draws whose inside shares underflow.
- Solve input-market (signed-rate) BLP level FOCs and level margins with the analytic share Jacobian from `elast(partial = TRUE)` instead of a numerical Jacobian of `calcShares()`.
- Select the input-market level-FOC solver and level margins by inheritance rather than exact class name, so subclasses of `LogitBLP` (and `Logit`) are solved like their parent classes.
- Fix nested shares in the `calcMeanval()` BLP contraction, which omitted the within-nest factor `H^(sigmaNest - 1)`; nested models now recover mean utilities consistent with `calcShares()`.
- Compute BLP nested shares and log inclusive values with one max-shifted kernel in `calcShares()`, `calcSlopes()`, `calcMeanval()` and `CV()`, so draws with utilities beyond the range of `exp()` stay finite and calibration no longer clips utilities at +/-700.

# antitrust 0.99.43

- Keep homogeneous Cournot calibration inside its finite objective domain and scale optimization across price and quantity units.
- Separate calibration, structural specification, updating, model translation and simulation while retaining established S4 result methods.
- Add reusable counterfactual paths and ownership-aware synthetic markets.
- Correct signed input BLP equilibria, fractional ownership weights, rare-share convergence, and legacy heterogeneous draw reconstruction.
- Price leadership and repeated-game models are provided by coordination; vertical bargaining is provided by vertical.
