# BLP calibration and integration audit

This report records the BLP-only work on `refactor`. The existing S4 demand and
conduct classes remain the economic implementation; the refactor adds a shared
node/weight boundary and an observed-data calibration adapter.

## Calibration contract

The price-only calibrated model is:

```text
alpha_r = alphaMean + sigma * z_r,
z_r ~ N(0, 1), sigma >= 0
```

The user supplies prices, unconditional product shares, proportional margins,
pre-merger ownership, and known `s0`. The implementation requires
`sum(shares) = 1 - s0` and does not estimate or rescale `s0`. For each candidate
`(alphaMean, sigma)`, BLP contraction recovers `delta`; the conduct-specific
margin routine then supplies the minimum-distance moments.

All four calibrated paths use the same fixed integration points and weights:

| conduct | class | demand object in the conduct calculation | observed/predicted margin units |
|---|---|---|---|
| Bertrand | `LogitBLP` | aggregate random-coefficient Jacobian `sum_r w_r alpha_r (diag(s_r) - s_r s_r')` | proportional |
| Cournot | `CournotBLP` | inverse aggregate demand Jacobian with the existing Cournot ownership FOC | proportional |
| second-score | `Auction2ndBLP` | draw-level firm winning probabilities and `log(1-S_Fr)/alpha_r`, integrated before the product margin is formed | proportional in the calibration API; level is available through `calcMargins(level = TRUE)` |
| bargaining | `BargainingBLP` | integrated draw-level bargaining kernel and surplus term, solved as one linear system | proportional |

The earlier effective-alpha shortcut was not retained for auction or bargaining.
The sigma-zero tests match the corresponding homogeneous Logit methods, and the
heterogeneous conduct regressions use independent draw-level evaluators in the
test files.

## Integration contract

`R/BLPIntegration.R` centralizes the representation as:

```text
points z_r
normalized weights w_r, sum(w_r) = 1
```

Supported rules are `auto`, `gauss-hermite`, `monte-carlo`, and `provided`.
Price-only `auto` selects Gauss-Hermite; multidimensional legacy BLP selects
Monte Carlo. Explicit Gauss-Hermite for multidimensional specifications fails
clearly. Legacy `sim()` calls retain fixed-draw Monte Carlo unless a new rule or
integration input is explicit. Provided points and weights do not consume the
caller RNG state; Monte Carlo draws are generated once and stored in the model.

The same weights are used for contraction, shares, derivatives, elasticities,
CV, and PriceLeadershipBLP mean-value recovery. CV trimming uses a weighted
empirical quantile and renormalizes surviving weights.

## Quadrature benchmark

The synthetic asymmetric market in `artifacts/blp_integration_benchmark.R`
was compared with high-accuracy `stats::integrate()` calculations. Demand-level
errors were:

| nodes | max share error | max derivative error |
|---:|---:|---:|
| 10 | 6.15e-08 | 3.88e-05 |
| 15 | 5.24e-08 | 1.14e-06 |
| 20 | 8.22e-10 | 6.17e-08 |
| 30 | 8.02e-13 | 6.17e-10 |
| 40 | 1.42e-13 | 9.54e-12 |

The full Bertrand calibration benchmark used direct-integration shares and
margins with a post-merger ownership change:

| nodes | alpha error | sigma error | post prices | CV | calibration seconds |
|---:|---:|---:|---|---:|---:|
| 10 | -5.07e-04 | 2.71e-04 | 1.558679, 2.828130, 2.600961 | 0.0561913 | 26.35 |
| 20 | 2.66e-05 | -2.46e-05 | 1.558661, 2.828050, 2.600961 | 0.0561781 | 29.23 |
| 30 | 2.58e-05 | -2.41e-05 | 1.558661, 2.828050, 2.600961 | 0.0561781 | 25.81 |
| 40 | 2.33e-05 | -2.18e-05 | 1.558662, 2.828050, 2.600961 | 0.0561781 | 28.51 |

The current default remains 31 nodes, which is between the 30- and 40-node
benchmarks and provides sub-`1e-9` derivative error in this test while keeping
runtime below the 40-node run. The benchmark is representative rather than a
claim of universal accuracy for every price scale or heterogeneity variance.

Fixed-draw Monte Carlo comparisons in the same benchmark produced the following
demand errors against direct integration:

| fixed draws | max share error | max derivative error |
|---:|---:|---:|
| 1,000 | 3.26e-04 | 2.54e-03 |
| 10,000 | 2.47e-03 | 6.21e-03 |

The script uses a fixed seed and is retained so the MC figures can be rerun on
the target machine without conflating machine timing with package behavior.

## Recovery checks

The non-circular recovery fixture uses `alphaMean = -5`, `sigma = 0.8`, an
asymmetric three-product market, fixed supplied nodes, and four conduct paths.
Representative recovered values were:

| conduct | true alphaMean | recovered alphaMean | true sigma | recovered sigma |
|---|---:|---:|---:|---:|
| Bertrand | -5.000000 | -4.999954 | 0.800000 | 0.799977 |
| Cournot | -5.000000 | -5.000034 | 0.800000 | 0.800018 |
| second-score | -5.000000 | -4.999836 | 0.800000 | 0.799906 |
| bargaining | -5.000000 | -4.999371 | 0.800000 | 0.799207 |

The same suite includes the `sigma = 0` boundary and diagnostics for the
contraction residual, margin RMSE, starts, boundary status, wrong-sign
probability, and the sigma profile. `preMergerFOCResidual` records the maximum
conduct-equation residual in the proportional-margin units used by calibration.
No optimizer Hessian standard errors are reported.

## Limitations preserved

* Observed-data calibration is limited to one normal random coefficient on
  price with known `s0`.
* Bargaining power is supplied, not estimated.
* Gauss-Hermite is limited to one-dimensional price heterogeneity.
* Existing multidimensional BLP paths retain Monte Carlo behavior.
* Legacy S4 model classes, equilibrium solvers, normalizations, and APIs remain
  in place.
