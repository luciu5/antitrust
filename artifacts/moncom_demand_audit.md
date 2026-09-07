# MonCom demand audit (antitrust/refactor)

This is an economics-first inventory of the demand families considered for
atomistic differentiated-product monopolistic competition.  The registry
continues to list complete implemented models; this table is not a proposal to
compose arbitrary demand and conduct modules.

The binding convention is

\[
q_j + (p_j-c_j)D^{MC}_j=0,
\]

where the relevant aggregate demand index, denominator, nest object, or
expenditure object is held fixed.  The direct derivative is derived from the
package's own demand equations rather than copied from a Bertrand ownership
system.

| Demand family | MonCom status | Direct perceived own derivative | Decision and reason |
| --- | --- | --- | --- |
| Flat Logit | Supported | \(D^{MC}_j=\alpha s_j\); output markup \(-1/\alpha\), input markup \(1/\alpha\) | The flat Logit numerator is \(\exp(\delta_j+\alpha p_j)\). Holding its denominator fixed gives the direct derivative. |
| Flat CES | Supported | \(d\log q_j/d\log p_j=-\gamma\); output markup \(p_j/\gamma\), input margin \((mc_j-p_j)/p_j=1/(-\gamma)\) | `calcShares()` implies \(q_j\propto p_j^{-\gamma}\) with the aggregate CES index fixed. This is deliberately not `diag(elast())`, whose full Marshallian value is \(-\gamma+(\gamma-1)s_j\). |
| Nested Logit | Unsupported for automatic lifecycle calibration | Holding the market denominator and nest inclusive value fixed gives \(D^{MC}_j=(\alpha/\sigma_g)s_j\) | The numerator slope is derivable, but observed MonCom margins identify only \(\sigma_g/\alpha\) within a nest. Without an explicit target alpha or nest parameter, automatic calibration would infer an unidentified primitive. No partial registry entry is advertised. |
| Nested CES | Unsupported for automatic lifecycle calibration | Holding the market and nest aggregates fixed gives \(d\log q_j/d\log p_j=-\sigma_g\) | The direct within-nest numerator slope is derivable, but the package's nested CES state has both within- and between-nest price indexes. A complete specified/calibrated MonCom path with explicit nest parameters is not yet validated, so the transition is rejected rather than guessed. |
| BLP / random-coefficients Logit | Supported for the validated price-random-coefficient integration path | \(D^{MC}_j=\sum_r w_r\alpha_r s_{jr}\) | The draw-level numerator derivative is integrated using the existing provided, Gauss-Hermite, and Monte Carlo integration engine. Ownership is intentionally absent from the FOC. The default Monte Carlo size remains 5,000. |
| Linear | Unsupported | Not adopted | Bertrand slope matrices cannot be reused as an atomistic FOC by identity, and the package has no separate validated MonCom Linear lifecycle. |
| LogLinear | Unsupported | Not adopted | The package's log-price slope/intercept representation has no validated atomistic conduct calibration and solver path. |
| AIDS | Unsupported | Not adopted | AIDS expenditure-share and market-expenditure responses are coupled; no fixed-expenditure-object MonCom calibration and solver is currently implemented. |
| PCAIDS / nested PCAIDS | Unsupported | Not adopted | PCAIDS is a calibrated Bertrand family with known-elasticity restrictions; it does not expose a validated atomistic own-product FOC. |
| LogitCap and other specialized descendants | Unsupported | Not adopted | Capacity, auction, bargaining, and vertical descendants have conduct-specific state and equations not audited for atomistic MonCom. |

## Implementation consequences

* `MonComLogit`, `MonComCES`, and `MonComBLP` are the only registered MonCom
  models on this head.
* `elast(CES)` is intentionally unchanged: it remains the full demand
  elasticity used by mature CES diagnostics and other conduct models.
  `MonComCES` uses a separate direct own derivative helper.
* Positive output-market CES calibration requires `gamma > 1`. Input-market
  CES uses the package's established orientation and requires `gamma < 0`,
  making the direct own derivative `-gamma` positive.
* BLP MonCom uses exactly the stored integration points and normalized weights;
  it does not create a second integration implementation or silently drop
  demographics/characteristic heterogeneity. The supported path is the
  existing validated price-random-coefficient engine.

