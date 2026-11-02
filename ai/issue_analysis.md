# LogitBLP Elasticity Issue Analysis

## Summary
The `elast` method for `LogitBLP` in `R/ElastMethods.R` incorrectly calculates product shares and elasticities when the nesting parameter `sigmaNest` is less than 1. The implementation fails to account for the inclusive value term correctly in the share equation, which propagates to incorrect elasticity calculations.

## Mathematical Derivation

### 1. Share Calculation Error

**Standard Nested Logit Theory:**
For a single nest of "inside" goods (nest $g$) and an outside good (nest 0), the share of product $j$ is:
$$s_j = s_{j|g} \cdot s_g$$
Where:
- $V_j = \delta_j + \alpha(p_j - p_0)$
- Inclusive Value of the nest: $D_g = \sum_{k \in g} \exp(V_k / \sigma)$
- Conditional share: $s_{j|g} = \frac{\exp(V_j / \sigma)}{D_g}$
- Nest share: $s_g = \frac{D_g^\sigma}{1 + D_g^\sigma}$

Combining these:
$$s_j = \frac{\exp(V_j / \sigma)}{D_g} \cdot \frac{D_g^\sigma}{1 + D_g^\sigma} = \frac{\exp(V_j / \sigma) \cdot D_g^{\sigma-1}}{1 + D_g^\sigma}$$

**Current Implementation:**
The code calculates:
```r
expUtil = exp(util / sigmaNest)  # exp(V_j / sigma)
sumExpUtil = sum(expUtil)        # D_g
insideIV = sumExpUtil^sigmaNest  # D_g^sigma
denom = 1 + insideIV             # 1 + D_g^sigma
sInd = expUtil / denom           # exp(V_j / sigma) / (1 + D_g^sigma)
```
Mathematically, the code computes:
$$s_j^{code} = \frac{\exp(V_j / \sigma)}{1 + D_g^\sigma}$$

**The Discrepancy:**
The code is missing the term $D_g^{\sigma-1}$ in the numerator.
$$s_j^{true} = s_j^{code} \cdot D_g^{\sigma-1}$$

When $\sigma = 1$, $D_g^0 = 1$, so the formulas match (Standard Logit).
When $\sigma < 1$, the calculated shares are incorrect.

### 2. Elasticity Calculation Error

Because the share formula is incorrect, the elasticity derivatives derived from it are also incorrect.

**Standard Nested Logit Elasticities:**
- **Own-Price:** $\epsilon_{jj} = \alpha p_j \left[ \frac{1}{\sigma} - s_{j|g}(\frac{1}{\sigma} - 1) - s_j \right]$
  (Often written as: $\frac{\alpha p_j}{\sigma} [ 1 - \sigma s_j - (1-\sigma) s_{j|g} ]$)

- **Cross-Price (same nest):** $\epsilon_{jk} = -\alpha p_k \left[ s_{k|g}(\frac{1}{\sigma} - 1) + s_k \right]$

**Current Implementation:**
The code attempts a vectorized calculation that resembles the standard logit formula but adjusted for nesting in a way that doesn't align with the standard derivation for the one-nest case.

## Impact
- **Shares:** The model predicts incorrect market shares for a given set of parameters.
- **Elasticities:** Price sensitivities are miscalculated, leading to incorrect merger simulation results (e.g., predicted price effects).
- **Calibration:** If this method is used in calibration (BLP contraction), the recovered mean valuations ($\delta$) will be wrong.

## Proposed Fix
Replace the custom logic with the standard Nested Logit formulas for shares and elasticities.
