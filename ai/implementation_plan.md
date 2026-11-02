# Fix LogitBLP Elasticity Calculations

The `elast` method for `LogitBLP` in `R/ElastMethods.R` contains incorrect formulas for share calculation and elasticity matrices when `sigmaNest` (nesting parameter) is not 1.

## User Review Required

> [!IMPORTANT]
> The current implementation assumes a specific (and likely incorrect) form for Nested Logit shares and elasticities when `sigmaNest < 1`. I will update it to use the standard Nested Logit formulas.

## Proposed Changes

### R

#### [MODIFY] [ElastMethods.R](file:///d:/Projects/antitrust/R/ElastMethods.R)
- Update the `elast` method for `LogitBLP`.
- Fix `sInd` calculation to correctly account for the inclusive value term: $s_{j} = s_{j|g} \cdot s_g$.
- Fix the elasticity matrix calculation to use the standard Nested Logit formulas.
- **[NEW]** Update `market=TRUE` logic to calculate market elasticity by aggregating the individual product elasticities (or using the consistent aggregate formula) rather than using the `mean(alphas)` approximation. This ensures consistency with the random coefficients and nesting.
    - Formula: $E_M = \sum_j s_j (\sum_k \epsilon_{jk})$ (share-weighted sum of row-sums of elasticity matrix).
- **[NEW]** Fix Aggregation Logic: The previous implementation averaged elasticities across draws ($\bar{E} = \frac{1}{R} \sum E_r$). This is incorrect. The correct approach is to average the share derivatives ($\frac{\partial s}{\partial p} = \frac{1}{R} \sum \frac{\partial s_r}{\partial p}$) and then compute elasticity ($E = \frac{p}{s} \frac{\partial s}{\partial p}$).

#### [MODIFY] [PricesMethods.R](file:///d:/Projects/antitrust/R/PricesMethods.R)
- Update `calcPrices` for `LogitBLP`.
- Fix the analytic Jacobian calculation `JAC`.
- The current implementation uses a simplified derivative $\frac{\alpha}{\sigma} s (1 - \sigma s)$ which is only correct for $\sigma=1$.
- Update to use the correct Nested Logit derivative: $\frac{\partial s_j}{\partial p_j} = \alpha s_j [ \frac{1}{\sigma} + s_{j|g} (1 - \frac{1}{\sigma} - s_g) ]$.
- This requires calculating $s_{j|g}$ (within-nest share) and $s_g$ (nest share) from the `shares_draw`.

#### [MODIFY] [CVMethods.R](file:///d:/Projects/antitrust/R/CVMethods.R)
- Audit `CV` method for `LogitBLP`.
- Remove incorrect `sigmaNest` multiplier from the Consumer Surplus calculation if confirmed to be a bug.
- Formula should be $CS = \frac{1}{\alpha} \ln(1 + (\sum e^{u/\sigma})^\sigma)$.

## Verification Plan

### Automated Tests
- [x] Run the reproduction script `reproduce_issue.R` (updated with the fix logic) to verify the new implementation matches the "True" logic.
- [x] Verify market elasticity is consistent with the product matrix.
- [x] **[NEW]** Verify `calcPrices` Jacobian by comparing numerical and analytic derivatives in the reproduction script. (Verified by successful price recovery)
- [x] **[NEW]** Verify `CV` calculation using `reproduce_cv.R`. (Verified by successful match with standard formula)

### Manual Verification
- None required beyond the script.
