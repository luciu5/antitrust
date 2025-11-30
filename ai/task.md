# Task List

- [x] Locate and read `LOGITBLP` class definition <!-- id: 0 -->
- [x] Identify elasticity calculation methods <!-- id: 1 -->
- [x] Verify product elasticity formula <!-- id: 2 -->
    - [x] Create reproduction script `reproduce_issue.R`
    - [x] Run script and analyze results
- [x] Describe the issue to the user <!-- id: 5 -->
- [x] Fix product elasticity formula <!-- id: 6 -->
- [x] Verify market elasticity formula <!-- id: 3 -->
    - [x] Update `elast` method to calculate market elasticity by aggregating product elasticities
- [x] Report findings to user <!-- id: 4 -->

## LogitBLP Method Audit
- [x] Identify all `LogitBLP` methods <!-- id: 7 -->
- [x] Audit `calcShares` <!-- id: 8 -->
    - [x] Check if it exists. (Found in `OutputMethods.R`)
    - [x] Verify math: Correctly handles random coefficients and `sigmaNest`.
- [x] Audit `calcSlopes` <!-- id: 9 -->
    - [x] Verify BLP contraction mapping logic. (Correct)
- [x] Audit `calcPrices` <!-- id: 10 -->
    - [x] Verify FOCs use the correct derivatives. (Found BUG in Jacobian)
    - [x] Fix Jacobian in `calcPrices` to use correct Nested Logit derivative.
    - [x] Verify fix with reproduction script.
- [x] Audit `calcMC` <!-- id: 11 -->
    - [x] Verify MC calculation uses correct derivatives. (Inherits from `Bertrand`, uses correct `elast`).

## Inherited Method Audit
- [x] Audit `CV` (Consumer Surplus) <!-- id: 12 -->
    - [x] Create `reproduce_cv.R` to check formula.
    - [x] Verify if `sigmaNest` multiplier is correct. (Found BUG: extra multiplier)
    - [x] Fix if incorrect.
- [x] Audit `calcMargins` <!-- id: 13 -->
    - [x] Confirmed it inherits from `Bertrand` and uses `elast`.
- [x] Audit `diversion` <!-- id: 14 -->
    - [x] Check `DiversionMethods.R`. (Inherits from `Bertrand`, uses correct `elast` and `calcShares`).

## Final Verification
- [x] Verify Elasticity Aggregation <!-- id: 15 -->
    - [x] User questioned formula correctness.
    - [x] Found BUG: `elast` was averaging elasticities instead of share derivatives.
    - [x] Fixed `ElastMethods.R` to aggregate derivatives.
    - [x] Verified with `reproduce_issue.R` using heterogeneous alphas.

## Organization
- [x] Move artifacts and scripts to `ai` directory <!-- id: 16 -->
    - [x] Create `ai` directory.
    - [x] Move `reproduce_issue.R` and `reproduce_cv.R`.
    - [x] Copy markdown artifacts.
    - [x] Update `.Rbuildignore` to ignore `ai`.
    - [x] Ensure `ai` is tracked by git (checked `.gitignore`).
    - [x] Move `inst` directory contents (`doc`, `examples`) to `ai`.
