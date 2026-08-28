## Submission Summary (v0.99.34)
This patch release fixes an algebraic sign discrepancy in the CES Cournot margin equations (`R/MarginsMethods.R` and `R/ParamsMethods.R`), ensuring market shares correctly scale Cournot equilibrium margins as specified in the vignette equations.

## Test Environments
* Local Ubuntu 22.04.5 LTS, R 4.6.1
* Local R-devel (`R Under development (unstable) 2026-06-06 r90114`)

## R CMD check Results
There were 0 ERRORS and 0 WARNINGS.

There was 1 NOTE:
* Standard local environment note regarding missing system 'tidy' utility.

## Downstream Dependencies
None affected.
