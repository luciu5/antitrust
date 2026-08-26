## Submission Summary (v0.99.33)
This patch release updates the validity method in `BertrandRUMClasses.R` to permit empty `weights` slots, ensuring 100% backward compatibility with reverse dependency `trade` (v0.8.3). It also includes R-devel documentation link anchor updates (`\code{\link[stats]{optim}}`, `\code{\link[stats]{constrOptim}}`).

## Test Environments
* Local Ubuntu 22.04.5 LTS, R 4.6.1
* Local R-devel (`R Under development (unstable) 2026-06-06 r90114`)

## R CMD check Results
There were 0 ERRORS and 0 WARNINGS.

There was 1 NOTE:
* Standard local environment note regarding missing system 'tidy' utility.

## Downstream Dependencies
Reverse dependency `trade` tested and verified compatible.
