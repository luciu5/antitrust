# antitrust package audit report

Audit date: 2026-06-01

## Formal QA expansion (2026-08-30)

The smoke harness has been supplemented with a reviewable `testthat` suite. The
suite uses synthetic fixtures only and currently records 221 successful
expectations in 41 tests. It covers every documented `sim()` demand/supply
combination, core S4 output methods, Cournot/Stackelberg, capacity auctions,
vertical bargaining, stochastic BLP paths, adversarial validation, and an
explicit API matrix (`api_coverage.csv`). Every matrix row is either linked to
a formal test or has a non-empty documented exclusion reason.

Run the complete evidence-producing workflow with:

```sh
Rscript ai/package_audit/run_qa.R . ai/package_audit/latest
```

The runner records the commit, package tarball hash (or the build disposition
when a local build times out), R/dependency/platform/BLAS information, RNG
kind, fixture checksums, expectation results, all warnings, and `sessionInfo()`.
Warnings are deliberately non-fatal, but are never discarded: the warning
classification policy is in `warning_dispositions.csv` and every run writes a
machine-readable `warnings.csv` plus a Markdown report.

### Independent numerical oracles

| Family or output | Independent check | Default tolerance |
| --- | --- | --- |
| Logit | `s_i = exp(delta_i + alpha p_i)/(1 + sum exp(...))`; Bertrand margin `-1/(alpha p_i(1-s_i))`; revenue/quantity round trip | `1e-7` to `1e-8` |
| CES | Revenue weights proportional to `meanval_i p_i^(1-gamma)` and quantity shares obtained by dividing by price | `1e-8` |
| Linear | Hand-coded diversion-to-slope map `B_ji/B_jj` and pre-merger margin recovery | `1e-8` |
| AIDS/PCAIDS | Pre-merger price-share transformation and finite calibrated parameter surface | `1e-7` |
| HHI | Common-ownership aggregation followed by squared percentage shares | `1e-10` |
| Equilibrium | Pre/post accounting, no-merger and zero-efficiency identities, permutation/label/unit scaling, subset exclusion, and solver parity | scale-specific (`1e-3` for scaled prices; otherwise `1e-4` or tighter) |

These checks establish software consistency on controlled fixtures. They do
not establish data quality, identification in an empirical market, economic
model validity, or any legal conclusion. Those remain independent economic and
legal review questions, consistent with the DOJ/FTC evidence and evaluation
guidance and Federal Rule of Evidence 702.

### Environment and version policy

`renv.lock` pins the canonical R 4.6.1 dependency set used for the recorded
fixture run. `.github/workflows/qa.yml` runs the formal tests on the minimum
compatibility line selected for this audit (`oldrel-1`), current release, and R-devel, and runs
`R CMD check --as-cran` on the canonical line. CI artifacts retain the same
metadata and warning evidence that a reviewer can rerun from a clean checkout.
The final local staged-source run also passed `R CMD build --no-build-vignettes
--no-manual --no-resave-data` and `R CMD check --no-manual --no-vignettes` with
status OK; the canonical CI job adds `--as-cran` where CRAN connectivity is
available.

## Scope

This pass inventories the exported R package surface and adds a smoke harness for representative constructors, helpers, and S4 methods. The generated API inventory is `ai/package_audit/api_coverage.csv`; the executable smoke script is `ai/examples/package_audit_smoke.R`.

## Findings

| Severity | Location | Affected API | Reproduction | Expected | Actual | Proposed fix | Validation |
| --- | --- | --- | --- | --- | --- | --- | --- |
| High | `R/CMCRCournotFunctions.R` | `upp.cournot()` | `upp.cournot(1, c(.2, .3), c("A", "B"))` | Return a two-product UPP vector by delegating to `upp.bertrand()` with a 2x2 diversion matrix. | Failed before this patch with `only matrix diagonals can be replaced`; after changing `rep()` to `matrix()`, positional argument matching passed `mcDelta` as `output`, causing a second error. | Fixed: build `diversions` with `matrix(1, ncol = 2, nrow = 2)` and call `upp.bertrand()` with explicit `ownerPost`, `mcDelta`, and `labels` names. | `Rscript -e "pkgload::load_all('.', quiet=TRUE); upp.cournot(1, c(.2,.3), c('A','B'))"` |
| High | `R/ParamsMethods.R` / `R/LinearFunctions.R` | `linear()` / `calcSlopes,Linear-method` | `timeout 5s Rscript -e "pkgload::load_all('.', quiet=TRUE); prices<-c(1,1.1,1.25); quantities<-c(22,18,12); margins<-c(.35,.3,.28); div<-matrix(c(-1,.55,.45,.5,-1,.5,.4,.6,-1),3,byrow=TRUE); linear(prices,quantities,margins,div, ownerPre=c('A','B','C'), ownerPost=c('A','A','C'), control.slopes=list(maxit=50,tol=1e-6))"` | Small linear examples should either calibrate quickly or fail with a clear validation error before entering the optimizer. | The 3-product smoke case exceeded a 5-second timeout; a 2-product case failed with `initial value is not in the interior of the feasible region`. | Audit optimizer starts and constraints in `calcSlopes,Linear-method`; validate or transform starts so invalid interior points fail immediately with a package-level message. | Add a focused regression under `ai/examples/` once the intended feasible parameterization is confirmed. |
| Fixed | `R/ParamsMethods.R`, `NAMESPACE` | `getNestsParms()` | `findMethods("getNestsParms", where=asNamespace("antitrust"))` | The helper should work for nested PCAIDS, nested Logit, nested auction Logit, nested CES, and vertical bargaining objects with nested downstream demand. | Previously the generic was exported, but only `PCAIDSNests` had a method. Calling it on `CESNests` errored with no inherited method. | Fixed: added methods for `LogitNests`, `Auction2ndLogitNests`, `CESNests`, `VertBargBertLogitNests`, and `VertBarg2ndLogitNests`. `LogitNestsALM` inherits the `LogitNests` method. | `Rscript ai/examples/package_audit_smoke.R` |
| Medium | `R/ParamsMethods.R` | `ces.nests()` / `calcSlopes,CESNests-method` | `Rscript ai/examples/package_audit_smoke.R` | A small valid nested CES example should calibrate cleanly or flag identification issues before optimization. | Smoke run completes but warns that singleton nests are normalized and that the optimizer may not have found a good solution. | Improve examples to avoid singleton nests; consider stricter validation or clearer warning text when nest parameters are unidentified. | `Rscript ai/examples/package_audit_smoke.R` |
| Low | `R/Auction2nd*`, `R/BargainingLogitFunctions.R` | `control.slopes` handling | Passing a shared control list with both `tol` and `reltol` to auction/bargaining constructors. | Optimizer controls should be documented per solver, or unsupported names should be ignored consistently. | `optimize()` requires `reltol`, while `BBoptim()` rejects `tol`/`reltol` combinations in some paths. | Document solver-specific control lists; optionally sanitize controls before passing to each optimizer. | Covered by the smoke script using solver-specific controls. |

## Fixes Applied

- Fixed `upp.cournot()` diversion matrix construction.
- Fixed `upp.cournot()` delegation to `upp.bertrand()` by using explicit argument names.
- Added `getNestsParms()` methods for nested Logit, nested auction Logit, nested CES, and nested vertical bargaining objects.
- Added `ai/examples/package_audit_smoke.R`, a representative smoke harness covering Logit, AIDS/PCAIDS, CES, nested/capacity Logit, auction, bargaining, standalone HHI/CMCR/UPP helpers, and `sim()`.
- Added `ai/package_audit/generate_api_inventory.R` to regenerate `api_coverage.csv`.

## Remaining Work

Prioritize the `linear()` calibration failure/timeout next. It is a correctness and usability risk because these are small examples and, per audit expectations, they should run quickly.
