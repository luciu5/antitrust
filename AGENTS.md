# Repository Guidelines

## Project Structure & Module Organization

This repository is an R package named `antitrust`. Core code lives in `R/`, with S4 classes, methods, and model-specific functions split by topic, for example `BertrandClasses.R`, `LogitFunctions.R`, and `HHIMethods.R`. Package metadata and load order are defined in `DESCRIPTION`; update `Collate` when adding files with class or method ordering requirements. Exported symbols are listed in `NAMESPACE`.

Generated help pages live in `man/` and should normally be refreshed from roxygen comments rather than edited by hand. Long-form material is in `vignettes/`. The `ai/` and `artifacts/` directories are excluded from package builds via `.Rbuildignore`; treat them as development notes and examples, not shipped source.

## Build, Test, and Development Commands

- `R CMD build .`: build the package tarball from the repository root.
- `R CMD check antitrust_*.tar.gz`: run the standard package check on a built tarball.
- `Rscript -e "devtools::document()"`: regenerate `NAMESPACE` and `man/` from roxygen comments, if `devtools` is installed.
- `Rscript -e "devtools::load_all()"`: load the package in place for interactive development.
- `Rscript ai/examples/quick_test.R`: run a local example script for touched model behavior.

## Coding Style & Naming Conventions

Use existing base-R/S4 style. Keep names consistent with public APIs such as `Logit*`, `CES*`, `calc*`, and `sim*` patterns. Use four-space indentation inside functions and avoid unrelated formatting churn. Prefer explicit argument names in numerical routines where silent positional mistakes would be costly. Keep roxygen documentation adjacent to exported functions.

## Testing Guidelines

There is no formal `tests/testthat/` suite in the current tree. Validate changes with `R CMD check` and targeted scripts under `ai/examples/`, especially for calibration, simulation, and price-leadership behavior. For bug fixes, add a small reproducible script or example when a package-level test is not practical. Check numerical changes against stable tolerances rather than exact floating-point equality.

## Commit & Pull Request Guidelines

Recent commits use short imperative or descriptive subjects, for example `Support input markets in bargaining logit` and `Fix sim BLP weights initialization`. Keep the first line concise and specific; mention affected models or methods when useful. Pull requests should describe the behavioral change, list validation commands, and note documentation updates. Include linked issues or reproduction scripts for bug fixes.

## Security & Configuration Tips

Do not commit local build outputs such as `antitrust_*.tar.gz`, `.Rproj.user/`, `.Rlib/`, `..Rcheck/`, or `antitrust.Rcheck/`. Keep package dependencies in `DESCRIPTION`; avoid relying on user-local library paths in committed code.
