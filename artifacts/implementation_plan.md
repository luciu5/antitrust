# Implementation Plan: Adding CES extensions to sim() and CRAN Compilation

## Objective
1. **Fix package CRAN checks**: Address the remaining `R CMD check` warnings related to missing documentation for the `lim` argument in `R/CVMethods.R`. 
2. **Extend `sim()` function**: Integrate `CESCournot` and `Auction2ndCES` models into the `sim()` wrapper function in `R/SimFunctions.R` to allow users to specify these models through the generalized simulation function.

## Steps
1. **Fix CVMethods Documentation**
   - Modify `R/CVMethods.R` to explicitly document the `lim` argument in the `@param` tag of the roxygen block. Because `setMethod("CV", "LogitBLP", function(object, lim = c(0, 1)))` has `lim`, roxygen requires an explicit `@param lim` tag, otherwise `R CMD check` complains about undocumented arguments.

2. **Verify Check Pass**
   - Run `devtools::check(vignettes=FALSE)` and ensure no ERRORs or WARNINGs persist.

3. **Extend `SimFunctions.R`**
   - Identify the `sim()` wrapper function and update its `demand` parameter to include `"CESCournot"` and `"Auction2ndCES"`.
   - Add parsing logic to handle these new demand models properly. `Auction2ndCES` expects parameters identical to `auction2nd.ces()` (e.g. `priceStart`), while `CESCournot` expects arguments for `ces.cournot()`.
   - Integrate them into the data generation logic (`simulate_data()`) or the `calibrate_and_simulate()` workflows to allow testing and usage via standard `sim()` calls.

4. **Verify Implementation**
   - Execute `R CMD check` once more to ensure adding the new methods hasn't disrupted structural integrity or introduced missing parameter warnings for `sim()`.
   - Provide summary to the user upon success.
