## Verify PLE Implementation
## Quick check that the implementation is working correctly

library(antitrust)

cat("\n=== PLE Implementation Verification ===\n\n")

## Simple test case
prices <- c(10, 12, 11, 9)
shares <- c(0.25, 0.20, 0.18, 0.17)
margins <- c(0.40, NA, 0.35, 0.25)  # Coalition: 1,2,3 with margins 0.40, NA, 0.35; Fringe: 4 with 0.25

ownerPre <- c("A", "A", "B", "C")
coalitionPre <- c(1, 2, 3)

cat("Test Case:\n")
cat("  Prices:    ", prices, "\n")
cat("  Shares:    ", shares, "\n")
cat("  Margins:   ", margins, "\n")
cat("  Coalition: ", coalitionPre, "\n")
cat("  Fringe:    ", 4, "\n\n")

## Test 1: Object creation
cat("Test 1: Creating PLE object...\n")
result <- tryCatch({
  ple(
    prices = prices,
    shares = shares,
    margins = margins,
    ownerPre = ownerPre,
    ownerPost = ownerPre,  # No merger
    coalitionPre = coalitionPre,
    coalitionPost = coalitionPre,  # Same coalition post-merger
    insideSize = 10000
  )
}, error = function(e) {
  cat("  ✗ Error:", e$message, "\n")
  traceback()
  NULL
})

if(is.null(result)){
  cat("\n✗ FAILED: Could not create PLE object\n")
  quit(status = 1)
}

cat("  ✓ PLE object created\n\n")

## Test 2: Check calibrated parameters
cat("Test 2: Checking calibrated parameters...\n")
cat("  Alpha (price coef):      ", round(result@slopes$alpha, 4), "\n")
cat("  Supermarkup Pre:         ", round(result@supermarkupPre, 4), "\n")
cat("  Supermarkup Post:        ", round(result@supermarkupPost, 4), "\n")
cat("  Timing parameter:        ", result@timingParam, "\n")
cat("  Binding firm:            ", result@bindingFirm, "\n\n")

## Test 3: Marginal costs
cat("Test 3: Marginal costs...\n")
cat("  MC Pre:  ", round(result@mcPre, 2), "\n")
cat("  MC Post: ", round(result@mcPost, 2), "\n\n")

## Test 4: Different pricing regimes
cat("Test 4: Testing different regimes...\n\n")

# Bertrand (no coordination)
cat("  Bertrand (no coordination):\n")
prices_bert <- tryCatch({
  calcPrices(result, preMerger = TRUE, regime = "bertrand")
}, error = function(e) {
  cat("    ✗ Error:", e$message, "\n")
  NULL
})
if(!is.null(prices_bert)){
  cat("    Prices: ", round(prices_bert, 2), "\n")
  cat("    ✓ Success\n\n")
}

# Coordination (coalition coordinates, no supermarkup)
cat("  Coordination (coalition coordinates, m=0):\n")
prices_coord <- tryCatch({
  calcPrices(result, preMerger = TRUE, regime = "coordination")
}, error = function(e) {
  cat("    ✗ Error:", e$message, "\n")
  NULL
})
if(!is.null(prices_coord)){
  cat("    Prices: ", round(prices_coord, 2), "\n")
  cat("    ✓ Success\n\n")
}

# Constrained (includes supermarkup)
cat("  Constrained (includes supermarkup):\n")
prices_const <- tryCatch({
  calcPrices(result, preMerger = TRUE, regime = "constrained")
}, error = function(e) {
  cat("    ✗ Error:", e$message, "\n")
  NULL
})
if(!is.null(prices_const)){
  cat("    Prices: ", round(prices_const, 2), "\n")
  cat("    ✓ Success\n\n")
}

## Test 5: Pre vs Post merger
cat("Test 5: Pre vs Post merger prices...\n")
cat("  Pre-merger:  ", round(result@pricePre, 2), "\n")
cat("  Post-merger: ", round(result@pricePost, 2), "\n\n")

## Test 6: Verify regime logic
if(!is.null(prices_bert) && !is.null(prices_coord) && !is.null(prices_const)){
  cat("Test 6: Verifying regime relationships...\n")
  
  # Coalition prices should satisfy: bertrand < coordination < constrained
  # (assuming positive supermarkup)
  
  coalition_diff_coord_bert <- mean(prices_coord[coalitionPre]) - mean(prices_bert[coalitionPre])
  coalition_diff_const_coord <- mean(prices_const[coalitionPre]) - mean(prices_coord[coalitionPre])
  
  cat("  Coalition avg price (Bertrand):     ", round(mean(prices_bert[coalitionPre]), 2), "\n")
  cat("  Coalition avg price (Coordination): ", round(mean(prices_coord[coalitionPre]), 2), "\n")
  cat("  Coalition avg price (Constrained):  ", round(mean(prices_const[coalitionPre]), 2), "\n\n")
  
  cat("  Coordination - Bertrand:  ", round(coalition_diff_coord_bert, 2), "\n")
  cat("  Constrained - Coordination: ", round(coalition_diff_const_coord, 2), "\n")
  
  if(result@supermarkupPre > 0){
    if(coalition_diff_const_coord > 0){
      cat("  ✓ Constrained > Coordination (supermarkup adds to price)\n")
    } else {
      cat("  ✗ Warning: Constrained should be > Coordination when supermarkup > 0\n")
    }
  }
  
  cat("\n")
}

cat("=== All Tests Complete ===\n\n")
