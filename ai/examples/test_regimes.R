## Test Regime Functionality
## Tests the 4 pricing regimes in calcMargins and calcPrices

library(antitrust)

set.seed(42)

## Simple 4-product market
prices <- c(10, 12, 11, 9)
shares <- c(0.25, 0.20, 0.18, 0.17)
margins_pct <- c(0.40, NA, 0.35, 0.25)

ownerPre <- c("A", "A", "B", "C")  
coalitionPre <- c(1, 2, 3)

cat("\n=== Testing Price Leadership Regimes ===\n\n")

## Create PLE object
result <- tryCatch({
  ple(
    prices = prices,
    shares = shares,
    margins = margins_pct,
    ownerPre = ownerPre,
    coalitionPre = coalitionPre,
    insideSize = 10000
  )
}, error = function(e) {
  cat("Error creating PLE object:", e$message, "\n")
  NULL
})

if(!is.null(result)){
  cat("✓ PLE object created successfully\n\n")
  
  ## Test Regime 1: Bertrand
  cat("Testing Regime 1: Bertrand (no coordination)...\n")
  margins_bert <- tryCatch({
    calcMargins(result, preMerger = TRUE, regime = "bertrand")
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    NULL
  })
  
  if(!is.null(margins_bert)){
    cat("  ✓ Bertrand margins calculated\n")
    cat("  Margins:", round(margins_bert, 4), "\n\n")
  }
  
  ## Test Regime 2: Deviation
  cat("Testing Regime 2: Deviation...\n")
  margins_dev <- tryCatch({
    calcMargins(result, preMerger = TRUE, regime = "deviation")
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    NULL
  })
  
  if(!is.null(margins_dev)){
    cat("  ✓ Deviation margins calculated\n")
    cat("  Margins:", round(margins_dev, 4), "\n\n")
  }
  
  ## Test Regime 3: Coordination
  cat("Testing Regime 3: Coordination...\n")
  margins_coord <- tryCatch({
    calcMargins(result, preMerger = TRUE, regime = "coordination")
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    NULL
  })
  
  if(!is.null(margins_coord)){
    cat("  ✓ Coordination margins calculated\n")
    cat("  Margins:", round(margins_coord, 4), "\n\n")
  }
  
  ## Test Regime 4: Constrained
  cat("Testing Regime 4: Constrained (with supermarkup)...\n")
  margins_const <- tryCatch({
    calcMargins(result, preMerger = TRUE, regime = "constrained")
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    NULL
  })
  
  if(!is.null(margins_const)){
    cat("  ✓ Constrained margins calculated\n")
    cat("  Margins:", round(margins_const, 4), "\n\n")
  }
  
  ## Test calcPrices with different regimes
  cat("Testing calcPrices with different regimes...\n\n")
  
  cat("  Bertrand prices:\n")
  prices_bert <- tryCatch({
    calcPrices(result, preMerger = TRUE, regime = "bertrand")
  }, error = function(e) {
    cat("    Error:", e$message, "\n")
    NULL
  })
  if(!is.null(prices_bert)){
    cat("    ✓ Success:", round(prices_bert, 2), "\n\n")
  }
  
  cat("  Coordination prices (default):\n")
  prices_coord <- tryCatch({
    calcPrices(result, preMerger = TRUE)
  }, error = function(e) {
    cat("    Error:", e$message, "\n")
    NULL
  })
  if(!is.null(prices_coord)){
    cat("    ✓ Success:", round(prices_coord, 2), "\n\n")
  }
  
  ## Comparison
  if(!is.null(margins_bert) && !is.null(margins_coord)){
    cat("Comparison of Regimes:\n")
    cat("  Bertrand vs Coordination margin difference:\n")
    diff <- margins_coord - margins_bert
    cat("    Products 1-3 (coalition):", round(diff[1:3], 4), "\n")
    cat("    Product 4 (fringe):      ", round(diff[4], 4), "\n")
    cat("  → Coalition margins should be higher with coordination\n\n")
  }
  
  if(!is.null(prices_bert) && !is.null(prices_coord)){
    cat("  Bertrand vs Coordination price difference:\n")
    diff <- prices_coord - prices_bert
    cat("    Products 1-3 (coalition):", round(diff[1:3], 2), "\n")
    cat("    Product 4 (fringe):      ", round(diff[4], 2), "\n")
    cat("  → Coalition prices should be higher with coordination\n\n")
  }
  
  cat("=== All Tests Complete ===\n")
  
} else {
  cat("✗ Failed to create PLE object\n")
}
