## Benchmark for LogitBLP Performance
## Tests speed improvements from vectorization

library(antitrust)

# Set seed for reproducibility
set.seed(123)

cat("\n=== LogitBLP Performance Benchmark ===\n\n")

# Test with increasing product counts
product_counts <- c(4, 10, 20)

for (nprods in product_counts) {
  cat("Testing with", nprods, "products...\n")
  
  # Generate test data
  prices <- runif(nprods, 8, 12)
  shares <- runif(nprods, 0.05, 0.15)
  shares <- shares / sum(shares) * 0.8  # Scale to sum to 0.8
  margins <- runif(nprods, 0.2, 0.4)
  
  # Random ownership
  ownerPre <- sample(paste0("Firm_", 1:ceiling(nprods/2)), nprods, replace=TRUE)
  ownerPost <- ownerPre
  ownerPost[1:2] <- ownerPre[1]  # Merge first two products
  
  demand.param <- list(
    alpha = -0.5,
    sigma = 0.3,
    sigmaNest = 0.85,
    nDraws = 500
  )
  
  # Run simulation and time it
  start_time <- Sys.time()
  
  result <- sim(
    prices = prices,
    shares = shares,
    margins = margins,
    demand = "LogitBLP",
    demand.param = demand.param,
    supply = "bertrand",
    ownerPre = ownerPre,
    ownerPost = ownerPost,
    insideSize = 10000,
    labels = paste0("Prod_", 1:nprods)
  )
  
  # Calculate elasticities (this is where vectorization helps most)
  elast_pre <- elast(result, preMerger = TRUE)
  
  end_time <- Sys.time()
  elapsed <- as.numeric(end_time - start_time, units = "secs")
  
  cat("  Total time:", round(elapsed, 2), "seconds\n")
  cat("  Price changes:", paste(round(100 * (result@pricePost - result@pricePre) / result@pricePre, 2), collapse=", "), "\n\n")
}

cat("=== Benchmark Complete ===\n")
cat("\nNote: Vectorization improvements are most significant for:\n")
cat("  - Large product counts (20+ products)\n")
cat("  - High draw counts (1000+ draws)\n")
cat("  - Multiple elasticity calculations\n")
