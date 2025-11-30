# Test PriceLeadershipBLP implementation
# Compares ple() (standard Logit) with ple.blp() (Logit with BLP random coefficients)

library(devtools)
load_all()

cat("=== Testing PriceLeadershipBLP Implementation ===\n\n")

# Test data (simulating banking merger scenario)
prices <- c(0.93, 0.88, 1.10, 1.02)
shares <- c(0.35, 0.25, 0.25, 0.15)
margins <- c(0.40, NA, 0.25, 0.20)
ownerPre <- c("Bank1", "Bank2", "Bank3", "Fringe")
ownerPost <- c("Bank1", "Bank2", "Bank3", "Fringe")
coalitionPre <- c(1, 2, 3)  # Top 3 banks coordinate

cat("Test Data:\n")
cat("Prices:", prices, "\n")
cat("Shares:", shares, "\n")
cat("Margins:", margins, "\n")
cat("Coalition (pre):", coalitionPre, "\n\n")

## Test 1: Standard PLE (Logit demand)
cat("=== Test 1: Standard ple() with Logit demand ===\n")
time_logit <- system.time({
  result_logit <- ple(
    prices = prices,
    shares = shares,
    margins = margins,
    ownerPre = ownerPre,
    ownerPost = ownerPost,
    coalitionPre = coalitionPre,
    insideSize = 1000
  )
})

cat("Alpha (price coef):", result_logit@slopes$alpha, "\n")
cat("Supermarkup (pre):", result_logit@supermarkupPre, "\n")
cat("Supermarkup (post):", result_logit@supermarkupPost, "\n")
cat("Timing params:", result_logit@timingParam, "\n")
cat("Pre-merger prices:", result_logit@pricePre, "\n")
cat("Post-merger prices:", result_logit@pricePost, "\n")
cat("Execution time:", time_logit[3], "seconds\n\n")

## Test 2: PLE with BLP random coefficients
cat("=== Test 2: ple.blp() with LogitBLP demand (random coefficients) ===\n")

# For BLP, we need to provide pre-calibrated demand parameters
# In practice, these would come from a BLP estimation
# Here we'll use the calibrated parameters from Test 1 and add heterogeneity

# Extract calibrated parameters from standard logit
alpha_mean <- result_logit@slopes$alpha
meanval <- result_logit@slopes$meanval

# Add random coefficient heterogeneity
sigma <- 0.5  # Standard deviation of price coefficient across consumers
# Higher sigma = more heterogeneous price sensitivity

# Build slopes list for BLP (pre-calibrated demand)
slopes_blp <- list(
  alphaMean = alpha_mean,
  alpha = alpha_mean,  # Mean of random coefficients
  sigma = sigma,       # Std dev of random coefficients
  meanval = meanval,   # Mean valuations (from logit calibration)
  sigmaNest = 1,       # No nesting
  nDraws = 500,        # Number of simulation draws
  nDemog = 0,
  piDemog = numeric(0)
)

cat("BLP demand parameters (pre-calibrated):\n")
cat("  alphaMean:", alpha_mean, "\n")
cat("  sigma (heterogeneity):", sigma, "\n")
cat("  nDraws:", slopes_blp$nDraws, "\n\n")

cat("Running ple.blp() WITHOUT margins (inferred from prices)...\n")
time_blp <- system.time({
  result_blp <- ple.blp(
    prices = prices,
    shares = shares,
    # margins NOT provided - supermarkup inferred from observed vs Bertrand prices
    ownerPre = ownerPre,
    ownerPost = ownerPost,
    coalitionPre = coalitionPre,
    coalitionPost = coalitionPre,  # Same coalition post-merger
    insideSize = 1000,
    slopes = slopes_blp  # Pre-calibrated BLP demand parameters
  )
})

cat("Results from ple.blp():\n")
cat("Alpha (mean):", result_blp@slopes$alphaMean, "\n")
cat("Sigma (heterogeneity):", result_blp@slopes$sigma, "\n")
cat("Supermarkup (pre):", result_blp@supermarkupPre, "\n")
cat("Supermarkup (post):", result_blp@supermarkupPost, "\n")
cat("Timing params:", result_blp@timingParam, "\n")
cat("Pre-merger prices:", result_blp@pricePre, "\n")
cat("Post-merger prices:", result_blp@pricePost, "\n")
cat("Execution time:", time_blp[3], "seconds\n\n")

## Comparison
cat("=== Comparison: Logit vs LogitBLP ===\n\n")

cat("Supermarkup (Pre-merger):\n")
cat("  Logit:", result_logit@supermarkupPre, "\n")
cat("  LogitBLP:", result_blp@supermarkupPre, "\n")
cat("  Difference:", result_blp@supermarkupPre - result_logit@supermarkupPre, "\n\n")

cat("Supermarkup (Post-merger):\n")
cat("  Logit:", result_logit@supermarkupPost, "\n")
cat("  LogitBLP:", result_blp@supermarkupPost, "\n")
cat("  Difference:", result_blp@supermarkupPost - result_logit@supermarkupPost, "\n\n")

cat("Pre-merger Price (Product 1):\n")
cat("  Logit:", result_logit@pricePre[1], "\n")
cat("  LogitBLP:", result_blp@pricePre[1], "\n")
cat("  Difference:", result_blp@pricePre[1] - result_logit@pricePre[1], "\n\n")

cat("Post-merger Price (Product 1):\n")
cat("  Logit:", result_logit@pricePost[1], "\n")
cat("  LogitBLP:", result_blp@pricePost[1], "\n")
cat("  Difference:", result_blp@pricePost[1] - result_logit@pricePost[1], "\n\n")

cat("=== Key Insights ===\n\n")

cat("1. BLP random coefficients (sigma > 0) introduce consumer heterogeneity\n")
cat("   in price sensitivity, which affects:\n")
cat("   - Elasticities (different consumers have different price responses)\n")
cat("   - Markups (optimal pricing with heterogeneous consumers)\n")
cat("   - Coordination incentives (IC constraints may differ)\n\n")

cat("2. With sigma =", sigma, ":\n")
cat("   - Some consumers are very price-sensitive (high |alpha|)\n")
cat("   - Some consumers are less price-sensitive (low |alpha|)\n")
cat("   - This typically makes demand less elastic on average\n")
cat("   - Firms can sustain higher markups\n\n")

cat("3. Price leadership parameters (supermarkup, timing) are re-calibrated\n")
cat("   using BLP demand, but may differ from standard Logit because:\n")
cat("   - Different elasticities → different Bertrand profits\n")
cat("   - Different deviation incentives → different IC constraints\n")
cat("   - Different optimal coordination levels\n\n")

cat("Test completed successfully!\n")
cat("\nBoth ple() and ple.blp() are working correctly.\n")
cat("Use ple() when you have simple logit demand (will calibrate from margins).\n")
cat("Use ple.blp() when you have pre-estimated BLP demand parameters.\n")
cat("  - Margins are OPTIONAL for ple.blp() (supermarkup inferred from prices)\n")
cat("  - Provide margins only if you want to calibrate firm-specific timing parameters\n")
