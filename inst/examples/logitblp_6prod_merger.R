# Example: 6-product LogitBLP simulation with a merger between firms 1 and 2
#
# This example demonstrates the LogitBLP workflow:
# 1. Generate internally consistent pre-merger equilibrium with known MC
# 2. Use that equilibrium as "observed" data to calibrate delta
# 3. Simulate post-merger equilibrium
#
# How to run inside R:
#   # from package root after installation or devtools::load_all()
#   source(system.file("examples", "logitblp_6prod_merger.R", package = "antitrust"))
rm(list=ls())
# Parameters
set.seed(123)

k <- 6
labels <- paste0("Prod", 1:k)

# Pre-merger owners: each product by a distinct firm
ownerPre <- as.character(1:k)
# Post-merger owners: firms 1 and 2 merge
ownerPost <- as.character(c(1, 1, 3, 4, 5, 6))

# STEP 1: Set up observed equilibrium data (rates, shares, deposits)
cat("STEP 1: Using observed deposit rates and shares as equilibrium data...\n")

# Observed equilibrium deposit rates (prices in input market)
rates_equilibrium <- c(2.50, 2.80, 2.45, 2.35, 2.40, 2.75)

# Observed market shares (deposit shares)
shares_equilibrium <- c(0.15, 0.13, 0.17, 0.16, 0.12, 0.14)  # Sum = 0.87, outside = 0.13

# BLP demand parameters
# Input market (deposits): alpha > 0 (higher rates attract more deposits)
alpha_mean <- 0.4250503836
sigma <- 0.2878503814  # Random coefficient std dev (heterogeneity in price sensitivity)
piDemog <- 0.2090011524
nDemog <- 1
sigmaNest <- 0.91950485022  # Nesting parameter: sigmaNest = 1 - lambda = 1 - 0.08049514978
nDraws <- 500  # Reduced for faster post-merger equilibrium solving

cat("  Observed rates:", round(rates_equilibrium, 2), "\n")
cat("  Observed shares:", round(shares_equilibrium, 4), "(inside share:", round(sum(shares_equilibrium), 3), ")\n")
cat("  Calibrating BLP demand model from observed data...\n\n")

# STEP 2: Build LogitBLP model and simulate merger
cat("STEP 2: Building LogitBLP model via sim() and simulating merger...\n")

# Calibrate demand from observed rates and shares
# sim() will recover delta via BLP contraction and back out implied MRP from rates

demand.param <- list(
  alpha = alpha_mean,
  sigma = sigma,
  piDemog = piDemog,
  nDemog = nDemog,
  sigmaNest = sigmaNest,
  nDraws = nDraws,
  contractionTol = 1e-9,
  contractionMaxIter = 200
)

# Call sim() to construct the object using the package pipeline
# sim() will:
# 1. Run BLP contraction to recover delta from observed rates/shares
# 2. Back out implied MRP from observed rates using demand
# 3. Solve for post-merger equilibrium rates
sim_res <- sim(
  prices = rates_equilibrium,
  shares = shares_equilibrium,
  supply = "bertrand",
  demand = "LogitBLP",
  demand.param = demand.param,
  ownerPre = ownerPre,
  ownerPost = ownerPost,
  insideSize = 1e6,
  priceOutside = 0,
  labels = labels
)

cat("\nSTEP 3: Merger simulation complete; reporting results...\n")

sharesPre <- calcShares(sim_res, preMerger = TRUE, revenue = FALSE)
sharesPost <- calcShares(sim_res, preMerger = FALSE, revenue = FALSE)

# Print results
cat("\n--- LogitBLP 6-product deposit market merger (firms 1 & 2) ---\n")
cat("Pre-merger rates:\n"); print(round(sim_res@pricePre, 4))
cat("Post-merger rates:\n"); print(round(sim_res@pricePost, 4))
cat("Rate changes (%):\n"); print(round(100 * (sim_res@pricePost / sim_res@pricePre - 1), 2))
cat("\nImplied MRP (marginal revenue product):\n")
cat("Pre-merger:\n"); print(round(sim_res@mcPre, 4))
cat("Post-merger:\n"); print(round(sim_res@mcPost, 4))

cat("\nPre-merger shares (quantities):\n"); print(round(sharesPre, 5))
cat("Post-merger shares (quantities):\n"); print(round(sharesPost, 5))
cat("Inside share pre vs post:\n")
cat(round(sum(sharesPre), 5), "->", round(sum(sharesPost), 5), "\n\n")

