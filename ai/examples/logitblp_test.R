## Simple Monte Carlo Example for LogitBLP Merger Simulation
## This example simulates a merger between two firms in a market with 4 products

library(antitrust)

# Set seed for reproducibility
set.seed(123)

# Market structure: 4 products from 3 firms
# Firm A: Products 1 and 2
# Firm B: Product 3
# Firm C: Product 4
# Merger: Firm A acquires Firm B

## 1. Define market parameters
nprods <- 4
prodNames <- c("Product_A1", "Product_A2", "Product_B", "Product_C")

# Pre-merger ownership
ownerPre <- c("Firm_A", "Firm_A", "Firm_B", "Firm_C")

# Post-merger ownership (Firm A acquires Firm B)
ownerPost <- c("Firm_A", "Firm_A", "Firm_A", "Firm_C")

## 2. Observed market data
# Prices
prices <- c(10, 12, 11, 9)

# Market shares (inside goods only, sum to ~0.8, outside share ~0.2)
shares <- c(0.25, 0.20, 0.18, 0.17)

# Margins (calibrated or observed)
margins <- c(0.30, 0.35, 0.28, 0.25)

# Total market size (units sold for inside products)
insideSize <- 10000

## 3. BLP demand parameters (typically estimated from data)
# These would normally come from BLP estimation
# For this example, we use assumed/calibrated values

# Mean price coefficient (negative for output market)
alphaMean <- -0.5

# Standard deviation of random coefficient on price
# Captures consumer heterogeneity in price sensitivity
sigma <- 0.3

# Nesting parameter (outside good in separate nest)
# sigmaNest = 1 is flat logit, sigmaNest < 1 increases correlation among inside products
sigmaNest <- 0.85

# Number of Monte Carlo draws for simulating consumer heterogeneity
nDraws <- 500

# Package parameters
demand.param <- list(
  alpha = alphaMean,        # Mean price coefficient
  sigma = sigma,            # Random coefficient std dev
  sigmaNest = sigmaNest,    # Nesting parameter
  nDraws = nDraws          # Number of simulation draws
)

## 4. Run merger simulation
cat("\n=== Running LogitBLP Merger Simulation ===\n\n")

result <- sim(
  prices = prices,
  shares = shares,
  margins = margins,
  demand = "LogitBLP",
  demand.param = demand.param,
  supply = "bertrand",
  ownerPre = ownerPre,
  ownerPost = ownerPost,
  insideSize = insideSize,
  labels = prodNames
)

## 5. Display results
cat("\n=== Simulation Results ===\n")
print(result)

cat("\n=== Detailed Summary ===\n")
summary(result)

## 6. Calculate merger effects
cat("\n=== Pre-Merger Prices ===\n")
print(result@pricePre)

cat("\n=== Post-Merger Prices ===\n")
print(result@pricePost)

cat("\n=== Price Changes (%) ===\n")
priceChange <- 100 * (result@pricePost - result@pricePre) / result@pricePre
names(priceChange) <- prodNames
print(round(priceChange, 2))

cat("\n=== Pre-Merger Shares ===\n")
sharesPre <- calcShares(result, preMerger = TRUE)
print(round(sharesPre, 4))

cat("\n=== Post-Merger Shares ===\n")
sharesPost <- calcShares(result, preMerger = FALSE)
print(round(sharesPost, 4))

cat("\n=== Pre-Merger Elasticities ===\n")
elastPre <- elast(result, preMerger = TRUE)
print(round(elastPre, 3))

cat("\n=== Post-Merger Elasticities ===\n")
elastPost <- elast(result, preMerger = FALSE)
print(round(elastPost, 3))

cat("\n=== Diversion Ratios (Pre-Merger) ===\n")
divPre <- diversion(result, preMerger = TRUE)
print(round(divPre, 3))

cat("\n=== Consumer Welfare Change ===\n")
cv <- CV(result)
cat("Compensating Variation: $", round(sum(cv), 2), "\n", sep="")
cat("Per-consumer CV: $", round(mean(cv), 2), "\n", sep="")

## 7. Test alternative specification with meanval provided
cat("\n\n=== Testing with Pre-Computed Meanval ===\n")

# Extract recovered delta and draws from first simulation
delta <- result@slopes$meanval
consDraws <- result@slopes$consDraws
demogDraws <- result@slopes$demogDraws

# Re-run with meanval and draws provided (should skip BLP contraction and use same draws)
demand.param2 <- list(
  alpha = alphaMean,
  meanval = delta,          # Provide pre-computed delta
  sigma = sigma,
  sigmaNest = sigmaNest,
  nDraws = nDraws,
  consDraws = consDraws,    # Reuse same random draws
  demogDraws = demogDraws   # Reuse same demographic draws
)

result2 <- sim(
  prices = prices,
  shares = shares,
  margins = margins,
  demand = "LogitBLP",
  demand.param = demand.param2,
  supply = "bertrand",
  ownerPre = ownerPre,
  ownerPost = ownerPost,
  insideSize = insideSize,
  labels = prodNames
)

cat("\nPrice changes should be identical:\n")
cat("Specification 1 (shares provided): ", round(priceChange[1], 4), "\n")
priceChange2 <- 100 * (result2@pricePost - result2@pricePre) / result2@pricePre
cat("Specification 2 (meanval provided): ", round(priceChange2[1], 4), "\n")

cat("\n=== Test Complete ===\n")
