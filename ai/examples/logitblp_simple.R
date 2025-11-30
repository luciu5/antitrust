## Minimal LogitBLP Example
## Tests basic functionality of BLP contraction and merger simulation

library(antitrust)
set.seed(42)

# Simple 3-product market
prices <- c(10, 15, 12)
shares <- c(0.3, 0.25, 0.25)  # Outside share = 0.2
margins <- c(0.3, 0.4, 0.35)

# Ownership: Firms 1 and 2 merge
ownerPre <- c("A", "B", "C")
ownerPost <- c("A", "A", "C")

# BLP parameters
demand.param <- list(
  alpha = -0.8,       # Price sensitivity
  sigma = 0.2,        # Random coefficient std dev
  sigmaNest = 0.9,    # Nesting (0.9 = weak nesting)
  nDraws = 200        # Monte Carlo draws
)

cat("Running BLP merger simulation...\n")

# Run simulation
result <- sim(
  prices = prices,
  shares = shares,
  margins = margins,
  demand = "LogitBLP",
  demand.param = demand.param,
  ownerPre = ownerPre,
  ownerPost = ownerPost,
  insideSize = 1000,
  labels = c("Prod_A", "Prod_B", "Prod_C")
)

# Show results
cat("\nPre-merger prices: ", round(result@pricePre, 2), "\n")
cat("Post-merger prices:", round(result@pricePost, 2), "\n")
cat("Price increase (%):", round(100 * (result@pricePost - result@pricePre) / result@pricePre, 1), "\n")

cat("\nRecovered mean utilities (delta):\n")
print(round(result@slopes$meanval, 3))

cat("\nSimulation complete!\n")
