## Example: Price Leadership with Merger Coalition Expansion
## Demonstrates how to explicitly specify coalition expansion when a merger brings
## coalition and non-coalition products under common ownership

library(antitrust)

## 4-product market setup
## Pre-merger:
##   - Products 1, 2: Coalition (Firm A owns 1, Firm B owns 2)
##   - Product 3: Fringe (Firm C)
##   - Product 4: Fringe (Firm D)
##
## Merger: Firm B acquires Firm C (products 2 and 3 merge)
##
## Expected Post-merger coalition:
##   - Products 1, 2, 3 (product 3 joins because it's now commonly owned with coalition product 2)
##   - Product 4: Fringe
##   - Note: Coalition expansion must now be specified explicitly

prices <- c(10, 11, 9, 8.5)
shares <- c(0.25, 0.22, 0.20, 0.13)
margins_pct <- c(0.38, 0.36, 0.22, NA)  # Coalition margins higher

# Pre-merger ownership
ownerPre <- c("Firm_A", "Firm_B", "Firm_C", "Firm_D")

# Post-merger: Firm B acquires Firm C
ownerPost <- c("Firm_A", "Firm_B", "Firm_B", "Firm_D")  # Products 2 & 3 now under Firm_B

# Pre-merger coalition: only products 1 and 2
coalitionPre <- c(1, 2)

cat("\n=== Coalition Expansion Example ===\n\n")

cat("Pre-merger Structure:\n")
cat("  Coalition: Products 1 (Firm A), 2 (Firm B)\n")
cat("  Fringe: Products 3 (Firm C), 4 (Firm D)\n\n")

cat("Merger: Firm B acquires Firm C\n")
cat("  Products 2 & 3 now commonly owned\n\n")

cat("Creating PriceLeadership model...\n")

result <- ple(
  prices = prices,
  shares = shares,
  margins = margins_pct,
  ownerPre = ownerPre,
  ownerPost = ownerPost,
  coalitionPre = coalitionPre,
  coalitionPost = c(1, 2, 3),  # Explicitly specify expanded coalition post-merger
  insideSize = 10000,
  labels = c("Prod_A", "Prod_B", "Prod_C", "Prod_D")
)

cat("\n=== Coalition Membership ===\n\n")
cat("Pre-merger coalition: ", result@coalitionPre, "\n")
cat("Post-merger coalition: ", result@coalitionPost, "\n")

if(length(result@coalitionPost) > length(result@coalitionPre)){
  cat("\n✓ Coalition explicitly expanded post-merger!\n")
  newMembers <- setdiff(result@coalitionPost, result@coalitionPre)
  cat("  New coalition members: Products", newMembers, "\n")
}

cat("\n=== Price Effects ===\n\n")
cat("Pre-merger prices:\n")
print(round(result@pricePre, 3))

cat("\nPost-merger prices:\n")
print(round(result@pricePost, 3))

cat("\nPrice changes (%):\n")
priceChange <- 100 * (result@pricePost - result@pricePre) / result@pricePre
print(round(priceChange, 2))

cat("\n=== Coalition Product Price Change ===\n")
cat("Product 3 (new coalition member) price change: ", 
    round(priceChange[3], 2), "%\n")
cat("This reflects coordination with products 1 & 2 post-merger.\n")
