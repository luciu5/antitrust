## Simple Test Example for Price Leadership Model
## Tests basic functionality of the PriceLeadership class

library(antitrust)

# Set seed for reproducibility
set.seed(123)

## Simple 4-product market
## Products 1-3: Coalition (coordinating)
## Product 4: Fringe (competitive)

## Market data
prices <- c(10, 12, 11, 9)
shares <- c(0.25, 0.20, 0.18, 0.17)  # Sum = 0.80, outside share = 0.20

# Margins: Coalition margins higher than fringe
# Provide 1 coalition and 1 fringe margin for calibration
margins_pct <- c(0.40, NA, 0.35, 0.25)  # Coalition at 40%, 35%; Fringe at 25%

# Ownership structure
ownerPre <- c("Firm_A", "Firm_A", "Firm_B", "Firm_C")  
ownerPost <- ownerPre  # No merger in this example

# Pre-merger coalition: Products 1-3
coalitionPre <- c(1, 2, 3)

cat("\n=== Price Leadership Model Test ===\n\n")

cat("Market Structure:\n")
cat("  Coalition Products: 1, 2, 3 (Firms A, A, B)\n")
cat("  Fringe Product: 4 (Firm C)\n\n")

cat("Observed Data:\n")
cat("  Prices:  ", prices, "\n")
cat("  Shares:  ", shares, "\n")
cat("  Margins: ", margins_pct, "(proportion)\n\n")

## Test the constructor
cat("Creating PriceLeadership object...\n")

result <- ple(
  prices = prices,
  shares = shares,
  margins = margins_pct,
  ownerPre = ownerPre,
  ownerPost = ownerPost,
  coalitionPre = coalitionPre,
  coalitionPost = coalitionPre,  # Same coalition post-merger (no merger)
  insideSize = 10000,
  labels = c("Prod_A1", "Prod_A2", "Prod_B", "Prod_C")
)

cat("\n=== Calibration Results ===\n\n")

cat("Price coefficient (alpha): ", result@slopes$alpha, "\n")
cat("Supermarkup (m): ", result@supermarkup, "\n")
cat("Timing parameter (delta): ", 
    ifelse(is.na(result@timingParam), "Not identified (unconstrained)", 
           as.character(round(result@timingParam, 4))), "\n")
cat("Binding firm: ", result@bindingFirm, "\n\n")

cat("Mean valuations:\n")
print(result@slopes$meanval)

cat("\n=== Incentive Compatibility Analysis ===\n")

if(!is.na(result@timingParam)){
  cat("\nSlack function values g_f(m):\n")
  print(round(result@slackValues, 4))
  cat("\nInterpretation: g_f(m) >= 0 means firm's IC constraint is satisfied.\n")
  cat("Binding firm has g_f(m) ≈ 0.\n")
} else {
  cat("\nImmediate payoff differences [π^PL - π^D]:\n")
  print(round(result@slackValues, 4))
  cat("\nAll values positive: Firms prefer coordination even without future considerations.\n")
}

cat("\n=== Pre-Merger Equilibrium ===\n")
cat("Prices:\n")
print(result@pricePre)

cat("\nMargins (calculated):\n")
marginsCalc <- calcMargins(result, preMerger = TRUE, level = FALSE)
print(marginsCalc)

cat("\nShares (calculated):\n")
sharesCalc <- calcShares(result, preMerger = TRUE)
print(sharesCalc)

cat("\n=== Validation ===\n")
cat("Share match: ", 
    isTRUE(all.equal(shares, sharesCalc, tolerance = 1e-3)), "\n")
cat("Margin match (fringe): ", 
    isTRUE(all.equal(margins_pct[4], marginsCalc[4], tolerance = 1e-3)), "\n")

cat("\n=== Additional Analysis ===\n")

# Test calcSlack function directly
cat("\nRecalculating slack functions:\n")
slackRecomputed <- calcSlack(result, preMerger = TRUE)
print(round(slackRecomputed, 4))

# Calculate Bertrand baseline
cat("\nBertrand equilibrium prices (punishment):\n")
bertrandPrices <- calcBertrandPrices(result, preMerger = TRUE)
print(round(bertrandPrices, 3))

cat("\nPrice leadership premium over Bertrand:\n")
premiums <- result@pricePre - bertrandPrices
names(premiums) <- c("Prod_A1", "Prod_A2", "Prod_B", "Prod_C")
print(round(premiums, 3))

cat("\n=== Test Complete ===\n")
