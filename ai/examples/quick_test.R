# Quick Performance Test
# Copy and paste into R console

# Load package
devtools::load_all()

# Test data
prices <- c(0.93, 0.88, 1.10, 1.02)
shares <- c(0.35, 0.25, 0.25, 0.15)
margins <- c(0.40, NA, 0.35, 0.25)  # Coalition margins and one fringe margin

# Time the execution
system.time({
  result <- ple(
    prices = prices,
    shares = shares,
    margins = margins,
    ownerPre = c("AB", "AB", "MC", "Fringe"),
    ownerPost = c("AB", "AB", "MC", "Fringe"),
    coalitionPre = c(1, 2, 3),
    coalitionPost = c(1, 2, 3),  # Same coalition post-merger
    insideSize = 1000
  )
})

# Check results
print(result@supermarkupPre)
print(result@timingParam)
print(result@pricePre)
print(result@pricePost)
