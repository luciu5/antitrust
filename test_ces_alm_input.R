## Verify ces.alm correctly calibrates model parameters when output=FALSE
devtools::load_all()

cat("\n==========================================\n")
cat("  CES ALM INPUT MARKET (output=FALSE)     \n")
cat("==========================================\n\n")

## --- Test 1: Cross-check ALM vs CES with known outside share ---
cat("--- Test 1: ALM vs CES cross-check ---\n")
# First, run CES with shares < 1 (known outside share)
prices1  <- c(8, 12)
shares1  <- c(0.25, 0.35)  # sum = 0.6, outside share = 0.4
margins1 <- c(0.4, 0.5)
ownerPre1  <- diag(2)
ownerPost1 <- matrix(1, 2, 2)

res_ces <- ces(prices1, shares1, margins1,
               ownerPre  = ownerPre1,
               ownerPost = ownerPost1,
               output    = FALSE)

cat("  CES gamma:       ", res_ces@slopes$gamma, "\n")
cat("  CES shareInside: ", res_ces@shareInside, "\n")
cat("  CES mktElast:    ",
    elast(res_ces, preMerger = TRUE, market = TRUE), "\n")

# Now run ALM with shares normalized to 1 + mktElast from CES
# But mktElast for input CES is positive, and validity rejects > -1.
# So we test ALM using only 2 margins (no mktElast).
shares1_alm <- shares1 / sum(shares1)

res_alm <- ces.alm(prices1, shares1_alm, margins1,
                   ownerPre  = ownerPre1,
                   ownerPost = ownerPost1,
                   output    = FALSE)

cat("  ALM gamma:       ", res_alm@slopes$gamma, "\n")
cat("  ALM shareInside: ", res_alm@shareInside, "\n")
cat("  ALM sOut:        ", 1 - res_alm@shareInside, "\n")

# gamma should be < 1 (can be negative for input markets)
stopifnot(res_alm@slopes$gamma < 1)
cat("  gamma < 1: OK\n")

# Implied margins should match supplied margins
implied1 <- calcMargins(res_alm, preMerger = TRUE, level = FALSE)
cat("  Supplied margins: ", margins1, "\n")
cat("  Implied  margins: ", implied1, "\n")
cat("  Max abs diff:     ",
    max(abs(margins1 - implied1)), "\n")
stopifnot(all.equal(margins1, as.numeric(implied1),
                    tolerance = 1e-3))
cat("  Margins match: OK\n")

# MC should be positive (input: MC = (1+margin)*price > price)
mc1 <- res_alm@mcPre
cat("  MC (pre):         ", mc1, "\n")
stopifnot(all(mc1 > 0))
cat("  MC > 0: OK\n")

# Pre-merger prices should reproduce observed prices
cat("  pricePre:         ", res_alm@pricePre, "\n")
stopifnot(all.equal(prices1, as.numeric(res_alm@pricePre),
                    tolerance = 1e-3))
cat("  Prices match: OK\n")

# Compare ALM recovered parameters to CES
cat("  Gamma diff:       ",
    abs(res_ces@slopes$gamma - res_alm@slopes$gamma), "\n")
cat("  shareIn diff:     ",
    abs(res_ces@shareInside - res_alm@shareInside), "\n")
stopifnot(all.equal(unname(res_ces@slopes$gamma),
                    unname(res_alm@slopes$gamma),
                    tolerance = 0.1))
cat("  Gamma approx match: OK\n")
stopifnot(all.equal(unname(res_ces@shareInside),
                    unname(res_alm@shareInside),
                    tolerance = 0.1))
cat("  shareInside approx match: OK\n")

cat("\n--- Test 1 PASSED ---\n")


## --- Test 2: Asymmetric 3-product (all margins provided) ---
cat("\n--- Test 2: 3-product asymmetric ---\n")
prices2  <- c(8, 12, 10)
shares2  <- c(0.30, 0.45, 0.25)
margins2 <- c(0.3, 0.5, 0.4)  # need 3 margins for 2 unknowns
ownerPre2  <- diag(3)
ownerPost2 <- matrix(c(1, 1, 0, 1, 1, 0, 0, 0, 1), 3, 3)

res2 <- ces.alm(prices2, shares2, margins2,
                ownerPre  = ownerPre2,
                ownerPost = ownerPost2,
                output    = FALSE)

cat("  gamma:", res2@slopes$gamma, "\n")
cat("  sOut: ", 1 - res2@shareInside, "\n")
stopifnot(res2@slopes$gamma < 1)
cat("  gamma < 1: OK\n")

implied2 <- calcMargins(res2, preMerger = TRUE, level = FALSE)
cat("  Supplied margins:", margins2, "\n")
cat("  Implied  margins:", implied2, "\n")
# Over-identified (3 margins, 2 params): approximate match
maxDiff2 <- max(abs(margins2 - implied2))
cat("  Max abs diff:    ", maxDiff2, "\n")
stopifnot(maxDiff2 < 0.15)
cat("  Margins approx match: OK\n")

stopifnot(all(res2@mcPre > 0))
cat("  MC > 0: OK\n")

cat("\n--- Test 2 PASSED ---\n")


## --- Test 3: Merger direction (monopsony) ---
cat("\n--- Test 3: Merger direction ---\n")
cat("  Pre-merger prices: ", res_alm@pricePre, "\n")
cat("  Post-merger prices:", res_alm@pricePost, "\n")
# Input market merger: prices to suppliers decrease
stopifnot(all(res_alm@pricePost < res_alm@pricePre))
cat("  Post < Pre: OK (monopsony effect)\n")

cat("\n--- Test 3 PASSED ---\n")


## --- Test 4: mktElast now accepted for input markets ---
cat("\n--- Test 4: mktElast with input market ---\n")
mktE <- elast(res_ces, preMerger = TRUE, market = TRUE)
cat("  CES input market elasticity:", mktE, "\n")
cat("  (Positive; should now be accepted)\n")

res4 <- ces.alm(prices1, shares1_alm, margins1,
                ownerPre  = ownerPre1,
                ownerPost = ownerPost1,
                mktElast  = mktE,
                output    = FALSE)

cat("  gamma:       ", res4@slopes$gamma, "\n")
cat("  shareInside: ", res4@shareInside, "\n")
stopifnot(res4@slopes$gamma < 1)
cat("  gamma < 1: OK\n")

implied4 <- calcMargins(res4, preMerger = TRUE, level = FALSE)
cat("  Supplied margins:", margins1, "\n")
cat("  Implied  margins:", implied4, "\n")
stopifnot(all.equal(margins1, as.numeric(implied4),
                    tolerance = 1e-2))
cat("  Margins match: OK\n")

# mktElast constraint should be satisfied
impliedMktE <- elast(res4, preMerger = TRUE, market = TRUE)
cat("  Supplied mktElast:", mktE, "\n")
cat("  Implied  mktElast:", impliedMktE, "\n")
stopifnot(all.equal(mktE, as.numeric(impliedMktE),
                    tolerance = 0.05))
cat("  mktElast match: OK\n")

# Cross-check: should closely match CES with known outside share
cat("  Gamma diff vs CES:",
    abs(res_ces@slopes$gamma - res4@slopes$gamma), "\n")
cat("  shareIn diff:     ",
    abs(res_ces@shareInside - res4@shareInside), "\n")
stopifnot(all.equal(unname(res_ces@slopes$gamma),
                    unname(res4@slopes$gamma),
                    tolerance = 0.05))
cat("  Gamma match: OK\n")
stopifnot(all.equal(unname(res_ces@shareInside),
                    unname(res4@shareInside),
                    tolerance = 0.05))
cat("  shareInside match: OK\n")

cat("\n--- Test 4 PASSED ---\n")


## --- Test 5: Output market still works ---
cat("\n--- Test 5: Output market (regression test) ---\n")
prices5  <- c(8, 12)
shares5  <- c(0.4, 0.6)
margins5 <- c(0.3, 0.5)
mktElast5 <- -2

res5 <- ces.alm(prices5, shares5, margins5,
                ownerPre  = diag(2),
                ownerPost = matrix(1, 2, 2),
                mktElast  = mktElast5,
                output    = TRUE)

stopifnot(res5@slopes$gamma > 1)
cat("  gamma:", res5@slopes$gamma, " (> 1): OK\n")

implied5 <- calcMargins(res5, preMerger = TRUE, level = FALSE)
cat("  Supplied margins:", margins5, "\n")
cat("  Implied  margins:", implied5, "\n")
# Over-identified: approximate match
stopifnot(max(abs(margins5 - implied5)) < 0.15)
cat("  Margins approx match: OK\n")

stopifnot(all(res5@mcPre > 0))
cat("  MC > 0: OK\n")

cat("\n--- Test 5 PASSED ---\n")


cat("\n==========================================\n")
cat("  ALL TESTS COMPLETED                     \n")
cat("==========================================\n")
