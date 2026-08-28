test_that("logit.alm recovers unobserved outside share s0 from multiple margins", {
  prices <- c(2.0, 2.5, 3.0)
  shares_inside <- c(0.40, 0.35, 0.25) # Inside shares sum to 1.0
  margins <- c(0.45, 0.40, 0.35)       # 3 observed margins
  ownerPre <- c("A", "B", "C")
  ownerPost <- c("A", "A", "C")

  fit_alm <- logit.alm(prices = prices, shares = shares_inside, margins = margins,
                       ownerPre = ownerPre, ownerPost = ownerPost)

  # ALM must estimate an outside share shareInside < 1
  expect_true(fit_alm@shareInside < 1.0)
  expect_true(fit_alm@shareInside > 0.0)

  # Pre-merger price diagnostic error must be near zero
  diag <- calcDiagnostics(fit_alm)
  expect_true(all(abs(diag$Prices) < 1.0))
})

test_that("ces.alm recovers unobserved outside share from multiple margins", {
  prices <- c(2.0, 2.5, 3.0)
  shares_inside <- c(0.40, 0.35, 0.25)
  margins <- c(0.45, 0.40, 0.35)
  ownerPre <- c("A", "B", "C")
  ownerPost <- c("A", "A", "C")

  fit_ces_alm <- ces.alm(prices = prices, shares = shares_inside, margins = margins,
                         ownerPre = ownerPre, ownerPost = ownerPost)

  expect_true(fit_ces_alm@shareInside <= 1.0)
  expect_true(fit_ces_alm@shareInside > 0.0)

  diag <- calcDiagnostics(fit_ces_alm)
  expect_true(all(abs(diag$Prices) < 1.0))
})

test_that("Partial margin calibration (N-1 missing margins) recovers implied margins for all products", {
  prices <- c(2.0, 2.5, 3.0, 2.8, 2.2)
  shares <- c(0.20, 0.15, 0.15, 0.10, 0.10) # Outside share s0 = 0.30
  margins <- c(0.40, NA, NA, NA, NA)         # Only 1 margin observed
  ownerPre <- c("A", "B", "C", "D", "E")
  ownerPost <- c("A", "A", "C", "D", "E")

  fit_partial <- logit(prices = prices, shares = shares, margins = margins,
                       ownerPre = ownerPre, ownerPost = ownerPost)

  # Parameter alpha must be calibrated from single observed margin
  expect_false(is.na(fit_partial@slopes$alpha))

  # Pre-merger margins for ALL 5 products must be non-NA and valid (between 0 and 1)
  pre_margins <- calcMargins(fit_partial, preMerger = TRUE)
  expect_equal(length(pre_margins), 5)
  expect_true(all(!is.na(pre_margins)))
  expect_true(all(pre_margins > 0 & pre_margins < 1))

  # First product margin must match observed margin 0.40 exactly
  expect_equal(unname(pre_margins[1]), 0.40, tolerance = 1e-5)
})

test_that("bargaining.ces.alm estimates unobserved outside share under Nash Bargaining", {
  prices <- c(1.5, 2.0, 2.5)
  shares_inside <- c(0.50, 0.30, 0.20)
  margins <- c(0.30, 0.25, 0.20)
  ownerPre <- c("A", "B", "C")
  ownerPost <- c("A", "A", "C")

  fit_barg_alm <- bargaining.ces.alm(prices = prices, shares = shares_inside, margins = margins,
                                     ownerPre = ownerPre, ownerPost = ownerPost)

  expect_false(is.na(fit_barg_alm@slopes$gamma))
  expect_true(fit_barg_alm@shareInside <= 1.0)
})
