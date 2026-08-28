test_that("Logit Bertrand pre-merger margin matches exact analytical formula", {
  prices <- c(2.0, 2.0, 2.0)
  shares <- c(0.25, 0.25, 0.25) # Outside share s0 = 0.25
  margins <- c(0.40, 0.40, 0.40)
  ownerPre <- c("A", "B", "C")
  ownerPost <- c("A", "A", "C")

  fit <- logit(prices = prices, shares = shares, margins = margins,
               ownerPre = ownerPre, ownerPost = ownerPost)

  # Analytical alpha for Logit: alpha = -1 / (m * p * (1 - s))
  alpha_true <- -1 / (margins[1] * prices[1] * (1 - shares[1]))
  expect_equal(unname(fit@slopes$alpha), alpha_true, tolerance = 1e-6)

  # Pre-merger FOC residual must be identically zero (< 1e-6)
  pre_margins <- calcMargins(fit, preMerger = TRUE)
  expect_equal(unname(pre_margins), margins, tolerance = 1e-6)

  # Round-trip check: pre-merger MC = price * (1 - margin)
  mc_true <- prices * (1 - margins)
  expect_equal(unname(fit@mcPre), mc_true, tolerance = 1e-6)
})

test_that("CES Bargaining pre-merger margin scales exactly by (1 - barg) from base CES margin", {
  prices <- c(1.5, 1.5)
  shares <- c(0.66, 0.34)
  margins <- c(0.30, 0.20)
  barg <- c(0.5, 0.5)
  ownerPre <- c("A", "B")
  ownerPost <- c("A", "A")

  fit_barg <- bargaining.ces(prices = prices, shares = shares, margins = margins,
                             bargpowerPre = barg, ownerPre = ownerPre, ownerPost = ownerPost)

  pre_margins_barg <- calcMargins(fit_barg, preMerger = TRUE)
  pre_margins_ces  <- calcMargins(as(fit_barg, "CES"), preMerger = TRUE)

  # Exact mathematical invariant: m_barg = (1 - b) * m_ces for the same model gamma
  expect_equal(unname(pre_margins_barg), unname((1 - barg) * pre_margins_ces), tolerance = 1e-6)
})

test_that("Symmetric firm merger produces identical post-merger price predictions", {
  prices <- c(2.0, 2.0, 2.0)
  shares <- c(0.25, 0.25, 0.25)
  margins <- c(0.40, 0.40, 0.40)
  ownerPre <- c("A", "B", "C")
  ownerPost <- c("A", "A", "C") # A and B merge symmetrically

  fit <- logit(prices = prices, shares = shares, margins = margins,
               ownerPre = ownerPre, ownerPost = ownerPost)

  delta <- calcPriceDelta(fit)
  # Merging symmetric firms 1 and 2 must have identical price deltas
  expect_equal(unname(delta[1]), unname(delta[2]), tolerance = 1e-6)
})
