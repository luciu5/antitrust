test_that("Pre- and post-merger equilibria match between nleqslv and AG solvers on Logit package examples", {
  prices <- c(2.0, 2.2, 2.5)
  shares <- c(0.35, 0.25, 0.20)
  margins <- c(0.40, 0.35, 0.30)
  ownerPre <- c("A", "B", "C")
  ownerPost <- c("A", "A", "C")

  fit_nleqslv <- logit(prices = prices, shares = shares, margins = margins,
                       ownerPre = ownerPre, ownerPost = ownerPost, solver = "nleqslv")

  fit_ag      <- logit(prices = prices, shares = shares, margins = margins,
                       ownerPre = ownerPre, ownerPost = ownerPost, solver = "ag")

  # Pre-merger price parity
  expect_equal(unname(fit_ag@pricePre), unname(fit_nleqslv@pricePre), tolerance = 1e-4)

  # Post-merger price parity
  expect_equal(unname(fit_ag@pricePost), unname(fit_nleqslv@pricePost), tolerance = 1e-4)

  # Price delta (% change) parity
  expect_equal(unname(calcPriceDelta(fit_ag)), unname(calcPriceDelta(fit_nleqslv)), tolerance = 1e-4)
})

test_that("Pre- and post-merger equilibria match between nleqslv and AG solvers on CES package examples", {
  prices <- c(1.5, 1.8, 2.0)
  shares <- c(0.50, 0.30, 0.20)
  margins <- c(0.30, 0.25, 0.20)
  ownerPre <- c("A", "B", "C")
  ownerPost <- c("A", "A", "C")

  fit_nleqslv <- ces(prices = prices, shares = shares, margins = margins,
                     ownerPre = ownerPre, ownerPost = ownerPost, solver = "nleqslv")

  fit_ag      <- ces(prices = prices, shares = shares, margins = margins,
                     ownerPre = ownerPre, ownerPost = ownerPost, solver = "ag")

  expect_equal(unname(fit_ag@pricePre), unname(fit_nleqslv@pricePre), tolerance = 1e-4)
  expect_equal(unname(fit_ag@pricePost), unname(fit_nleqslv@pricePost), tolerance = 1e-4)
  expect_equal(unname(calcPriceDelta(fit_ag)), unname(calcPriceDelta(fit_nleqslv)), tolerance = 1e-4)
})

test_that("Pre- and post-merger equilibria match on simulated 5-product asymmetric markets with efficiencies", {
  nprods <- 5
  prices <- c(1.0, 1.2, 1.5, 1.8, 2.0)
  shares <- c(0.25, 0.20, 0.15, 0.10, 0.10) # Outside share s0 = 0.20
  margins <- c(0.40, 0.35, 0.30, 0.25, 0.20)
  ownerPre <- c("Firm1", "Firm2", "Firm3", "Firm4", "Firm5")
  ownerPost <- c("Firm1", "Firm1", "Firm3", "Firm4", "Firm5") # Firm1 and Firm2 merge
  mcDelta <- c(-0.05, 0, 0, 0, 0) # 5% efficiency on Product 1

  fit_nleqslv <- logit(prices = prices, shares = shares, margins = margins,
                       ownerPre = ownerPre, ownerPost = ownerPost, mcDelta = mcDelta,
                       solver = "nleqslv")

  fit_ag      <- logit(prices = prices, shares = shares, margins = margins,
                       ownerPre = ownerPre, ownerPost = ownerPost, mcDelta = mcDelta,
                       solver = "ag")

  expect_equal(unname(fit_ag@pricePre), unname(fit_nleqslv@pricePre), tolerance = 1e-4)
  expect_equal(unname(fit_ag@pricePost), unname(fit_nleqslv@pricePost), tolerance = 1e-4)
  expect_equal(unname(calcPriceDelta(fit_ag)), unname(calcPriceDelta(fit_nleqslv)), tolerance = 1e-4)
})

test_that("Pre- and post-merger equilibria match on Bargaining Logit models", {
  prices <- c(2.0, 2.2, 2.5)
  shares <- c(0.35, 0.25, 0.20)
  margins <- c(0.40, 0.35, 0.30)
  ownerPre <- c("A", "B", "C")
  ownerPost <- c("A", "A", "C")

  fit_nleqslv <- bargaining.logit(prices = prices, shares = shares, margins = margins,
                                  ownerPre = ownerPre, ownerPost = ownerPost, solver = "nleqslv")

  fit_ag      <- bargaining.logit(prices = prices, shares = shares, margins = margins,
                                  ownerPre = ownerPre, ownerPost = ownerPost, solver = "ag")

  expect_equal(unname(fit_ag@pricePre), unname(fit_nleqslv@pricePre), tolerance = 1e-4)
  expect_equal(unname(fit_ag@pricePost), unname(fit_nleqslv@pricePost), tolerance = 1e-4)
  expect_equal(unname(calcPriceDelta(fit_ag)), unname(calcPriceDelta(fit_nleqslv)), tolerance = 1e-4)
})
