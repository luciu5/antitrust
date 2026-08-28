test_that("solver = 'ag' argument in constructor functions uses Aggregative Games solver", {
  prices <- c(1.0, 1.0)
  shares <- c(0.40, 0.30)
  margins <- c(0.40, 0.30)
  ownerPre <- c("A", "B")
  ownerPost <- c("A", "A")

  fit_std <- logit(prices = prices, shares = shares, margins = margins,
                   ownerPre = ownerPre, ownerPost = ownerPost, solver = "nleqslv")

  fit_ag  <- logit(prices = prices, shares = shares, margins = margins,
                   ownerPre = ownerPre, ownerPost = ownerPost, solver = "ag")

  expect_equal(unname(calcPriceDelta(fit_ag)), unname(calcPriceDelta(fit_std)), tolerance = 1e-4)
})

test_that("solver = 'ag' in bargaining.ces constructor works cleanly", {
  prices <- c(1.0, 1.0)
  shares <- c(0.66, 0.34)
  margins <- c(0.30, 0.20)
  ownerPre <- c("A", "B")
  ownerPost <- c("A", "A")

  fit_barg_ag <- bargaining.ces(prices = prices, shares = shares, margins = margins,
                                ownerPre = ownerPre, ownerPost = ownerPost, solver = "ag")

  expect_true(all(calcPriceDelta(fit_barg_ag) > 0))
})
