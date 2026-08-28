test_that("cmcr.bertrand calculates compensating marginal cost reductions", {
  prices <- c(1.0, 1.0)
  margins <- c(0.40, 0.30)
  ownerPre <- c(1, 2)
  ownerPost <- matrix(1, 2, 2)
  diversions <- matrix(c(-1, 0.4, 0.6, -1), nrow = 2)

  res <- cmcr.bertrand(prices = prices, margins = margins,
                       ownerPre = ownerPre, ownerPost = ownerPost,
                       diversions = diversions)

  expect_true(all(res > 0))
})

test_that("upp.bertrand calculates upward pricing pressure", {
  prices <- c(1.0, 1.0)
  margins <- c(0.40, 0.30)
  ownerPre <- c(1, 2)
  ownerPost <- matrix(1, 2, 2)
  diversions <- matrix(c(-1, 0.4, 0.6, -1), nrow = 2)

  res <- upp.bertrand(prices = prices, margins = margins,
                      ownerPre = ownerPre, ownerPost = ownerPost,
                      diversions = diversions)

  expect_true(all(res > 0))
})
