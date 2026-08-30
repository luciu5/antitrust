test_that("Output market elasticities have negative own-price and positive cross-price signs", {
  prices <- c(1.0, 1.0)
  shares <- c(0.40, 0.30)
  margins <- c(0.40, 0.30)
  ownerPre <- c("A", "B")
  ownerPost <- c("A", "A")

  # Logit Output Market (output = TRUE)
  fit_logit <- logit(prices = prices, shares = shares, margins = margins,
                     ownerPre = ownerPre, ownerPost = ownerPost, output = TRUE)
  E_logit <- elast(fit_logit)

  # Own-price elasticities (diagonal) must be strictly negative
  expect_true(all(diag(E_logit) < 0))

  # Cross-price elasticities (off-diagonal) must be strictly positive
  expect_true(E_logit[1, 2] > 0)
  expect_true(E_logit[2, 1] > 0)

  # CES Output Market (output = TRUE)
  fit_ces <- ces(prices = prices, shares = c(0.45, 0.30),
                 margins = c(0.6451613, 0.5882353),
                 ownerPre = ownerPre, ownerPost = ownerPost, output = TRUE)
  E_ces <- elast(fit_ces)

  expect_true(all(diag(E_ces) < 0))
  expect_true(E_ces[1, 2] > 0)
  expect_true(E_ces[2, 1] > 0)
})

test_that("Input market elasticities have positive own-price and negative cross-price signs", {
  prices <- c(1.0, 1.0)
  shares <- c(0.40, 0.30)
  margins <- c(0.40, 0.30)
  ownerPre <- c("A", "B")
  ownerPost <- c("A", "A")

  # Logit Input Market (output = FALSE)
  fit_logit_in <- logit(prices = prices, shares = shares, margins = margins,
                        ownerPre = ownerPre, ownerPost = ownerPost, output = FALSE)
  E_logit_in <- elast(fit_logit_in)

  # Own-price elasticities (diagonal) must be strictly positive in input markets
  expect_true(all(diag(E_logit_in) > 0))

  # Cross-price elasticities (off-diagonal) must be strictly negative in input markets
  expect_true(E_logit_in[1, 2] < 0)
  expect_true(E_logit_in[2, 1] < 0)

  # Bargaining CES Input Market (output = FALSE)
  fit_ces_in <- bargaining.ces(prices = prices, shares = c(0.66, 0.34),
                              margins = c(0.3623188, 0.1381224),
                              ownerPre = ownerPre, ownerPost = ownerPost, output = FALSE)
  E_ces_in <- elast(fit_ces_in)

  expect_true(all(diag(E_ces_in) > 0))
  expect_true(E_ces_in[1, 2] < 0)
  expect_true(E_ces_in[2, 1] < 0)
})
