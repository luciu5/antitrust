test_that("LogitNests assigns zero shares to an entirely inactive nest", {
  fit <- sim(
    prices = c(2, 2.2, 2.5),
    supply = "bertrand", demand = "LogitNests",
    demand.param = list(
      alpha = -1, meanval = c(.4, .2, .1), sigma = .7
    ),
    nests = c("N1", "N1", "N2"),
    ownerPre = c("A", "B", "C"),
    ownerPost = c("A", "A", "C"),
    subset = c(FALSE, FALSE, TRUE)
  )

  post_quantity_shares <- calcShares(fit, preMerger = FALSE)
  post_revenue_shares <- calcShares(fit, preMerger = FALSE, revenue = TRUE)

  expect_equal(unname(post_quantity_shares[1:2]), c(0, 0))
  expect_equal(unname(post_revenue_shares[1:2]), c(0, 0))
  expect_true(is.finite(post_quantity_shares[3]))
  expect_true(is.finite(post_revenue_shares[3]))
})
