input_blp_fit <- function(sigma_nest = 1) {
  suppressWarnings(specify(
    "blp", "bertrand", prices = c(.2, .1, .3, .25),
    parameters = list(alphaMean = .8, sigma = .3, sigmaNest = sigma_nest,
                      meanval = c(0, .1, .2, -.1)),
    ownerPre = c("A", "B", "C", "D"), insideSize = 100, output = FALSE))
}

numeric_share_jacobian <- function(model, pre) {
  active <- if (pre) rep(TRUE, length(model@labels)) else model@subset
  prices <- if (pre) model@pricePre else model@pricePost
  numDeriv::jacobian(function(p) {
    m <- model
    if (pre) m@pricePre[active] <- p else m@pricePost[active] <- p
    calcShares(m, pre, revenue = FALSE)[active]
  }, prices[active])
}

test_that("input LogitBLP solves its level FOCs with the analytic share Jacobian", {
  ## Guards replacing numDeriv by elast(partial = TRUE) in calcPrices and
  ## calcMargins: the analytic Jacobian must be the derivative of calcShares
  ## (flat and nested, pre and post, with a removed product), and the solved
  ## rates must satisfy the level FOCs evaluated with a numerical Jacobian.
  for (sigma_nest in c(1, .6)) {
    result <- simulate(input_blp_fit(sigma_nest), c("A", "A", "C", "D"), price_domain = "real")
    gone <- result
    gone@subset <- c(TRUE, TRUE, FALSE, TRUE)
    for (m in list(result, gone)) for (pre in c(TRUE, FALSE)) {
      active <- if (pre) rep(TRUE, 4) else m@subset
      expect_equal(unname(elast(m, pre, partial = TRUE)[active, active]),
                   numeric_share_jacobian(m, pre), tolerance = 1e-7)
    }
    s <- calcShares(result, FALSE, revenue = FALSE)
    owner <- result@ownerPost
    residual <- (s * diag(owner) - as.vector((t(numeric_share_jacobian(result, FALSE)) * owner) %*%
      (result@mcPost - result@pricePost))) / s
    expect_lt(max(abs(residual)), 1e-6)
    expect_equal(unname(calcMargins(result, FALSE, level = TRUE)),
                 unname(as.vector(solve(t(numeric_share_jacobian(result, FALSE)) * owner, s * diag(owner)))),
                 tolerance = 1e-6)
  }
})

test_that("LogitBLP subclasses take the same input-market solver and margins as LogitBLP", {
  ## Guards exact class-name dispatch (class(object)[[1L]] %in% ...), which sent
  ## subclasses to the generic fallback solver and the revenue-elasticity margins.
  methods::setClass("TestChildLogitBLP", contains = "LogitBLP", where = globalenv())
  on.exit(methods::removeClass("TestChildLogitBLP", where = globalenv()), add = TRUE)
  for (sigma_nest in c(1, .6)) {
    parent <- simulate(input_blp_fit(sigma_nest), c("A", "A", "C", "D"), price_domain = "real")
    child <- methods::new("TestChildLogitBLP", parent)
    for (pre in c(TRUE, FALSE)) {
      expect_identical(calcMargins(child, pre, level = TRUE), calcMargins(parent, pre, level = TRUE))
    }
    expect_identical(calcPrices(child, FALSE), calcPrices(parent, FALSE))
  }
})
