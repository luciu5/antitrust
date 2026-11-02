# Category A: permanent economic/API regression; input level margins survive zero aggregate revenue.
test_that("input level margins survive zero aggregate revenue", {
  market <- suppressWarnings(sim(
    prices = c(.1, .1), supply = "bertrand", demand = "Logit",
    demand.param = list(alpha = .8, meanval = c(-.08, .08)),
    ownerPre = c("A", "B"), ownerPost = c("A", "A"),
    priceOutside = 0, insideSize = 100, price_domain = "real"
  ))
  market@pricePost <- c(.1, -.1)
  share <- calcShares(market, FALSE, revenue = FALSE)
  expect_equal(as.numeric(share), rep(1 / 3, 2), tolerance = 1e-12)
  level <- calcMargins(market, FALSE, level = TRUE)
  derivative <- .8 * (diag(share) - tcrossprod(share))
  expect_true(all(is.finite(level)))
  expect_lt(max(abs(as.vector((t(derivative) * market@ownerPost) %*%
                               level) - share)), 1e-10)
})

# Category A: permanent economic/API regression; input BLP solves its level FOCs instead of accepting a tail root.
test_that("input BLP solves its level FOCs instead of accepting a tail root", {
  args <- list(
    prices = c(.2, .1, .3), supply = "bertrand", demand = "BLP",
    demand.param = list(alpha = .8, sigma = .1,
                        meanval = c(0, .1, .2),
                        integration = "gauss-hermite", nNodes = 7L),
    ownerPre = c("A", "B", "C"), ownerPost = c("A", "A", "C"),
    insideSize = 100, mcDelta = rep(-.2, 3)
  )
  expect_s3_class(tryCatch(suppressMessages(do.call(sim, args)),
                           error = identity), "antitrust_price_domain_error")
  market <- suppressMessages(do.call(sim,
    c(args, list(price_domain = "real"))))
  expect_true(all(is.finite(market@pricePost)))
  expect_true(any(market@pricePost < 0))
  share_at <- function(prices) {
    changed <- market
    changed@pricePost <- prices
    as.numeric(calcShares(changed, FALSE, revenue = FALSE))
  }
  shares <- share_at(market@pricePost)
  derivative <- numDeriv::jacobian(share_at,
                                    as.numeric(market@pricePost))
  foc <- shares - as.vector((t(derivative) * market@ownerPost) %*%
    (market@mcPost - market@pricePost))
  expect_lt(max(abs(foc)), 1e-6)
  expect_lt(max(abs(market@mcPost - market@pricePost -
                    calcMargins(market, FALSE, level = TRUE))), 1e-6)
})

# Category A: permanent economic/API regression; signed input BLP uses both quadrature factors in its FOCs.
test_that("signed input BLP uses both quadrature factors in its FOCs", {
  market <- suppressMessages(suppressWarnings(sim(
    prices = c(.2, .1, .3), supply = "bertrand", demand = "BLP",
    demand.param = list(alpha = .8, sigma = .05, piDemog = .04,
      meanval = c(0, .1, .2), integration = "gauss-hermite",
      nNodes = c(5L, 7L)),
    ownerPre = c("A", "B", "C"), ownerPost = c("A", "A", "C"),
    insideSize = 100, mcDelta = rep(-.2, 3), price_domain = "real"
  )))
  expect_equal(dim(market@slopes$integrationPoints), c(35L, 2L))
  expect_true(any(market@pricePost < 0))
  share_at <- function(prices) {
    changed <- market
    changed@pricePost <- prices
    as.numeric(calcShares(changed, FALSE, revenue = FALSE))
  }
  shares <- share_at(market@pricePost)
  derivative <- numDeriv::jacobian(share_at,
                                    as.numeric(market@pricePost))
  foc <- shares - as.vector((t(derivative) * market@ownerPost) %*%
    (market@mcPost - market@pricePost))
  expect_lt(max(abs(foc)), 1e-6)
})

# Category A: permanent economic/API regression; signed input Cournot BLP satisfies its inverse-demand FOCs.
test_that("signed input Cournot BLP satisfies its inverse-demand FOCs", {
  market <- suppressMessages(suppressWarnings(sim(
    prices = c(.2, .1, .3), shares = c(.405, .315, .18),
    supply = "cournot", demand = "BLP",
    demand.param = list(alpha = .8, sigma = .05, piDemog = .04,
      meanval = c(.5, .1, .2), demogMean = .4,
      demogCov = matrix(.25, 1L, 1L),
      integration = "gauss-hermite", nNodes = c(5L, 7L)),
    ownerPre = c("A", "B", "C"), ownerPost = c("A", "A", "C"),
    insideSize = 100, mcDelta = rep(-.2, 3), price_domain = "real"
  )))
  expect_equal(dim(market@slopes$integrationPoints), c(35L, 2L))
  expect_equal(as.numeric(market@pricePre), c(.2, .1, .3),
               tolerance = 1e-8)
  expect_true(all(is.finite(market@pricePost)))
  expect_true(any(market@pricePost < 0))

  share_at <- function(prices) {
    changed <- market
    changed@pricePost <- prices
    as.numeric(calcShares(changed, FALSE, revenue = FALSE))
  }
  prices <- as.numeric(market@pricePost)
  shares <- share_at(prices)
  derivative <- numDeriv::jacobian(share_at, prices)
  owner <- unname(market@ownerPost)
  cournot_foc <- as.numeric(market@mcPost) - prices -
    as.vector((owner * t(solve(derivative))) %*% shares)
  bertrand_foc <- shares - as.vector((t(derivative) * owner) %*%
    (as.numeric(market@mcPost) - prices))
  expect_lt(max(abs(cournot_foc)), 1e-6)
  expect_gt(max(abs(bertrand_foc)), 1e-3)
})

# Category A: permanent economic/API regression; legacy vector price draws keep independent demographic draws.
test_that("legacy vector price draws keep independent demographic draws", {
  price_draws <- seq(-1, 1, length.out = 5L)
  set.seed(276)
  market <- suppressMessages(suppressWarnings(sim(
    prices = c(2, 2.2), shares = c(.4, .3),
    demand = "BLP", supply = "bertrand",
    demand.param = list(alpha = -1.5, sigma = .1,
                        piDemog = .05, consDraws = price_draws),
    ownerPre = c("A", "B"), ownerPost = c("A", "A"),
    insideSize = 100
  )))
  expect_identical(market@slopes$consDraws, price_draws)
  expect_equal(dim(market@slopes$demogDraws), c(5L, 1L))
  expect_equal(market@slopes$alphas,
               -1.5 + .1 * price_draws +
                 .05 * market@slopes$demogDraws[, 1L],
               tolerance = 1e-13)
  expect_true(is.finite(CV(market)))
})

# Category A: permanent economic/API regression; recalibrating a legacy BLP object keeps random characteristics.
test_that("recalibrating a legacy BLP object keeps random characteristics", {
  characteristics <- matrix(c(1, 2, 3), ncol = 1L)
  market <- suppressMessages(suppressWarnings(sim(
    prices = c(2, 2.2, 2.5), shares = c(.35, .25, .20),
    supply = "bertrand", demand = "BLP",
    demand.param = list(alpha = -1.5, sigma = .1,
      prodChar = characteristics, beta = 0, sigmaChar = .2,
      integration = "gauss-hermite", nNodes = c(5L, 5L)),
    ownerPre = c("A", "B", "C"), ownerPost = c("A", "A", "C"),
    insideSize = 100
  )))
  legacy <- market
  legacy@slopes$integrationPoints <- NULL
  legacy@slopes$factorOrder <- NULL
  legacy@slopes$nodesPerAxis <- NULL
  legacy@slopes$charDraws <- NULL
  legacy@slopes$legacyVectorIntegration <- NULL
  set.seed(123)
  legacy <- suppressMessages(suppressWarnings(calcSlopes(legacy)))
  expect_gt(stats::sd(legacy@slopes$char_random[, 1L]), .01)
  expect_equal(legacy@slopes$char_random,
               sweep(legacy@slopes$charDraws, 2L, .2, "*") %*%
                 t(characteristics), tolerance = 1e-13)
})

# Category A: fractional self-ownership scales both sides of a firm's FOC;
# unit-diagonal ownership tests cannot detect a missing own-quantity weight.
test_that("input level FOCs preserve fractional self-ownership weights", {
  alpha <- .8
  meanval <- c(.2, .1, .3)
  owners <- diag(rep(.5, 3))
  for (demand in c("logit", "blp")) {
    params <- if (demand == "logit") list(alpha = alpha, meanval = meanval) else
      list(alphaMean = alpha, sigma = .05, meanval = meanval)
    fit <- suppressWarnings(specify(demand, "bertrand", prices = c(.2, .1, .3),
      parameters = params, ownerPre = owners, insideSize = 100, output = FALSE))
    result <- suppressWarnings(simulate(fit, ownerPost = owners,
      mcDelta = rep(-.1, 3), price_domain = "real"))
    share_at <- function(prices) {
      model <- result
      model@pricePost <- prices
      as.numeric(calcShares(model, FALSE, revenue = FALSE))
    }
    shares <- share_at(result@pricePost)
    derivative <- numDeriv::jacobian(share_at, as.numeric(result@pricePost))
    margins <- calcMargins(result, FALSE, level = TRUE)
    expect_lt(max(abs(shares * diag(owners) -
      as.vector((t(derivative) * owners) %*% margins))), 1e-7)
    expect_lt(max(abs(shares * diag(owners) -
      as.vector((t(derivative) * owners) %*%
        (result@mcPost - result@pricePost)))), 1e-6)
  }
})

# Category A: an absolute share-level residual can falsely certify tiny shares
# while the input price equation remains materially wrong.
test_that("rare-share input BLP checks relative FOCs before accepting prices", {
  result <- suppressMessages(suppressWarnings(sim(
    prices = c(.2, .1, .3), demand = "BLP", supply = "bertrand",
    demand.param = list(alpha = .8, sigma = .05, meanval = rep(-16, 3),
      integration = "gauss-hermite", nNodes = 7L),
    ownerPre = c("A", "B", "C"), ownerPost = c("A", "B", "C"),
    insideSize = 100, mcDelta = rep(-.2, 3), price_domain = "real")))
  share_at <- function(prices) {
    model <- result
    model@pricePost <- prices
    as.numeric(calcShares(model, FALSE, revenue = FALSE))
  }
  shares <- share_at(result@pricePost)
  derivative <- numDeriv::jacobian(share_at, as.numeric(result@pricePost))
  # Distinct owners remove cross-product profit effects.
  relative_foc <- (shares - diag(derivative) *
    (result@mcPost - result@pricePost)) / shares
  expect_lt(max(abs(relative_foc)), 1e-6)
  expect_lt(max(abs(result@mcPost - result@pricePost -
    calcMargins(result, FALSE, level = TRUE))), 1e-6)
})
