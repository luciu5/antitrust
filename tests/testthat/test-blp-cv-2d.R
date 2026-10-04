test_that("joint price and demographic quadrature gives the BLP CV integral", {
  market <- suppressMessages(suppressWarnings(sim(
    prices = c(2, 2.2, 2.5), shares = c(.35, .25, .20),
    supply = "bertrand", demand = "BLP",
    demand.param = list(alpha = -1.8, sigma = .12, piDemog = .09,
                        demogMean = .4, demogCov = matrix(.25, 1L, 1L),
                        integration = "gauss-hermite", nNodes = c(5L, 7L)),
    ownerPre = c("A", "B", "C"), ownerPost = c("A", "A", "C"),
    insideSize = 100
  )))

  first <- antitrust:::.blp_normal_nodes(5L)
  second <- antitrust:::.blp_normal_nodes(7L)
  price_z <- rep(first$nodes, 7L)
  demog_z <- rep(second$nodes, each = 5L)
  weights <- as.vector(outer(first$weights, second$weights))
  alpha <- -1.8 + .12 * price_z + .09 * .5 * demog_z
  expect_identical(market@slopes$factorOrder, c("price", "demog1"))
  expect_equal(unname(market@slopes$integrationPoints),
               unname(cbind(price_z, demog_z)), tolerance = 1e-14)
  expect_equal(market@slopes$drawWeights, weights, tolerance = 1e-14)
  expect_equal(market@slopes$alphas, alpha, tolerance = 1e-13)

  utility <- function(prices) {
    sweep(outer(alpha, prices - market@priceOutside), 2L,
          as.numeric(market@slopes$meanval), "+")
  }
  pre <- utility(market@pricePre)
  post <- utility(market@pricePost)
  draw_shares <- exp(pre) / (1 + rowSums(exp(pre)))
  expect_equal(as.numeric(market@shares),
               as.vector(crossprod(weights, draw_shares)),
               tolerance = 1e-7)
  manual_cv <- market@mktSize * sum(weights *
    (log1p(rowSums(exp(post))) - log1p(rowSums(exp(pre)))) / alpha)
  expect_equal(CV(market), manual_cv, tolerance = 1e-9)
})

test_that("supplied joint points and weights are not flattened", {
  points <- matrix(c(-1, 0, 1, 1, -1, .5), ncol = 2L)
  weights <- c(.2, .3, .5)
  rule <- calcBLPintegration(list(
    sigma = .1, piDemog = .2, nDemog = 1L,
    integrationPoints = points, integrationWeights = weights
  ))
  expect_identical(rule$integrationPoints, points)
  expect_identical(rule$factorOrder, c("price", "demog1"))
  expect_equal(rule$weights, weights, tolerance = 0)
  market <- suppressMessages(suppressWarnings(sim(
    prices = c(2, 2.2, 2.5), shares = c(.35, .25, .20),
    supply = "bertrand", demand = "BLP",
    demand.param = list(alpha = -1.8, sigma = .1, piDemog = .2,
                        integrationPoints = points,
                        integrationWeights = weights),
    ownerPre = c("A", "B", "C"), ownerPost = c("A", "A", "C"),
    insideSize = 100
  )))
  expect_identical(market@nDraws, 3L)
  expect_identical(market@slopes$integrationPoints, points)
  expect_equal(market@slopes$drawWeights, weights, tolerance = 0)
  expect_true(is.finite(CV(market)))
  expect_error(calcBLPintegration(list(
    sigma = .1, piDemog = .2, nDemog = 1L, draws = points
  )), "use 'integrationPoints'")
  expect_error(calcBLPintegration(list(
    sigma = .1, piDemog = .2, nDemog = 1L,
    integrationPoints = points[, 1, drop = FALSE]
  )), "at least two columns")
})

test_that("two demographic factors reproduce a correlated covariance", {
  covariance <- matrix(c(4, 1.2, 1.2, 2), 2L, 2L)
  mean <- c(.5, -1)
  rule <- calcBLPintegration(list(
    sigma = 0, nDemog = 2L, piDemog = c(.2, 0),
    integration = "gauss-hermite", nNodes = c(5L, 7L)
  ))
  draws <- antitrust:::.blp_materialize_draws(
    rule, alphaMean = -2, nDemog = 2L, piDemog = c(.2, 0),
    demogMean = mean, demogCov = covariance
  )$demogDraws
  centered <- sweep(draws, 2L, mean, "-")
  expect_identical(rule$factorOrder, c("demog1", "demog2"))
  expect_equal(as.vector(crossprod(rule$weights, draws)), mean,
               tolerance = 1e-12)
  expect_equal(unname(crossprod(sweep(centered, 1L,
                                      sqrt(rule$weights), "*"))),
               covariance, tolerance = 1e-12)
  expect_error(calcBLPintegration(list(
    sigma = .1, nDemog = 2L, piDemog = c(.2, .1),
    integration = "gauss-hermite", nNodes = 5L
  )), "at most two active dimensions")
})

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

test_that("price leadership retains legacy multi-factor Monte Carlo", {
  shares <- c(.35, .25, .25, .15)
  args <- list(
    prices = c(.93, .88, 1.10, 1.02), shares = shares,
    ownerPre = c("Bank1", "Bank2", "Bank3", "Fringe"),
    ownerPost = c("Bank1", "Bank2", "Bank3", "Fringe"),
    coalitionPre = 1:3, coalitionPost = 1:3, insideSize = 1000,
    slopes = list(alphaMean = -5.767013, alpha = -5.767013,
      sigma = .5, piDemog = .05, nDemog = 1L,
      meanval = c(0, log(shares[-1] / shares[1]) -
        (-5.767013) * (c(.88, 1.10, 1.02) - .93)),
      sigmaNest = 1),
    integration = "monte-carlo", nDraws = 12L
  )
  set.seed(2718)
  market <- suppressMessages(suppressWarnings(do.call(ple.blp, args)))
  expect_identical(market@slopes$integration, "monte-carlo")
  expect_length(market@slopes$consDraws, 12L)
  expect_equal(dim(market@slopes$demogDraws), c(12L, 1L))
  expect_true(is.finite(CV(market)))
  args$integration <- "gauss-hermite"
  args$nDraws <- NULL
  args$nNodes <- 5L
  expect_error(do.call(ple.blp, args),
               "Two-dimensional BLP quadrature is not supported")
})

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
