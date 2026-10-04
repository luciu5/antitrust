test_that("input Logit signed rates require an explicit real domain", {
  prices <- c(.2, .1, .3)
  shares <- c(.405, .315, .18)
  args <- list(
    prices = prices, supply = "bertrand", demand = "Logit",
    demand.param = list(alpha = .8,
      meanval = log(shares / .1) - .8 * (prices - 1)),
    ownerPre = c("A", "B", "C"), ownerPost = c("A", "A", "C"),
    priceOutside = 1, insideSize = 100, mcDelta = rep(-.2, 3)
  )
  domain_error <- tryCatch(do.call(sim, args), error = identity)
  expect_s3_class(domain_error, "antitrust_price_domain_error")
  expect_identical(domain_error$category, "positive_domain_violation")
  market <- do.call(sim, c(args, list(price_domain = "real")))
  expect_true(any(market@pricePost < 0))
  post_share <- calcShares(market, FALSE, revenue = FALSE)
  derivative <- .8 * (diag(post_share) - tcrossprod(post_share))
  foc <- post_share - as.vector((t(derivative) * market@ownerPost) %*%
    (market@mcPost - market@pricePost))
  expect_lt(max(abs(foc)), 1e-7)

  zero <- market
  zero@pricePre[[1L]] <- 0
  expect_true(all(is.finite(calcMargins(zero, TRUE, level = TRUE))))
  expect_true(is.na(calcMargins(zero, TRUE)[[1L]]))
  expect_error(do.call(sim, c(args, list(price_domain = "invalid"))),
               "arg")
  output_args <- args
  output_args$demand.param$alpha <- -.8
  expect_error(do.call(sim, c(output_args, list(price_domain = "real"))),
               "supports input")
})

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

test_that("nested input Logit retains input orientation under real rates", {
  market <- suppressWarnings(suppressMessages(sim(
    prices = c(.2, .1, .3, .15), shares = c(.32, .25, .18, .15),
    supply = "bertrand", demand = "LogitNests",
    demand.param = list(alpha = .8, meanval = c(0, -.1, .2, .05),
                        sigma = c(.7, .8)),
    ownerPre = c("A", "B", "C", "D"),
    ownerPost = c("A", "A", "C", "D"),
    nests = c("N1", "N1", "N2", "N2"), insideSize = 100,
    mcDelta = rep(-.2, 4), price_domain = "real"
  )))
  expect_false(market@output)
  expect_true(any(market@pricePost < 0))
  expect_lt(max(abs(market@mcPost - market@pricePost -
                    calcMargins(market, FALSE, level = TRUE))), 1e-6)
})

test_that("signed input rates reach the supported legacy conduct solvers", {
  prices <- c(.2, .1, .3)
  shares <- c(.405, .315, .18)
  for (conduct in c("cournot", "auction2nd", "bargaining",
                    "bargaining2nd")) {
    market <- suppressWarnings(suppressMessages(sim(
      prices = prices, supply = conduct, demand = "Logit",
      demand.param = list(alpha = .8,
        meanval = log(shares / .1) - .8 * (prices - 1)),
      ownerPre = c("A", "B", "C"),
      ownerPost = c("A", "A", "C"),
      priceOutside = 1, insideSize = 100,
      mcDelta = rep(-.2, 3), price_domain = "real"
    )))
    expect_false(market@output)
    expect_true(all(is.finite(market@pricePost)))
    expect_true(any(market@pricePost < 0))
  }
})

test_that("capacity input FOCs remain regular at signed rates", {
  prices <- c(.2, .1, .3)
  shares <- c(.405, .315, .18)
  args <- list(
    prices = prices, shares = shares, supply = "bertrand",
    demand = "LogitCap",
    demand.param = list(alpha = .8,
      meanval = log(shares / .1) - .8 * prices, mktSize = 100),
    ownerPre = c("A", "B", "C"), ownerPost = c("A", "A", "C"),
    capacities = c(45, 35, 20), insideSize = 100,
    mcDelta = rep(-.2, 3)
  )
  expect_s3_class(tryCatch(suppressWarnings(do.call(sim, args)),
                           error = identity), "antitrust_price_domain_error")
  market <- suppressWarnings(do.call(sim,
    c(args, list(price_domain = "real"))))
  expect_true(any(market@pricePost < 0))
  quantities <- calcQuantities(market, FALSE)
  expect_true(all(quantities <= market@capacitiesPost + 1e-6))
  post_share <- calcShares(market, FALSE, revenue = FALSE)
  derivative <- market@mktSize * .8 *
    (diag(post_share) - tcrossprod(post_share))
  foc <- quantities - as.vector((t(derivative) * market@ownerPost) %*%
    (market@mcPost - market@pricePost))
  slack <- quantities < market@capacitiesPost - 1e-5
  expect_true(any(slack))
  expect_lt(max(abs(foc[slack])), 1e-6)
  expect_true(all(foc[!slack] <= 1e-6))
})

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
