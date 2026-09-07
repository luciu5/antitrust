moncom_logit_fixture <- function(output = TRUE) {
  prices <- c(2, 2.5, 3)
  shares <- c(.30, .25, .25)
  alpha <- if (output) -1.5 else 1.5
  margins <- 1 / (abs(alpha) * prices)
  list(
    prices = prices,
    shares = shares,
    margins = margins,
    ownerPre = c("A", "A", "B"),
    alpha = alpha
  )
}

moncom_ces_fixture <- function() {
  prices <- c(2, 2.5, 3)
  shares <- c(.30, .25, .25)
  gamma <- 2
  list(
    prices = prices,
    shares = shares,
    margins = 1 / (gamma - (gamma - 1) * shares),
    ownerPre = c("A", "A", "B"),
    gamma = gamma
  )
}

test_that("MonCom Logit calibration uses own-product FOCs", {
  x <- moncom_logit_fixture()
  fit <- calibrate("logit", "moncom", prices = x$prices,
                   shares = x$shares, margins = x$margins,
                   ownerPre = x$ownerPre, insideSize = 100)

  expect_s4_class(fit@model, "MonComLogit")
  expect_equal(unname(fit@parameters$alpha), x$alpha, tolerance = 1e-10)
  expect_equal(unname(calcShares(fit@model, TRUE)), x$shares,
               tolerance = 1e-10)
  expect_equal(unname(calcMargins(fit@model, TRUE, level = TRUE)),
               rep(-1 / x$alpha, length(x$prices)), tolerance = 1e-10)
  expect_equal(unname(calcMC(fit@model, TRUE)), x$prices + 1 / x$alpha,
               tolerance = 1e-10)
  shares <- calcShares(fit@model, TRUE)
  foc <- shares + (x$prices - calcMC(fit@model, TRUE)) *
    unname(fit@parameters$alpha) * shares
  expect_lt(max(abs(foc)), 1e-10)
  expect_lt(fit@diagnostics$foc_residual_pre, 1e-10)
  expect_lt(fit@diagnostics$foc_residual_post, 1e-10)
  expect_true(isTRUE(fit@diagnostics$ownership_irrelevant))

  diagnostics <- calcDiagnostics(fit@model)
  expect_lt(attr(diagnostics, "moncom")$foc_residual_pre, 1e-10)
})

test_that("MonCom CES calibration uses the CES own elasticity", {
  x <- moncom_ces_fixture()
  fit <- calibrate("ces", "moncom", prices = x$prices,
                   shares = x$shares, margins = x$margins,
                   ownerPre = x$ownerPre, insideSize = 100,
                   priceOutside = 1)

  expect_s4_class(fit@model, "MonComCES")
  expect_equal(unname(fit@parameters$gamma), x$gamma, tolerance = 1e-10)
  expect_equal(unname(calcShares(fit@model, TRUE, revenue = TRUE)), x$shares,
               tolerance = 1e-10)
  expected_margin <- 1 / (x$gamma - (x$gamma - 1) * x$shares)
  expect_equal(unname(calcMargins(fit@model, TRUE)), expected_margin,
               tolerance = 1e-10)
  quantity_shares <- calcShares(fit@model, TRUE)
  revenue_shares <- calcShares(fit@model, TRUE, revenue = TRUE)
  own_elast <- -x$gamma + (x$gamma - 1) * revenue_shares
  foc <- quantity_shares +
    (x$prices - calcMC(fit@model, TRUE)) * quantity_shares *
      own_elast / x$prices
  expect_lt(max(abs(foc)), 1e-10)
  expect_lt(fit@diagnostics$foc_residual_pre, 1e-10)
})

test_that("MonCom CES preserves the confirmed input-market sign convention", {
  x <- moncom_ces_fixture()
  gamma <- -2
  margins <- 1 / (-gamma + (gamma - 1) * x$shares)
  fit <- calibrate("ces", "moncom", prices = x$prices,
                   shares = x$shares, margins = margins,
                   ownerPre = x$ownerPre, output = FALSE,
                   insideSize = 100, priceOutside = 1)

  expect_equal(unname(fit@parameters$gamma), gamma, tolerance = 1e-10)
  expect_equal(unname(calcMargins(fit@model, TRUE)), margins,
               tolerance = 1e-10)
})

test_that("MonCom supports input-market Logit sign conventions", {
  x <- moncom_logit_fixture(output = FALSE)
  fit <- calibrate("logit", "moncom", prices = x$prices,
                   shares = x$shares, margins = x$margins,
                   ownerPre = x$ownerPre, output = FALSE,
                   insideSize = 100)

  expect_equal(unname(fit@parameters$alpha), x$alpha, tolerance = 1e-10)
  expect_true(isFALSE(fit@model@output))
  expect_equal(unname(calcMC(fit@model, TRUE)), x$prices + 1 / x$alpha,
               tolerance = 1e-10)
  expect_equal(unname(calcMargins(fit@model, TRUE, level = TRUE)),
               rep(1 / x$alpha, length(x$prices)), tolerance = 1e-10)
})

test_that("MonCom known primitives use the ordinary specify and sim APIs", {
  x <- moncom_logit_fixture()
  parameters <- list(
    alpha = x$alpha,
    meanval = log(x$shares / .2) - x$alpha * x$prices
  )
  fit <- specify("logit", "moncom", prices = x$prices,
                 parameters = parameters, ownerPre = x$ownerPre,
                 shares = x$shares, insideSize = 100)
  legacy <- sim(x$prices, shares = x$shares, ownerPre = x$ownerPre,
                ownerPost = x$ownerPre, supply = "moncom", demand = "Logit",
                demand.param = parameters, insideSize = 100)

  expect_s4_class(fit@model, "MonComLogit")
  expect_equal(fit@model@pricePost, legacy@pricePost, tolerance = 1e-10)
  expect_equal(fit@model@mcPre, legacy@mcPre, tolerance = 1e-10)
  expect_equal(fit@diagnostics$source, "specified")
})

test_that("MonCom is a first-class conduct transition", {
  x <- moncom_logit_fixture()
  moncom <- calibrate("logit", "moncom", prices = x$prices,
                      shares = x$shares, margins = x$margins,
                      ownerPre = x$ownerPre, insideSize = 100)

  bertrand <- respecify(moncom, conduct = "bertrand")
  expect_s4_class(bertrand@model, "Logit")
  expect_equal(bertrand@model@pricePre, x$prices, tolerance = 1e-10)
  expect_equal(bertrand@parameters$alpha, moncom@parameters$alpha,
               tolerance = 0)
  expect_equal(bertrand@diagnostics$transition$kind, "conduct_change")
  expect_true("source conduct-specific supply state" %in%
              bertrand@diagnostics$transition$discarded)
  expect_null(bertrand@diagnostics$calibration_args)
  expect_true(is.list(bertrand@diagnostics$source_calibration_args))

  back <- respecify(bertrand, conduct = "moncom")
  expect_s4_class(back@model, "MonComLogit")
  expect_equal(back@model@pricePre, x$prices, tolerance = 1e-10)
  expect_equal(back@diagnostics$transition$kind, "conduct_change")

  cournot <- respecify(moncom, conduct = "cournot")
  expect_s4_class(cournot@model, "LogitCournot")
  expect_equal(cournot@diagnostics$transition$kind, "conduct_change")

  expect_error(respecify(moncom, conduct = "bargaining"),
               "requires explicit target primitive")
  bargaining <- respecify(moncom, conduct = "bargaining",
                          bargpowerPre = rep(.5, length(x$prices)))
  expect_s4_class(bargaining@model, "BargainingLogit")

  updated <- update(moncom, margins = x$margins)
  expect_s4_class(updated@model, "MonComLogit")
  expect_equal(unname(updated@parameters$alpha),
               unname(moncom@parameters$alpha), tolerance = 1e-10)
})

test_that("MonCom Cournot transitions work for CES when both targets are complete", {
  x <- moncom_ces_fixture()
  moncom <- calibrate("ces", "moncom", prices = x$prices,
                      shares = x$shares, margins = x$margins,
                      ownerPre = x$ownerPre, insideSize = 100,
                      priceOutside = 1)
  cournot <- respecify(moncom, conduct = "cournot")
  expect_s4_class(cournot@model, "CESCournot")
  back <- respecify(cournot, conduct = "moncom")
  expect_s4_class(back@model, "MonComCES")
  expect_equal(back@model@pricePre, x$prices, tolerance = 1e-8)
})

test_that("ownership, costs, quality, and product-set changes use MonCom economics", {
  x <- moncom_logit_fixture()
  fit <- calibrate("logit", "moncom", prices = x$prices,
                   shares = x$shares, margins = x$margins,
                   ownerPre = x$ownerPre, insideSize = 100)
  owner_post <- c("A", "A", "A")

  ownership <- simulate(fit, counterfactual(ownership = owner_post))
  expect_equal(unname(ownership@pricePost), x$prices, tolerance = 1e-10)

  costs <- simulate(fit, counterfactual(costs = c(-.1, 0, 0)))
  expected_mc <- (x$prices + 1 / x$alpha) * c(.9, 1, 1)
  expect_equal(unname(costs@mcPost), expected_mc, tolerance = 1e-10)
  expect_equal(unname(costs@pricePost), expected_mc - 1 / x$alpha,
               tolerance = 1e-10)
  
  quality <- simulate(fit, counterfactual(quality = c(Prod1 = .1)))
  expect_gt(calcShares(quality, preMerger = FALSE)[1],
            calcShares(fit@model, preMerger = TRUE)[1])
  ## Flat Logit MonCom has a constant absolute markup, so a quality shock
  ## changes shares without changing the own-product price equation.
  expect_equal(unname(quality@pricePost), x$prices, tolerance = 1e-10)

  exited <- simulate(fit, counterfactual(exit = "Prod1"))
  expect_false(exited@subset[1])
  expect_true(all(is.finite(exited@pricePost[exited@subset])))

  entered <- simulate(
    fit,
    counterfactual(entry = entrant("New", meanval = 3, cost = 1.8,
                                   priceStart = 2.5))
  )
  expect_s4_class(entered, "MonComLogit")
  expect_equal(length(entered@labels), 4L)
  expect_true(all(is.finite(entered@pricePost[entered@subset])))
})

test_that("Bertrand and MonCom differ on multiproduct ownership but both respond to costs", {
  x <- moncom_logit_fixture()
  moncom <- calibrate("logit", "moncom", prices = x$prices,
                      shares = x$shares, margins = x$margins,
                      ownerPre = x$ownerPre, insideSize = 100)
  bertrand <- specify(
    "logit", "bertrand", prices = x$prices,
    parameters = list(alpha = x$alpha,
                      meanval = log(x$shares / .2) - x$alpha * x$prices),
    ownerPre = x$ownerPre, shares = x$shares, insideSize = 100
  )

  moncom_merger <- simulate(moncom,
                            counterfactual(ownership = c("A", "A", "A")))
  bertrand_merger <- simulate(bertrand,
                              counterfactual(ownership = c("A", "A", "A")))
  expect_equal(unname(moncom_merger@pricePost), x$prices, tolerance = 1e-10)
  expect_true(any(abs(bertrand_merger@pricePost - x$prices) > 1e-6))

  moncom_cost <- simulate(moncom, counterfactual(costs = c(-.1, 0, 0)))
  bertrand_cost <- simulate(bertrand, counterfactual(costs = c(-.1, 0, 0)))
  expect_true(abs(moncom_cost@pricePost[1] - x$prices[1]) > 1e-6)
  expect_true(abs(bertrand_cost@pricePost[1] - x$prices[1]) > 1e-6)
})

test_that("unsupported MonCom demand families fail explicitly", {
  for (demand in c("logit_nests", "ces_nests", "blp", "linear",
                   "loglin", "aids", "pcaids", "pcaids_nests")) {
    expect_error(model_spec(demand, "moncom"), "currently not supported")
  }
})

test_that("flat Logit MonCom has the atomistic markup, not Bertrand's share wedge", {
  x <- moncom_logit_fixture()
  moncom <- calibrate("logit", "moncom", prices = x$prices,
                      shares = x$shares, margins = x$margins,
                      ownerPre = x$ownerPre, insideSize = 100)
  bertrand <- calibrate("logit", "bertrand", prices = x$prices,
                        shares = x$shares,
                        margins = c(1 / (1.5 * x$prices[1] * (1 - x$shares[1])),
                                    NA, NA),
                        ownerPre = c("A", "B", "C"), insideSize = 100)

  expect_equal(unname(calcMargins(moncom@model, TRUE, level = TRUE)),
               rep(-1 / x$alpha, length(x$prices)), tolerance = 1e-10)
  bertrand_alpha <- unname(bertrand@parameters$alpha)
  expect_equal(unname(calcMargins(bertrand@model, TRUE, level = TRUE)),
               -1 / (bertrand_alpha * (1 - x$shares)), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(
    unname(calcMargins(moncom@model, TRUE, level = TRUE)),
    unname(calcMargins(bertrand@model, TRUE, level = TRUE)))))
})
