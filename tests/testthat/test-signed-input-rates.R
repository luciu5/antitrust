test_that("input Logit simulation requires an explicit real domain for signed rates", {
  p <- c(.2, .1, .3)
  q <- c(.405, .315, .18)
  s0 <- .1
  for (conduct in c("bertrand", "moncom", "cournot", "auction2nd",
                    "bargaining", "bargaining2nd")) {
    args <- list("logit", conduct, prices = p,
      parameters = list(alpha = .8,
        meanval = log(q / s0) - .8 * (p - 1)),
      ownerPre = c("a", "b", "c"),
      priceOutside = 1, insideSize = 100, output = FALSE,
      baseline = "observed")
    if (conduct %in% c("bargaining", "bargaining2nd"))
      args$bargpowerPre <- rep(.5, 3)
    fit <- do.call(specify, args)
    domain_error <- tryCatch(
      simulate(fit, c("a", "a", "c"), mcDelta = rep(-.2, 3)),
      error = identity)
    expect_s3_class(domain_error, "antitrust_price_domain_error")
    expect_identical(domain_error$category, "positive_domain_violation")
    expect_lt(domain_error$minimum_rate, 0)
    result <- simulate(fit, c("a", "a", "c"), mcDelta = rep(-.2, 3),
      price_domain = "real")
    expect_true(any(result@pricePost < 0))
    expect_lt(max(abs(result@mcPost - result@pricePost -
      calcMargins(result, FALSE, level = TRUE))), 1e-6)
    if (conduct == "moncom") {
      zero <- fit@model
      zero@pricePre[1] <- 0
      expect_true(all(is.finite(calcMargins(zero, TRUE, level = TRUE))))
    }
  }
  direct_args <- list(prices = p, supply = "bertrand", demand = "Logit",
    demand.param = list(alpha = .8, meanval = log(q / s0) - .8 * (p - 1)),
    ownerPre = c("a", "b", "c"), ownerPost = c("a", "a", "c"),
    priceOutside = 1, insideSize = 100, output = FALSE,
    mcDelta = rep(-.2, 3))
  expect_s3_class(tryCatch(do.call(sim, direct_args), error = identity),
                  "antitrust_price_domain_error")
  direct <- do.call(sim, c(direct_args, list(price_domain = "real")))
  expect_true(any(direct@pricePost < 0))
})

test_that("signed input rates persist through sequential and resumed policies", {
  p <- c(.2, .1, .3)
  q <- c(.405, .315, .18)
  fit <- specify(
    "logit", "bertrand", prices = p,
    parameters = list(alpha = .8,
      meanval = log(q / .1) - .8 * (p - 1)),
    ownerPre = c("a", "b", "c"), priceOutside = 1,
    insideSize = 100, output = FALSE, baseline = "observed")
  first <- counterfactual(ownership = c("a", "a", "c"),
                          costs = rep(-.2, 3))
  policy <- add_step(first, costs = rep(.1, 3))
  path <- simulate(fit, policy, price_domain = "real")
  expect_s4_class(path, "CounterfactualPath")
  expect_true(any(result_at(path, 1)@pricePost < 0))
  expect_true(any(final_result(path)@pricePost < 0))
  expect_equal(result_at(path, 2)@pricePre,
               result_at(path, 1)@pricePost)
  resumed <- simulate(path, counterfactual(costs = rep(0, 3)),
                      price_domain = "real")
  expect_true(any(final_result(resumed)@pricePost < 0))
  expect_s3_class(tryCatch(simulate(path,
    counterfactual(costs = rep(0, 3))), error = identity),
    "antitrust_price_domain_error")
})

test_that("input Logit dollar markdowns stay finite at an exact zero rate", {
  p <- c(.2, .1, .3)
  q <- c(.405, .315, .18)
  fit <- specify(
    "logit", "bertrand", prices = p,
    parameters = list(alpha = .8,
      meanval = log(q / .1) - .8 * (p - 1)),
    ownerPre = c("a", "b", "c"), priceOutside = 1,
    insideSize = 100, output = FALSE, baseline = "observed")
  model <- fit@model
  model@pricePre[1] <- 0
  level <- calcMargins(model, TRUE, level = TRUE)
  share <- calcShares(model, TRUE, revenue = FALSE)
  derivative <- .8 * (diag(share) - tcrossprod(share))
  owner <- ownerToMatrix(model, TRUE)
  expect_true(all(is.finite(level)))
  expect_lt(max(abs(as.vector((t(derivative) * owner) %*% level) - share)),
            1e-10)
  expect_true(is.na(calcMargins(model, TRUE)[1]))
})

test_that("nested Logit and BLP preserve input orientation and signed rates", {
  nested_args <- list(
    demand = "logit_nests", conduct = "bertrand",
    prices = c(.2, .1, .3, .15),
    shares = c(.32, .25, .18, .15),
    ownerPre = c("A", "B", "C", "D"),
    parameters = list(alpha = .8, meanval = c(0, -.1, .2, .05),
                      sigma = c(.7, .8)),
    nests = c("N1", "N1", "N2", "N2"),
    insideSize = 100)
  nested <- suppressWarnings(do.call(specify,
    c(nested_args, list(output = FALSE))))
  inferred <- suppressWarnings(do.call(specify, nested_args))
  expect_false(inferred@model@output)
  blp <- suppressWarnings(specify(
    "blp", "bertrand", prices = c(.2, .1, .3),
    parameters = list(alphaMean = .8, sigma = .1,
                      meanval = c(0, .1, .2)),
    ownerPre = c("A", "B", "C"), insideSize = 100,
    output = FALSE))
  for (fit in list(nested, blp)) {
    expect_false(fit@model@output)
    ownership <- if (fit@spec$demand == "blp")
      c("A", "A", "C") else c("A", "A", "C", "D")
    shock <- rep(-.2, length(ownership))
    positive <- tryCatch(simulate(fit, ownership, mcDelta = shock),
                         error = identity)
    expect_s3_class(positive, "antitrust_price_domain_error")
    result <- simulate(fit, ownership, mcDelta = shock,
                       price_domain = "real")
    expect_true(any(result@pricePost < 0))
    expect_lt(max(abs(result@mcPost - result@pricePost -
                      calcMargins(result, FALSE, level = TRUE))), 1e-6)
    zero <- fit@model
    zero@pricePre[1] <- 0
    level <- calcMargins(zero, TRUE, level = TRUE)
    share <- calcShares(zero, TRUE, revenue = FALSE)
    expect_true(all(is.finite(level)))
    h <- 1e-5
    derivative <- vapply(seq_along(share), function(j) {
      above <- below <- zero
      above@pricePre[j] <- above@pricePre[j] + h
      below@pricePre[j] <- below@pricePre[j] - h
      (calcShares(above, TRUE, revenue = FALSE) -
         calcShares(below, TRUE, revenue = FALSE)) / (2 * h)
    }, numeric(length(share)))
    expect_lt(max(abs(as.vector((t(derivative) *
      ownerToMatrix(zero, TRUE)) %*% level) - share)), 1e-5)
  }
})

test_that("capacity Logit input rates satisfy level FOCs and capacity bounds", {
  prices <- c(.2, .1, .3)
  shares <- c(.405, .315, .18)
  fit <- suppressWarnings(specify(
    "logit_cap", "bertrand", prices = prices, shares = shares,
    parameters = list(alpha = .8,
                      meanval = log(shares / .1) - .8 * prices,
                      mktSize = 100),
    ownerPre = c("A", "B", "C"), insideSize = 100,
    capacities = c(45, 35, 20), output = FALSE))
  neutral <- simulate(fit, c("A", "B", "C"), mcDelta = rep(0, 3))
  expect_equal(unname(neutral@pricePost), prices, tolerance = 1e-7)
  ownership <- c("A", "A", "C")
  shock <- rep(-.2, 3)
  expect_s3_class(tryCatch(simulate(fit, ownership, mcDelta = shock),
                           error = identity), "antitrust_price_domain_error")
  result <- simulate(fit, ownership, mcDelta = shock,
                     price_domain = "real")
  expect_true(any(result@pricePost < 0))
  quantities <- calcQuantities(result, FALSE)
  capacities <- result@capacitiesPost
  expect_true(all(quantities <= capacities + 1e-6))
  share <- calcShares(result, FALSE, revenue = FALSE)
  derivative <- result@mktSize * result@slopes$alpha *
    (diag(share) - tcrossprod(share))
  foc <- quantities - as.vector((t(derivative) *
    ownerToMatrix(result, FALSE)) %*% (result@mcPost - result@pricePost))
  slack <- quantities < capacities - 1e-5
  expect_true(any(slack))
  expect_lt(max(abs(foc[slack])), 1e-6)
  expect_true(all(foc[!slack] <= 1e-6))
  binding <- suppressWarnings(simulate(fit, c("A", "B", "C"),
                      mcDelta = c(.1, -.2, -.2),
                      capacitiesPost = c(41, 35, 20),
                      price_domain = "real"))
  expect_true(any(binding@pricePost < 0))
  binding_quantities <- calcQuantities(binding, FALSE)
  binding_shares <- calcShares(binding, FALSE, revenue = FALSE)
  binding_derivative <- binding@mktSize * binding@slopes$alpha *
    (diag(binding_shares) - tcrossprod(binding_shares))
  binding_foc <- binding_quantities - as.vector((t(binding_derivative) *
    ownerToMatrix(binding, FALSE)) %*%
    (binding@mcPost - binding@pricePost))
  at_capacity <- abs(binding_quantities - binding@capacitiesPost) < 1e-5
  expect_true(any(at_capacity))
  expect_true(any(!at_capacity))
  expect_true(all(binding_quantities <= binding@capacitiesPost + 1e-6))
  expect_lt(max(abs(binding_foc[!at_capacity])), 1e-6)
  expect_true(all(binding_foc[at_capacity] <= 1e-6))
  initial_path <- simulate(fit,
    add_step(counterfactual(ownership = c("A", "B", "C")),
             costs = rep(0, 3)))
  resumed <- suppressWarnings(simulate(initial_path,
    counterfactual(costs = rep(-.2, 3)), price_domain = "real"))
  expect_true(any(final_result(resumed)@pricePost < 0))
  zero <- fit@model
  zero@pricePre[1] <- 0
  expect_true(all(is.finite(calcMargins(zero, TRUE, level = TRUE))))
})

test_that("BLP conduct models retain finite negative input rates", {
  for (conduct in c("moncom", "cournot", "auction2nd", "bargaining")) {
    args <- list("blp", conduct, prices = c(.2, .1, .3),
      parameters = list(alphaMean = .8, sigma = .1,
                        meanval = c(0, .1, .2)),
      ownerPre = c("A", "B", "C"), insideSize = 100,
      output = FALSE)
    if (conduct == "bargaining") args$bargpowerPre <- rep(.5, 3)
    fit <- suppressWarnings(do.call(specify, args))
    owner <- c("A", "A", "C")
    shock <- rep(-.2, 3)
    expect_s3_class(tryCatch(simulate(fit, owner, mcDelta = shock),
                             error = identity), "antitrust_price_domain_error")
    result <- suppressWarnings(simulate(fit, owner, mcDelta = shock,
                                        price_domain = "real"))
    expect_true(any(result@pricePost < 0))
    expect_lt(max(abs(result@mcPost - result@pricePost -
                      calcMargins(result, FALSE, level = TRUE))), 1e-6)
    zero <- fit@model
    zero@pricePre[1] <- 0
    expect_true(all(is.finite(calcMargins(zero, TRUE, level = TRUE))))
  }
})
