test_that("mixed product retention satisfies physical Bertrand profit FOCs", {
    prices <- c(2, 2.2, 2.5)
    retention <- c(.70, .85, 1.10)
    fit <- suppressWarnings(specify(
        "logit", "bertrand", prices = prices,
        parameters = list(alpha = -1.2, meanval = c(.3, .1, -.2)),
        ownerPre = c("A", "A", "B"), insideSize = 100,
        baseline = "observed", revenueRetentionPre = retention
    ))
    model <- fit@model
    owner <- model@ownerPre
    kappa <- model@mcPre
    physical_cost <- retention * kappa
    profit <- function(i, price_i) {
        changed <- model
        changed@pricePre[i] <- price_i
        p <- prices
        p[i] <- price_i
        q <- calcShares(changed, preMerger = TRUE)
        sum(owner[i, ] * (retention * p - physical_cost) * q)
    }
    for (i in seq_along(prices)) {
        gradient <- (profit(i, prices[i] + 1e-5) -
                     profit(i, prices[i] - 1e-5)) / 2e-5
        expect_lt(abs(gradient), 1e-6)
    }
    expect_equal(model@ownerPre, owner)
    expect_equal(getRetention(fit), retention)
    expect_equal(getRetention(fit, FALSE), retention)
    expect_equal(unname(calcProducerSurplus(model)),
                 unname(retention * (prices - kappa) *
                        calcQuantities(model)), tolerance = 1e-10)

    post <- c(.8, .9, 1.2)
    result <- simulate(fit, ownerPost = c("A", "A", "B"),
                       revenueRetentionPost = post)
    expect_equal(getRetention(result, FALSE), post)
    expect_equal(unname(result@mcPost * post),
                 unname(physical_cost), tolerance = 1e-10)
    expect_equal(result@ownerPost, owner)
})

test_that("retention API validates vectors and preserves omitted states", {
    fit <- suppressWarnings(specify(
        "logit", "bertrand", prices = c(2, 2.2),
        parameters = list(alpha = -1.2, meanval = c(.3, .1)),
        ownerPre = c("A", "A"), insideSize = 100,
        baseline = "observed"
    ))
    expect_equal(getRetention(fit), c(1, 1))
    changed <- setRetention(fit, retentionPre = c(Prod2 = .9, Prod1 = .8))
    expect_equal(getRetention(changed), c(.8, .9))
    expect_equal(getRetention(changed, FALSE), c(.8, .9))
    changed <- setRetention(changed, retentionPost = 1.2)
    expect_equal(getRetention(changed), c(.8, .9))
    expect_equal(getRetention(changed, FALSE), c(1.2, 1.2))
    expect_error(setRetention(fit, retentionPre = 0), "positive finite")
    expect_error(setRetention(fit, retentionPre = c(.8, NA)), "positive finite")
    expect_error(setRetention(fit, retentionPre = c(.8, .9, 1)), "one value per product")
})

test_that("mixed-retention Logit bargaining satisfies product Nash gradients", {
    prices <- c(2, 2.2, 2.5)
    retention <- c(.70, .85, 1.10)
    fit <- suppressWarnings(specify(
        "logit", "bargaining", prices = prices,
        parameters = list(alpha = -1.2, meanval = c(.3, .1, -.2)),
        ownerPre = c("A", "A", "B"), insideSize = 100,
        baseline = "observed", bargpowerPre = c(.2, .3, .1),
        revenueRetentionPre = retention
    ))
    model <- fit@model
    owner <- model@ownerPre
    kappa <- model@mcPre
    value <- function(i, price_i) {
        changed <- model
        changed@pricePre[i] <- price_i
        q <- calcShares(changed)
        delta <- q - q / (1 - q[i])
        delta[i] <- q[i]
        p <- prices
        p[i] <- price_i
        seller_gain <- sum(owner[i, ] * retention *
                           (p - kappa) * delta)
        buyer_gain <- log1p(-q[i]) / model@slopes$alpha
        b <- model@bargpowerPre[i]
        b * log(buyer_gain) + (1 - b) * log(seller_gain)
    }
    for (i in seq_along(prices)) {
        gradient <- (value(i, prices[i] + 1e-5) -
                     value(i, prices[i] - 1e-5)) / 2e-5
        expect_lt(abs(gradient), 1e-6)
    }
    expect_equal(owner, model@ownerPre)
})

test_that("mixed-retention CES Cournot satisfies physical quantity FOCs", {
    prices <- c(2, 2.2, 2.4, 2.6)
    retention <- c(1, .9, .7, .6)
    fit <- suppressWarnings(specify(
        "ces", "cournot", prices = prices,
        parameters = list(gamma = 5, meanval = c(1, .8, .6, .4)),
        ownerPre = c("A", "A", "B", "C"), insideSize = 100,
        baseline = "observed", revenueRetentionPre = retention
    ))
    model <- fit@model
    q <- calcQuantities(model)
    c <- retention * model@mcPre
    owner <- model@ownerPre
    profit_at_quantity <- function(i, change) {
        target <- q
        target[i] <- target[i] + change
        demand_error <- function(log_prices) {
            candidate <- model
            candidate@pricePre <- exp(log_prices)
            (calcQuantities(candidate) - target) / pmax(q, 1)
        }
        root <- nleqslv::nleqslv(log(prices), demand_error,
                                 control = list(ftol = 1e-12))
        expect_lt(max(abs(demand_error(root$x))), 1e-9)
        candidate <- model
        candidate@pricePre <- exp(root$x)
        sum(owner[i, ] * (retention * candidate@pricePre - c) *
            calcQuantities(candidate))
    }
    for (i in seq_along(prices)) {
        gradient <- (profit_at_quantity(i, 1e-4) -
                     profit_at_quantity(i, -1e-4)) / 2e-4
        expect_lt(abs(gradient), 1e-6)
    }
    expect_equal(unname(calcMC(model, TRUE)), unname(model@mcPre),
                 tolerance = 1e-12)
    expect_equal(model@ownerPre, owner)
})

test_that("CES Cournot rejects a failed positive interior equilibrium", {
    fit <- suppressWarnings(specify(
        "ces", "cournot", prices = c(2, 2.2, 2.4, 2.6),
        parameters = list(gamma = 5, meanval = c(1, .8, .6, .4)),
        ownerPre = c("A", "A", "B", "C"), insideSize = 100,
        baseline = "observed"
    ))
    model <- setRetention(fit@model, retentionPost = c(1, .8, .7, .6))
    model@mcPost <- model@mcPre / getRetention(model, FALSE)
    expect_error(calcPrices(model, FALSE),
                 "failed to find a valid positive interior equilibrium")
})

test_that("mixed retention within an auction firm is rejected", {
    auction <- suppressWarnings(auction2nd.logit(
        prices = c(2, 2.2, 2.5), shares = c(.35, .25, .2),
        margins = c(.4, .35, .3), ownerPre = c("A", "A", "B"),
        ownerPost = c("A", "A", "B")
    ))
    mixed <- setRetention(auction, c(.8, .9, 1))
    expect_error(calcMargins(mixed, level = TRUE),
                 "mixed revenue retention within an auction firm")
    expect_error(calcShares(mixed),
                 "mixed revenue retention within an auction firm")
    tiny <- setRetention(auction, c(1e-12, 1e-11, 1))
    expect_error(calcShares(tiny),
                 "mixed revenue retention within an auction firm")
    expect_error(calcPrices(tiny),
                 "mixed revenue retention within an auction firm")
    uniform <- setRetention(auction, c(.8, .8, 1))
    expect_true(all(is.finite(calcMargins(uniform, level = TRUE))))
})
