## Migration parity A: permanent independent economic checks.

test_that("Logit Bertrand preserves heterogeneous costs and all-product FOCs", {
    s <- c(.18, .12, .21, .19, .30)
    c <- c(40, 70, 55, 80, 100)
    fit <- synthetic_market(n_firms = 2, n_products = 2,
        shares = s, costs = c, reference_margin = .25, seed = 7)
    market <- fit@diagnostics$synthetic_market
    D <- diag(s) - tcrossprod(s)
    z <- solve(t(market$ownership * D), s)
    pref <- c[5] / .75
    alpha <- -z[5] / (pref - c[5])
    expect_equal(market$prices, unname(c - z / alpha), tolerance = 1e-8)
    expect_equal(market$costs, c)
    expect_equal(market$observed$costs, c)
    expect_equal(unname(fit@model@mcPre), c, tolerance = 1e-7)
    expect_equal(unname(calcShares(fit@model, TRUE)), s, tolerance = 1e-10)
    expect_equal(unname(fit@model@slopes$alpha), unname(alpha), tolerance = 1e-8)
    expect_lt(max(abs(s + t(market$ownership * (alpha * D)) %*%
                      (market$prices - c))), 1e-12)
    expect_lt(abs((market$prices[5] - c[5]) / market$prices[5] - .25), 1e-12)
    expect_equal(market$products$firm_id, c(1, 1, 2, 2, 3))
})

test_that("single-product Logit margin changes demand scale and prices", {
    s <- c(.2, .3, .5); costs <- c(20, 70, 100)
    fits <- lapply(c(.1, .3), function(m) synthetic_market(
        n_firms = 2, shares = s, costs = costs,
        reference_margin = m, seed = 12))
    for (j in seq_along(fits)) {
        m <- c(.1, .3)[j]
        markup <- costs[3] * m / (1 - m) * (1 - s[3]) / (1 - s)
        expect_equal(unname(fits[[j]]@model@pricePre),
                     unname(costs + markup), tolerance = 1e-8)
    }
    expect_true(all(fits[[2]]@model@pricePre > fits[[1]]@model@pricePre))
    expect_lt(abs(fits[[2]]@model@slopes$alpha),
              abs(fits[[1]]@model@slopes$alpha))
})

test_that("heterogeneous ownership and cost scales satisfy independent Logit FOCs", {
    s <- c(.08, .11, .12, .14, .18, .17, .20)
    for (scale in c(.5, 2)) {
        costs <- scale * c(20, 30, 40, 50, 60, 70, 80)
        fit <- synthetic_market(n_firms = 3, n_products = c(1, 2, 3),
            shares = s, costs = costs, reference_margin = .22, seed = 13)
        market <- fit@diagnostics$synthetic_market
        D <- fit@model@slopes$alpha * (diag(s) - tcrossprod(s))
        expect_lt(max(abs(s + t(market$ownership * D) %*%
                          (market$prices - costs))), 1e-10)
        expect_equal(market$costs, costs)
        expect_true(all(is.finite(market$prices) & market$prices > costs))
    }
})

test_that("each implemented conduct reproduces shares, costs, and native FOCs", {
    s <- c(.2, .2, .2, .3, .1)
    costs <- c(50, 60, 70, 80, 90)
    grid <- expand.grid(demand = c("logit", "ces"),
        supply = c("bertrand", "cournot", "moncom", "auction2nd",
                   "bargaining", "bargaining2nd"),
        stringsAsFactors = FALSE)
    for (k in seq_len(nrow(grid))) {
        d <- grid$demand[k]; conduct <- grid$supply[k]
        margin <- if (conduct == "cournot" && d == "ces") .3 else .2
        fit <- synthetic_market(demand = d, supply = conduct,
            n_firms = 2, n_products = 2, shares = s, costs = costs,
            reference_margin = margin, seed = 20)
        market <- fit@diagnostics$synthetic_market
        expect_equal(market$costs, costs, info = paste(d, conduct))
        expect_equal(unname(fit@model@mcPre), costs, tolerance = 1e-7,
                     info = paste(d, conduct))
        expect_lt(fit@diagnostics$synthetic$share_residual, 1e-9)
        expect_lt(fit@diagnostics$synthetic$foc_residual, 1e-8)
        expect_lt(abs(fit@diagnostics$synthetic$reference_margin_residual),
                  1e-10)
        expect_true(all(is.finite(market$prices) & market$prices > costs))
        if (grepl("bargaining", conduct)) {
            expect_equal(fit@model@bargpowerPre, rep(.5, 5))
        }
    }
})

test_that("conduct alternatives imply different prices with common design", {
    args <- list(n_firms = 2, n_products = 2,
        shares = c(.2, .2, .2, .3, .1), costs = c(50, 60, 70, 80, 90),
        reference_margin = .2, seed = 7)
    b <- do.call(synthetic_market, c(list(supply = "bertrand"), args))
    m <- do.call(synthetic_market, c(list(supply = "moncom"), args))
    a <- do.call(synthetic_market, c(list(supply = "auction2nd"), args))
    expect_false(isTRUE(all.equal(b@model@pricePre, m@model@pricePre)))
    expect_false(isTRUE(all.equal(b@model@pricePre, a@model@pricePre)))
})

test_that("Cournot, auction and bargaining invert their distinct margin equations", {
    s <- c(.2, .2, .2, .3, .1)
    costs <- c(50, 60, 70, 80, 90)
    owners <- c(1, 1, 2, 2, 3)
    firm_share <- as.numeric(tapply(s, owners, sum)[owners])
    args <- list(n_firms = 2, n_products = 2, shares = s,
                 costs = costs, seed = 20)
    logit_c <- do.call(synthetic_market, c(list(demand = "logit",
        supply = "cournot", reference_margin = .2), args))
    markup <- logit_c@model@pricePre - costs
    expect_equal(unname(markup / markup[5]),
                 (1 + firm_share / s[5]) / 2, tolerance = 1e-9)

    logit_a <- do.call(synthetic_market, c(list(demand = "logit",
        supply = "auction2nd", reference_margin = .2), args))
    auction_coef <- -log1p(-firm_share) / firm_share
    markup <- logit_a@model@pricePre - costs
    expect_equal(unname(markup / markup[5]),
                 auction_coef / auction_coef[5], tolerance = 1e-9)

    power <- c(.3, .4, .5, .6, .2)
    logit_b <- do.call(synthetic_market, c(list(demand = "logit",
        supply = "bargaining2nd", reference_margin = .2,
        bargpowerPre = power), args))
    markup <- logit_b@model@pricePre - costs
    expected <- (1 - power) * auction_coef
    expect_equal(unname(markup / markup[5]),
                 expected / expected[5], tolerance = 1e-9)

    ces_c <- do.call(synthetic_market, c(list(demand = "ces",
        supply = "cournot", reference_margin = .3), args))
    gamma <- (1 - firm_share[5]) / (.3 - firm_share[5])
    expected <- firm_share + (1 - firm_share) / gamma
    expect_equal(unname((ces_c@model@pricePre - costs) /
                 ces_c@model@pricePre), expected, tolerance = 1e-9)

    ces_a <- do.call(synthetic_market, c(list(demand = "ces",
        supply = "auction2nd", reference_margin = .2), args))
    gamma <- 1 + log1p(-firm_share[5]) / log1p(-.2)
    expected <- 1 - (1 - firm_share)^(1 / (gamma - 1))
    expect_equal(unname((ces_a@model@pricePre - costs) /
                 ces_a@model@pricePre), expected, tolerance = 1e-9)
})

test_that("nested Logit uses supplied nesting curvature and its own Jacobian", {
    shares <- c(.2, .1, .15, .25, .3)
    costs <- c(40, 50, 60, 70, 80)
    nests <- c("a", "a", "b", "b", "c")
    fit <- synthetic_market(demand = "logit_nests", supply = "bertrand",
        n_firms = 2, n_products = 2, shares = shares, costs = costs,
        reference_margin = .2, nests = nests,
        parameters = list(sigma = c(a = .7, b = .8, c = 1)), seed = 7)
    market <- fit@diagnostics$synthetic_market
    D <- numDeriv::jacobian(function(prices) {
        changed <- fit@model
        changed@pricePre <- prices
        as.numeric(calcShares(changed, TRUE))
    }, market$prices)
    expect_lt(max(abs(shares + t(market$ownership * D) %*%
                      (market$prices - costs))), 1e-8)
    expect_equal(market$costs, costs)
    expect_equal(unname(calcShares(fit@model, TRUE)), shares,
                 tolerance = 1e-10)
    expect_error(synthetic_market(demand = "logit_nests",
        n_firms = 2, n_products = 2, shares = shares, costs = costs,
        reference_margin = .2, nests = nests, seed = 7),
        "parameters\\$sigma")
})

test_that("nested CES identifies cross-nest curvature from a revenue FOC", {
    shares <- c(.2, .1, .15, .25, .3)
    costs <- c(40, 50, 60, 70, 80)
    nests <- c("a", "a", "b", "b", "b")
    sigma <- c(a = 5, b = 6)
    fit <- synthetic_market(demand = "ces_nests", supply = "bertrand",
        n_firms = 2, n_products = 2, shares = shares, costs = costs,
        reference_margin = .25, nests = nests,
        parameters = list(sigma = sigma), seed = 7)
    market <- fit@diagnostics$synthetic_market
    prices <- market$prices
    E <- numDeriv::jacobian(function(log_prices) {
        changed <- fit@model
        changed@pricePre <- exp(log_prices)
        log(as.numeric(calcShares(changed, TRUE, revenue = TRUE)))
    }, log(prices)) - diag(length(shares))
    foc <- shares + (t(E) * market$ownership) %*%
        (shares * (prices - costs) / prices)
    expect_lt(max(abs(foc)), 1e-7)
    expect_equal(unname(calcShares(fit@model, TRUE, revenue = TRUE)),
                 shares, tolerance = 1e-10)
    expect_equal(market$costs, costs)
    expect_gt(fit@model@slopes$gamma, 1)
    expect_lt(fit@model@slopes$gamma, min(sigma))
    ## A valid curvature can land exactly on the bounded solver grid. The
    ## absence of a strict sign crossing must not discard that solution.
    node_fit <- synthetic_market(demand = "ces_nests", supply = "bertrand",
        n_firms = 2, n_products = 2, shares = shares, costs = costs,
        reference_margin = 35 / 144, nests = nests,
        parameters = list(sigma = sigma), seed = 7)
    expect_equal(unname(node_fit@model@slopes$gamma), 3, tolerance = 1e-9)
    expect_error(synthetic_market(demand = "ces_nests",
        n_firms = 2, n_products = 2, shares = shares, costs = costs,
        reference_margin = .25, nests = nests, seed = 7),
        "parameters\\$sigma")
})

test_that("Linear Bertrand identifies only slope scale from the reference margin", {
    shares <- c(.2, .1, .15, .25, .3)
    costs <- c(40, 50, 60, 70, 80)
    shape <- matrix(.1, 5, 5); diag(shape) <- -2
    fit <- synthetic_market(demand = "linear", supply = "bertrand",
        n_firms = 2, n_products = 2, shares = shares, costs = costs,
        reference_margin = .2, parameters = list(slopes = shape), seed = 7)
    market <- fit@diagnostics$synthetic_market
    B <- fit@model@slopes
    expect_equal(B / B[1, 1], shape / shape[1, 1], tolerance = 1e-10)
    expect_lt(max(abs(shares + t(market$ownership * B) %*%
                      (market$prices - costs))), 1e-10)
    expect_equal(market$costs, costs)
    expect_equal(unname(calcQuantities(fit@model, TRUE)), shares,
                 tolerance = 1e-10)
    expect_error(synthetic_market(demand = "linear", supply = "bertrand",
        n_firms = 2, n_products = 2, shares = shares, costs = costs,
        reference_margin = .2, seed = 7), "parameters\\$slopes")
})

test_that("Log-linear Bertrand identifies elasticity scale from a singleton reference FOC", {
    shares <- c(.2, .1, .15, .25, .3)
    costs <- c(40, 50, 60, 70, 80)
    shape <- matrix(.1, 5, 5); diag(shape) <- -2
    fit <- synthetic_market(demand = "loglin", supply = "bertrand",
        n_firms = 2, n_products = 2, shares = shares, costs = costs,
        reference_margin = .2, parameters = list(slopes = shape), seed = 7)
    market <- fit@diagnostics$synthetic_market
    D <- numDeriv::jacobian(function(prices) {
        changed <- fit@model
        changed@pricePre <- prices
        as.numeric(calcQuantities(changed, TRUE))
    }, market$prices)
    foc <- shares + t(market$ownership * D) %*%
        (market$prices - costs)
    expect_lt(max(abs(foc)), 1e-7)
    expect_equal(unname(calcQuantities(fit@model, TRUE)), shares,
                 tolerance = 1e-10)
    expect_equal(market$costs, costs)
    expect_equal(fit@model@slopes / fit@model@slopes[1, 1],
                 shape / shape[1, 1], tolerance = 1e-10)
    expect_error(synthetic_market(demand = "loglin", supply = "bertrand",
        n_firms = 2, n_products = 2, shares = shares, costs = costs,
        reference_margin = .2, seed = 7), "parameters\\$slopes")
})

test_that("PCAIDS defaults identify demand from the reference elasticity", {
    shares <- c(.2, .1, .15, .25, .3)
    costs <- c(40, 50, 60, 70, 80)
    fit <- synthetic_market(demand = "pcaids", supply = "bertrand",
        n_firms = 2, n_products = 2, shares = shares, costs = costs,
        reference_margin = .2, seed = 7)
    market <- fit@diagnostics$synthetic_market
    scale <- (1 / .2 - 1) / (1 - shares[5])
    B <- -scale * (diag(shares) - tcrossprod(shares))
    expect_equal(unname(fit@model@slopes), B, tolerance = 1e-7)
    expect_equal(fit@model@mktElast, -1)
    E <- sweep(B, 1, shares, "/") - diag(length(shares))
    proportional_margin <- (market$prices - costs) / market$prices
    foc <- shares + (t(E) * market$ownership) %*%
        (shares * proportional_margin)
    expect_lt(max(abs(foc)), 1e-8)
    expect_equal(unname(calcShares(fit@model, TRUE, revenue = TRUE)),
                 shares, tolerance = 1e-10)
    expect_equal(market$costs, costs)
})

test_that("ALM Bertrand retains a passive outside share and active reference firm", {
    conditional <- c(.2, .1, .15, .25, .3)
    costs <- c(40, 50, 60, 70, 80)
    s0 <- .3
    actual <- (1 - s0) * conditional
    for (demand in c("LogitALM", "CESALM")) {
        fit <- synthetic_market(demand = demand, supply = "bertrand",
            n_firms = 2, n_products = 2, shares = conditional,
            costs = costs, reference_margin = .2,
            passive_outside_share = s0, seed = 7)
        market <- fit@diagnostics$synthetic_market
        prices <- market$prices
        expect_equal(market$shares, conditional)
        expect_equal(market$observed$unconditional_shares, actual)
        expect_equal(1 - fit@model@shareInside, s0, tolerance = 1e-9)
        expect_equal(market$costs, costs)
        if (fit@spec$demand == "logit") {
            D <- diag(actual) - tcrossprod(actual)
            foc <- actual + fit@model@slopes$alpha *
                t(market$ownership * D) %*% (prices - costs)
        } else {
            gamma <- fit@model@slopes$gamma
            E <- matrix(rep((gamma - 1) * actual,
                            each = length(actual)), length(actual))
            diag(E) <- diag(E) - gamma
            foc <- actual + (t(E) * market$ownership) %*%
                (actual * (prices - costs) / prices)
        }
        expect_lt(max(abs(foc)), 1e-8)
        expect_lt(fit@diagnostics$synthetic$share_residual, 1e-9)
    }
    drawn <- synthetic_market(demand = "LogitALM", supply = "bertrand",
        n_firms = 2, n_products = 2, shares = conditional,
        costs = costs, reference_margin = .2, seed = 9)
    expect_gt(drawn@diagnostics$synthetic$passive_outside_share, .1)
    expect_lt(drawn@diagnostics$synthetic$passive_outside_share, .5)
    expect_error(synthetic_market(demand = "LogitALM",
        n_firms = 2, n_products = 2, shares = conditional,
        costs = costs, reference_margin = .2,
        passive_outside_share = 0, seed = 9), "passive outside share")
})

test_that("ALM Cournot auction and bargaining obey their distinct margin equations", {
    conditional <- c(.2, .1, .15, .25, .3)
    costs <- c(40, 50, 60, 70, 80)
    s0 <- .3
    q <- conditional * (1 - s0)
    for (demand in c("LogitALM", "CESALM")) {
        for (conduct in c("cournot", "auction2nd", "bargaining")) {
            fit <- synthetic_market(demand = demand, supply = conduct,
                n_firms = 2, n_products = 2, shares = conditional,
                costs = costs, reference_margin = .2,
                passive_outside_share = s0, seed = 7)
            market <- fit@diagnostics$synthetic_market
            firm <- as.vector(market$ownership %*% q)
            price <- market$prices
            observed <- (price - costs) / price
            if (fit@spec$demand == "logit") {
                alpha <- fit@model@slopes$alpha
                if (conduct == "cournot") {
                    implied <- -(1 + firm / s0) / alpha / price
                } else if (conduct == "auction2nd") {
                    implied <- log1p(-firm) / (alpha * firm) / price
                } else {
                    M <- -market$ownership * q
                    diag(M) <- diag(market$ownership) + diag(M)
                    power <- fit@model@bargpowerPre
                    direct <- -log1p(-q) / ((power / (1 - power)) *
                        q / (1 - q) - log1p(-q))
                    implied <- as.vector(solve(t(M), direct)) / (-alpha * price)
                }
            } else {
                gamma <- fit@model@slopes$gamma
                if (conduct == "cournot") {
                    a <- s0 / (1 - s0)
                    implied <- 1 / gamma + (gamma - 1) * (1 + a) /
                        (gamma * (1 + gamma * a)) * firm
                } else if (conduct == "auction2nd") {
                    implied <- 1 - (1 - firm)^(1 / (gamma - 1))
                } else {
                    implied <- (1 - fit@model@bargpowerPre) /
                        (gamma * (1 - firm) + firm)
                }
            }
            expect_equal(observed, implied, tolerance = 1e-8,
                         info = paste(demand, conduct))
            expect_equal(market$costs, costs)
            expect_equal(unname(market$ownership),
                         unname(fit@model@ownerPre))
            expect_equal(unname(calcShares(fit@model, TRUE,
                revenue = fit@spec$demand == "ces")), q,
                tolerance = 1e-8)
            expect_equal(observed[length(observed)], .2, tolerance = 1e-8)
            expect_lt(fit@diagnostics$synthetic$foc_residual, 1e-6)
            if (conduct == "bargaining") {
                expect_equal(fit@model@bargpowerPre, rep(.5, length(q)))
            }
        }
    }
})

test_that("nested Logit ALM uses a passive outside option and supplied nests", {
    conditional <- c(.2, .1, .15, .25, .3)
    costs <- c(40, 50, 60, 70, 80)
    nests <- c("a", "a", "b", "b", "ref")
    sigma <- c(a = .7, b = .8, ref = 1)
    fit <- synthetic_market(demand = "LogitNestsALM",
        supply = "bertrand", n_firms = 2, n_products = 2,
        shares = conditional, costs = costs, reference_margin = .2,
        passive_outside_share = .3, nests = nests,
        parameters = list(sigma = sigma), seed = 7)
    market <- fit@diagnostics$synthetic_market
    prices <- market$prices
    actual <- conditional * .7
    D <- numDeriv::jacobian(function(candidate) {
        changed <- fit@model
        changed@pricePre <- candidate
        as.numeric(calcShares(changed, TRUE))
    }, prices)
    foc <- actual + t(market$ownership * D) %*% (prices - costs)
    expect_lt(max(abs(foc)), 1e-6)
    expect_equal(unname(calcShares(fit@model, TRUE)), actual,
                 tolerance = 1e-9)
    expect_equal(unname(fit@model@slopes$sigma[names(sigma)]),
                 unname(sigma), tolerance = 1e-7)
    expect_equal(market$costs, costs)
    expect_equal((prices[5] - costs[5]) / prices[5], .2,
                 tolerance = 1e-9)
    expect_error(synthetic_market(demand = "LogitNestsALM",
        n_firms = 2, n_products = 2, shares = conditional,
        costs = costs, reference_margin = .2,
        passive_outside_share = .3, seed = 7),
        "requires all-product 'nests'")
})

test_that("inadmissible and underidentified observed designs fail explicitly", {
    expect_error(synthetic_market(costs = c(1, -1), seed = 1), "costs")
    expect_error(synthetic_market(reference_margin = 1, seed = 1),
                 "reference_margin")
    expect_error(synthetic_market(reference_price = 100, seed = 1),
                 "observed mode uses")
    expect_error(synthetic_market(demand = "blp", supply = "bertrand", seed = 1),
                 "heterogeneity dispersion and integration draws")
    bad_shape <- matrix(1, 5, 5); diag(bad_shape) <- -1
    expect_error(synthetic_market(demand = "linear", supply = "bertrand",
        n_firms = 2, n_products = 2, shares = c(.2, .1, .15, .25, .3),
        costs = c(40, 50, 60, 70, 80), reference_margin = .2,
        parameters = list(slopes = bad_shape), seed = 1),
        "singular or ill-conditioned")
    expect_error(synthetic_market(demand = "ces", supply = "cournot",
        n_firms = 2, shares = c(.2, .2, .6), costs = c(50, 60, 70),
        reference_margin = .2, seed = 1), "margin above")
})

test_that("primitives mode still passes supplied parameters", {
    fit <- synthetic_market(mode = "primitives", n_firms = 2,
        parameters = list(alpha = -.05), reference_price = 100, seed = 7)
    expect_equal(fit@diagnostics$route, "specify")
    expect_equal(unname(fit@parameters$alpha), -.05)
    expect_equal(unname(calcShares(fit@model, TRUE)),
                 fit@diagnostics$synthetic_market$shares, tolerance = 1e-10)
})
