test_that("observed synthetic Logit markets use the full multi-product FOC", {
    fit <- synthetic_market(
        demand = "logit", supply = "bertrand", n_firms = 2,
        n_products = 2, reference_price = 100, outside_margin = 20,
        seed = 7
    )
    market <- fit@diagnostics$synthetic_market

    expect_s4_class(fit, "AntitrustFit")
    expect_equal(fit@spec$demand, "logit")
    expect_equal(fit@spec$conduct, "bertrand")
    expect_equal(market$design$reference_product, 5)
    expect_equal(unname(calcShares(fit@model, TRUE)), market$shares,
                 tolerance = 1e-10)
    expect_lt(abs(fit@diagnostics$synthetic$reference_markup_error), 1e-6)
    expect_lt(fit@diagnostics$synthetic$foc_residual, 1e-8)
    expect_equal(fit@diagnostics$synthetic$foc_rank, 5)
    expect_lt(fit@diagnostics$synthetic$foc_condition_number,
              fit@diagnostics$synthetic$foc_condition_limit)
    d_unit <- diag(market$shares) - tcrossprod(market$shares)
    z <- solve(t(market$ownership * d_unit), market$shares)
    alpha <- -z[5] / 20
    expect_equal(fit@diagnostics$synthetic$derivative, alpha * d_unit,
                 tolerance = 1e-12)
    expect_equal(fit@diagnostics$synthetic$retained_cost,
                 unname(market$prices + z / alpha), tolerance = 1e-12)
    expect_equal(unname(fit@model@pricePre), market$prices, tolerance = 1e-10)
    expect_true(all(is.finite(fit@model@mcPre)))
    expect_equal(market$products$firm_id, c(1, 1, 2, 2, 3))
    expect_equal(unname(market$ownership[1:2, 1:2]), matrix(1, 2, 2))
    expect_equal(unname(market$ownership[5, 1:4]), rep(0, 4))
})

test_that("one-product synthetic Logit reduces to the analytic special case", {
    fit <- synthetic_market(
        demand = "logit", supply = "bertrand", n_firms = 3,
        reference_price = 100, outside_margin = 20, seed = 12
    )
    share_ref <- fit@model@shares[length(fit@model@shares)]
    expected <- -1 / (20 * (1 - share_ref))
    expect_equal(unname(fit@model@slopes$alpha), expected, tolerance = 1e-8)
    expect_lt(fit@diagnostics$synthetic$foc_residual, 1e-8)
})

test_that("MonCom Logit diagnostics use the direct own-product FOC", {
    fit <- synthetic_market(
        demand = "logit", supply = "moncom", n_firms = 2,
        reference_price = 100, outside_margin = 10, seed = 20
    )

    expect_equal(fit@diagnostics$synthetic$foc_status, "verified")
    expect_lt(fit@diagnostics$synthetic$foc_residual, 1e-8)
    expect_true(isTRUE(fit@diagnostics$synthetic$equilibrium_check))
    expected <- diag(-.1 * fit@diagnostics$synthetic$quantity)
    expect_equal(fit@diagnostics$synthetic$derivative, expected,
                 tolerance = 1e-12)
})

test_that("MonCom CES diagnostics use the direct own-product FOC", {
    fit <- synthetic_market(
        demand = "ces", supply = "moncom", n_firms = 2,
        reference_price = 100, outside_margin = 10, seed = 20
    )

    expect_equal(fit@diagnostics$synthetic$foc_status, "verified")
    expect_lt(fit@diagnostics$synthetic$foc_residual, 1e-8)
    expect_true(isTRUE(fit@diagnostics$synthetic$equilibrium_check))
    expected <- diag(-10 * fit@diagnostics$synthetic$quantity /
                     fit@model@pricePre)
    expect_equal(fit@diagnostics$synthetic$derivative, expected,
                 tolerance = 1e-12)
})

test_that("unsupported synthetic diagnostics remain unavailable", {
    fit <- synthetic_market(
        demand = "ces", supply = "bertrand", n_firms = 2,
        reference_price = 100, outside_margin = 10, seed = 20
    )

    expect_equal(fit@diagnostics$synthetic$status, "unavailable")
    expect_equal(fit@diagnostics$synthetic$foc_status, "unavailable")
    expect_true(is.na(fit@diagnostics$synthetic$foc_residual))
    expect_true(is.na(fit@diagnostics$synthetic$equilibrium_check))
})

test_that("MonCom BLP diagnostics integrate supplied points independently", {
    fit <- synthetic_market(
        demand = "blp", supply = "moncom", mode = "primitives",
        n_firms = 2, n_products = 1, reference_price = 100,
        parameters = list(alphaMean = -0.05, sigma = 0.01,
                          draws = c(-1, 0, 1),
                          drawWeights = c(.2, .6, .2)), seed = 21
    )
    model <- fit@model
    prices <- unname(model@pricePre)
    delta <- unname(model@slopes$meanval)
    alphas <- -.05 + .01 * c(-1, 0, 1)
    weights <- c(.2, .6, .2)
    draw_shares <- matrix(NA_real_, nrow = length(prices), ncol = 3L)
    for (r in seq_len(3L)) {
        utility <- delta + alphas[r] *
            (prices - model@priceOutside)
        expu <- exp(utility - max(utility))
        draw_shares[, r] <- expu / sum(expu)
    }
    direct <- as.vector(draw_shares %*% (weights * alphas))
    expected <- diag(model@insideSize * direct)

    expect_equal(model@slopes$consDraws, c(-1, 0, 1), tolerance = 0)
    expect_equal(model@slopes$drawWeights, weights, tolerance = 0)
    expect_equal(fit@diagnostics$synthetic$derivative, expected,
                 tolerance = 1e-12)
    expect_lt(fit@diagnostics$synthetic$foc_residual, 1e-10)
    expect_true(isTRUE(fit@diagnostics$synthetic$equilibrium_check))
})

test_that("synthetic oracle fails perturbed price, cost, and Bertrand ownership", {
    fit <- synthetic_market(
        demand = "logit", supply = "bertrand", n_firms = 2,
        n_products = 2, reference_price = 100, outside_margin = 20, seed = 27
    )
    market <- fit@diagnostics$synthetic_market
    ref <- market$design$reference_product
    check <- function(candidate, candidate_market = market) {
        antitrust:::.antitrust_synthetic_diagnostics(
            candidate, candidate_market, market$observed$reference_markup,
            "observed", parameter_truth = list()
        )
    }

    price_bad <- fit
    price_bad@model@pricePre[1] <- price_bad@model@pricePre[1] + 1
    price_diag <- check(price_bad)
    expect_gt(price_diag$price_residual, 1e-8)
    expect_false(isTRUE(price_diag$equilibrium_check))

    market_price_bad <- market
    market_price_bad$prices[1] <- market_price_bad$prices[1] + 1
    market_price_diag <- check(fit, market_price_bad)
    expect_gt(market_price_diag$price_residual, 1e-8)
    expect_false(isTRUE(market_price_diag$equilibrium_check))

    cost_bad_market <- market
    cost_bad_market$products$cost[1] <- cost_bad_market$products$cost[1] + 1
    cost_diag <- check(fit, cost_bad_market)
    expect_gt(cost_diag$foc_residual, 1e-8)
    expect_false(isTRUE(cost_diag$equilibrium_check))

    owner_bad <- fit
    owner_bad@model@ownerPre[1, ref] <- 1
    owner_bad@model@ownerPre[ref, 1] <- 1
    owner_diag <- check(owner_bad)
    expect_gt(owner_diag$ownership_residual, 1e-8)
    expect_false(isTRUE(owner_diag$equilibrium_check))

    market_owner_bad <- market
    market_owner_bad$ownership[1, ref] <- 1
    market_owner_bad$ownership[ref, 1] <- 1
    market_owner_diag <- check(fit, market_owner_bad)
    expect_gt(market_owner_diag$market_ownership_residual, 1e-8)
    expect_false(isTRUE(market_owner_diag$equilibrium_check))
})

test_that("MonCom oracle remains valid when ownership changes", {
    fit <- synthetic_market(
        demand = "logit", supply = "moncom", n_firms = 2,
        n_products = 2, reference_price = 100, outside_margin = 20, seed = 28
    )
    changed <- fit
    changed@model@ownerPre[,] <- 1
    diagnostics <- antitrust:::.antitrust_synthetic_diagnostics(
        changed, fit@diagnostics$synthetic_market, 20, "observed"
    )
    expect_true(isTRUE(diagnostics$equilibrium_check))
    expect_true(is.na(diagnostics$ownership_residual))
})

test_that("synthetic Bertrand Logit reacts to zero and positive cost shocks", {
    fit <- synthetic_market(
        demand = "logit", supply = "bertrand", n_firms = 2,
        reference_price = 100, outside_margin = 20, seed = 29
    )
    n <- length(fit@model@pricePre)
    unchanged <- simulate(fit, counterfactual(costs = rep(0, n)))
    increased <- simulate(fit, counterfactual(costs = rep(.10, n)))

    expect_equal(unname(unchanged@pricePost), unname(fit@model@pricePre),
                 tolerance = 1e-8)
    expect_true(any(increased@pricePost > fit@model@pricePre + 1e-8))
})

test_that("synthetic antitrust markets preserve heterogeneous ownership", {
    fit <- synthetic_market(
        demand = "logit", supply = "bertrand", n_firms = 3,
        n_products = c(1, 2, 3), reference_price = 100,
        outside_margin = 20, seed = 13
    )
    market <- fit@diagnostics$synthetic_market
    expect_equal(market$design$products_per_firm, c(1, 2, 3))
    expect_equal(market$products$firm_id, c(1, 2, 2, 3, 3, 3, 4))
    expect_lt(fit@diagnostics$synthetic$foc_residual, 1e-8)
})

test_that("primitives mode uses specify and retains truth", {
    fit <- synthetic_market(
        demand = "logit", supply = "bertrand", mode = "primitives",
        n_firms = 2, n_products = 2, parameters = list(alpha = -0.05),
        reference_price = 100, seed = 7
    )
    market <- fit@diagnostics$synthetic_market

    expect_s4_class(fit, "AntitrustFit")
    expect_equal(fit@diagnostics$route, "specify")
    expect_equal(unname(fit@parameters$alpha), -0.05)
    expect_equal(unname(fit@diagnostics$synthetic_truth$alpha), -0.05)
    expect_equal(unname(calcShares(fit@model, TRUE)), market$shares,
                 tolerance = 1e-10)
    expect_equal(fit@diagnostics$synthetic_parameter_error$alpha, 0,
                 tolerance = 1e-12)
})

test_that("known primitives can be hidden and recovered from the reference markup", {
    qa_skip_unless_tier("extended")
    known <- synthetic_market(
        demand = "logit", supply = "bertrand", mode = "primitives",
        n_firms = 2, parameters = list(alpha = -0.05),
        outside_beta = c(2, 20), reference_price = 100, seed = 91
    )
    ref <- length(known@model@shares)
    reference_markup <- calcMargins(known@model, TRUE, level = TRUE)[ref]
    recovered <- synthetic_market(
        demand = "logit", supply = "bertrand", mode = "observed",
        n_firms = 2, outside_beta = c(2, 20), reference_price = 100,
        outside_margin = reference_markup, seed = 91
    )

    expect_equal(unname(recovered@parameters$alpha), -0.05, tolerance = 1e-7)
    expect_equal(unname(recovered@diagnostics$synthetic_recovered$alpha),
                 -0.05, tolerance = 1e-7)
    expect_lt(abs(recovered@diagnostics$synthetic$foc_residual), 1e-8)
})

test_that("synthetic fits plug into the counterfactual workflow", {
    fit <- synthetic_market(
        demand = "logit", supply = "bertrand", n_firms = 2,
        reference_price = 100, outside_margin = 20, seed = 7
    )
    before <- fit@diagnostics$synthetic_market$shares
    owner_post <- fit@model@ownerPre
    owner_post[1, 2] <- owner_post[2, 1] <- 1

    result <- simulate(fit, counterfactual(ownership = owner_post))
    expect_s4_class(result, "Logit")
    expect_equal(fit@diagnostics$synthetic_market$shares, before)
})
