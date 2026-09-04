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
