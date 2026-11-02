## Migration parity A: standalone Logit realization has independent cost-first FOCs.

test_that("realize_market solves prices from observed costs and margin", {
    shares <- c(.2, .1, .15, .25, .3)
    costs <- c(40, 50, 60, 70, 80)
    market <- fake_market(n_firms = 2, n_products = 2,
        shares = shares, costs = costs, reference_margin = .2, seed = 7)
    realized <- realize_market(market, model_spec("logit", "bertrand"))
    D <- diag(shares) - tcrossprod(shares)
    z <- solve(t(market$ownership * D), shares)
    alpha <- -z[5] / (costs[5] * .2 / .8)
    expect_equal(realized$prices, unname(costs - z / alpha),
                 tolerance = 1e-10)
    expect_equal(realized$costs, costs)
    expect_equal(realized$observed$costs, costs)
    expect_lt(realized$diagnostics$foc_residual, 1e-10)
    expect_lt(realized$diagnostics$share_residual, 1e-10)
    expect_lt(abs(realized$diagnostics$reference_margin_residual), 1e-10)
})

test_that("realize_market rejects unidentified models and conflicting slope", {
    market <- fake_market(n_firms = 2, costs = c(50, 60, 70),
        reference_margin = .2, seed = 7)
    expect_error(realize_market(market, model_spec("ces", "bertrand")),
                 "standard Bertrand Logit only")
    expect_error(realize_market(market, model_spec("logit", "bertrand"),
                                alpha = -.1), "omit 'alpha'")
})
