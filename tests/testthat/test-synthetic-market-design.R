## Migration parity A: design inputs must not become residual outputs.

test_that("observed design retains complete heterogeneous costs and ownership", {
    costs <- c(10, 20, 30, 40, 50, 60, 70)
    shares <- c(.08, .11, .12, .14, .18, .17, .20)
    market <- fake_market(n_firms = 3, n_products = c(1, 2, 3),
        costs = costs, shares = shares, reference_margin = .25, seed = 7)
    expect_equal(market$observed$costs, costs)
    expect_equal(market$products$cost, costs)
    expect_equal(market$shares, shares)
    expect_true(all(is.na(market$prices)))
    expect_equal(market$products$firm_id, c(1, 2, 2, 3, 3, 3, 4))
    expect_equal(market$ownership[2, 3], 1)
    expect_equal(market$ownership[2, 4], 0)
    expect_equal(market$observed$reference_margin, .25)
})

test_that("common and bounded heterogeneous cost rules are reproducible", {
    common <- fake_market(n_firms = 2, cost_level = 42, seed = 11)
    a <- fake_market(n_firms = 2, cost_rule = "uniform",
        cost_range = c(15, 25), seed = 19)
    b <- fake_market(n_firms = 2, cost_rule = "uniform",
        cost_range = c(15, 25), seed = 19)
    expect_equal(common$costs, rep(42, 3))
    expect_true(all(a$costs >= 15 & a$costs <= 25))
    expect_equal(a$costs, b$costs)
    expect_equal(a$shares, b$shares)
})

test_that("invalid cost, share, and margin designs are rejected", {
    expect_error(fake_market(n_firms = 2, costs = c(1, 2)), "costs")
    expect_error(fake_market(n_firms = 2, costs = c(1, -2, 3)), "costs")
    expect_error(fake_market(n_firms = 2, costs = c(1, Inf, 3)), "costs")
    expect_error(fake_market(n_firms = 2, shares = c(.2, .3, .4)), "shares")
    expect_error(fake_market(reference_margin = 0), "reference_margin")
    expect_error(fake_market(reference_margin = 1), "reference_margin")
    expect_error(fake_market(prices = c(1, 2, 3, 4)), "price-first")
})

test_that("primitives design still uses its complete price contract", {
    market <- fake_market(mode = "primitives", n_firms = 2,
        parameters = list(alpha = -.1), prices = c(10, 20, 30),
        reference_price = 30, seed = 7)
    expect_equal(market$prices, c(10, 20, 30))
    expect_true(all(is.na(market$costs)))
    expect_equal(market$truth$alpha, -.1)
})

test_that("batch rejection records failed economic designs", {
    batch <- simulate_markets(2, seed = 7, max_attempts = 2,
        generator = function(seed) fake_market(n_firms = 2,
            costs = c(1, -1, 3), seed = seed))
    expect_equal(batch$diagnostics$n_rejected, 2)
    expect_true(all(vapply(batch$diagnostics$rejection_reasons,
        length, integer(1)) == 2))
})
