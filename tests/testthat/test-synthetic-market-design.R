test_that("Dirichlet shares are product-level and preserve totals", {
    expect_error(
        fake_market(n_firms = 3, n_products = 2, dirichlet_alpha = c(1, 2)),
        "dirichlet_alpha.*length 6"
    )
    expect_error(
        fake_market(n_firms = 3, n_products = 2,
                    dirichlet_alpha = c(1, 2, 3, 4, 5, 0)),
        "dirichlet_alpha"
    )
    expect_error(
        fake_market(n_firms = 3, n_products = 2,
                    dirichlet_alpha = c(1, 2, 3, 4, 5, Inf)),
        "dirichlet_alpha"
    )

    market <- fake_market(
        n_firms = 3, n_products = 2,
        dirichlet_alpha = c(1, 2, 3, 4, 5, 8), outside_beta = c(2, 3),
        seed = 10
    )
    expect_equal(length(market$design$dirichlet_alpha), 6)
    expect_equal(sum(market$design$relative_product_shares), 1,
                 tolerance = 1e-14)
    expect_equal(sum(market$design$firm_shares),
                 1 - market$design$outside_share, tolerance = 1e-14)
    expect_equal(sum(market$shares), 1, tolerance = 1e-14)
    expect_equal(sum(market$products$product_share), 1, tolerance = 1e-14)
    expect_true(market$design$outside_share > 0)
    expect_true(market$design$outside_share < 1)
    for (firm in seq_len(market$design$n_firms)) {
        idx <- market$products$firm_id == firm
        expect_equal(sum(market$products$product_share[idx]),
                     market$firms$firm_share[firm], tolerance = 1e-14)
    }
})

test_that("symmetric and asymmetric product-level designs are reproducible", {
    symmetric_a <- fake_market(
        n_firms = 2, n_products = 2, dirichlet_alpha = rep(2, 4), seed = 123
    )
    symmetric_b <- fake_market(
        n_firms = 2, n_products = 2, dirichlet_alpha = rep(2, 4), seed = 123
    )
    asymmetric <- fake_market(
        n_firms = 2, n_products = 2, dirichlet_alpha = c(1, 2, 3, 8), seed = 123
    )
    different_seed <- fake_market(
        n_firms = 2, n_products = 2, dirichlet_alpha = rep(2, 4), seed = 124
    )
    expect_equal(symmetric_a$shares, symmetric_b$shares)
    expect_equal(symmetric_a$prices, symmetric_b$prices)
    expect_equal(symmetric_a$observed$reference_markup,
                 symmetric_b$observed$reference_markup)
    expect_false(isTRUE(all.equal(symmetric_a$shares, different_seed$shares)))
    expect_false(isTRUE(all.equal(symmetric_a$shares, asymmetric$shares)))
})

test_that("n_products supports common and heterogeneous product counts", {
    expect_equal(fake_market(n_firms = 2, seed = 1)$design$n_products, 1)
    expect_error(fake_market(n_products = 0), "n_products")
    expect_error(fake_market(n_products = 1.5), "n_products")
    expect_error(fake_market(n_firms = 3, n_products = c(1, 2)), "n_products")
    market <- fake_market(n_firms = 3, n_products = 2, seed = 1)
    expect_equal(market$design$n_products, 2)
    expect_equal(market$design$products_per_firm, c(2, 2, 2))
    expect_equal(market$design$n_inside_products, 6)
    expect_equal(market$design$n_total_products, 7)
    expect_equal(market$firms$n_products, c(2, 2, 2, 1))
    expect_equal(nrow(market$products), 7)
    heterogeneous <- fake_market(n_firms = 3, n_products = c(1, 2, 3), seed = 1)
    expect_equal(heterogeneous$design$n_products, c(1, 2, 3))
    expect_equal(heterogeneous$design$products_per_firm, c(1, 2, 3))
    expect_equal(heterogeneous$design$n_inside_products, 6)
    expect_equal(heterogeneous$firms$n_products, c(1, 2, 3, 1))
    expect_equal(heterogeneous$products$firm_id, c(1, 2, 2, 3, 3, 3, 4))
    expect_equal(dim(heterogeneous$ownership), c(7, 7))
})

test_that("outside share Beta parameters are explicit and validated", {
    expect_error(fake_market(outside_beta = c(1)), "outside_beta.*length 2")
    expect_error(fake_market(outside_beta = c(1, 0)), "outside_beta")
    expect_error(fake_market(outside_beta = c(1, NA)), "outside_beta")
    a <- fake_market(outside_beta = c(1, 9), seed = 33)
    b <- fake_market(outside_beta = c(1, 9), seed = 33)
    expect_equal(a$design$outside_share, b$design$outside_share)
})

test_that("product-level shares aggregate to firms and define ownership", {
    market <- fake_market(n_firms = 2, n_products = 2, seed = 4)
    expect_equal(market$products$firm_id, c(1, 1, 2, 2, 3))
    expect_equal(market$design$ownership_map$firm_id, market$products$firm_id)
    expect_equal(unname(market$ownership[1:2, 1:2]), matrix(1, 2, 2))
    expect_equal(unname(market$ownership[1:2, 3:4]), matrix(0, 2, 2))
    expect_equal(unname(market$ownership[5, 1:4]), rep(0, 4))
    for (firm in 1:2) {
        idx <- market$products$firm_id == firm
        expect_equal(sum(market$products$product_share[idx]),
                     market$firms$firm_share[firm], tolerance = 1e-14)
        expect_equal(market$products$firm_share[idx],
                     rep(market$firms$firm_share[firm], 2))
    }
})

test_that("reference product is active, priced, and separately owned", {
    market <- fake_market(n_firms = 2, n_products = 2, price_level = 25, seed = 4)
    ref <- market$design$reference_product
    expect_equal(market$products$reference_product,
                 seq_len(nrow(market$products)) == ref)
    expect_equal(market$products$price[ref], 25)
    expect_equal(market$products$firm_id[ref], market$design$reference_firm)
    expect_equal(market$reference_product, ref)
    expect_equal(market$reference_share, market$shares[ref])
    expect_equal(market$reference_price, market$prices[ref])
    expect_equal(unname(market$ownership[ref, -ref]),
                 rep(0, nrow(market$products) - 1))
    expect_equal(market$ownership[ref, ref], 1)
    expect_equal(market$observed$reference_share, market$shares[ref])
    expect_equal(market$metadata$reference_normalization,
                 "mean utility only; reference price is real")
})

test_that("prices and open-boundary outside-margin conventions are explicit", {
    supplied <- fake_market(n_firms = 2, n_products = 2,
                            prices = c(10, 20, 30, 40), reference_price = 50,
                            seed = 5)
    expect_equal(supplied$prices, c(10, 20, 30, 40, 50))
    expect_equal(supplied$design$price_rule, "user-supplied-inside")
    all_prices <- fake_market(n_firms = 2, n_products = 2,
                              prices = c(10, 20, 30, 40, 50), seed = 5)
    expect_equal(all_prices$design$price_rule, "user-supplied")
    expect_error(fake_market(n_firms = 2, n_products = 2,
                             prices = c(0, 2, 3, 4, 5)), "prices")
    market <- fake_market(n_firms = 2, price_level = 100, seed = 5)
    markup <- market$observed$reference_markup
    expect_true(markup > 0 && markup < 100)
    expect_equal(market$metadata$units$markup, "price level")
    expect_equal(market$design$markup_rule, "uniform-open-U(0,100)")
    expect_equal(market$observed$outside_margin, markup)
    expect_equal(market$design$outside_margin, markup)
    expect_equal(market$products$observed_markup[
        market$design$reference_product], markup)
    expect_equal(market$products$outside_margin[
        market$design$reference_product], markup)
    supplied_margin <- fake_market(n_firms = 2, outside_margin = 20, seed = 5)
    expect_equal(supplied_margin$observed$outside_margin, 20)
    expect_equal(supplied_margin$observed$reference_markup, 20)
    expect_error(
        fake_market(mode = "primitives", parameters = list(alpha = -1),
                    outside_margin = 20), "outside margin.*only valid"
    )
    expect_equal(fake_market(mode = "observed_information", seed = 5)$design$mode,
                 "observed")
    expect_equal(fake_market(mode = "known_primitives",
                             parameters = list(alpha = -1), seed = 5)$design$mode,
                 "primitives")
})

test_that("RNG state is preserved and Monte Carlo seeds are deterministic", {
    set.seed(99)
    expected <- runif(1)
    set.seed(99)
    invisible(fake_market(seed = 101))
    observed <- runif(1)
    expect_equal(observed, expected)
    first <- simulate_markets(3, fake_market, seed = 8,
                              n_firms = 2, dirichlet_alpha = c(1, 3))
    second <- simulate_markets(3, fake_market, seed = 8,
                               n_firms = 2, dirichlet_alpha = c(1, 3))
    expect_equal(first$seeds, second$seeds)
    expect_equal(first$diagnostics$rejection_rate, 0)
    expect_equal(lapply(first$markets, `[[`, "shares"),
                 lapply(second$markets, `[[`, "shares"))
    expect_equal(length(first$markets), 3)
    realized <- simulate_markets(
        2, fake_market, seed = 12, realizer = identity,
        n_firms = 2, dirichlet_alpha = c(1, 3)
    )
    expect_true(realized$diagnostics$realized)
    expect_equal(realized$diagnostics$n_rejected, 0)
})
