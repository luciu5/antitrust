test_that("multi-product Bertrand Logit realizes a known-primitives market", {
  skip_if_not_installed("iopolicy")
  market <- iopolicy::fake_market(
    mode = "primitives", n_firms = 2, dirichlet_alpha = c(2, 3),
    outside_beta = c(3, 5), products_per_firm = c(2, 1),
    price_level = 100, parameters = list(alpha = -1), seed = 101
  )
  source_products <- market$products
  realized <- iopolicy::realize_market(market, model_spec("logit", "bertrand"))

  expect_s3_class(realized, "SyntheticMarket")
  expect_equal(realized$diagnostics$equilibrium_status, "realized")
  expect_lt(realized$diagnostics$foc_residual, 1e-8)
  expect_equal(realized$diagnostics$recovered_parameters$alpha, -1,
               tolerance = 1e-12)
  expect_equal(realized$products$cost,
               realized$prices - realized$products$markup, tolerance = 1e-12)
  expect_equal(realized$products$mean_value[realized$design$reference_product], 0,
               tolerance = 1e-14)
  expect_equal(market$products, source_products)
  expect_true(all(is.na(market$products$cost)))
})

test_that("observed reference markup identifies alpha with the full ownership system", {
  skip_if_not_installed("iopolicy")
  market <- iopolicy::fake_market(
    mode = "observed", n_firms = 3, dirichlet_alpha = rep(1, 3),
    outside_beta = c(4, 4), products_per_firm = 1, price_level = 100,
    observed_markup = 20, seed = 202
  )
  realized <- iopolicy::realize_market(market, model_spec("logit", "bertrand"))
  ref <- market$design$reference_product
  expected_alpha <- -1 / (market$observed$reference_markup * (1 - market$shares[ref]))

  expect_equal(realized$diagnostics$recovered_parameters$alpha, expected_alpha,
               tolerance = 1e-12)
  expect_lt(realized$diagnostics$foc_residual, 1e-8)
  expect_equal(realized$products$markup[ref], market$observed$reference_markup,
               tolerance = 1e-12)
})

test_that("full multi-product FOCs differ from the one-product shortcut", {
  skip_if_not_installed("iopolicy")
  market <- iopolicy::fake_market(
    mode = "observed", n_firms = 2, dirichlet_alpha = c(3, 2),
    outside_beta = c(3, 7), products_per_firm = c(2, 2),
    price_level = 100, observed_markup = 20, seed = 303
  )
  realized <- iopolicy::realize_market(market, model_spec("logit", "bertrand"))
  shares <- market$shares
  ownership <- market$ownership
  d_unit <- diag(shares) - tcrossprod(shares, shares)
  z <- solve(t(ownership * d_unit), shares)
  expected_markup <- unname(-z / realized$diagnostics$recovered_parameters$alpha)

  expect_equal(realized$products$markup, expected_markup, tolerance = 1e-12)
  expect_lt(realized$diagnostics$foc_residual, 1e-8)
  shortcut <- -1 / (realized$products$markup[1] * (1 - shares[1]))
  expect_false(isTRUE(all.equal(shortcut,
                                realized$diagnostics$recovered_parameters$alpha,
                                tolerance = 1e-6)))
})

test_that("known-primitive QA records parameter error and mean-share recovery", {
  skip_if_not_installed("iopolicy")
  market <- iopolicy::fake_market(
    mode = "primitives", n_firms = 2, dirichlet_alpha = c(1, 4),
    products_per_firm = c(2, 1), price_level = 100,
    parameters = list(alpha = -0.75), seed = 404
  )
  realized <- iopolicy::realize_market(market, model_spec("logit", "bertrand"))
  expect_equal(realized$diagnostics$parameter_error$alpha, 0,
               tolerance = 1e-12)
  expect_lt(realized$diagnostics$share_residual, 1e-12)
  expect_equal(sum(realized$shares), 1, tolerance = 1e-14)
})
