moncom_blp_specify <- function(integration, sigma = .3, draws = NULL,
                               drawWeights = NULL, nNodes = NULL,
                               nDraws = NULL) {
    parameters <- list(
        alphaMean = -1.2,
        sigma = sigma,
        meanval = c(.5, .2, -.1),
        integration = integration
    )
    if (!is.null(draws)) parameters$draws <- draws
    if (!is.null(drawWeights)) parameters$drawWeights <- drawWeights
    if (!is.null(nNodes)) parameters$nNodes <- nNodes
    if (!is.null(nDraws)) parameters$nDraws <- nDraws
    specify(
        "blp", "moncom",
        prices = c(1.5, 1.8, 2.1),
        shares = c(.30, .25, .15),
        parameters = parameters,
        ownerPre = c("A", "A", "B"),
        insideSize = 100,
        s0 = .30
    )
}


test_that("MonCom BLP uses the integrated draw-level own derivative", {
    fit <- moncom_blp_specify(
        integration = "provided",
        draws = c(-1, 0, 1),
        drawWeights = c(.2, .6, .2)
    )
    model <- fit@model
    shares_draw <- calcShares(model, aggregate = FALSE)
    weights <- model@slopes$drawWeights
    alpha <- model@slopes$alphas
    direct <- as.vector(shares_draw %*% (weights * alpha))
    aggregate <- calcShares(model)
    expected <- -aggregate / direct

    expect_s4_class(model, "MonComBLP")
    expect_equal(unname(calcMargins(model, level = TRUE)),
                 unname(expected), tolerance = 1e-12)
    expect_lt(fit@diagnostics$foc_residual_pre, 1e-12)
    expect_true(isTRUE(fit@diagnostics$ownership_irrelevant))
})


test_that("MonCom BLP calibration uses the integrated own derivative", {
    known <- moncom_blp_specify(
        integration = "provided", sigma = .2,
        draws = c(-1, 0, 1), drawWeights = c(.2, .6, .2)
    )
    shares <- calcShares(known@model)
    margins <- calcMargins(known@model)
    fit <- suppressWarnings(calibrate(
        "blp", "moncom", prices = known@model@pricePre,
        shares = shares, margins = margins,
        ownerPre = c("A", "A", "B"), s0 = 1 - sum(shares),
        integration = "provided", draws = c(-1, 0, 1),
        drawWeights = c(.2, .6, .2)
    ))
    expect_s4_class(fit@model, "MonComBLP")
    expect_equal(unname(fit@parameters$alphaMean), -1.2,
                 tolerance = 1e-4)
    expect_equal(unname(fit@parameters$sigma), .2,
                 tolerance = 1e-4)
    expect_lt(fit@diagnostics$foc_residual_pre, 1e-8)
})


test_that("MonCom BLP preserves provided, Gauss-Hermite, and Monte Carlo integration", {
    set.seed(20260907)
    fits <- list(
        provided = moncom_blp_specify(
            integration = "provided", draws = c(-1, 0, 1),
            drawWeights = c(.2, .6, .2)
        ),
        gauss_hermite = moncom_blp_specify(
            integration = "gauss-hermite", nNodes = 11L
        ),
        monte_carlo = moncom_blp_specify(
            integration = "monte-carlo", nDraws = 17L
        )
    )

    expect_identical(fits$provided@model@slopes$integration, "provided")
    expect_identical(fits$gauss_hermite@model@slopes$integration, "gauss-hermite")
    expect_identical(fits$monte_carlo@model@slopes$integration, "monte-carlo")
    expect_length(fits$gauss_hermite@model@slopes$consDraws, 11L)
    expect_length(fits$monte_carlo@model@slopes$consDraws, 17L)
    for (fit in fits) {
        expect_true(all(is.finite(fit@model@pricePost)))
        expect_lt(fit@diagnostics$foc_residual_pre, 1e-10)
    }
})


test_that("MonCom BLP sigma zero has the homogeneous Logit limit", {
    fit <- moncom_blp_specify(
        integration = "provided", sigma = 0,
        draws = c(-1, 0, 1), drawWeights = c(.2, .6, .2)
    )
    expect_equal(unname(calcMargins(fit@model, level = TRUE)),
                 rep(1 / 1.2, 3), tolerance = 1e-12)
})


test_that("MonCom BLP ownership is irrelevant but cost shocks are not", {
    fit <- moncom_blp_specify(
        integration = "provided",
        draws = c(-1, 0, 1), drawWeights = c(.2, .6, .2)
    )
    merger <- simulate(fit, counterfactual(ownership = c("A", "A", "A")))
    cost <- simulate(fit, counterfactual(costs = c(-.1, 0, 0)))
    expect_equal(unname(merger@pricePost), unname(fit@model@pricePre),
                 tolerance = 1e-9)
    expect_true(abs(cost@pricePost[1] - fit@model@pricePre[1]) > 1e-6)
})
