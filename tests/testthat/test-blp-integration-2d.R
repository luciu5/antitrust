# Focused economic tests for two-dimensional BLP integration.  Expected
# moments and shares are formed directly from the normal-product rule rather
# than by calling a second package implementation.

test_that("two-dimensional Gauss-Hermite is a correctly ordered product rule", {
    set.seed(20260909)
    before <- .Random.seed
    rule <- antitrust:::.blp_integration(list(
        integration = "gauss-hermite", sigma = .2,
        nDemog = 1L, piDemog = .3
    ))
    after <- .Random.seed
    one <- antitrust:::.blp_normal_nodes(15L)

    expect_identical(after, before)
    expect_identical(rule$factorOrder, c("price", "demog1"))
    expect_identical(rule$nodesPerAxis, c(15L, 15L))
    expect_equal(dim(rule$integrationPoints), c(225L, 2L))
    expect_equal(rule$integrationPoints[, 1], rep(one$nodes, 15L), tolerance = 0)
    expect_equal(rule$integrationPoints[, 2], rep(one$nodes, each = 15L), tolerance = 0)
    expect_equal(rule$weights, as.vector(outer(one$weights, one$weights)),
                 tolerance = 0)
    expect_equal(sum(rule$weights), 1, tolerance = 1e-14)
    expect_equal(as.vector(crossprod(rule$weights, rule$integrationPoints)),
                 c(0, 0), tolerance = 1e-13)
    second <- crossprod(
        sweep(rule$integrationPoints, 1L, sqrt(rule$weights), "*")
    )
    expect_equal(unname(second), diag(2), tolerance = 1e-12)
})


test_that("two-dimensional node counts are per axis", {
    square <- antitrust:::.blp_integration(list(
        integration = "gauss-hermite", nNodes = 5L,
        sigma = .2, nDemog = 1L, piDemog = .3
    ))
    rectangle <- antitrust:::.blp_integration(list(
        integration = "gauss-hermite", nNodes = c(5L, 7L),
        sigma = .2, nDemog = 1L, piDemog = .3
    ))

    expect_identical(square$nodesPerAxis, c(5L, 5L))
    expect_equal(dim(square$integrationPoints), c(25L, 2L))
    expect_identical(rectangle$nodesPerAxis, c(5L, 7L))
    expect_equal(dim(rectangle$integrationPoints), c(35L, 2L))
})


test_that("provided two-dimensional points have an unambiguous API", {
    points <- matrix(c(-1, 0, 1, 2, -2, .5), ncol = 2L)
    weights <- c(1, 2, 3)
    rule <- antitrust:::.blp_integration(list(
        integration = "provided", integrationPoints = points,
        integrationWeights = weights,
        sigma = .2, nDemog = 1L, piDemog = .3
    ))

    expect_identical(rule$integrationPoints, points)
    expect_equal(rule$weights, weights / sum(weights), tolerance = 0)
    public_rule <- calcBLPintegration(list(
        integration = "provided", integrationPoints = points,
        integrationWeights = weights,
        sigma = .2, nDemog = 1L, piDemog = .3
    ))
    expect_identical(public_rule$integrationPoints, points)
    expect_equal(public_rule$weights, weights / sum(weights), tolerance = 0)
    expect_error(
        antitrust:::.blp_integration(list(
            draws = points, sigma = .2, nDemog = 1L, piDemog = .3
        )),
        "use 'integrationPoints'"
    )
    expect_error(
        antitrust:::.blp_integration(list(
            integrationPoints = points[, 1, drop = FALSE], sigma = .2
        )),
        "at least two columns"
    )
    expect_error(
        antitrust:::.blp_integration(list(
            integrationPoints = points, nDraws = 4L,
            sigma = .2, nDemog = 1L, piDemog = .3
        )),
        "must equal the number"
    )
})


test_that("price and demographic draws materialize from one tensor state", {
    integration <- antitrust:::.blp_integration(list(
        integration = "gauss-hermite", nNodes = c(5L, 7L),
        sigma = .25, nDemog = 1L, piDemog = .4
    ))
    materialized <- antitrust:::.blp_materialize_draws(
        integration, alphaMean = -2, sigma = .25,
        nDemog = 1L, piDemog = .4,
        demogMean = .3, demogCov = matrix(2.25, 1L, 1L)
    )
    points <- integration$integrationPoints
    expected_demog <- .3 + 1.5 * points[, 2]
    expected_alpha <- -2 + .25 * points[, 1] +
        .4 * (expected_demog - .3)

    expect_equal(materialized$priceDraws, points[, 1], tolerance = 0)
    expect_equal(materialized$demogDraws[, 1], expected_demog, tolerance = 1e-14)
    expect_equal(materialized$alphas, expected_alpha, tolerance = 1e-14)
})


test_that("correlated demographic quadrature uses the full Cholesky block", {
    covariance <- matrix(c(4, 1.2, 1.2, 2), 2L, 2L)
    means <- c(.5, -1)
    integration <- antitrust:::.blp_integration(list(
        integration = "gauss-hermite", nNodes = c(7L, 7L),
        sigma = 0, nDemog = 2L, piDemog = c(.2, 0)
    ))
    materialized <- antitrust:::.blp_materialize_draws(
        integration, alphaMean = -2, sigma = 0,
        nDemog = 2L, piDemog = c(.2, 0),
        demogMean = means, demogCov = covariance
    )
    centered <- sweep(materialized$demogDraws, 2L, means, "-")
    weighted_mean <- as.vector(crossprod(integration$weights,
                                         materialized$demogDraws))
    weighted_cov <- crossprod(sweep(centered, 1L,
                                    sqrt(integration$weights), "*"))

    expect_identical(integration$factorOrder, c("demog1", "demog2"))
    expect_equal(weighted_mean, means, tolerance = 1e-12)
    expect_equal(unname(weighted_cov), covariance, tolerance = 1e-11)
})


test_that("price and characteristic heterogeneity enter utility independently", {
    prod_char <- matrix(c(1, 2, -1), ncol = 1L)
    integration <- antitrust:::.blp_integration(list(
        integration = "gauss-hermite", nNodes = c(5L, 7L),
        sigma = .2, sigmaChar = .4, prodChar = prod_char,
        nDemog = 0L
    ))
    materialized <- antitrust:::.blp_materialize_draws(
        integration, alphaMean = -1.5, sigma = .2,
        prodChar = prod_char, sigmaChar = .4
    )
    points <- integration$integrationPoints

    expect_identical(integration$factorOrder, c("price", "char1"))
    expect_equal(materialized$alphas, -1.5 + .2 * points[, 1],
                 tolerance = 1e-14)
    expect_equal(materialized$char_random,
                 outer(.4 * points[, 2], prod_char[, 1]),
                 tolerance = 1e-14)
})


test_that("higher dimensions fall back to Monte Carlo but explicit GH fails", {
    set.seed(20260909)
    rule <- antitrust:::.blp_integration(list(
        integration = "auto", nDraws = 7L,
        sigma = .2, nDemog = 2L, piDemog = c(.1, .2)
    ))
    expect_identical(rule$rule, "monte-carlo")
    expect_equal(dim(rule$integrationPoints), c(7L, 3L))
    set.seed(20260909)
    repeated <- antitrust:::.blp_integration(list(
        integration = "auto", nDraws = 7L,
        sigma = .2, nDemog = 2L, piDemog = c(.1, .2)
    ))
    expect_identical(repeated$integrationPoints, rule$integrationPoints)
    reused <- antitrust:::.blp_integration(list(
        integration = "provided", integrationPoints = rule$integrationPoints,
        integrationWeights = rule$weights,
        sigma = .2, nDemog = 2L, piDemog = c(.1, .2)
    ))
    expect_identical(reused$integrationPoints, rule$integrationPoints)
    expect_error(
        antitrust:::.blp_integration(list(
            integration = "gauss-hermite", nNodes = 5L,
            sigma = .2, nDemog = 2L, piDemog = c(.1, .2)
        )),
        "at most two active dimensions"
    )
})


test_that("a specified two-dimensional BLP stores and reuses exact state", {
    fit <- suppressMessages(suppressWarnings(specify(
        demand = "blp", conduct = "bertrand",
        prices = c(1.5, 1.8, 2.1), shares = c(.30, .25, .15),
        ownerPre = c("A", "B", "C"), insideSize = 100,
        parameters = list(
            alpha = -1.5, sigma = .2, meanval = c(.5, .2, -.1),
            piDemog = .3, demogMean = .2,
            demogCov = matrix(.64, 1L, 1L),
            integration = "gauss-hermite", nNodes = 5L
        )
    )))
    model <- fit@model
    points <- model@slopes$integrationPoints
    weights <- model@slopes$drawWeights
    utility <- outer(model@slopes$alphas,
                     model@pricePre - model@priceOutside, "*")
    utility <- sweep(utility, 2L, model@slopes$meanval, "+")
    exp_utility <- exp(utility)
    expected_draw_shares <- t(sweep(exp_utility, 1L,
                                    1 + rowSums(exp_utility), "/"))

    expect_equal(dim(points), c(25L, 2L))
    expect_identical(model@nDraws, 25)
    expect_identical(model@slopes$nNodes, c(5L, 5L))
    expect_equal(unname(calcShares(model, aggregate = FALSE)), expected_draw_shares,
                 tolerance = 1e-13)
    expect_equal(as.vector(calcShares(model)),
                 as.vector(expected_draw_shares %*% weights),
                 tolerance = 1e-13)
    expected_derivative <- matrix(0, 3L, 3L)
    for (r in seq_len(ncol(expected_draw_shares))) {
        share_r <- expected_draw_shares[, r]
        expected_derivative <- expected_derivative +
            weights[r] * model@slopes$alphas[r] *
            (diag(share_r) - tcrossprod(share_r))
    }
    expect_equal(unname(elast(model, partial = TRUE)), expected_derivative,
                 tolerance = 1e-12)

    set.seed(917)
    before <- .Random.seed
    reused <- antitrust:::.blp_object_integration(model)
    after <- .Random.seed
    expect_identical(after, before)
    expect_identical(reused$integrationPoints, points)
    expect_identical(reused$weights, weights)
    expect_identical(reused$factorOrder, c("price", "demog1"))
    expect_identical(reused$nodesPerAxis, c(5L, 5L))

    result <- suppressWarnings(simulate(fit, counterfactual()))
    expect_identical(result@slopes$integrationPoints, points)
    expect_identical(result@slopes$drawWeights, weights)
    expect_identical(result@slopes$nNodes, c(5L, 5L))
})


test_that("calcMeanval contracts over the stored two-dimensional rule", {
    fit <- suppressMessages(suppressWarnings(specify(
        demand = "blp", conduct = "bertrand",
        prices = c(1.5, 1.8, 2.1), shares = c(.30, .25, .15),
        ownerPre = c("A", "B", "C"), insideSize = 100,
        parameters = list(
            alpha = -1.5, sigma = .2, meanval = c(.5, .2, -.1),
            piDemog = .3, demogMean = .2,
            demogCov = matrix(.64, 1L, 1L),
            integration = "gauss-hermite", nNodes = 5L
        )
    )))
    model <- fit@model
    target <- as.vector(calcShares(model))
    points <- model@slopes$integrationPoints
    weights <- model@slopes$drawWeights
    model@shares <- target
    model@shareInside <- sum(target)
    model@slopes$meanval <- NULL

    set.seed(919)
    before <- .Random.seed
    contracted <- suppressMessages(suppressWarnings(calcMeanval(model)))
    after <- .Random.seed

    expect_identical(after, before)
    expect_equal(as.vector(calcShares(contracted)), target, tolerance = 1e-8)
    expect_identical(contracted@slopes$integrationPoints, points)
    expect_identical(contracted@slopes$drawWeights, weights)
    expect_identical(contracted@slopes$factorOrder, c("price", "demog1"))
    expect_identical(contracted@slopes$nNodes, c(5L, 5L))
    expect_identical(contracted@nDraws, 25)
})


test_that("legacy materialized vector states retain every heterogeneous component", {
    fit <- suppressMessages(suppressWarnings(specify(
        demand = "blp", conduct = "bertrand",
        prices = c(1.5, 1.8, 2.1), shares = c(.30, .25, .15),
        ownerPre = c("A", "B", "C"), insideSize = 100,
        parameters = list(
            alpha = -1.5, sigma = .2, meanval = c(.5, .2, -.1),
            piDemog = .3, demogMean = .2,
            demogCov = matrix(.64, 1L, 1L),
            integration = "gauss-hermite", nNodes = 5L
        )
    )))
    legacy <- fit@model
    expected_alpha <- legacy@slopes$alphas
    expected_demog <- legacy@slopes$demogDraws
    legacy@slopes$integrationPoints <- NULL
    legacy@slopes$factorOrder <- NULL
    legacy@slopes$nodesPerAxis <- NULL

    set.seed(918)
    before <- .Random.seed
    rebuilt <- suppressMessages(suppressWarnings(calcSlopes(legacy)))
    after <- .Random.seed
    expect_identical(after, before)
    expect_equal(rebuilt@slopes$alphas, expected_alpha, tolerance = 0)
    expect_equal(rebuilt@slopes$demogDraws, expected_demog, tolerance = 0)
})


test_that("multidimensional BLP calibration is explicitly unsupported", {
    expect_error(
        calibrate(
            demand = "blp", conduct = "bertrand",
            prices = c(2, 2.3, 2.7), shares = c(.3, .25, .15),
            margins = c(.4, .35, .3), ownerPre = 1:3, s0 = .3,
            integrationPoints = matrix(c(-1, -1, 1, 1), ncol = 2L)
        ),
        "estimates only price random-coefficient heterogeneity"
    )
})


test_that("advanced BLP conduct routing is explicit and truthful", {
    parameters <- list(
        alpha = -1.5, sigma = .2, meanval = c(.5, .2, -.1),
        piDemog = .3, demogMean = .2,
        demogCov = matrix(.64, 1L, 1L),
        integration = "gauss-hermite", nNodes = 3L
    )
    expected_classes <- c(cournot = "CournotBLP", moncom = "MonComBLP")
    for (conduct in names(expected_classes)) {
        fit <- suppressMessages(suppressWarnings(specify(
            demand = "blp", conduct = conduct,
            prices = c(1.5, 1.8, 2.1), shares = c(.30, .25, .15),
            ownerPre = c("A", "B", "C"), insideSize = 100,
            parameters = parameters
        )))
        expect_s4_class(fit@model, expected_classes[[conduct]])
        expect_equal(dim(fit@model@slopes$integrationPoints), c(9L, 2L))
    }
    for (conduct in c("auction2nd", "bargaining")) {
        expect_error(
            specify(
                demand = "blp", conduct = conduct,
                prices = c(1.5, 1.8, 2.1), shares = c(.30, .25, .15),
                ownerPre = c("A", "B", "C"), insideSize = 100,
                parameters = parameters
            ),
            "remain price-random-coefficient-only"
        )
    }
})
