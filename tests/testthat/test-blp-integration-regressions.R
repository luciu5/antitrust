# Regression tests for BLP integration and aggregation semantics.
#
# These tests deliberately use supplied one-dimensional integration points so
# that the expected aggregates can be formed independently of any optimizer
# or random-number generation.

blp_integration_test_model <- function(nodes, weights,
                                       alphaMean = -1.2, sigma = 0.3,
                                       meanval = c(.5, .2, -.1)) {
    ## Use the internal model builder so these tests isolate the integration
    ## contract from the legacy constructor's calibration/cost plumbing.
    ## The internal builder deliberately initializes ownerPost = ownerPre;
    ## suppress that unrelated no-merger warning in this demand-only fixture.
    suppressWarnings(antitrust:::.blp_model(
        conduct = "bertrand",
        prices = c(1.5, 1.8, 2.1),
        shares = c(.30, .25, .15),
        margins = c(.20, .25, .30),
        ownerPre = c("A", "B", "C"),
        alphaMean = alphaMean, sigma = sigma, meanval = meanval,
        draws = nodes, drawWeights = weights, s0 = .30
    ))
}


test_that("Gauss-Hermite normal rule reproduces standard moments", {
    for (n in c(10L, 15L, 20L)) {
        rule <- antitrust:::.blp_normal_nodes(n)
        expect_equal(sum(rule$weights), 1, tolerance = 1e-14)
        moments <- vapply(c(0, 1, 2, 3, 4, 6), function(power) {
            sum(rule$weights * rule$nodes^power)
        }, numeric(1))
        expect_equal(moments, c(1, 0, 1, 0, 3, 15), tolerance = 1e-10,
                     info = paste("n =", n))
    }
})


test_that("Gauss-Hermite shares agree with direct one-dimensional integration", {
    prices <- c(1.5, 1.8, 2.1)
    delta <- c(.5, .2, -.1)
    alpha <- -1.2
    sigma <- .3
    share_at <- function(z) {
        utility <- delta + (alpha + sigma * z) * prices
        scale <- max(c(0, utility))
        exp(utility - scale) /
            (exp(-scale) + sum(exp(utility - scale)))
    }
    reference <- vapply(seq_along(prices), function(j) {
        stats::integrate(function(z) {
            vapply(z, function(zz) share_at(zz)[j] * stats::dnorm(zz),
                   numeric(1))
        }, -Inf, Inf,
                         subdivisions = 200L, rel.tol = 1e-11)$value
    }, numeric(1))
    rule <- antitrust:::.blp_normal_nodes(30L)
    quadrature <- vapply(seq_along(prices), function(j) {
        sum(rule$weights * vapply(rule$nodes,
                                  function(z) share_at(z)[j], numeric(1)))
    }, numeric(1))
    expect_equal(quadrature, reference, tolerance = 1e-10)
})


test_that("BLP aggregate derivatives pass a finite-difference check", {
    nodes <- c(-1.5, -.25, .75, 1.75)
    weights <- c(.05, .15, .30, .50)
    model <- blp_integration_test_model(nodes, weights)
    analytic <- unname(elast(model, preMerger = TRUE, partial = TRUE))
    numerical <- matrix(0, nrow = length(model@pricePre),
                        ncol = length(model@pricePre))
    h <- 1e-6
    for (j in seq_along(model@pricePre)) {
        plus <- model
        minus <- model
        plus@pricePre[j] <- plus@pricePre[j] + h
        minus@pricePre[j] <- minus@pricePre[j] - h
        numerical[, j] <- (calcShares(plus) - calcShares(minus)) / (2 * h)
    }
    expect_equal(analytic, numerical, tolerance = 2e-9)
})


test_that("integration node and weight validation is strict", {
    expect_error(antitrust:::.blp_normal_nodes(10.5), "positive integer")
    expect_error(antitrust:::.blp_integration(list(
        integration = "provided", draws = c(-1, 0),
        integrationWeights = c(1, -1)
    )), "non-negative")
    expect_error(antitrust:::.blp_integration(list(
        integration = "gauss-hermite", nNodes = 5,
        prodChar = matrix(1, nrow = 3, ncol = 1)
    )), "one-dimensional")
})


test_that("BLP aggregate shares honor supplied integration weights", {
    nodes <- c(-1, 0, 2)
    weights <- c(.10, .20, .70)
    model <- blp_integration_test_model(nodes, weights / sum(weights))

    prices <- model@pricePre
    alphas <- model@slopes$alphas
    meanval <- model@slopes$meanval
    utility <- outer(alphas, prices - model@priceOutside, "*") +
        matrix(meanval, nrow = length(alphas), ncol = length(prices),
               byrow = TRUE)
    exp_utility <- exp(utility)
    denominator <- 1 + rowSums(exp_utility)
    expected_draws <- t(sweep(exp_utility, 1, denominator, "/"))
    expected <- as.vector(expected_draws %*% weights)

    observed_draws <- calcShares(model, preMerger = TRUE, aggregate = FALSE)
    observed <- calcShares(model, preMerger = TRUE, aggregate = TRUE)

    expect_equal(unname(observed_draws), expected_draws, tolerance = 1e-14)
    expect_equal(unname(observed), expected, tolerance = 1e-14)
})


test_that("BLP aggregate derivatives honor supplied integration weights", {
    nodes <- c(-1, 0, 2)
    weights <- c(.10, .20, .70)
    model <- blp_integration_test_model(nodes, weights)
    draw_shares <- unname(calcShares(model, aggregate = FALSE))
    alphas <- model@slopes$alphas

    expected <- matrix(0, nrow = nrow(draw_shares), ncol = nrow(draw_shares))
    for (r in seq_along(alphas)) {
        s <- draw_shares[, r]
        expected <- expected + weights[r] * alphas[r] *
            (diag(s) - tcrossprod(s))
    }

    observed <- unname(elast(model, preMerger = TRUE, partial = TRUE))
    expect_equal(observed, expected, tolerance = 1e-13)
})


test_that("provided BLP points and weights are normalized without RNG use", {
    nodes <- c(-2, -.5, .5, 2)
    weights <- c(1, 2, 4, 3)
    set.seed(20260903)
    before <- .Random.seed
    rule <- antitrust:::.blp_integration(list(
        draws = nodes, drawWeights = weights
    ))
    after <- .Random.seed

    expect_identical(after, before)
    expect_equal(rule$draws, nodes, tolerance = 0)
    expect_equal(rule$weights, weights / sum(weights), tolerance = 0)
    expect_identical(rule$rule, "provided")

    set.seed(20260904)
    before_fit <- .Random.seed
    model <- blp_integration_test_model(nodes, weights / sum(weights))
    after_fit <- .Random.seed
    expect_identical(after_fit, before_fit)
    expect_equal(model@slopes$consDraws, nodes, tolerance = 0)
    expect_equal(model@slopes$drawWeights, weights / sum(weights),
                 tolerance = 0)
})


test_that("BLP repeated share and derivative evaluations are deterministic", {
    nodes <- c(-1.5, -.25, .75, 1.75)
    weights <- c(.05, .15, .30, .50)
    model <- blp_integration_test_model(nodes, weights)

    shares_one <- calcShares(model, preMerger = TRUE)
    shares_two <- calcShares(model, preMerger = TRUE)
    derivatives_one <- elast(model, preMerger = TRUE, partial = TRUE)
    derivatives_two <- elast(model, preMerger = TRUE, partial = TRUE)

    expect_identical(shares_one, shares_two)
    expect_identical(derivatives_one, derivatives_two)
})


test_that("BLP fits reuse their integration rule across counterfactual simulations", {
    nodes <- c(-1.5, -.25, .75, 1.75)
    weights <- c(.05, .15, .30, .50)
    fit <- specify(
        demand = "blp", conduct = "bertrand",
        prices = c(1.5, 1.8, 2.1),
        shares = c(.30, .25, .15),
        ownerPre = c("A", "B", "C"),
        parameters = list(
            alphaMean = -1.2, sigma = .3,
            meanval = c(.5, .2, -.1),
            draws = nodes, drawWeights = weights
        )
    )
    original_prices <- fit@model@pricePre
    original_draws <- fit@model@slopes$consDraws
    original_weights <- fit@model@slopes$drawWeights
    first <- simulate(fit, ownerPost = c("A", "A", "C"))
    second <- simulate(fit, ownerPost = c("A", "A", "C"))

    expect_identical(fit@model@pricePre, original_prices)
    expect_identical(fit@model@slopes$consDraws, original_draws)
    expect_identical(fit@model@slopes$drawWeights, original_weights)
    expect_equal(first@pricePost, second@pricePost, tolerance = 1e-12)
    expect_equal(first@slopes$drawWeights, weights / sum(weights),
                 tolerance = 0)
})


test_that("BLP CV trimming uses integration-weighted quantiles and means", {
    nodes <- c(-2, -.5, .5, 2)
    weights <- c(.05, .55, .20, .20)
    model <- blp_integration_test_model(
        nodes, weights, alphaMean = -1, sigma = .25,
        meanval = c(.2, .1, 0)
    )
    model@pricePost <- model@pricePre + c(.2, .1, .3)

    alphas <- model@slopes$alphas
    meanval <- model@slopes$meanval
    pre_utility <- outer(alphas, model@pricePre - model@priceOutside) +
        matrix(meanval, nrow = length(alphas), ncol = length(meanval),
               byrow = TRUE)
    post_utility <- outer(alphas, model@pricePost - model@priceOutside) +
        matrix(meanval, nrow = length(alphas), ncol = length(meanval),
               byrow = TRUE)
    v_pre <- log(1 + rowSums(exp(pre_utility)))
    v_post <- log(1 + rowSums(exp(post_utility)))
    cv_by_draw <- (v_post - v_pre) / alphas

    ## The weighted median is alpha[2] here (cumulative weight .60), so
    ## draws 2:4 survive lim = c(.5, 1).  Ordinary unweighted quantiles
    ## would instead start at the midpoint between alpha[2] and alpha[3].
    keep <- 2:4
    expected <- sum(weights[keep] * cv_by_draw[keep]) / sum(weights[keep])
    expect_equal(CV(model, lim = c(.5, 1)), expected, tolerance = 1e-13)
})
