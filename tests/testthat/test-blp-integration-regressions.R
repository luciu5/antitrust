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


test_that("flat BLP draw-level slopes reduce to the Logit derivative", {
    nodes <- c(-1, 0, 1)
    weights <- c(.2, .6, .2)
    model <- blp_integration_test_model(
        nodes, weights, alphaMean = -1.2, sigma = .3
    )
    shares_draw <- calcShares(model, preMerger = TRUE, aggregate = FALSE)
    alpha <- model@slopes$alphas
    expected <- matrix(0, nrow = nrow(shares_draw), ncol = nrow(shares_draw))
    for (r in seq_along(alpha)) {
        s <- shares_draw[, r]
        expected <- expected + weights[r] * alpha[r] *
            (diag(s) - tcrossprod(s))
    }
    expect_equal(unname(elast(model, preMerger = TRUE, partial = TRUE)),
                 unname(expected), tolerance = 1e-12)
})


test_that("BLP dimension selection ignores empty demographics and mean-only characteristics", {
    expect_false(antitrust:::.blp_multidimensional(list(
        sigma = .1, piDemog = numeric(0), nDemog = 0
    )))
    expect_identical(
        antitrust:::.blp_integration(list(
            sigma = .1, piDemog = numeric(0), nDemog = 0
        ))$rule,
        "gauss-hermite"
    )
    expect_false(antitrust:::.blp_multidimensional(list(
        sigma = .1, prodChar = matrix(1, nrow = 3, ncol = 1),
        beta = 1, nDemog = 0
    )))
})


test_that("Monte Carlo BLP integration defaults to 5000 draws", {
    set.seed(20260906)
    rule <- antitrust:::.blp_integration(list(integration = "monte-carlo"))
    expect_identical(rule$rule, "monte-carlo")
    expect_length(rule$draws, 5000L)
    expect_equal(rule$weights, rep(1 / 5000, 5000), tolerance = 0)
})


test_that("BLP specify forwards market metadata into the new builder", {
    fit <- suppressWarnings(specify(
        demand = "blp", conduct = "bertrand",
        prices = c(2, 2.5, 3), shares = c(.2, .3, .25),
        ownerPre = 1:3, insideSize = 75, priceOutside = 1,
        labels = LETTERS[1:3],
        parameters = list(alpha = -1.5, sigma = .1, integration = "auto")
    ))
    expect_equal(fit@model@insideSize, 75, tolerance = 0)
    expect_equal(fit@model@mktSize, 100, tolerance = 0)
    expect_equal(fit@model@priceOutside, 1, tolerance = 0)
    expect_identical(fit@model@labels, LETTERS[1:3])
})


test_that("BLP specify derives no-outside shares from supplied meanval", {
    prices <- c(2, 2.2, 2.5)
    alpha <- -1.5
    delta <- c(.2, .1, -.1)
    expected <- exp(delta + alpha * prices)
    expected <- expected / sum(expected)
    fit <- suppressWarnings(specify(
        demand = "blp", conduct = "bertrand", prices = prices,
        ownerPre = c("A", "B", "C"), s0 = 0,
        parameters = list(alpha = alpha, sigma = 0, meanval = delta,
                          draws = 0, drawWeights = 1)
    ))
    expect_equal(unname(fit@model@shares), expected, tolerance = 1e-12)
    expect_equal(sum(fit@model@shares), 1, tolerance = 1e-12)
    expect_equal(fit@diagnostics$s0, 0, tolerance = 0)
    expect_equal(fit@model@normIndex, 1, tolerance = 0)
})


test_that("default BLP specify dispatches registered auction conduct", {
    fit <- suppressWarnings(specify(
        demand = "blp", conduct = "auction2nd",
        prices = c(2, 2.5, 3), shares = c(.2, .3, .1),
        ownerPre = 1:3,
        parameters = list(
            alpha = -1.5, sigma = 0,
            meanval = log(c(.2, .3, .1) / .4) + 1.5 * c(2, 2.5, 3)
        )
    ))
    expect_s4_class(fit@model, "Auction2ndBLP")
    expect_identical(fit@model@slopes$integration, "gauss-hermite")
})


test_that("default BLP specify dispatches registered bargaining conduct", {
    fit <- suppressWarnings(specify(
        demand = "blp", conduct = "bargaining",
        prices = c(2, 2.5, 3), shares = c(.2, .3, .1),
        ownerPre = 1:3,
        parameters = list(
            alpha = -1.5, sigma = 0,
            meanval = log(c(.2, .3, .1) / .4) + 1.5 * c(2, 2.5, 3)
        )
    ))
    expect_s4_class(fit@model, "BargainingBLP")
    expect_identical(fit@model@slopes$integration, "gauss-hermite")
})


test_that("one demographic with zero price sigma defaults to Gauss-Hermite", {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(20260904)
    before <- .Random.seed
    first <- qa_value(sim(
        prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .20),
        supply = "bertrand", demand = "BLP",
        demand.param = list(
            alpha = -1, sigma = 0, piDemog = .1,
            demogMean = .3, demogCov = matrix(.25, nrow = 1, ncol = 1),
            meanval = c(.5, .3, .1)
        ),
        ownerPre = c("A", "B", "C"),
        ownerPost = c("A", "A", "C"),
        insideSize = 100
    ), "one-dimensional demographic BLP")
    after <- .Random.seed

    expect_identical(before, after)
    expect_identical(first@slopes$integration, "gauss-hermite")
    expect_length(first@slopes$consDraws, 31L)
    expect_equal(
        first@slopes$demogDraws[, 1],
        .3 + .5 * first@slopes$consDraws,
        tolerance = 1e-14
    )
})


test_that("price-only BLP simulation selects the requested integration rule", {
    args <- list(
        prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .20),
        supply = "bertrand", demand = "BLP",
        demand.param = list(
            alpha = -1, sigma = .1, meanval = c(.5, .3, .1),
            integration = "auto"
        ),
        ownerPre = c("A", "B", "C"),
        ownerPost = c("A", "A", "C"),
        insideSize = 100
    )
    quadrature <- suppressWarnings(do.call(sim, c(args, list(nNodes = 15L))))
    expect_identical(quadrature@slopes$integration, "gauss-hermite")
    expect_length(quadrature@slopes$consDraws, 15L)
    expect_equal(sum(quadrature@slopes$drawWeights), 1, tolerance = 1e-14)

    set.seed(918273)
    mc_args <- args
    mc_args$demand.param$integration <- "monte-carlo"
    monte_carlo <- suppressWarnings(do.call(sim, c(mc_args, list(nDraws = 15L))))
    expect_identical(monte_carlo@slopes$integration, "monte-carlo")
    expect_length(monte_carlo@slopes$consDraws, 15L)
    expect_equal(monte_carlo@slopes$drawWeights, rep(1 / 15, 15),
                 tolerance = 0)
})


test_that("PriceLeadershipBLP uses the shared integration rules", {
    qa_skip_unless_tier("extended")
    shares <- c(.35, .25, .25, .15)
    prices <- c(.93, .88, 1.10, 1.02)
    alpha <- -5.767013
    common <- list(
        prices = prices, shares = shares,
        ownerPre = c("Bank1", "Bank2", "Bank3", "Fringe"),
        ownerPost = c("Bank1", "Bank2", "Bank3", "Fringe"),
        coalitionPre = 1:3, coalitionPost = 1:3,
        insideSize = 1000,
        slopes = list(
            alphaMean = alpha, alpha = alpha, sigma = .5,
            meanval = c(0, log(shares[-1] / shares[1]) -
                alpha * (prices[-1] - prices[1])), sigmaNest = 1
        )
    )
    default_rule <- suppressWarnings(do.call(ple.blp, common))
    expect_identical(default_rule@slopes$integration, "gauss-hermite")
    expect_length(default_rule@slopes$consDraws, 31L)

    quadrature <- suppressWarnings(do.call(
        ple.blp, c(common, list(integration = "gauss-hermite", nNodes = 15L))
    ))
    expect_identical(quadrature@slopes$integration, "gauss-hermite")
    expect_length(quadrature@slopes$consDraws, 15L)
    expect_equal(sum(quadrature@slopes$drawWeights), 1, tolerance = 1e-14)

    monte_carlo <- suppressWarnings(do.call(
        ple.blp, c(common, list(integration = "monte-carlo", nDraws = 15L))
    ))
    expect_identical(monte_carlo@slopes$integration, "monte-carlo")
    expect_length(monte_carlo@slopes$consDraws, 15L)
    expect_equal(monte_carlo@slopes$drawWeights, rep(1 / 15, 15),
                 tolerance = 0)
})


test_that("PriceLeadershipBLP uses shared Gauss-Hermite nodes for one demographic", {
    qa_skip_unless_tier("extended")
    shares <- c(.35, .25, .25, .15)
    prices <- c(.93, .88, 1.10, 1.02)
    common <- list(
        prices = prices, shares = shares,
        ownerPre = c("Bank1", "Bank2", "Bank3", "Fringe"),
        ownerPost = c("Bank1", "Bank2", "Bank3", "Fringe"),
        coalitionPre = 1:3, coalitionPost = 1:3,
        insideSize = 1000,
        slopes = list(
            alphaMean = -5.767013, alpha = -5.767013, sigma = 0,
            piDemog = .4, demogMean = .2,
            demogCov = matrix(4, nrow = 1L, ncol = 1L),
            meanval = c(0, log(shares[-1] / shares[1]) -
                (-5.767013) * (prices[-1] - prices[1])), sigmaNest = 1
        )
    )
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(20260907)
    before <- .Random.seed
    fit <- suppressWarnings(do.call(ple.blp, common))
    after <- .Random.seed
    nodes <- antitrust:::.blp_normal_nodes(31L)$nodes

    expect_identical(before, after)
    expect_identical(fit@slopes$nDemog, 1L)
    expect_identical(fit@slopes$integration, "gauss-hermite")
    expect_equal(fit@slopes$demogDraws, matrix(.2 + 2 * nodes, ncol = 1L),
                 tolerance = 0)
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
        sigmaChar = .2
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


test_that("BLP consDraws and integrationWeights aliases are exact", {
    nodes <- c(-1, 0, 1)
    weights <- c(.2, .5, .3)
    rule <- antitrust:::.blp_integration(list(
        consDraws = nodes, integrationWeights = weights
    ))
    expect_identical(rule$rule, "provided")
    expect_equal(rule$draws, nodes, tolerance = 0)
    expect_equal(rule$weights, weights / sum(weights), tolerance = 0)

    fit <- suppressWarnings(specify(
        demand = "blp", conduct = "bertrand",
        prices = c(1.5, 1.8, 2.1), shares = c(.30, .25, .15),
        ownerPre = c("A", "B", "C"),
        parameters = list(
            alpha = -1.2, sigma = .3, meanval = c(.5, .2, -.1),
            consDraws = nodes, integrationWeights = weights
        )
    ))
    expect_identical(fit@model@slopes$integration, "provided")
    expect_equal(fit@model@slopes$consDraws, nodes, tolerance = 0)
    expect_equal(fit@model@slopes$drawWeights,
                 weights / sum(weights), tolerance = 0)
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


test_that("legacy BLP retains wrong-sign draws under the explicit domain contract", {
    model <- blp_integration_test_model(
        nodes = c(-1, 0, 1), weights = rep(1 / 3, 3),
        alphaMean = -1, sigma = 1.5
    )
    expect_warning({ model <- calcSlopes(model) }, "retained")
    expect_equal(model@slopes$alphas, c(-2.5, -1, .5), tolerance = 0)
})


test_that("BLP fits reuse their integration rule across counterfactual simulations", {
    qa_skip_unless_tier("extended")
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
    qa_skip_unless_tier("extended")
    nodes <- c(-2, -.5, .5, 2)
    weights <- c(.05, .55, .20, .20)
    model <- blp_integration_test_model(
        nodes, weights, alphaMean = -1, sigma = .25,
        meanval = c(.2, .1, 0)
    )
    expect_equal(model@mktSize, 1 / sum(model@shares), tolerance = 1e-14)
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
    expected <- model@mktSize *
        sum(weights[keep] * cv_by_draw[keep]) / sum(weights[keep])
    expect_equal(CV(model, lim = c(.5, 1)), expected, tolerance = 1e-13)
})
