# Observed-data BLP calibration recovery tests.
#
# The synthetic margin moments are formed directly from demand primitives and
# conduct equations in this test file.  They are intentionally not generated
# by calcMargins(), so the recovery tests do not merely reproduce the method
# under test by construction.

# Recovery is split by cost: deterministic integration contracts remain in
# fast tests elsewhere, representative calibrated boundaries run in
# extended, and the all-conduct/multi-rule recovery grids run nightly.

.blp_recovery_fixture <- function(conduct, alpha = -5, sigma = .8) {
    prices <- c(1.5, 2, 2.6)
    owner <- c("A", "A", "B")
    nodes <- c(-1.5, -.5, .5, 1.5)
    weights <- c(.10, .20, .30, .40)
    target_s0 <- .20
    base_delta <- c(.7, .35, .05)

    aggregate_share <- function(shift) {
        antitrust:::.blp_stable_shares(
            base_delta + shift, prices, alpha + sigma * nodes,
            nodes, weights, outside = TRUE
        )$aggregate
    }
    shift <- uniroot(
        function(value) sum(aggregate_share(value)) - (1 - target_s0),
        c(-10, 10), tol = 1e-12
    )$root
    delta <- base_delta + shift
    shares <- aggregate_share(shift)
    s0 <- 1 - sum(shares)

    object <- suppressWarnings(antitrust:::.blp_model(
        conduct = conduct, prices = prices, shares = shares,
        margins = rep(.2, length(prices)), ownerPre = owner,
        alphaMean = alpha, sigma = sigma, meanval = delta,
        draws = nodes, drawWeights = weights, s0 = s0,
        output = TRUE,
        bargpowerPre = if (conduct == "bargaining") rep(.4, 3) else NULL
    ))
    draw_shares <- calcShares(object, aggregate = FALSE)
    draw_weights <- object@slopes$drawWeights
    alphas <- object@slopes$alphas
    owner_matrix <- object@ownerPre

    if (conduct == "bertrand") {
        derivative <- matrix(0, nrow = 3, ncol = 3)
        for (r in seq_along(alphas)) {
            s <- draw_shares[, r]
            derivative <- derivative + draw_weights[r] * alphas[r] *
                (diag(s) - tcrossprod(s))
        }
        revenue <- prices * shares / sum(prices * shares)
        elasticity <- derivative * outer(1 / shares, prices)
        margins <- -as.vector(
            solve(t(elasticity) * owner_matrix) %*%
                (revenue * diag(owner_matrix))
        ) / revenue
    } else if (conduct == "cournot") {
        derivative <- matrix(0, nrow = 3, ncol = 3)
        for (r in seq_along(alphas)) {
            s <- draw_shares[, r]
            derivative <- derivative + draw_weights[r] * alphas[r] *
                (diag(s) - tcrossprod(s))
        }
        margins <- -as.vector(
            (owner_matrix * t(solve(derivative))) %*% shares
        ) / prices
    } else if (conduct == "auction2nd") {
        firm_draw_shares <- owner_matrix %*% draw_shares
        firm_shares <- as.vector(firm_draw_shares %*% draw_weights)
        numerator <- as.vector(
            (log(1 - firm_draw_shares) /
                 matrix(alphas, nrow = 3, ncol = length(alphas), byrow = TRUE)) %*%
                draw_weights
        )
        margins <- numerator / firm_shares / prices
    } else if (conduct == "bargaining") {
        bargaining <- rep(.4, 3) / (1 - rep(.4, 3))
        aggregate_shares <- as.vector(draw_shares %*% draw_weights)
        derivative <- matrix(0, nrow = 3, ncol = 3)
        buyer_surplus <- numeric(3)
        for (r in seq_along(alphas)) {
            s <- draw_shares[, r]
            derivative <- derivative + draw_weights[r] * alphas[r] *
                (diag(s) - tcrossprod(s))
            buyer_surplus <- buyer_surplus + draw_weights[r] *
                log1p(-s) / alphas[r]
        }
        aggregate_elast <- derivative * outer(1 / aggregate_shares, prices)
        revenue <- prices * aggregate_shares
        margin_system <- t(
            diag(1 / revenue) %*%
                (t(aggregate_elast * owner_matrix) %*%
                     diag(aggregate_shares))
        )
        own_normalized <- diag(derivative) / aggregate_shares
        right_hand_side <- own_normalized /
            (-1 * (own_normalized - bargaining * aggregate_shares /
                   buyer_surplus))
        right_hand_side <- diag(owner_matrix) * right_hand_side
        margins <- as.vector(solve(t(margin_system), right_hand_side)) / prices
    } else {
        stop("unknown recovery conduct")
    }

    list(
        prices = prices, shares = shares, margins = margins,
        ownerPre = owner, s0 = s0, nodes = nodes, weights = weights,
        alphaMean = alpha, sigma = sigma, delta = delta,
        bargpowerPre = if (conduct == "bargaining") rep(.4, 3) else NULL
    )
}


test_that("BLP calibration recovers price heterogeneity under all supported conducts", {
    qa_skip_if_not_nightly()
    for (conduct in c("bertrand", "cournot", "auction2nd", "bargaining")) {
        fixture <- .blp_recovery_fixture(conduct)
        fit <- calibrate(
            demand = "blp", conduct = conduct,
            prices = fixture$prices, shares = fixture$shares,
            margins = fixture$margins, ownerPre = fixture$ownerPre,
            s0 = fixture$s0, output = TRUE,
            bargpowerPre = fixture$bargpowerPre,
            integration = "provided", draws = fixture$nodes,
            integrationWeights = fixture$weights,
            optimizer_control = list(maxit = 150, factr = 1e3, pgtol = 1e-8)
        )

        expect_equal(fit@parameters$alphaMean, fixture$alphaMean,
                     tolerance = 2e-3, info = conduct)
        expect_equal(fit@parameters$sigma, fixture$sigma,
                     tolerance = 2e-3, info = conduct)
        expect_equal(unname(calcShares(fit@model, preMerger = TRUE)),
                     fixture$shares, tolerance = 1e-10, info = conduct)
        expect_equal(fit@diagnostics$maxAbsResidual, 0, tolerance = 2e-6,
                     info = conduct)
        expect_equal(fit@diagnostics$preMergerFOCResidual,
                     fit@diagnostics$maxAbsResidual, tolerance = 0,
                     info = conduct)
        expect_equal(fit@diagnostics$integration$weights,
                     fixture$weights / sum(fixture$weights), tolerance = 0,
                     info = conduct)
    }
})


test_that("BLP calibration retains the sigma-zero boundary", {
    qa_skip_unless_tier("extended")
    for (conduct in c("bertrand", "cournot", "auction2nd", "bargaining")) {
        fixture <- .blp_recovery_fixture(conduct, sigma = 0)
        fit <- calibrate(
            demand = "blp", conduct = conduct,
            prices = fixture$prices, shares = fixture$shares,
            margins = fixture$margins, ownerPre = fixture$ownerPre,
            s0 = fixture$s0, output = TRUE,
            bargpowerPre = fixture$bargpowerPre,
            integration = "provided", draws = fixture$nodes,
            integrationWeights = fixture$weights,
            optimizer_control = list(maxit = 150, factr = 1e3, pgtol = 1e-8)
        )
        expect_equal(fit@parameters$alphaMean, fixture$alphaMean,
                     tolerance = 2e-3, info = conduct)
        expect_equal(fit@parameters$sigma, 0, tolerance = 2e-6,
                     info = conduct)
        expect_true(isTRUE(fit@diagnostics$sigmaOnBoundary), info = conduct)
    }
})


test_that("BLP diagnostics expose the fixed outside share and contraction state", {
    qa_skip_unless_tier("extended")
    fixture <- .blp_recovery_fixture("bertrand")
    fit <- calibrate(
        demand = "blp", conduct = "bertrand",
        prices = fixture$prices, shares = fixture$shares,
        margins = fixture$margins, ownerPre = fixture$ownerPre,
        s0 = fixture$s0, output = TRUE,
        integration = "provided", draws = fixture$nodes,
        integrationWeights = fixture$weights,
        optimizer_control = list(maxit = 150, factr = 1e3, pgtol = 1e-8)
    )
    expect_equal(fit@diagnostics$s0, fixture$s0, tolerance = 1e-12)
    expect_equal(fit@diagnostics$marginMoments, 3)
    expect_true(fit@diagnostics$contraction$converged)
    expect_lt(fit@diagnostics$contraction$maxError, 1e-9)
    expect_true(is.finite(fit@diagnostics$wrongSignProbability))
    expect_true(is.function(fit@diagnostics$profile_sigma))
})


test_that("adaptive BLP multistart is economically identical to exhaustive mode", {
    qa_skip_if_not_nightly()
    fixture <- .blp_recovery_fixture("bertrand")
    args <- list(
        demand = "blp", conduct = "bertrand",
        prices = fixture$prices, shares = fixture$shares,
        margins = fixture$margins, ownerPre = fixture$ownerPre,
        s0 = fixture$s0, output = TRUE,
        integration = "provided", draws = fixture$nodes,
        integrationWeights = fixture$weights,
        optimizer_control = list(maxit = 100, factr = 1e3, pgtol = 1e-8)
    )
    adaptive <- do.call(calibrate, c(args, list(multistart = "adaptive")))
    exhaustive <- do.call(calibrate, c(args, list(multistart = "exhaustive")))

    expect_equal(adaptive@parameters$alphaMean,
                 exhaustive@parameters$alphaMean, tolerance = 2e-6)
    expect_equal(adaptive@parameters$sigma, exhaustive@parameters$sigma,
                 tolerance = 2e-6)
    expect_equal(adaptive@diagnostics$objective,
                 exhaustive@diagnostics$objective, tolerance = 1e-10)
    expect_equal(unname(calcShares(adaptive@model)),
                 unname(calcShares(exhaustive@model)), tolerance = 1e-10)
    expect_equal(adaptive@diagnostics$preMergerFOCResidual,
                 exhaustive@diagnostics$preMergerFOCResidual, tolerance = 1e-10)

    adaptive_post <- simulate(adaptive, ownerPost = c("A", "A", "C"))
    exhaustive_post <- simulate(exhaustive, ownerPost = c("A", "A", "C"))
    expect_equal(adaptive_post@pricePost, exhaustive_post@pricePost,
                 tolerance = 2e-7)
    expect_equal(CV(adaptive_post), CV(exhaustive_post), tolerance = 2e-7)
    expect_lt(adaptive@diagnostics$multistart$evaluatedStarts,
              exhaustive@diagnostics$multistart$evaluatedStarts)
    expect_identical(exhaustive@diagnostics$multistart$strategy, "exhaustive")
    expect_equal(exhaustive@diagnostics$multistart$evaluatedStarts, 12L)
    expect_null(adaptive@diagnostics$profile_sigma_grid)
    expect_null(adaptive@diagnostics$profile_sigma_values)

    profiled <- profileBLP(adaptive, c(0, fixture$sigma))
    expect_equal(nrow(profiled), 2L)
    expect_true(is.finite(profiled$objective[[1L]]))
    expect_true(is.list(attr(profiled, "performance")))
})


test_that("BLP contraction reuses a price utility component without changing shares", {
    nodes <- c(-1, 0, 1)
    weights <- c(.2, .5, .3)
    prices <- c(1.5, 2, 2.6)
    alpha <- c(-4.8, -5, -5.2)
    delta <- c(.4, .2, -.1)
    expected <- antitrust:::.blp_stable_shares(
        delta, prices, alpha, nodes, weights, outside = TRUE
    )
    reused <- antitrust:::.blp_stable_shares(
        delta, prices, alpha, nodes, weights, outside = TRUE,
        priceUtility = outer(alpha, prices)
    )
    expect_equal(reused$draw, expected$draw, tolerance = 0)
    expect_equal(reused$aggregate, expected$aggregate, tolerance = 0)
})


test_that("no-demographics price-random-coefficient calibration works under both integration rules", {
    qa_skip_if_not_nightly()
    make_fixture <- function(nodes, weights) {
        prices <- c(1.5, 2, 2.6)
        owner <- c("A", "A", "B")
        alpha <- -5
        sigma <- .8
        base_delta <- c(.7, .35, .05)
        target_s0 <- .2
        aggregate_share <- function(delta) {
            antitrust:::.blp_stable_shares(
                delta, prices, alpha + sigma * nodes,
                nodes, weights, outside = TRUE
            )$aggregate
        }
        shift <- uniroot(
            function(value) {
                sum(aggregate_share(base_delta + value)) - (1 - target_s0)
            },
            c(-10, 10), tol = 1e-12
        )$root
        delta <- base_delta + shift
        shares <- aggregate_share(delta)
        model <- suppressWarnings(antitrust:::.blp_model(
            conduct = "bertrand", prices = prices, shares = shares,
            margins = rep(.2, length(prices)), ownerPre = owner,
            alphaMean = alpha, sigma = sigma, meanval = delta,
            draws = nodes, drawWeights = weights, s0 = 1 - sum(shares),
            output = TRUE
        ))
        draw_shares <- calcShares(model, aggregate = FALSE)
        alphas <- model@slopes$alphas
        derivative <- matrix(0, nrow = 3, ncol = 3)
        for (r in seq_along(alphas)) {
            s <- draw_shares[, r]
            derivative <- derivative + weights[r] * alphas[r] *
                (diag(s) - tcrossprod(s))
        }
        elasticity <- derivative * outer(1 / shares, prices)
        revenue <- prices * shares / sum(prices * shares)
        margins <- -as.vector(
            solve(t(elasticity) * model@ownerPre) %*%
                (revenue * diag(model@ownerPre))
        ) / revenue
        list(prices = prices, shares = shares, margins = margins,
             ownerPre = owner, s0 = 1 - sum(shares), alpha = alpha,
             sigma = sigma)
    }

    for (rule in c("gauss-hermite", "monte-carlo")) {
        if (identical(rule, "gauss-hermite")) {
            quadrature <- antitrust:::.blp_normal_nodes(15L)
            fixture <- make_fixture(quadrature$nodes, quadrature$weights)
            integration_args <- list(integration = rule, nNodes = 15L)
        } else {
            set.seed(20260906)
            nodes <- rnorm(15L)
            fixture <- make_fixture(nodes, rep(1 / 15, 15))
            set.seed(20260906)
            integration_args <- list(integration = rule, nDraws = 15L)
        }
        fit <- do.call(calibrate, c(list(
            demand = "blp", conduct = "bertrand",
            prices = fixture$prices, shares = fixture$shares,
            margins = fixture$margins, ownerPre = fixture$ownerPre,
            s0 = fixture$s0, output = TRUE,
            optimizer_control = list(maxit = 150, factr = 1e3, pgtol = 1e-8)
        ), integration_args))

        expect_equal(fit@parameters$alphaMean, fixture$alpha, tolerance = 2e-3,
                     info = rule)
        expect_equal(fit@parameters$sigma, fixture$sigma, tolerance = 2e-3,
                     info = rule)
        expect_identical(fit@diagnostics$integration$rule, rule)
        expect_length(fit@diagnostics$integration$nodes, 15L)
        expect_equal(unname(calcShares(fit@model, preMerger = TRUE)),
                     fixture$shares, tolerance = 2e-10, info = rule)
        expect_lt(fit@diagnostics$maxAbsResidual, 2e-5)
    }
})
