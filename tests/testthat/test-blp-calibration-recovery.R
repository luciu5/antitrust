# Observed-data BLP calibration recovery tests.
#
# The synthetic margin moments are formed directly from demand primitives and
# conduct equations in this test file.  They are intentionally not generated
# by calcMargins(), so the recovery tests do not merely reproduce the method
# under test by construction.

# Multi-start BLP recovery is intentionally extended: the fast integration
# tests still protect nodes, weights, derivatives, and deterministic reuse.
qa_skip_if_not_extended()

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
        margin_system <- matrix(0, nrow = 3, ncol = 3)
        right_hand_side <- numeric(3)
        for (r in seq_along(alphas)) {
            s <- draw_shares[, r]
            kernel <- -owner_matrix * rep(s, times = 3)
            diag(kernel) <- diag(owner_matrix) + diag(kernel)
            margin_system <- margin_system + draw_weights[r] * kernel
            div <- s / (1 - s)
            term <- log(1 - s) /
                (alphas[r] * (bargaining * div - log(1 - s)))
            right_hand_side <- right_hand_side + draw_weights[r] *
                diag(owner_matrix) * term
        }
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
