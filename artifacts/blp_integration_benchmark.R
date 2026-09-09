## Development benchmark for the one-dimensional price-random-coefficient
## BLP path.  This file is not part of the package build.

devtools::load_all(quiet = TRUE)

prices <- c(1.5, 2.0, 2.6)
base_delta <- c(.7, .35, .05)
alpha_mean <- -5
sigma <- .8
owner <- c("A", "B", "C")
s0_target <- .2

delta <- base_delta
share_sum <- function(shift) {
    sum(vapply(seq_along(prices), function(j) {
        integrate(function(z) {
            vapply(z, function(zz) {
                utility <- base_delta + shift +
                    (alpha_mean + sigma * zz) * prices
                scale <- max(c(0, utility))
                (exp(utility[j] - scale) /
                     (exp(-scale) + sum(exp(utility - scale)))) * dnorm(zz)
            }, numeric(1))
        }, -Inf, Inf, rel.tol = 1e-10)$value
    }, numeric(1)))
}
delta <- base_delta + uniroot(
    function(shift) share_sum(shift) - (1 - s0_target),
    c(-10, 10), tol = 1e-11
)$root

share_at <- function(z) {
    utility <- delta + (alpha_mean + sigma * z) * prices
    scale <- max(c(0, utility))
    exp(utility - scale) /
        (exp(-scale) + sum(exp(utility - scale)))
}

share_reference <- vapply(seq_along(prices), function(j) {
    integrate(function(z) {
        vapply(z, function(zz) share_at(zz)[j] * dnorm(zz), numeric(1))
    }, -Inf, Inf, rel.tol = 1e-11)$value
}, numeric(1))
s0 <- 1 - sum(share_reference)

derivative_reference <- matrix(0, nrow = length(prices), ncol = length(prices))
for (i in seq_along(prices)) for (j in seq_along(prices)) {
    derivative_reference[i, j] <- integrate(function(z) {
        vapply(z, function(zz) {
            s <- share_at(zz)
            (alpha_mean + sigma * zz) *
                (s[i] * (as.numeric(i == j) - s[j])) * dnorm(zz)
        }, numeric(1))
    }, -Inf, Inf, rel.tol = 1e-10)$value
}

revenue <- prices * share_reference / sum(prices * share_reference)
owner_matrix <- ownerToMatrix(
    suppressWarnings(.blp_model(
        conduct = "bertrand", prices = prices, shares = share_reference,
        margins = rep(.2, length(prices)), ownerPre = owner,
        alphaMean = alpha_mean, sigma = sigma, meanval = delta,
        draws = 0, drawWeights = 1, s0 = s0
    )), preMerger = TRUE
)
elasticity_reference <- derivative_reference * outer(1 / share_reference, prices)
margin_reference <- -as.vector(
    solve(t(elasticity_reference) * owner_matrix) %*%
        (revenue * diag(owner_matrix))
) / revenue

## A calibration benchmark reports phase costs and contraction work rather than
## only elapsed time.  The retained profile is invoked explicitly here to
## reproduce the former eager-grid cost while keeping calibrate() on its
## on-demand path.
benchmark_profile_grid <- function(fit) {
    alpha_scale <- max(abs(alpha_mean), 1e-4)
    sigma_fit <- fit@parameters$sigma
    sigma_starts <- c(0, .10, .25, .50) * alpha_scale
    unique(c(0, sigma_fit, sigma_starts,
             seq(0, max(2 * sigma_fit, alpha_scale), length.out = 7)))
}

benchmark_calibration <- function(n, strategy = c("adaptive", "exhaustive")) {
    strategy <- match.arg(strategy)
    started <- proc.time()[["elapsed"]]
    fit <- calibrate(
        demand = "blp", conduct = "bertrand", prices = prices,
        shares = share_reference, margins = margin_reference,
        ownerPre = owner, s0 = s0, integration = "gauss-hermite",
        nNodes = n, multistart = strategy,
        optimizer_control = list(maxit = 100, factr = 1e3, pgtol = 1e-8)
    )
    calibration_elapsed <- proc.time()[["elapsed"]] - started
    before <- fit@diagnostics$performance$metrics
    profile_started <- proc.time()[["elapsed"]]
    profile <- profileBLP(fit, benchmark_profile_grid(fit))
    profile_elapsed <- proc.time()[["elapsed"]] - profile_started
    after <- attr(profile, "performance")$metrics
    metric_delta <- function(name) {
        before_value <- if (is.null(before[[name]])) 0 else before[[name]]
        after_value <- if (is.null(after[[name]])) 0 else after[[name]]
        after_value - before_value
    }
    post <- simulate(fit, ownerPost = c("A", "A", "C"))
    data.frame(
        nNodes = n, strategy = strategy,
        runtime_seconds = calibration_elapsed,
        multistart_seconds = fit@diagnostics$performance$multistartSeconds,
        identification_seconds = fit@diagnostics$performance$identificationSeconds,
        profile_seconds = profile_elapsed,
        evaluated_starts = fit@diagnostics$multistart$evaluatedStarts,
        convergence_count = fit@diagnostics$multistart$convergedStarts,
        fallback = fit@diagnostics$multistart$fallback,
        objective_evaluations_multistart = fit@diagnostics$performance$objectiveEvaluationsMultistart,
        objective_evaluations_identification = fit@diagnostics$performance$objectiveEvaluationsIdentification,
        objective_evaluations_profile = metric_delta("objective_evaluations_profile"),
        contraction_calls_multistart = fit@diagnostics$performance$metrics$contraction_calls_multistart,
        contraction_calls_profile = metric_delta("contraction_calls_profile"),
        contraction_iterations_multistart = fit@diagnostics$performance$metrics$contraction_iterations_multistart,
        contraction_iterations_identification = fit@diagnostics$performance$metrics$contraction_iterations_identification,
        contraction_iterations_profile = metric_delta("contraction_iterations_profile"),
        alpha_error = fit@parameters$alphaMean - alpha_mean,
        sigma_error = fit@parameters$sigma - sigma,
        objective = fit@diagnostics$objective,
        post_price_1 = post@pricePost[[1L]],
        post_price_2 = post@pricePost[[2L]],
        post_price_3 = post@pricePost[[3L]],
        cv = CV(post)
    )
}

measure_rule <- function(n) {
    rule <- .blp_normal_nodes(n)
    shares <- vapply(seq_along(prices), function(j) {
        sum(rule$weights * vapply(rule$nodes,
                                  function(z) share_at(z)[j], numeric(1)))
    }, numeric(1))
    derivative <- matrix(0, nrow = length(prices), ncol = length(prices))
    for (r in seq_along(rule$nodes)) {
        s <- share_at(rule$nodes[r])
        derivative <- derivative + rule$weights[r] *
            (alpha_mean + sigma * rule$nodes[r]) *
            (diag(s) - tcrossprod(s))
    }
    c(
        nNodes = n,
        share_error = max(abs(shares - share_reference)),
        derivative_error = max(abs(derivative - derivative_reference))
    )
}

rule_results <- do.call(rbind, lapply(c(10L, 15L, 20L, 30L, 40L), measure_rule))
print(rule_results)

calibration_result <- function(n) {
    started <- proc.time()[["elapsed"]]
    fit <- calibrate(
        demand = "blp", conduct = "bertrand", prices = prices,
        shares = share_reference, margins = margin_reference,
        ownerPre = owner, s0 = s0, integration = "gauss-hermite",
        nNodes = n,
        optimizer_control = list(maxit = 100, factr = 1e3, pgtol = 1e-8)
    )
    elapsed <- proc.time()[["elapsed"]] - started
    post <- simulate(fit, ownerPost = c("A", "A", "C"))
    c(
        nNodes = n,
        alpha_error = fit@parameters$alphaMean - alpha_mean,
        sigma_error = fit@parameters$sigma - sigma,
        post_price_1 = post@pricePost[[1L]],
        post_price_2 = post@pricePost[[2L]],
        post_price_3 = post@pricePost[[3L]],
        cv = CV(post),
        objective = fit@diagnostics$objective,
        runtime_seconds = elapsed
    )
}

calibration_results <- do.call(rbind, lapply(c(10L, 20L, 30L, 40L),
                                             calibration_result))
print(calibration_results)

## Exhaustive mode is intentionally opt-in for this expensive comparison.
## Set the tier variable when producing the before/after report; ordinary
## development runs retain the smaller integration-rule benchmark above.
if (identical(Sys.getenv("ANTITRUST_BLP_BENCHMARK_TIER"), "nightly")) {
    benchmark_results <- do.call(rbind, list(
        benchmark_calibration(10L, "adaptive"),
        benchmark_calibration(10L, "exhaustive")
    ))
    print(benchmark_results)
}

mc_measure <- function(n, seed = 20260903) {
    started <- proc.time()[["elapsed"]]
    set.seed(seed)
    draws <- rnorm(n)
    weights <- rep(1 / n, n)
    shares <- vapply(seq_along(prices), function(j) {
        sum(weights * vapply(draws, function(z) share_at(z)[j], numeric(1)))
    }, numeric(1))
    derivative <- matrix(0, nrow = length(prices), ncol = length(prices))
    for (r in seq_along(draws)) {
        s <- share_at(draws[r])
        derivative <- derivative + weights[r] *
            (alpha_mean + sigma * draws[r]) *
            (diag(s) - tcrossprod(s))
    }
    elapsed <- proc.time()[["elapsed"]] - started
    c(
        nDraws = n,
        share_error = max(abs(shares - share_reference)),
        derivative_error = max(abs(derivative - derivative_reference)),
        runtime_seconds = elapsed
    )
}

mc_results <- rbind(mc_measure(1000L), mc_measure(10000L))
print(mc_results)
