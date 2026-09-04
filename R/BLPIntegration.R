# Shared integration helpers for random-coefficients Logit models.
#
# The economic code consumes only integration points and normalized weights.
# Deterministic quadrature and fixed-draw Monte Carlo therefore follow the
# same aggregation path.

.blp_single_demographic_dimension <- function(dots) {
    sigma <- dots$sigma
    n_demog <- dots$nDemog
    if (is.null(n_demog)) n_demog <- length(dots$piDemog)
    is_single_demog <- length(n_demog) == 1L && is.numeric(n_demog) &&
        is.finite(n_demog) && n_demog == 1L
    is_zero_sigma <- length(sigma) == 1L && is.numeric(sigma) &&
        is.finite(sigma) && as.numeric(sigma) == 0
    no_random_characteristics <- is.null(dots$sigmaChar) ||
        length(dots$sigmaChar) == 0L
    is_single_demog && is_zero_sigma && no_random_characteristics
}


.blp_multidimensional <- function(dots) {
    ## A single demographic with sigma = 0 is one normal integration
    ## dimension, even when demographic interactions load product
    ## characteristics.  Keep the existing conservative Monte Carlo behavior
    ## for all other characteristic/demographic specifications.
    if (.blp_single_demographic_dimension(dots)) return(FALSE)
    indicators <- c("prodChar", "sigmaChar", "pi", "piDemog",
                    "demogMean", "demogCov")
    any(vapply(indicators, function(name) !is.null(dots[[name]]), logical(1))) ||
        (!is.null(dots$nDemog) && isTRUE(as.numeric(dots$nDemog) > 0))
}


.blp_validate_points <- function(points, weights = NULL) {
    if (!is.numeric(points) || length(points) < 1L ||
        any(!is.finite(points))) {
        stop("BLP integration points must be a non-empty finite numeric vector.")
    }
    points <- as.numeric(points)
    if (is.null(weights)) {
        weights <- rep(1 / length(points), length(points))
    }
    if (!is.numeric(weights) || length(weights) != length(points) ||
        any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
        stop("BLP integration weights must be finite, non-negative, and match the integration points.")
    }
    list(points = points, weights = as.numeric(weights / sum(weights)))
}


.blp_draw_weights <- function(object, nDraws = length(object@slopes$alphas)) {
    weights <- object@slopes$drawWeights
    if (is.null(weights)) weights <- object@slopes$integrationWeights
    if (is.null(weights)) weights <- rep(1 / nDraws, nDraws)
    .blp_validate_points(seq_len(nDraws), weights)$weights
}


.blp_normal_nodes <- function(n) {
    if (length(n) != 1L || !is.numeric(n) || !is.finite(n) ||
        n < 1 || n != as.integer(n)) {
        stop("'nNodes' must be a positive integer.")
    }
    n <- as.integer(n)
    if (n == 1L) return(list(nodes = 0, weights = 1))

    ## Golub-Welsch nodes and weights for a standard normal integral.
    off_diag <- sqrt(seq_len(n - 1L) / 2)
    jacobi <- matrix(0, nrow = n, ncol = n)
    jacobi[cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)] <- off_diag
    jacobi[cbind(seq_len(n - 1L) + 1L, seq_len(n - 1L))] <- off_diag
    eig <- eigen(jacobi, symmetric = TRUE)
    order_nodes <- order(eig$values)
    list(
        nodes = sqrt(2) * eig$values[order_nodes],
        weights = eig$vectors[1, order_nodes]^2
    )
}


.blp_quadrature_demog_draws <- function(points, nDemog, demogMean = NULL,
                                        demogCov = NULL) {
    if (length(nDemog) != 1L || !is.numeric(nDemog) || nDemog != 1L) {
        stop("Gauss-Hermite demographic integration requires exactly one demographic dimension.")
    }
    if (is.null(demogMean)) demogMean <- 0
    if (is.null(demogCov)) demogCov <- matrix(1, nrow = 1L, ncol = 1L)
    if (!is.numeric(demogMean) || length(demogMean) != 1L ||
        !is.finite(demogMean) || !is.matrix(demogCov) ||
        !identical(dim(demogCov), c(1L, 1L)) || !is.finite(demogCov[1, 1]) ||
        demogCov[1, 1] <= 0) {
        stop("demogCov must be positive definite for Gauss-Hermite demographic integration.")
    }
    matrix(as.numeric(demogMean) + sqrt(demogCov[1, 1]) * as.numeric(points),
           ncol = 1L)
}


.blp_integration <- function(dots) {
    supplied_draws <- dots$draws
    if (is.null(supplied_draws)) supplied_draws <- dots$consDraws
    supplied_weights <- dots$integrationWeights
    if (is.null(supplied_weights)) supplied_weights <- dots$drawWeights

    requested <- dots$integration
    if (is.null(requested)) requested <- "auto"
    requested <- match.arg(requested,
                           c("auto", "gauss-hermite", "monte-carlo", "provided"))

    if (!is.null(supplied_draws)) {
        if (!identical(requested, "auto") && !identical(requested, "provided")) {
            stop("supplied BLP integration points conflict with integration = '",
                 requested, "'.")
        }
        validated <- .blp_validate_points(supplied_draws, supplied_weights)
        return(list(draws = validated$points, weights = validated$weights,
                    rule = "provided"))
    }

    multidimensional <- .blp_multidimensional(dots)
    method <- if (identical(requested, "auto")) {
        if (multidimensional) "monte-carlo" else "gauss-hermite"
    } else {
        requested
    }

    if (identical(method, "gauss-hermite") && multidimensional) {
        stop("Gauss-Hermite integration is currently supported only for one-dimensional BLP heterogeneity; use integration = 'monte-carlo'.")
    }

    if (identical(method, "provided")) {
        stop("integration = 'provided' requires supplied BLP integration points.")
    }

    if (identical(method, "gauss-hermite")) {
        n <- if (!is.null(dots$nNodes)) dots$nNodes else 31L
        rule <- .blp_normal_nodes(n)
        return(list(draws = rule$nodes, weights = rule$weights,
                    rule = "gauss-hermite"))
    }

    n <- if (is.null(dots$nDraws)) 1000L else dots$nDraws
    if (length(n) != 1L || !is.finite(n) || n < 1 || n != as.integer(n)) {
        stop("'nDraws' must be a positive integer for Monte Carlo integration.")
    }
    ## Draw exactly once, outside every contraction/objective evaluation.
    draws <- rnorm(as.integer(n))
    list(draws = draws, weights = rep(1 / n, n), rule = "monte-carlo")
}


.blp_object_integration <- function(object, legacy_default = "monte-carlo") {
    slopes <- object@slopes
    has_points <- !is.null(slopes$draws) || !is.null(slopes$consDraws)
    requested <- slopes$integration
    if (has_points) requested <- "provided"
    if (is.null(requested)) {
        ## Legacy objects retain their Monte Carlo default, except when the
        ## only heterogeneous dimension is one demographic and sigma = 0.
        ## That case is exactly one normal dimension and can be integrated
        ## deterministically without changing other legacy defaults.
        requested <- if (identical(legacy_default, "monte-carlo") &&
                         .blp_single_demographic_dimension(slopes)) {
            "auto"
        } else {
            legacy_default
        }
    }
    result <- .blp_integration(list(
        integration = requested,
        nNodes = slopes$nNodes,
        nDraws = object@nDraws,
        draws = slopes$draws,
        consDraws = slopes$consDraws,
        integrationWeights = slopes$integrationWeights,
        drawWeights = slopes$drawWeights,
        prodChar = slopes$prodChar,
        sigmaChar = slopes$sigmaChar,
        sigma = slopes$sigma,
        pi = slopes$pi,
        piDemog = slopes$piDemog,
        nDemog = slopes$nDemog,
        demogMean = slopes$demogMean,
        demogCov = slopes$demogCov
    ))
    ## A fitted object records the rule used to create its stored points.  The
    ## points/weights remain authoritative when re-used by a downstream path.
    if (!is.null(slopes$integration) &&
        slopes$integration %in% c("gauss-hermite", "monte-carlo", "provided")) {
        result$rule <- slopes$integration
    }
    result
}


.blp_weighted_quantile <- function(x, probabilities, weights) {
    stopifnot(length(x) == length(weights), length(probabilities) >= 1L)
    validated <- .blp_validate_points(x, weights)
    x <- validated$points
    weights <- validated$weights
    ## Retain the historical interpolation behavior for equal-weight MC.
    if (max(abs(weights - rep(1 / length(weights), length(weights)))) < 1e-14) {
        return(as.numeric(stats::quantile(x, probabilities,
                                          names = FALSE, type = 7)))
    }
    order_x <- order(x)
    x <- x[order_x]
    weights <- weights[order_x]
    cweights <- cumsum(weights)
    vapply(probabilities, function(probability) {
        if (probability <= 0) return(x[[1L]])
        if (probability >= 1) return(x[[length(x)]])
        x[[which(cweights >= probability)[1L]]]
    }, numeric(1))
}
