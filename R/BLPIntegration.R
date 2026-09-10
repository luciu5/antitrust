# Shared integration helpers for random-coefficients Logit models.
#
# The economic code consumes only integration points and normalized weights.
# Deterministic quadrature and fixed-draw Monte Carlo therefore follow the
# same aggregation path.

#' @title Select BLP Integration Points and Weights
#' @description Public entry point that selects the integration points and
#'   normalized weights for a BLP (random-coefficients Logit) demand
#'   specification, given a list of demand parameters. It chooses between
#'   Gauss-Hermite quadrature, fixed-draw Monte Carlo, or caller-supplied
#'   points according to the \code{integration} rule and the dimensionality of
#'   the heterogeneity. This is a thin, stable wrapper over the internal
#'   integration selection logic; downstream packages building models on the
#'   BLP demand system (e.g. price leadership with BLP demand) can call it
#'   instead of re-implementing integration-rule selection.
#' @param slopes A list of BLP demand parameters. Recognized elements include
#'   \code{integration} (\code{"auto"}, \code{"gauss-hermite"},
#'   \code{"monte-carlo"}, or \code{"provided"}), \code{nNodes}, \code{nDraws},
#'   \code{integrationPoints}, \code{consDraws}/\code{draws},
#'   \code{integrationWeights}/\code{drawWeights},
#'   \code{sigma}, \code{nDemog}, \code{piDemog}, \code{sigmaChar}, and
#'   \code{pi}.
#' @return A list containing \code{draws}, \code{weights} (normalized to sum
#'   to one), and \code{rule}. Two-dimensional rules additionally return an
#'   \code{integrationPoints} matrix, canonical \code{factorOrder}, and
#'   \code{nodesPerAxis} for Gauss-Hermite quadrature.
#' @export
calcBLPintegration <- function(slopes) {
  .blp_integration(slopes)
}

.blp_single_demographic_dimension <- function(dots) {
    sigma <- dots[["sigma"]]
    n_demog <- dots[["nDemog"]]
    if (is.null(n_demog)) n_demog <- length(dots[["piDemog"]])
    is_single_demog <- length(n_demog) == 1L && is.numeric(n_demog) &&
        is.finite(n_demog) && n_demog == 1L
    is_zero_sigma <- length(sigma) == 1L && is.numeric(sigma) &&
        is.finite(sigma) && as.numeric(sigma) == 0
    no_random_characteristics <- is.null(dots[["sigmaChar"]]) ||
        length(dots[["sigmaChar"]]) == 0L
    is_single_demog && is_zero_sigma && no_random_characteristics
}


.blp_active_factors <- function(dots) {
    sigma <- dots[["sigma"]]
    price_active <- !is.null(sigma) && length(sigma) == 1L &&
        is.numeric(sigma) && is.finite(sigma) && as.numeric(sigma) != 0

    n_demog <- dots[["nDemog"]]
    if (is.null(n_demog)) n_demog <- length(dots[["piDemog"]])
    if (length(n_demog) != 1L || !is.numeric(n_demog) ||
        !is.finite(n_demog) || n_demog < 0 || n_demog != as.integer(n_demog)) {
        return(NULL)
    }
    n_demog <- as.integer(n_demog)

    ## Demographics form one correlated latent block.  If any coordinate is
    ## loaded, retain every latent coordinate so that the supplied covariance
    ## matrix is represented with its complete Cholesky factor.
    pi_demog <- dots[["piDemog"]]
    if (is.null(pi_demog)) pi_demog <- numeric(0)
    pi <- dots[["pi"]]
    demog_active <- logical(n_demog)
    if (n_demog > 0L) {
        if (length(pi_demog) == n_demog && is.numeric(pi_demog)) {
            demog_active <- demog_active | abs(pi_demog) > 0
        }
        if (is.matrix(pi) && nrow(pi) == n_demog) {
            demog_active <- demog_active | apply(abs(pi) > 0, 1L, any)
        }
    }

    sigma_char <- dots[["sigmaChar"]]
    char_active <- if (is.null(sigma_char)) logical(0) else {
        if (!is.numeric(sigma_char) || any(!is.finite(sigma_char))) return(NULL)
        abs(sigma_char) > 0
    }
    if (any(demog_active)) demog_active[] <- TRUE
    factor_names <- c(
        if (price_active) "price",
        if (any(demog_active)) paste0("demog", which(demog_active)),
        if (any(char_active)) paste0("char", which(char_active))
    )
    list(
        names = factor_names,
        price = price_active,
        demog = which(demog_active),
        char = which(char_active),
        nDemog = n_demog,
        dimension = length(factor_names)
    )
}


.blp_integration_dimensions <- function(dots) {
    factors <- .blp_active_factors(dots)
    if (is.null(factors)) return(Inf)
    factors$dimension
}


.blp_multidimensional <- function(dots) {
    .blp_integration_dimensions(dots) > 1L
}


.blp_n_points <- function(points) {
    if (is.matrix(points)) return(nrow(points))
    length(points)
}


.blp_validate_nnodes <- function(nNodes, dimension) {
    if (is.null(nNodes)) {
        return(if (dimension > 1L) c(15L, 15L) else 31L)
    }
    if (!is.numeric(nNodes) || any(!is.finite(nNodes)) ||
        any(nNodes < 1) || any(nNodes != as.integer(nNodes)) ||
        length(nNodes) > 2L || length(nNodes) < 1L) {
        stop("'nNodes' must be a positive integer, or a length-two vector for two-dimensional quadrature.")
    }
    nNodes <- as.integer(nNodes)
    if (dimension <= 1L && length(nNodes) != 1L) {
        stop("'nNodes' must be scalar for one-dimensional BLP integration.")
    }
    if (dimension > 2L) {
        stop("Gauss-Hermite integration supports at most two active dimensions.")
    }
    if (dimension == 2L && length(nNodes) == 1L) nNodes <- rep(nNodes, 2L)
    nNodes
}


.blp_tensor_normal_rule <- function(nNodes) {
    nNodes <- as.integer(nNodes)
    rules <- lapply(nNodes, .blp_normal_nodes)
    if (length(rules) == 1L) {
        return(list(points = rules[[1L]]$nodes, weights = rules[[1L]]$weights,
                    nodes = rules[[1L]]$nodes, nodesPerAxis = nNodes))
    }
    grid <- expand.grid(lapply(rules, `[[`, "nodes"), KEEP.OUT.ATTRS = FALSE)
    points <- as.matrix(grid)
    weights <- as.vector(outer(rules[[1L]]$weights, rules[[2L]]$weights))
    list(points = points, weights = weights, nodes = points,
         nodesPerAxis = nNodes)
}


.blp_validate_points <- function(points, weights = NULL) {
    if (is.matrix(points)) {
        stop("BLP 'draws'/'consDraws' integration points must be a numeric vector; use 'integrationPoints' for two-dimensional points.")
    }
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


.blp_validate_integration_points <- function(points, weights = NULL,
                                             dimension = NULL) {
    if (!is.matrix(points) || !is.numeric(points) || nrow(points) < 1L ||
        ncol(points) < 2L || any(!is.finite(points))) {
        stop("'integrationPoints' must be a non-empty finite numeric matrix with at least two columns.")
    }
    if (!is.null(dimension) && (is.null(dimension) ||
                                ncol(points) != dimension)) {
        stop("'integrationPoints' columns must match the active BLP factor dimension.")
    }
    if (is.null(weights)) weights <- rep(1 / nrow(points), nrow(points))
    if (!is.numeric(weights) || length(weights) != nrow(points) ||
        any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
        stop("BLP integration weights must be finite, non-negative, and match the integration point rows.")
    }
    list(points = points, weights = as.numeric(weights / sum(weights)))
}


.blp_draw_weights <- function(object, nDraws = length(object@slopes[["alphas"]])) {
    weights <- object@slopes[["drawWeights"]]
    if (is.null(weights)) weights <- object@slopes[["integrationWeights"]]
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
    ## Use exact list lookup throughout this adapter.  In particular,
    ## `$integration` partially matches `integrationWeights` when callers
    ## omit the optional integration rule.
    factors <- .blp_active_factors(dots)
    dimension <- if (is.null(factors)) Inf else factors$dimension
    supplied_points <- dots[["integrationPoints"]]
    supplied_draws <- dots[["draws"]]
    if (is.null(supplied_draws)) supplied_draws <- dots[["consDraws"]]
    supplied_weights <- dots[["integrationWeights"]]
    if (is.null(supplied_weights)) supplied_weights <- dots[["drawWeights"]]

    ## Validate nDraws even when caller-supplied points take precedence.  The
    ## legacy BLP boundary treats nDraws as a required positive scalar, and
    ## silently accepting zero here leaves an invalid parameter object in the
    ## migrated path.
    supplied_n <- dots[["nDraws"]]
    if (!is.null(supplied_n) &&
        (!is.numeric(supplied_n) || length(supplied_n) != 1L ||
         !is.finite(supplied_n) || supplied_n < 1)) {
        stop("'nDraws' must be a positive scalar.")
    }

    requested <- dots[["integration"]]
    if (is.null(requested)) requested <- "auto"
    requested <- match.arg(requested,
                           c("auto", "gauss-hermite", "monte-carlo", "provided"))

    if (!is.null(supplied_points) && !is.null(supplied_draws)) {
        stop("supply BLP points through either 'integrationPoints' or 'draws'/'consDraws', not both.")
    }
    supplied_count <- if (!is.null(supplied_points)) {
        nrow(supplied_points)
    } else if (!is.null(supplied_draws)) {
        .blp_n_points(supplied_draws)
    } else NULL
    if (!is.null(supplied_n) && !is.null(supplied_count) &&
        supplied_n != supplied_count) {
        stop("'nDraws' must equal the number of supplied BLP integration points.")
    }

    if (!is.null(supplied_points)) {
        if (!identical(requested, "auto") && !identical(requested, "provided")) {
            stop("supplied BLP integration points conflict with integration = '",
                 requested, "'.")
        }
        if (!is.finite(dimension)) {
            stop("the active BLP factor dimension must be valid when 'integrationPoints' is supplied.")
        }
        validated <- .blp_validate_integration_points(
            supplied_points, supplied_weights, dimension
        )
        return(list(
            draws = validated$points,
            points = validated$points,
            integrationPoints = validated$points,
            weights = validated$weights,
            rule = "provided", factorOrder = factors$names,
            factorNames = factors$names, factors = factors, dimension = dimension,
            nDraws = nrow(validated$points), nodesPerAxis = NULL
        ))
    }

    if (!is.null(supplied_draws)) {
        if (!identical(requested, "auto") && !identical(requested, "provided")) {
            stop("supplied BLP integration points conflict with integration = '",
                 requested, "'.")
        }
        if (is.matrix(supplied_draws)) {
            ## .blp_validate_points emits the public alias guidance.
            .blp_validate_points(supplied_draws, supplied_weights)
        }
        if (is.finite(dimension) && dimension > 1L) {
            stop("vector 'draws'/'consDraws' can represent only one BLP factor; use 'integrationPoints' for two-dimensional points.")
        }
        validated <- .blp_validate_points(supplied_draws, supplied_weights)
        return(list(draws = validated$points, points = validated$points,
                    weights = validated$weights,
                    integrationPoints = NULL, rule = "provided",
                    factorOrder = factors$names, factorNames = factors$names,
                    factors = factors, dimension = dimension,
                    nDraws = length(validated$points),
                    nodesPerAxis = NULL))
    }

    method <- if (identical(requested, "auto")) {
        if (is.finite(dimension) && dimension > 2L) "monte-carlo" else "gauss-hermite"
    } else {
        requested
    }

    if (identical(method, "gauss-hermite") &&
        (!is.finite(dimension) || dimension > 2L)) {
        stop("Gauss-Hermite integration supports at most two active dimensions; use integration = 'monte-carlo'.")
    }

    if (identical(method, "provided")) {
        stop("integration = 'provided' requires supplied BLP integration points.")
    }

    if (identical(method, "gauss-hermite")) {
        n <- .blp_validate_nnodes(dots[["nNodes"]], dimension)
        rule <- .blp_tensor_normal_rule(n)
        points <- rule$points
        if (is.matrix(points) && ncol(points) == 1L) points <- as.numeric(points[, 1L])
        return(list(
            draws = points, points = points,
            integrationPoints = if (dimension > 1L) points else NULL,
            weights = rule$weights, rule = "gauss-hermite",
            factorOrder = factors$names, factorNames = factors$names,
            factors = factors, dimension = dimension,
            nDraws = length(rule$weights),
            nodesPerAxis = rule$nodesPerAxis
        ))
    }

    n <- if (is.null(dots[["nDraws"]])) 5000L else dots[["nDraws"]]
    if (length(n) != 1L || !is.finite(n) || n < 1 || n != as.integer(n)) {
        stop("'nDraws' must be a positive integer for Monte Carlo integration.")
    }
    ## Draw exactly once, outside every contraction/objective evaluation.
    n <- as.integer(n)
    draws <- if (is.finite(dimension) && dimension > 1L) {
        matrix(rnorm(n * dimension), nrow = n, ncol = dimension)
    } else {
        rnorm(n)
    }
    list(
        draws = draws, points = draws,
        integrationPoints = if (is.matrix(draws)) draws else NULL,
        weights = rep(1 / n, n), rule = "monte-carlo",
        factorOrder = factors$names, factorNames = factors$names,
        factors = factors, dimension = dimension, nDraws = n,
        nodesPerAxis = NULL
    )
}


## Materialize every heterogeneous component from one canonical standardized
## draw matrix.  The first column is price, followed by the complete
## demographic block, followed by random characteristic factors.
.blp_materialize_draws <- function(integration, alphaMean, sigma = 0,
                                   nDemog = 0L, piDemog = numeric(0),
                                   demogMean = NULL, demogCov = NULL,
                                   prodChar = NULL, sigmaChar = NULL, pi = NULL,
                                   output = TRUE, storedDemogDraws = NULL,
                                   storedCharDraws = NULL) {
    weights <- integration$weights
    points <- integration$draws
    n <- length(weights)
    z <- if (is.matrix(points)) points else matrix(as.numeric(points), ncol = 1L)
    if (nrow(z) != n) stop("BLP integration points and weights have incompatible sizes.")
    factors <- integration$factorOrder
    if (is.null(factors)) factors <- character(0)
    if (length(factors) != ncol(z) && ncol(z) > 1L) {
        stop("BLP integration factor metadata does not match the point matrix.")
    }
    factor_col <- function(name) {
        if (!length(factors)) return(NA_integer_)
        match(name, factors)
    }
    price_col <- factor_col("price")
    price_z <- if (is.na(price_col)) rep(0, n) else z[, price_col]
    ## Legacy one-dimensional consDraws are the original standardized node
    ## vector even when that sole factor is demographic or characteristic.
    cons_draws <- if (is.matrix(points)) price_z else as.numeric(points)
    if (!length(factors) && !is.matrix(points)) price_z <- cons_draws

    if (is.null(demogMean)) demogMean <- rep(0, nDemog)
    if (nDemog > 0L &&
        (!is.numeric(demogMean) || length(demogMean) != nDemog ||
         any(!is.finite(demogMean)))) {
        stop("'demogMean' must be finite and match 'nDemog'.")
    }
    if (nDemog > 0L) {
        if (is.null(demogCov)) demogCov <- diag(nDemog)
        if (!is.matrix(demogCov) || !identical(dim(demogCov), c(nDemog, nDemog)) ||
            any(!is.finite(demogCov))) {
            stop("'demogCov' must be a finite nDemog by nDemog matrix.")
        }
        demog_chol <- tryCatch(chol(demogCov), error = function(e)
            stop("demogCov must be positive definite: ", e$message))
        demog_cols <- vapply(seq_len(nDemog), function(j) factor_col(paste0("demog", j)), integer(1))
        if (any(!is.na(demog_cols))) {
            z_demog <- matrix(0, nrow = n, ncol = nDemog)
            z_demog[, !is.na(demog_cols)] <- z[, demog_cols[!is.na(demog_cols)], drop = FALSE]
            demogDraws <- sweep(z_demog %*% demog_chol, 2L, demogMean, "+")
        } else if (!is.null(storedDemogDraws) && is.matrix(storedDemogDraws) &&
                   identical(dim(storedDemogDraws), c(n, nDemog))) {
            demogDraws <- storedDemogDraws
        } else {
            demogDraws <- matrix(0, nrow = n, ncol = nDemog)
            demogDraws <- sweep(demogDraws, 2L, demogMean, "+")
        }
    } else {
        demogDraws <- matrix(0, nrow = n, ncol = 0L)
        demogMean <- numeric(0)
    }

    alphas <- as.numeric(alphaMean) + as.numeric(sigma) * price_z
    if (nDemog > 0L && length(piDemog)) {
        if (!is.numeric(piDemog) || length(piDemog) != nDemog || any(!is.finite(piDemog))) {
            stop("'piDemog' must be finite and match 'nDemog'.")
        }
        alphas <- alphas + as.vector(sweep(demogDraws, 2L, demogMean, "-") %*% piDemog)
    }

    nChar <- if (is.matrix(prodChar)) ncol(prodChar) else 0L
    charDraws <- NULL
    char_random <- matrix(0, nrow = n, ncol = if (is.matrix(prodChar)) nrow(prodChar) else 0L)
    if (nChar > 0L && !is.null(sigmaChar)) {
        if (!is.numeric(sigmaChar) || length(sigmaChar) != nChar || any(!is.finite(sigmaChar))) {
            stop("'sigmaChar' must be finite and match the characteristic count.")
        }
        char_cols <- vapply(seq_len(nChar), function(j) factor_col(paste0("char", j)), integer(1))
        if (any(!is.na(char_cols))) {
            charDraws <- matrix(0, nrow = n, ncol = nChar)
            charDraws[, !is.na(char_cols)] <- z[, char_cols[!is.na(char_cols)], drop = FALSE]
        } else if (!is.null(storedCharDraws) && is.matrix(storedCharDraws) &&
                   identical(dim(storedCharDraws), c(n, nChar))) {
            charDraws <- storedCharDraws
        } else {
            charDraws <- matrix(0, nrow = n, ncol = nChar)
        }
        char_random <- char_random + sweep(charDraws, 2L, sigmaChar, "*") %*% t(prodChar)
    }
    if (nChar > 0L && nDemog > 0L && !is.null(pi)) {
        if (!is.matrix(pi) || !identical(dim(pi), c(nDemog, nChar)) || any(!is.finite(pi))) {
            stop("'pi' must have dimensions nDemog by the characteristic count.")
        }
        char_random <- char_random + demogDraws %*% pi %*% t(prodChar)
    }
    expected_sign <- if (isTRUE(output)) -1 else 1
    wrong_sign <- if (expected_sign > 0) alphas <= 0 else alphas >= 0
    list(
        standardizedDraws = z, consDraws = cons_draws, priceDraws = price_z,
        demogDraws = demogDraws, charDraws = charDraws, alphas = alphas,
        char_random = char_random, weights = weights,
        wrongSignMass = sum(weights[wrong_sign])
    )
}


.blp_object_integration <- function(object, legacy_default = "monte-carlo") {
    slopes <- object@slopes
    has_points <- !is.null(slopes[["integrationPoints"]]) ||
        !is.null(slopes[["draws"]]) || !is.null(slopes[["consDraws"]])
    requested <- slopes[["integration"]]
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
    ## A materialized matrix is authoritative; do not pass a legacy alias
    ## alongside it as an ambiguous second specification.
    stored_points <- slopes[["integrationPoints"]]
    stored_draws <- if (is.null(stored_points)) slopes[["draws"]] else NULL
    stored_cons_draws <- if (is.null(stored_points)) slopes[["consDraws"]] else NULL
    legacy_points <- if (!is.null(stored_draws)) stored_draws else stored_cons_draws
    legacy_vector_state <- is.null(stored_points) && !is.null(legacy_points) &&
        !is.matrix(legacy_points) &&
        (!is.null(slopes[["demogDraws"]]) || !is.null(slopes[["charDraws"]]))
    integration_sigma <- if (legacy_vector_state) 0 else slopes[["sigma"]]
    integration_n_demog <- if (legacy_vector_state) 0L else slopes[["nDemog"]]
    integration_pi_demog <- if (legacy_vector_state) numeric(0) else slopes[["piDemog"]]
    integration_sigma_char <- if (legacy_vector_state) NULL else slopes[["sigmaChar"]]
    integration_pi <- if (legacy_vector_state) NULL else slopes[["pi"]]
    result <- .blp_integration(list(
        integration = requested,
        nNodes = slopes[["nNodes"]],
        nDraws = if (has_points) NULL else object@nDraws,
        integrationPoints = stored_points,
        draws = stored_draws,
        consDraws = stored_cons_draws,
        integrationWeights = slopes[["integrationWeights"]],
        drawWeights = slopes[["drawWeights"]],
        prodChar = if (legacy_vector_state) NULL else slopes[["prodChar"]],
        sigmaChar = integration_sigma_char,
        sigma = integration_sigma,
        pi = integration_pi,
        piDemog = integration_pi_demog,
        nDemog = integration_n_demog,
        demogMean = slopes[["demogMean"]],
        demogCov = slopes[["demogCov"]]
    ))
    ## A fitted object records the rule used to create its stored points.  The
    ## points/weights remain authoritative when re-used by a downstream path.
    if (!is.null(slopes[["integration"]]) &&
        slopes[["integration"]] %in% c("gauss-hermite", "monte-carlo", "provided")) {
        result$rule <- slopes[["integration"]]
    }
    if (!is.null(slopes[["factorOrder"]])) {
        result$factorOrder <- slopes[["factorOrder"]]
        result$factorNames <- slopes[["factorOrder"]]
    }
    if (!is.null(slopes[["nodesPerAxis"]])) {
        result$nodesPerAxis <- slopes[["nodesPerAxis"]]
    }
    if (isTRUE(slopes[["integrationWeightsNormalized"]])) {
        stored_weights <- slopes[["integrationWeights"]]
        if (is.null(stored_weights)) stored_weights <- slopes[["drawWeights"]]
        if (!is.null(stored_weights)) result$weights <- as.numeric(stored_weights)
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
