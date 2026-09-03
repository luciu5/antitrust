# BLP calibration through the refactor architecture.
#
# This file intentionally contains only the observed-data calibration adapter
# and the two small conduct-specific BLP subclasses needed by the existing
# auction and bargaining methods.  The legacy S4 demand and conduct equations
# remain the source of truth.

#' @rdname BertrandRUM-Classes
#' @export
setClass(
    Class = "Auction2ndBLP",
    contains = "Auction2ndLogit",
    slots = list(nDraws = "numeric"),
    prototype = prototype(nDraws = 31)
)

#' @rdname Bargaining-Classes
#' @export
setClass(
    Class = "BargainingBLP",
    contains = "BargainingLogit",
    slots = list(nDraws = "numeric"),
    prototype = prototype(nDraws = 31)
)


.blp_draw_weights <- function(object, nDraws = length(object@slopes$alphas)) {
    weights <- object@slopes$drawWeights
    if (is.null(weights)) weights <- rep(1 / nDraws, nDraws)
    if (length(weights) != nDraws || any(!is.finite(weights)) ||
        any(weights < 0) || sum(weights) <= 0) {
        stop("BLP draw weights must be finite, non-negative, and match the number of draws.")
    }
    as.numeric(weights / sum(weights))
}


.blp_normal_nodes <- function(n) {
    n <- as.integer(n)
    if (length(n) != 1L || is.na(n) || n < 1L) {
        stop("'nDraws' must be a positive integer.")
    }
    if (n == 1L) return(list(nodes = 0, weights = 1))

    ## Golub-Welsch nodes and weights for a standard normal integral.  This
    ## avoids adding an integration dependency while keeping the rule fixed
    ## throughout the outer calibration optimizer.
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


.blp_integration <- function(dots) {
    supplied_draws <- dots$draws
    supplied_weights <- dots$drawWeights
    integration <- dots$integration
    if (is.null(integration)) {
        integration <- if (is.null(supplied_draws)) "gauss-hermite" else "provided"
    }
    integration <- match.arg(integration,
                             c("gauss-hermite", "monte-carlo", "provided"))

    if (!is.null(supplied_draws)) {
        if (!is.numeric(supplied_draws) || length(supplied_draws) < 1L ||
            any(!is.finite(supplied_draws))) {
            stop("'draws' must be a finite numeric vector when supplied.")
        }
        draws <- as.numeric(supplied_draws)
        weights <- if (is.null(supplied_weights)) {
            rep(1 / length(draws), length(draws))
        } else {
            supplied_weights
        }
        integration <- "provided"
    } else if (identical(integration, "gauss-hermite")) {
        n <- if (is.null(dots$nDraws)) 31L else dots$nDraws
        rule <- .blp_normal_nodes(n)
        draws <- rule$nodes
        weights <- rule$weights
    } else {
        n <- if (is.null(dots$nDraws)) 1000L else dots$nDraws
        if (length(n) != 1L || !is.finite(n) || n < 1 || n != as.integer(n)) {
            stop("'nDraws' must be a positive integer.")
        }
        ## Monte Carlo draws are made once here, outside every objective
        ## evaluation.  They are therefore fixed conditional on the caller's
        ## RNG state and cannot introduce objective noise.
        draws <- rnorm(as.integer(n))
        weights <- rep(1 / n, n)
    }

    if (!is.numeric(weights) || length(weights) != length(draws) ||
        any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
        stop("'drawWeights' must be finite, non-negative, and match 'draws'.")
    }
    list(
        draws = draws,
        weights = as.numeric(weights / sum(weights)),
        rule = integration
    )
}


.blp_stable_shares <- function(delta, prices, alpha, draws, weights,
                               priceOutside = 0, outside = TRUE) {
    utility <- outer(alpha, prices - priceOutside, "*")
    utility <- sweep(utility, 2, delta, "+")
    max_utility <- apply(utility, 1, max)
    max_utility <- if (outside) pmax(0, max_utility) else max_utility
    exp_inside <- exp(utility - max_utility)
    denominator <- rowSums(exp_inside)
    if (outside) denominator <- denominator + exp(-max_utility)
    draw_shares <- exp_inside / denominator
    weighted <- as.vector(crossprod(weights, draw_shares))
    list(draw = draw_shares, aggregate = weighted)
}


.blp_contract <- function(prices, shares, alphaMean, sigma, draws, weights,
                          s0, priceOutside = 0, tol = 1e-10,
                          maxIter = 2000L, initial = NULL) {
    alpha <- alphaMean + sigma * draws
    outside <- s0 > 0
    delta <- if (is.null(initial)) log(shares) else as.numeric(initial)
    if (length(delta) != length(shares) || any(!is.finite(delta))) {
        delta <- log(shares)
    }

    converged <- FALSE
    max_error <- Inf
    for (iter in seq_len(as.integer(maxIter))) {
        predicted <- .blp_stable_shares(
            delta, prices, alpha, draws, weights, priceOutside, outside
        )$aggregate
        if (any(!is.finite(predicted)) || any(predicted <= 0)) break
        error <- log(shares) - log(predicted)
        max_error <- max(abs(error))
        delta_new <- delta + error
        if (max_error < tol) {
            delta <- delta_new
            converged <- TRUE
            break
        }
        delta <- delta_new
    }

    if (!outside) delta <- delta - delta[1]
    final <- .blp_stable_shares(
        delta, prices, alpha, draws, weights, priceOutside, outside
    )
    list(
        delta = delta,
        predicted = final$aggregate,
        drawShares = final$draw,
        alpha = alpha,
        converged = converged,
        iterations = if (exists("iter")) iter else 0L,
        maxError = max_error
    )
}


.blp_validate_inputs <- function(prices, shares, margins, ownerPre, s0,
                                  output = TRUE) {
    n <- length(prices)
    if (!is.numeric(prices) || n < 2L || any(!is.finite(prices)) ||
        any(prices <= 0)) stop("'prices' must be a finite, positive numeric vector.")
    if (!is.numeric(shares) || length(shares) != n || any(!is.finite(shares)) ||
        any(shares <= 0) || any(shares >= 1)) {
        stop("'shares' must be a positive finite vector strictly below one.")
    }
    if (!is.numeric(margins) || length(margins) != n) {
        stop("'margins' must be a numeric vector with the same length as 'prices'.")
    }
    if (sum(!is.na(margins)) < 2L) {
        stop("BLP calibration requires at least two observed margin moments.")
    }
    if (any(!is.finite(margins[!is.na(margins)])) ||
        any(margins[!is.na(margins)] <= 0) ||
        (isTRUE(output) && any(margins[!is.na(margins)] > 1))) {
        stop("observed BLP margins must be positive proportional margins; output-market margins must not exceed one.")
    }
    if (!is.numeric(s0) || length(s0) != 1L || !is.finite(s0) ||
        s0 < 0 || s0 >= 1) stop("'s0' must be a single number in [0, 1).")
    if (!isTRUE(all.equal(sum(shares), 1 - s0, tolerance = 1e-8))) {
        stop("'shares' must sum to 1 - s0; BLP calibration does not infer or estimate s0.")
    }
    if (!is.logical(output) || length(output) != 1L || is.na(output)) {
        stop("'output' must be a single logical value.")
    }
    if (is.matrix(ownerPre)) {
        if (nrow(ownerPre) != n || ncol(ownerPre) != n) {
            stop("matrix 'ownerPre' must be square with one row and column per product.")
        }
    } else if (length(ownerPre) != n) {
        stop("'ownerPre' must have one element per product or be an n-by-n matrix.")
    }
    invisible(TRUE)
}


.blp_alpha_domain <- function(alpha, output) {
    if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha)) {
        stop("'alphaMean' must be a finite scalar.")
    }
    if (isTRUE(output) && alpha >= 0) {
        stop("output-market BLP 'alphaMean' must be negative.")
    }
    if (!isTRUE(output) && alpha <= 0) {
        stop("input-market BLP 'alphaMean' must be positive.")
    }
    invisible(TRUE)
}


.blp_parameters_validate <- function(parameters, output = TRUE) {
    alpha <- if (!is.null(parameters$alphaMean)) parameters$alphaMean else
        if (!is.null(parameters$alpha)) parameters$alpha else parameters$alpha_mean
    if (is.null(alpha)) stop("BLP parameters must include 'alphaMean' (or 'alpha').")
    .blp_alpha_domain(alpha, output)
    sigma <- parameters$sigma
    if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma < 0) {
        stop("BLP 'sigma' must be a finite scalar greater than or equal to zero.")
    }
    list(alpha = as.numeric(alpha), sigma = as.numeric(sigma))
}


.blp_model <- function(conduct, prices, shares, margins, ownerPre,
                       alphaMean, sigma, meanval, draws, drawWeights,
                       s0, output = TRUE, priceOutside = 0,
                       insideSize = 1, labels = NULL, bargpowerPre = NULL,
                       bargpowerPost = NULL, weights = NULL) {
    n <- length(prices)
    labels <- if (is.null(labels)) paste0("Prod", seq_len(n)) else labels
    if (is.null(weights)) weights <- rep(1, n)
    normIndex <- if (s0 == 0) 1L else NA_integer_
    slopes <- list(
        alpha = as.numeric(alphaMean), alphaMean = as.numeric(alphaMean),
        sigma = as.numeric(sigma), meanval = as.numeric(meanval),
        sigmaNest = 1, piDemog = numeric(0), nDemog = 0,
        consDraws = as.numeric(draws), demogDraws = matrix(numeric(0), nrow = length(draws), ncol = 0),
        alphas = as.numeric(alphaMean + sigma * draws),
        drawWeights = as.numeric(drawWeights)
    )
    class_name <- switch(conduct,
        bertrand = "LogitBLP",
        cournot = "CournotBLP",
        auction2nd = "Auction2ndBLP",
        bargaining = "BargainingBLP",
        stop("unsupported BLP conduct: ", conduct)
    )
    args <- list(
        prices = as.numeric(prices), shares = as.numeric(shares),
        margins = as.numeric(margins), normIndex = normIndex,
        shareInside = 1 - s0, weights = as.numeric(weights),
        priceOutside = as.numeric(priceOutside), insideSize = as.numeric(insideSize),
        mktSize = 1, priceStart = as.numeric(prices),
        mcDelta = rep(0, n), subset = rep(TRUE, n),
        ownerPre = ownerPre, ownerPost = ownerPre,
        pricePre = as.numeric(prices), pricePost = as.numeric(prices),
        mcPre = rep(NA_real_, n), mcPost = rep(NA_real_, n),
        output = output, labels = as.character(labels), slopes = slopes,
        nDraws = length(draws)
    )
    if (conduct == "bargaining") {
        args$bargpowerPre <- as.numeric(bargpowerPre)
        args$bargpowerPost <- as.numeric(if (is.null(bargpowerPost)) bargpowerPre else bargpowerPost)
    }
    result <- do.call(methods::new, c(list(Class = class_name), args))
    result@ownerPre <- ownerToMatrix(result, preMerger = TRUE)
    result@ownerPost <- ownerToMatrix(result, preMerger = FALSE)
    validObject(result)
    result
}


.blp_effective_alpha <- function(object, preMerger = TRUE) {
    shares <- calcShares(object, preMerger = preMerger, revenue = FALSE)
    derivative <- elast(object, preMerger = preMerger, partial = TRUE)
    denominator <- shares * (1 - shares)
    result <- diag(derivative) / denominator
    if (any(!is.finite(result))) {
        stop("BLP demand does not provide finite local price slopes for the conduct equations.")
    }
    result
}


# The existing second-score equation is retained; the scalar homogeneous
# slope is replaced by the BLP product-local slope.  This is exactly the
# homogeneous equation at sigma = 0 and makes heterogeneity enter through the
# same demand derivatives used by the rest of the BLP implementation.
setMethod(
    f = "calcMargins", signature = "Auction2ndBLP",
    definition = function(object, preMerger = TRUE, exAnte = FALSE, level = TRUE) {
        output <- ifelse(object@output, 1, -1)
        n <- length(object@shares)
        if (preMerger) {
            owner <- object@ownerPre
            subset <- rep(TRUE, n)
            prices <- object@pricePre
        } else {
            owner <- object@ownerPost
            subset <- object@subset
            prices <- object@pricePost
        }
        owner <- owner[subset, subset]
        shares <- calcShares(object, preMerger = preMerger, revenue = FALSE)[subset]
        firmShares <- drop(owner %*% shares)
        alpha <- .blp_effective_alpha(object, preMerger)[subset]
        margins <- rep(NA_real_, n)
        margins[subset] <- output * log(1 - firmShares) / (alpha * firmShares)
        if (exAnte) margins[subset] <- margins[subset] * shares
        if (level) margins <- margins * prices
        names(margins) <- object@labels
        as.vector(margins)
    }
)


setMethod(
    f = "calcMargins", signature = "BargainingBLP",
    definition = function(object, preMerger = TRUE, level = FALSE) {
        output <- ifelse(object@output, -1, 1)
        if (preMerger) {
            prices <- object@pricePre
            owner <- object@ownerPre
            barg <- object@bargpowerPre
        } else {
            prices <- object@pricePost
            owner <- object@ownerPost
            barg <- object@bargpowerPost
        }
        if (any(barg >= 1)) stop("Bargaining BLP requires bargaining power strictly below one.")
        barg <- barg / (1 - barg)
        shares <- calcShares(object, preMerger, revenue = FALSE)
        div <- shares / (1 - shares)
        alpha <- .blp_effective_alpha(object, preMerger)
        margin_matrix <- -owner * shares
        diag(margin_matrix) <- diag(owner) + diag(margin_matrix)
        inverse_matrix <- try(solve(t(margin_matrix)), silent = TRUE)
        if (inherits(inverse_matrix, "try-error")) inverse_matrix <- MASS::ginv(t(margin_matrix))
        term <- (log(1 - shares) * diag(owner))
        term <- sweep(term, 2, -1 * output * alpha *
                          (barg * div - log(1 - shares)), "/")
        margins <- as.vector(inverse_matrix %*% term)
        if (!level) margins <- margins / prices
        names(margins) <- object@labels
        as.vector(margins)
    }
)


setMethod(
    f = "calcPrices", signature = "BargainingBLP",
    definition = function(object, preMerger = TRUE, isMax = FALSE, subset, ...) {
        n <- length(object@shares)
        if (missing(subset)) subset <- rep(TRUE, n)
        if (!is.logical(subset) || length(subset) != n || !any(subset)) {
            stop("'subset' must be a logical vector with at least one TRUE value.")
        }
        start <- object@priceStart
        if (!preMerger && all(is.finite(object@pricePre))) start <- object@pricePre
        start <- start[subset]
        mc <- if (preMerger) object@mcPre else object@mcPost
        output <- object@output
        foc <- function(price) {
            if (preMerger) object@pricePre[subset] <- price else object@pricePost[subset] <- price
            predicted <- calcMargins(object, preMerger = preMerger, level = TRUE)[subset]
            actual <- if (output) price - mc[subset] else mc[subset] - price
            actual - predicted
        }
        maxit <- as.integer(object@control.equ$maxit)
        if (length(maxit) == 0L || is.na(maxit) || maxit < 1L) maxit <- 150L
        tol <- object@control.equ$tol
        if (length(tol) == 0L || !is.finite(tol)) tol <- 1e-10
        solution <- try(nleqslv::nleqslv(
            start, foc, control = list(ftol = tol, maxit = maxit)
        ), silent = TRUE)
        if (inherits(solution, "try-error") || solution$termcd != 1) {
            solution <- BB::BBsolve(start, foc,
                                    control = list(tol = tol, maxit = maxit), quiet = TRUE)
            prices <- if (!inherits(solution, "try-error") && solution$convergence == 0) {
                solution$par
            } else {
                warning("BargainingBLP price solver may not have fully converged.")
                if (inherits(solution, "try-error")) start else solution$par
            }
        } else {
            prices <- solution$x
        }
        result <- rep(NA_real_, n)
        result[subset] <- prices
        names(result) <- object@labels
        result
    }
)


setMethod(
    f = "calcShares", signature = "Auction2ndBLP",
    definition = function(object, preMerger = TRUE, revenue = FALSE, aggregate = TRUE) {
        methods::selectMethod("calcShares", "LogitBLP")(
            object, preMerger = preMerger, revenue = revenue, aggregate = aggregate
        )
    }
)
setMethod(
    f = "calcShares", signature = "BargainingBLP",
    definition = function(object, preMerger = TRUE, revenue = FALSE, aggregate = TRUE) {
        methods::selectMethod("calcShares", "LogitBLP")(
            object, preMerger = preMerger, revenue = revenue, aggregate = aggregate
        )
    }
)
setMethod(
    f = "elast", signature = "Auction2ndBLP",
    definition = function(object, preMerger = TRUE, market = FALSE, partial = FALSE) {
        methods::selectMethod("elast", "LogitBLP")(
            object, preMerger = preMerger, market = market, partial = partial
        )
    }
)
setMethod(
    f = "elast", signature = "BargainingBLP",
    definition = function(object, preMerger = TRUE, market = FALSE, partial = FALSE) {
        methods::selectMethod("elast", "LogitBLP")(
            object, preMerger = preMerger, market = market, partial = partial
        )
    }
)


.blp_new_model <- function(conduct, prices, shares, margins, ownerPre,
                           alphaMean, sigma, delta, integration, s0,
                           output, dots, bargpowerPre = NULL,
                           bargpowerPost = NULL, weights = NULL) {
    .blp_model(
        conduct = conduct, prices = prices, shares = shares, margins = margins,
        ownerPre = ownerPre, alphaMean = alphaMean, sigma = sigma,
        meanval = delta, draws = integration$draws,
        drawWeights = integration$weights, s0 = s0, output = output,
        priceOutside = if (is.null(dots$priceOutside)) 0 else dots$priceOutside,
        insideSize = if (is.null(dots$insideSize)) 1 else dots$insideSize,
        labels = dots$labels, bargpowerPre = bargpowerPre,
        bargpowerPost = bargpowerPost, weights = weights
    )
}


.blp_objective <- function(par, context, details = FALSE) {
    alpha <- par[1]
    sigma <- par[2]
    if (!is.finite(alpha) || !is.finite(sigma) || sigma < 0) return(1e100)
    if ((context$output && alpha >= 0) || (!context$output && alpha <= 0)) return(1e100)
    contracted <- try(.blp_contract(
        prices = context$prices, shares = context$shares,
        alphaMean = alpha, sigma = sigma, draws = context$integration$draws,
        weights = context$integration$weights, s0 = context$s0,
        priceOutside = context$priceOutside, tol = context$contractionTol,
        maxIter = context$contractionMaxIter
    ), silent = TRUE)
    if (inherits(contracted, "try-error") || !contracted$converged ||
        any(!is.finite(contracted$delta))) return(1e100)
    model <- try(.blp_new_model(
        conduct = context$conduct, prices = context$prices,
        shares = context$shares, margins = context$margins,
        ownerPre = context$ownerPre, alphaMean = alpha, sigma = sigma,
        delta = contracted$delta, integration = context$integration,
        s0 = context$s0, output = context$output, dots = context$dots,
        bargpowerPre = context$bargpowerPre,
        bargpowerPost = context$bargpowerPost, weights = context$weights
    ), silent = TRUE)
    if (inherits(model, "try-error")) return(1e100)
    predicted <- try(calcMargins(model, preMerger = TRUE, level = FALSE), silent = TRUE)
    if (inherits(predicted, "try-error") || any(!is.finite(predicted[context$moment_index]))) return(1e100)
    residuals <- predicted - context$margins
    residuals[is.na(residuals)] <- 0
    objective <- sum(context$weights[context$moment_index] *
                         residuals[context$moment_index]^2)
    if (!is.finite(objective)) objective <- 1e100
    if (!details) return(objective)
    list(objective = objective, model = model, delta = contracted$delta,
         predicted = predicted, residuals = residuals,
         contraction = contracted)
}


.blp_logit_start <- function(context) {
    constructor <- switch(context$conduct,
        bertrand = "logit", cournot = "logit.cournot",
        auction2nd = "auction2nd.logit", bargaining = "bargaining.logit"
    )
    args <- list(
        prices = context$prices, shares = context$shares,
        margins = context$margins, ownerPre = context$ownerPre,
        ownerPost = context$ownerPre
    )
    if (context$conduct == "bargaining") {
        args$bargpowerPre <- context$bargpowerPre
        args$bargpowerPost <- context$bargpowerPre
    }
    candidate <- try(do.call(.legacy_constructor(constructor), args), silent = TRUE)
    if (!inherits(candidate, "try-error") && methods::is(candidate, "Bertrand") &&
        is.finite(candidate@slopes$alpha)) {
        value <- as.numeric(candidate@slopes$alpha)
        if ((context$output && value < 0) || (!context$output && value > 0)) return(value)
    }
    scale <- median(context$prices[is.finite(context$prices)])
    sign_value <- if (context$output) -1 else 1
    sign_value / max(scale, 1e-6)
}


.blp_profile_function <- function(context, alpha_bounds) {
    function(sigma_grid) {
        if (!is.numeric(sigma_grid) || any(!is.finite(sigma_grid)) ||
            any(sigma_grid < 0)) stop("'sigma_grid' must be non-negative and finite.")
        result <- lapply(as.numeric(sigma_grid), function(sigma) {
            start <- .blp_logit_start(context)
            fit <- optim(
                par = start,
                fn = function(alpha) .blp_objective(c(alpha, sigma), context),
                method = "L-BFGS-B", lower = alpha_bounds[1], upper = alpha_bounds[2]
            )
            c(sigma = sigma, alphaMean = fit$par, objective = fit$value,
              convergence = fit$convergence)
        })
        as.data.frame(do.call(rbind, result), row.names = NULL)
    }
}


.calibrate_blp_fit <- function(spec, prices, shares, margins, ownerPre, s0,
                               dots, calibration_args) {
    if (is.null(shares) || is.null(margins)) {
        stop("BLP calibration requires observed 'shares' and 'margins'.")
    }
    output <- if (is.null(dots$output)) TRUE else dots$output
    .blp_validate_inputs(prices, shares, margins, ownerPre, s0, output)
    alpha_start <- .blp_logit_start(list(
        conduct = spec$conduct, prices = prices, shares = shares,
        margins = margins, ownerPre = ownerPre, output = output,
        bargpowerPre = dots$bargpowerPre
    ))
    if (spec$conduct == "bargaining" && is.null(dots$bargpowerPre)) {
        stop("'bargpowerPre' must be supplied for BLP bargaining calibration.")
    }
    barg_pre <- dots$bargpowerPre
    if (spec$conduct == "bargaining") {
        if (!is.numeric(barg_pre) || length(barg_pre) != length(prices) ||
            any(!is.finite(barg_pre)) || any(barg_pre < 0) || any(barg_pre >= 1)) {
            stop("'bargpowerPre' must be finite, length-k, and in [0, 1) for BLP bargaining calibration.")
        }
    }
    integration <- .blp_integration(dots)
    context <- list(
        conduct = spec$conduct, output = output, prices = as.numeric(prices),
        shares = as.numeric(shares), margins = as.numeric(margins),
        ownerPre = ownerPre, s0 = as.numeric(s0),
        priceOutside = if (is.null(dots$priceOutside)) 0 else dots$priceOutside,
        integration = integration, dots = dots,
        weights = if (is.null(dots$weights)) rep(1, length(prices)) else dots$weights,
        moment_index = which(!is.na(margins)),
        contractionTol = if (is.null(dots$contractionTol)) 1e-10 else dots$contractionTol,
        contractionMaxIter = if (is.null(dots$contractionMaxIter)) 2000L else dots$contractionMaxIter,
        bargpowerPre = barg_pre,
        bargpowerPost = if (is.null(dots$bargpowerPost)) barg_pre else dots$bargpowerPost
    )
    if (!is.numeric(context$weights) || length(context$weights) != length(prices) ||
        any(!is.finite(context$weights)) || any(context$weights < 0) ||
        any(context$weights[context$moment_index] <= 0)) {
        stop("'weights' must be finite and positive for observed BLP margin moments.")
    }
    alpha_scale <- max(abs(alpha_start), 1e-4)
    alpha_bounds <- if (output) c(-max(1e6, 1000 * alpha_scale), -1e-8) else
        c(1e-8, max(1e6, 1000 * alpha_scale))
    sigma_scale <- max(alpha_scale, 1e-4)
    sigma_upper <- max(1e6, 1000 * sigma_scale)
    sigma_starts <- c(0, .10, .25, .50) * alpha_scale
    starts <- expand.grid(
        alpha = alpha_start * c(.5, 1, 2),
        sigma = sigma_starts,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    starts$alpha <- pmin(pmax(starts$alpha, alpha_bounds[1]), alpha_bounds[2])
    starts$sigma <- pmin(starts$sigma, sigma_upper)
    fits <- lapply(seq_len(nrow(starts)), function(i) {
        opt <- try(optim(
            par = c(starts$alpha[i], starts$sigma[i]),
            fn = function(par) .blp_objective(par, context),
            method = "L-BFGS-B", lower = c(alpha_bounds[1], 0),
            upper = c(alpha_bounds[2], sigma_upper),
            control = dots$optimizer_control
        ), silent = TRUE)
        if (inherits(opt, "try-error")) return(NULL)
        list(optim = opt, converged = isTRUE(opt$convergence == 0), value = opt$value)
    })
    valid <- vapply(fits, function(x) !is.null(x) && x$converged && is.finite(x$value), logical(1))
    if (!any(valid)) stop("BLP calibration did not converge from any deterministic starting value.")
    values <- vapply(fits[valid], `[[`, numeric(1), "value")
    best_index <- which(valid)[which.min(values)]
    best <- fits[[best_index]]$optim
    details <- .blp_objective(best$par, context, details = TRUE)
    model <- details$model
    model@mcPre <- calcMC(model, preMerger = TRUE)
    model@mcPost <- calcMC(model, preMerger = FALSE)
    model@pricePre <- as.numeric(prices)
    model@pricePost <- as.numeric(prices)

    starts_report <- data.frame(
        alphaMean = starts$alpha, sigma = starts$sigma,
        objective = vapply(fits, function(x) if (is.null(x)) NA_real_ else x$value, numeric(1)),
        convergence = vapply(fits, function(x) if (is.null(x)) NA_integer_ else x$optim$convergence, integer(1))
    )
    residuals <- details$residuals
    rmse <- sqrt(mean((residuals[context$moment_index])^2))
    wrong_sign <- if (best$par[2] == 0) 0 else if (output) {
        1 - pnorm((0 - best$par[1]) / best$par[2])
    } else {
        pnorm((0 - best$par[1]) / best$par[2])
    }
    profile <- .blp_profile_function(context, alpha_bounds)
    profile_grid <- unique(c(0, best$par[2], sigma_starts,
                             seq(0, max(2 * best$par[2], alpha_scale), length.out = 7)))
    profile_values <- profile(profile_grid)
    diagnostics <- list(
        status = "completed", source = "calibrate", route = "calibrate",
        model_class = class(model)[[1]], calibration_args = calibration_args,
        integration = list(rule = integration$rule, nodes = integration$draws,
                           weights = integration$weights),
        contraction = details$contraction,
        objective = best$value, residuals = residuals,
        weightedRMSE = sqrt(sum(context$weights[context$moment_index] *
                                    residuals[context$moment_index]^2) /
                                sum(context$weights[context$moment_index])),
        maxAbsResidual = max(abs(residuals[context$moment_index])),
        marginMoments = length(context$moment_index), weights = context$weights,
        s0 = s0, alphaMean = best$par[1], sigma = best$par[2],
        wrongSignProbability = wrong_sign,
        sigmaOnBoundary = isTRUE(all.equal(best$par[2], 0)),
        starts = starts_report,
        profile_sigma_grid = profile_values$sigma,
        profile_sigma_values = profile_values$objective,
        profile_sigma = profile,
        optimizer = list(method = "L-BFGS-B", convergence = best$convergence,
                         message = best$message)
    )
    params <- model@slopes
    if (spec$conduct == "bargaining") {
        params$bargpowerPre <- model@bargpowerPre
        params$bargpowerPost <- model@bargpowerPost
    }
    new(
        "AntitrustFit", spec = spec, model = model, parameters = params,
        observed = list(prices = prices, shares = shares, margins = margins,
                        ownerPre = ownerPre, s0 = s0), diagnostics = diagnostics
    )
}


.specify_blp_conduct_fit <- function(spec, prices, parameters, ownerPre,
                                     shares, margins, insideSize, output,
                                     dots, specification_args) {
    if (is.null(shares)) stop("'shares' must be supplied for BLP parameter loading.")
    if (is.null(output)) output <- TRUE
    .blp_parameters_validate(parameters, output)
    alpha <- if (!is.null(parameters$alphaMean)) parameters$alphaMean else
        if (!is.null(parameters$alpha)) parameters$alpha else parameters$alpha_mean
    sigma <- parameters$sigma
    integration <- .blp_integration(c(dots, parameters[intersect(names(parameters), c("draws", "drawWeights", "nDraws"))]))
    s0 <- if (is.null(dots$s0)) 1 - sum(shares) else dots$s0
    .blp_validate_inputs(prices, shares, if (is.null(margins)) rep(1 / length(prices), length(prices)) else margins,
                         ownerPre, s0, output)
    delta <- parameters$meanval
    if (is.null(delta)) {
        contracted <- .blp_contract(
            prices, shares, alpha, sigma, integration$draws,
            integration$weights, s0,
            priceOutside = if (is.null(dots$priceOutside)) 0 else dots$priceOutside
        )
        if (!contracted$converged) stop("BLP supplied-parameter contraction did not converge.")
        delta <- contracted$delta
    } else if (!is.numeric(delta) || length(delta) != length(prices) || any(!is.finite(delta))) {
        stop("BLP 'meanval' must be a finite length-k vector when supplied.")
    }
    barg_pre <- dots$bargpowerPre
    if (spec$conduct == "bargaining" && is.null(barg_pre)) {
        barg_pre <- rep(.5, length(prices))
    }
    model <- .blp_new_model(
        conduct = spec$conduct, prices = prices, shares = shares,
        margins = if (is.null(margins)) rep(1 / length(prices), length(prices)) else margins,
        ownerPre = ownerPre, alphaMean = alpha, sigma = sigma, delta = delta,
        integration = integration, s0 = s0, output = output,
        dots = c(dots, list(insideSize = insideSize)),
        bargpowerPre = barg_pre, bargpowerPost = dots$bargpowerPost,
        weights = dots$weights
    )
    model@mcPre <- calcMC(model, TRUE)
    model@mcPost <- calcMC(model, FALSE)
    model@pricePre <- prices
    model@pricePost <- prices
    params <- model@slopes
    if (spec$conduct == "bargaining") {
        params$bargpowerPre <- model@bargpowerPre
        params$bargpowerPost <- model@bargpowerPost
    }
    new("AntitrustFit", spec = spec, model = model, parameters = params,
        observed = list(prices = prices, shares = shares, margins = margins,
                        ownerPre = ownerPre, s0 = s0),
        diagnostics = list(status = "completed", source = "specified", route = "specify",
                           model_class = class(model)[[1]], specification_args = specification_args,
                           integration = list(rule = integration$rule, nodes = integration$draws,
                                              weights = integration$weights)))
}


profileBLP <- function(fit, sigma_grid) {
    if (!methods::is(fit, "AntitrustFit") || fit@spec$demand != "blp") {
        stop("'fit' must be an AntitrustFit for BLP demand.")
    }
    if (is.null(fit@diagnostics$profile_sigma)) {
        stop("the BLP fit does not retain profile diagnostics.")
    }
    fit@diagnostics$profile_sigma(sigma_grid)
}
