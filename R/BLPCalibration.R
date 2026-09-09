# BLP calibration through the refactor architecture.
#
# This file intentionally contains only the observed-data calibration adapter
# and conduct-specific BLP subclasses needed by the existing methods.  The
# legacy S4 demand equations remain the source of truth; MonCom uses the
# explicit atomistic integrated own derivative in MonComMethods.R.

#' @rdname BertrandRUM-Classes
#' @export
setClass(
    Class = "Auction2ndBLP",
    contains = "Auction2ndLogit",
    slots = list(nDraws = "numeric"),
    prototype = prototype(nDraws = 5000)
)

#' @rdname Bargaining-Classes
#' @export
setClass(
    Class = "BargainingBLP",
    contains = "BargainingLogit",
    slots = list(nDraws = "numeric"),
    prototype = prototype(nDraws = 5000)
)


.blp_stable_shares <- function(delta, prices, alpha, draws, weights,
                               priceOutside = 0, outside = TRUE,
                               priceUtility = NULL) {
    ## `priceUtility` is an optional precomputed alpha-by-price matrix.  The
    ## public/internal call retains the historical calculation when it is not
    ## supplied; contraction can reuse this fixed component on every draw
    ## iteration.
    utility <- if (is.null(priceUtility)) {
        outer(alpha, prices - priceOutside, "*")
    } else {
        priceUtility
    }
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
                          maxIter = 2000L, initial = NULL, dampFactor = 1,
                          metrics = NULL, phase = NULL) {
    alpha <- alphaMean + sigma * draws
    price_utility <- outer(alpha, prices - priceOutside, "*")
    if (is.environment(metrics)) {
        metrics$contraction_calls <- if (is.null(metrics$contraction_calls)) 1L else
            metrics$contraction_calls + 1L
        if (!is.null(phase)) {
            call_name <- paste0("contraction_calls_", phase)
            metrics[[call_name]] <- if (is.null(metrics[[call_name]])) 1L else
                metrics[[call_name]] + 1L
        }
        if (is.null(metrics$contraction_iterations)) metrics$contraction_iterations <- 0L
        metric_name <- if (is.null(phase)) "contraction_iterations" else
            paste0("contraction_iterations_", phase)
        if (is.null(metrics[[metric_name]])) metrics[[metric_name]] <- 0L
    }
    outside <- s0 > 0
    delta <- if (is.null(initial)) log(shares) else as.numeric(initial)
    if (length(delta) != length(shares) || any(!is.finite(delta))) {
        delta <- log(shares)
    }

    converged <- FALSE
    max_error <- Inf
    for (iter in seq_len(as.integer(maxIter))) {
        predicted <- .blp_stable_shares(
            delta, prices, alpha, draws, weights, priceOutside, outside,
            priceUtility = price_utility
        )$aggregate
        if (any(!is.finite(predicted)) || any(predicted <= 0)) break
        error <- log(shares) - log(predicted)
        max_error <- max(abs(error))
        delta_new <- delta + error
        if (isTRUE(dampFactor == 1)) {
            if (max_error < tol) {
                delta <- delta_new
                converged <- TRUE
                break
            }
            delta <- delta_new
        } else {
            ## The parameterized sim() compatibility path retains the legacy
            ## damped fixed-point stopping rule and iteration count.  The
            ## calibration path keeps the direct contraction above.
            change <- delta_new - delta
            abs_diff <- max(abs(change))
            rel_diff <- max(abs(change / (abs(delta) + 1e-8)))
            if (rel_diff < tol || abs_diff < tol * 1e-2) {
                delta <- delta_new
                converged <- TRUE
                break
            }
            delta <- delta + dampFactor * change
        }
    }

    if (is.environment(metrics)) {
        metric_name <- if (is.null(phase)) "contraction_iterations" else
            paste0("contraction_iterations_", phase)
        iterations <- if (exists("iter")) as.integer(iter) else 0L
        metrics$contraction_iterations <- metrics$contraction_iterations + iterations
        if (!identical(metric_name, "contraction_iterations")) {
            metrics[[metric_name]] <- metrics[[metric_name]] + iterations
        }
        if (!is.null(initial)) {
            warm_name <- if (is.null(phase)) "warm_contractions" else
                paste0("warm_contractions_", phase)
            if (is.null(metrics[[warm_name]])) metrics[[warm_name]] <- 0L
            metrics[[warm_name]] <- metrics[[warm_name]] + 1L
        }
    }

    if (!outside) delta <- delta - delta[1]
    final <- .blp_stable_shares(
        delta, prices, alpha, draws, weights, priceOutside, outside,
        priceUtility = price_utility
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
                       bargpowerPost = NULL, weights = NULL,
                       integrationRule = "provided", validate = TRUE) {
    n <- length(prices)
    labels <- if (is.null(labels)) paste0("Prod", seq_len(n)) else labels
    if (is.null(weights)) weights <- rep(1, n)
    normIndex <- if (s0 == 0) 1L else NA_integer_
    slopes <- list(
        ## Match the legacy LogitBLP slope representation so parameterized
        ## sim() remains interchangeable with .sim_legacy().
        alpha = as.numeric(alphaMean),
        alphaMean = structure(as.numeric(alphaMean), names = "alphaMean"),
        meanval = structure(as.numeric(meanval), names = as.character(labels)),
        sigma = structure(as.numeric(sigma), names = "sigma"),
        sigmaNest = structure(1, names = "sigmaNest"),
        piDemog = numeric(0),
        nDemog = 0,
        alphas = as.numeric(alphaMean + sigma * draws),
        consDraws = as.numeric(draws),
        demogDraws = NULL,
        drawWeights = as.numeric(drawWeights),
        integrationWeights = as.numeric(drawWeights),
        integration = integrationRule,
        nNodes = if (identical(integrationRule, "gauss-hermite")) length(draws) else NULL
    )
    class_name <- switch(conduct,
        bertrand = "LogitBLP",
        moncom = "MonComBLP",
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
        mktSize = as.numeric(insideSize / (1 - s0)), priceStart = as.numeric(prices),
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
    if (isTRUE(validate)) validObject(result)
    result
}


# The second-score margin is an expectation over the heterogeneous buyer
# types.  For a draw r, the legacy Gumbel expression is
#   log(1 - S_Fr) / (alpha_r S_Fr),
# where S_Fr is the probability that the relevant firm wins.  A product's
# ex-ante conditional margin is therefore the integrated numerator divided by
# the integrated winning probability.  This reduces exactly to the legacy
# Auction2ndLogit expression when sigma = 0; it is not an average-alpha
# substitution.
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
        shares_draw <- calcShares(object, preMerger = preMerger,
                                  revenue = FALSE, aggregate = FALSE)
        shares_draw <- shares_draw[subset, , drop = FALSE]
        draw_weights <- .blp_draw_weights(object, ncol(shares_draw))
        alpha <- object@slopes$alphas
        if (length(alpha) != ncol(shares_draw) || any(!is.finite(alpha))) {
            stop("BLP auction demand has invalid price-coefficient integration points.")
        }
        firm_shares_draw <- owner %*% shares_draw
        firm_shares <- drop(owner %*% as.vector(shares_draw %*% draw_weights))
        if (any(firm_shares <= 0 | firm_shares >= 1) ||
            any(firm_shares_draw <= 0 | firm_shares_draw >= 1)) {
            stop("BLP auction winning probabilities must be strictly between zero and one.")
        }
        numerator <- as.vector(
            (log(1 - firm_shares_draw) /
                 matrix(alpha, nrow = nrow(firm_shares_draw),
                        ncol = ncol(firm_shares_draw), byrow = TRUE)) %*%
                draw_weights
        )
        margins <- rep(NA_real_, n)
        margins[subset] <- output * numerator / firm_shares
        if (exAnte) margins[subset] <- margins[subset] *
            as.vector(shares_draw %*% draw_weights)
        if (!level) margins[subset] <- margins[subset] / prices[subset]
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
        active <- if (preMerger) rep(TRUE, length(object@shares)) else object@subset
        prices <- prices[active]
        owner <- owner[active, active, drop = FALSE]
        if (any(barg >= 1)) stop("Bargaining BLP requires bargaining power strictly below one.")
        barg <- barg[active] / (1 - barg[active])
        shares_draw <- calcShares(object, preMerger, revenue = FALSE,
                                  aggregate = FALSE)
        shares_draw <- shares_draw[active, , drop = FALSE]
        draw_weights <- .blp_draw_weights(object, ncol(shares_draw))
        shares <- as.vector(shares_draw %*% draw_weights)
        alpha <- object@slopes$alphas
        if (length(alpha) != ncol(shares_draw) || any(!is.finite(alpha)) ||
            any(alpha == 0)) {
            stop("BLP bargaining demand has invalid price-coefficient integration points.")
        }

        ## Build the aggregate demand Jacobian from the consumer draws.  A
        ## draw-wise inverse followed by averaging is not the Nash bargaining
        ## FOC for an aggregate market.  With S indexed by product x draw,
        ## a = w * alpha gives the same aggregate as the legacy draw loop:
        ## D = diag(S %*% a) - (S * a) %*% t(S).
        a <- draw_weights * alpha
        derivative <- diag(drop(shares_draw %*% a),
                           nrow = nrow(shares_draw), ncol = nrow(shares_draw)) -
            sweep(shares_draw, 2L, a, "*") %*% t(shares_draw)
        buyer_surplus <- drop(log1p(-shares_draw) %*%
                              (draw_weights / alpha))
        if (any(!is.finite(buyer_surplus)) || any(buyer_surplus == 0)) {
            stop("BLP bargaining buyer surplus is not finite under the supplied price-coefficient draws.")
        }

        ## Normalize each price FOC by aggregate demand and revenue.  In the
        ## homogeneous case the own derivative and buyer-surplus terms reduce
        ## exactly to the legacy BargainingLogit formula.  With zero buyer
        ## bargaining power it is the aggregate Bertrand FOC in level-margin
        ## units.
        ## Express the system in level-margin units using the same elasticity
        ## normalization as the legacy Bertrand method.  Let E_ij be the
        ## aggregate elasticity of share i with respect to price j.  The
        ## legacy price FOC is
        ##   diag(1/(p*s)) %*% t(E * owner) %*% diag(s) %*% margin
        ##       = output * diag(owner).
        ## This normalization makes bargpower = 0 exactly the Bertrand
        ## boundary for ownership vectors and fractional ownership matrices.
        aggregate_elast <- derivative * outer(1 / shares, prices)
        revenue <- prices * shares
        margin_matrix <- t(
            diag(1 / revenue) %*%
                (t(aggregate_elast * owner) %*% diag(shares))
        )
        own_normalized <- diag(derivative) / shares
        rhs <- own_normalized /
            (output * (own_normalized - barg * shares / buyer_surplus))
        rhs <- diag(owner) * rhs

        inverse_matrix <- try(solve(t(margin_matrix)), silent = TRUE)
        if (inherits(inverse_matrix, "try-error")) {
            inverse_matrix <- MASS::ginv(t(margin_matrix))
        }
        margins_active <- as.vector(inverse_matrix %*% rhs)
        margins <- rep(NA_real_, length(object@shares))
        margins[active] <- margins_active
        ## `prices` is already restricted to the active products above.
        if (!level) margins[active] <- margins[active] / prices
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
        nprods <- length(object@shares)
        active <- if (preMerger) rep(TRUE, nprods) else object@subset
        alpha <- as.numeric(object@slopes$alphas)
        meanval <- as.numeric(object@slopes$meanval)
        if (length(alpha) < 1L || length(meanval) != nprods ||
            any(!is.finite(alpha)) || any(!is.finite(meanval))) {
            stop("BLP auction demand requires finite mean utilities and price coefficients.")
        }
        draw_weights <- .blp_draw_weights(object, length(alpha))

        ## BLP meanval is the utility index at zero price.  A second-score
        ## allocation depends on values net of marginal cost, rather than on
        ## the posted price.  Recover the auction value at the observed
        ## pre-merger prices, then apply only cost shocks in the counterfactual.
        baseline <- meanval + tcrossprod(
            object@pricePre - object@priceOutside, alpha
        )
        if (!preMerger) {
            mc_delta_out <- if (is.na(object@normIndex)) {
                object@priceOutside
            } else {
                object@mcDelta[object@normIndex]
            }
            baseline <- baseline + tcrossprod(
                object@mcDelta - mc_delta_out, alpha
            )
        }
        baseline[!active, ] <- -Inf

        ## Stable softmax over the auction allocation values.  The outside
        ## option remains available whenever the BLP object has an outside
        ## share (normIndex = NA).
        max_value <- apply(baseline, 2L, max)
        outside <- is.na(object@normIndex)
        if (outside) max_value <- pmax(0, max_value)
        exp_value <- exp(sweep(baseline, 2L, max_value, "-"))
        denominator <- colSums(exp_value)
        if (outside) denominator <- denominator + exp(-max_value)
        shares_draw <- sweep(exp_value, 2L, denominator, "/")
        shares_draw[!active, ] <- NA_real_

        if (aggregate) {
            shares <- as.vector(shares_draw %*% draw_weights)
            if (revenue) {
                prices <- if (preMerger) object@pricePre else object@pricePost
                total_inside <- sum(prices * shares, na.rm = TRUE)
                total <- if (outside) {
                    total_inside + object@priceOutside *
                        (1 - sum(shares, na.rm = TRUE))
                } else {
                    total_inside
                }
                shares <- prices * shares / total
            }
            shares[!active] <- NA_real_
            names(shares) <- object@labels
            return(shares)
        }
        if (revenue) {
            prices <- if (preMerger) object@pricePre else object@pricePost
            total_inside <- colSums(prices * shares_draw, na.rm = TRUE)
            total <- if (outside) {
                total_inside + object@priceOutside *
                    (1 - colSums(shares_draw, na.rm = TRUE))
            } else {
                total_inside
            }
            shares_draw <- sweep(prices * shares_draw, 2L, total, "/")
        }
        rownames(shares_draw) <- object@labels
        shares_draw
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
                           bargpowerPost = NULL, weights = NULL,
                           validate = TRUE) {
    .blp_model(
        conduct = conduct, prices = prices, shares = shares, margins = margins,
        ownerPre = ownerPre, alphaMean = alphaMean, sigma = sigma,
        meanval = delta, draws = integration$draws,
        drawWeights = integration$weights, s0 = s0, output = output,
        priceOutside = if (is.null(dots$priceOutside)) 0 else dots$priceOutside,
        insideSize = if (is.null(dots$insideSize)) 1 else dots$insideSize,
        labels = dots$labels, bargpowerPre = bargpowerPre,
        bargpowerPost = bargpowerPost, weights = weights,
        integrationRule = integration$rule, validate = validate
    )
}


.blp_objective <- function(par, context, details = FALSE, initial = NULL) {
    if (is.environment(context$metrics)) {
        context$metrics$objective_evaluations <- if (is.null(context$metrics$objective_evaluations)) 1L else
            context$metrics$objective_evaluations + 1L
        metric_name <- if (is.null(context$phase)) "objective_evaluations" else
            paste0("objective_evaluations_", context$phase)
        if (is.null(context$metrics[[metric_name]])) context$metrics[[metric_name]] <- 0L
        context$metrics[[metric_name]] <- context$metrics[[metric_name]] + 1L
    }
    alpha <- par[1]
    sigma <- par[2]
    if (!is.finite(alpha) || !is.finite(sigma) || sigma < 0) return(1e100)
    if ((context$output && alpha >= 0) || (!context$output && alpha <= 0)) return(1e100)
    contracted <- try(.blp_contract(
        prices = context$prices, shares = context$shares,
        alphaMean = alpha, sigma = sigma, draws = context$integration$draws,
        weights = context$integration$weights, s0 = context$s0,
        priceOutside = context$priceOutside, tol = context$contractionTol,
        maxIter = context$contractionMaxIter, initial = initial,
        metrics = context$metrics, phase = context$phase
    ), silent = TRUE)
    if (inherits(contracted, "try-error") || !contracted$converged ||
        any(!is.finite(contracted$delta))) return(1e100)
    model <- try(withCallingHandlers(
        .blp_new_model(
            conduct = context$conduct, prices = context$prices,
            shares = context$shares, margins = context$margins,
            ownerPre = context$ownerPre, alphaMean = alpha, sigma = sigma,
            delta = contracted$delta, integration = context$integration,
            s0 = context$s0, output = context$output, dots = context$dots,
            bargpowerPre = context$bargpowerPre,
            bargpowerPost = context$bargpowerPost,
            weights = context$weights, validate = FALSE
        ),
        warning = function(condition) {
            ## The optimizer constructs a pre-merger placeholder with the
            ## same ownership in both slots.  Muffle only that expected
            ## internal validity warning; all economic/numerical warnings
            ## remain observable.
            if (grepl("'ownerPost' and 'ownerPre' are the same",
                      conditionMessage(condition), fixed = TRUE)) {
                invokeRestart("muffleWarning")
            }
        }
    ), silent = TRUE)
    if (inherits(model, "try-error")) return(1e100)
    predicted <- try(calcMargins(model, preMerger = TRUE, level = FALSE), silent = TRUE)
    if (inherits(predicted, "try-error") || any(!is.finite(predicted[context$moment_index]))) return(1e100)
    residuals <- predicted - context$margins
    residuals[is.na(residuals)] <- 0
    objective <- sum(context$weights[context$moment_index] *
                         residuals[context$moment_index]^2)
    if (!is.finite(objective)) objective <- 1e100
    if (!details) {
        ## Preserve a numeric scalar for optim() while allowing the owning
        ## optimizer closure to carry the converged delta to its next call.
        return(structure(objective, blp_delta = contracted$delta))
    }
    list(objective = objective, model = model, delta = contracted$delta,
         predicted = predicted, residuals = residuals,
         contraction = contracted)
}


## An optimizer run owns its contraction warm start.  A warm contraction is
## only an initial value: if it fails to meet the regular stopping rule, retry
## from log(shares), and never return the unconverged warm result.
.blp_objective_runner <- function(context) {
    warm <- NULL
    function(par, details = FALSE) {
        first <- .blp_objective(par, context, details = details, initial = warm)
        invalid <- if (details) {
            !is.list(first) || !is.list(first$contraction) ||
                !isTRUE(first$contraction$converged)
        } else {
            !is.finite(first) || isTRUE(first >= 1e100)
        }
        if (invalid && !is.null(warm)) {
            if (is.environment(context$metrics)) {
                name <- if (is.null(context$phase)) "cold_retries" else
                    paste0("cold_retries_", context$phase)
                if (is.null(context$metrics[[name]])) context$metrics[[name]] <- 0L
                context$metrics[[name]] <- context$metrics[[name]] + 1L
            }
            first <- .blp_objective(par, context, details = details, initial = NULL)
        }
        usable <- if (details) {
            is.list(first) && is.list(first$contraction) &&
                isTRUE(first$contraction$converged) &&
                all(is.finite(first$delta))
        } else {
            is.finite(first) && isTRUE(first < 1e100)
        }
        if (usable) {
            warm <<- if (details) first$delta else attr(first, "blp_delta")
        } else if (is.environment(context$metrics)) {
            name <- if (is.null(context$phase)) "invalid_evaluations" else
                paste0("invalid_evaluations_", context$phase)
            if (is.null(context$metrics[[name]])) context$metrics[[name]] <- 0L
            context$metrics[[name]] <- context$metrics[[name]] + 1L
        }
        first
    }
}


.blp_multistart_decision <- function(convergence_count, objective_values,
                                     artificial_boundary = FALSE,
                                     invalid_evaluations = 0L,
                                     strategy = c("adaptive", "exhaustive")) {
    strategy <- match.arg(strategy)
    if (identical(strategy, "exhaustive")) {
        return(list(fallback = FALSE, reasons = character(),
                    convergence_count = as.integer(convergence_count),
                    objective_agreement = NA,
                    artificial_boundary = isTRUE(artificial_boundary),
                    invalid_evaluations = as.integer(invalid_evaluations)))
    }
    finite_values <- objective_values[is.finite(objective_values)]
    objective_agreement <- if (length(finite_values) >= 2L) {
        max(finite_values) - min(finite_values) <= 1e-10 +
            1e-4 * max(abs(finite_values))
    } else {
        NA
    }
    reasons <- character()
    if (convergence_count < 2L) {
        reasons <- c(reasons, "fewer than two converged pilot starts")
    }
    if (!isTRUE(objective_agreement)) {
        reasons <- c(reasons, "pilot objectives do not agree")
    }
    if (isTRUE(artificial_boundary)) {
        reasons <- c(reasons,
                     "pilot solution contacts an artificial parameter bound")
    }
    if (invalid_evaluations > 0L) {
        reasons <- c(reasons,
                     "pilot encountered an invalid contraction or diagnostic evaluation")
    }
    list(fallback = length(reasons) > 0L, reasons = reasons,
         convergence_count = as.integer(convergence_count),
         objective_agreement = objective_agreement,
         artificial_boundary = isTRUE(artificial_boundary),
         invalid_evaluations = as.integer(invalid_evaluations))
}


.blp_logit_start <- function(context) {
    constructor <- switch(context$conduct,
        bertrand = "logit", moncom = "logit", cournot = "logit.cournot",
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
    candidate <- try(withCallingHandlers(
        do.call(.legacy_constructor(constructor), args),
        warning = function(condition) {
            ## The legacy start-value constructor receives the pre-merger
            ## ownership in both slots.  Its validity warning is expected for
            ## this internal placeholder and would otherwise be repeated for
            ## every profile/objective evaluation.  Other warnings remain
            ## visible to the caller.
            if (grepl("'ownerPost' and 'ownerPre' are the same",
                      conditionMessage(condition), fixed = TRUE)) {
                invokeRestart("muffleWarning")
            }
        }
    ), silent = TRUE)
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
            profile_context <- context
            profile_context$phase <- "profile"
            runner <- .blp_objective_runner(profile_context)
            fit <- try(optim(
                par = start,
                fn = function(alpha) runner(c(alpha, sigma)),
                method = "L-BFGS-B", lower = alpha_bounds[1], upper = alpha_bounds[2]
            ), silent = TRUE)
            if (inherits(fit, "try-error") || !isTRUE(fit$convergence == 0) ||
                !is.finite(fit$value) || fit$value >= 1e100) {
                stop("BLP profile optimization did not converge at sigma = ", sigma, ".")
            }
            c(sigma = sigma, alphaMean = fit$par, objective = fit$value,
              convergence = fit$convergence)
        })
        result <- as.data.frame(do.call(rbind, result), row.names = NULL)
        ## Keep profiling's lightweight measurement available to the
        ## development benchmark without changing its data-frame contract.
        attr(result, "performance") <- list(metrics = as.list(context$metrics))
        result
    }
}


.blp_identification <- function(context, parameters, moment_index) {
    ## A two-parameter BLP fit needs locally independent margin moments for
    ## alphaMean and sigma.  This is a numerical rank diagnostic for the
    ## observed-margin map, not a covariance or standard-error calculation.
    parameters <- as.numeric(parameters)
    if (length(parameters) != 2L || any(!is.finite(parameters))) {
        return(list(status = "unavailable", identified = NA,
                    rank = NA_integer_, singularValues = numeric(0)))
    }
    identification_context <- context
    identification_context$phase <- "identification"
    steps <- pmax(1e-6, abs(parameters) * 1e-5)
    jacobian <- matrix(NA_real_, nrow = length(moment_index), ncol = 2L,
                       dimnames = list(NULL, c("alphaMean", "sigma")))
    for (column in seq_len(2L)) {
        plus <- parameters
        minus <- parameters
        plus[column] <- plus[column] + steps[column]
        minus[column] <- minus[column] - steps[column]
        if (column == 2L && minus[column] < 0) {
            minus[column] <- parameters[column]
        }
        plus_details <- try(.blp_objective(plus, identification_context, details = TRUE),
                            silent = TRUE)
        minus_details <- try(.blp_objective(minus, identification_context, details = TRUE),
                             silent = TRUE)
        plus_ok <- !inherits(plus_details, "try-error") &&
            is.list(plus_details) && is.finite(plus_details$objective)
        minus_ok <- !inherits(minus_details, "try-error") &&
            is.list(minus_details) && is.finite(minus_details$objective)
        if (plus_ok && minus_ok) {
            jacobian[, column] <-
                (plus_details$predicted[moment_index] -
                     minus_details$predicted[moment_index]) /
                (plus[column] - minus[column])
        } else if (plus_ok) {
            jacobian[, column] <-
                (plus_details$predicted[moment_index] -
                     context$margins[moment_index]) / steps[column]
        } else if (minus_ok) {
            jacobian[, column] <-
                (context$margins[moment_index] -
                     minus_details$predicted[moment_index]) /
                steps[column]
        }
    }
    if (any(!is.finite(jacobian))) {
        return(list(status = "unavailable", identified = NA,
                    rank = NA_integer_, singularValues = numeric(0),
                    jacobian = jacobian))
    }
    singular_values <- svd(jacobian, nu = 0L, nv = 0L)$d
    scale <- if (length(singular_values)) max(singular_values) else 0
    tolerance <- max(dim(jacobian)) * sqrt(.Machine$double.eps) *
        max(1, scale)
    rank <- sum(singular_values > tolerance)
    list(status = if (rank < 2L) "unidentified" else "identified",
         identified = isTRUE(rank >= 2L), rank = as.integer(rank),
         singularValues = singular_values, tolerance = tolerance,
         jacobian = jacobian)
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
        bargpowerPost = if (is.null(dots$bargpowerPost)) barg_pre else dots$bargpowerPost,
        metrics = new.env(parent = emptyenv()), phase = "multistart"
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
    optimizer_control <- dots$optimizer_control
    if (is.null(optimizer_control)) {
        ## The calibration objective is deterministic but can be very flat at
        ## a well-fitting point. A tighter function tolerance avoids treating
        ## a nearby local plateau as the best multi-start solution.
        optimizer_control <- list(factr = 1e3, pgtol = 1e-8)
    }
    ## The original grid is retained verbatim. Adaptive mode starts from a
    ## deterministic pilot including the sigma-zero boundary and two moderate
    ## heterogeneous starts. The first argument to `expand.grid()` varies
    ## fastest, so these are retained rows 2, 5, and 8. The larger-heterogeneity
    ## starts remain available whenever the pilot is inconclusive.
    if (!is.null(dots$blp_multistart) || !is.null(dots$exhaustive)) {
        stop("use the documented 'multistart' option; aliases are unsupported.")
    }
    multistart_option <- dots$multistart
    if (is.null(multistart_option)) multistart_option <- "adaptive"
    if (!is.character(multistart_option) || length(multistart_option) != 1L ||
        !multistart_option %in% c("adaptive", "exhaustive")) {
        stop("'multistart' must be either 'adaptive' or 'exhaustive'.")
    }
    strategy <- multistart_option
    n_starts <- nrow(starts)
    pilot_indices <- c(2L, 5L, 8L)
    pilot_indices <- pilot_indices[pilot_indices <= n_starts]
    evaluated_indices <- if (identical(strategy, "exhaustive")) {
        seq_len(n_starts)
    } else {
        pilot_indices
    }
    fits <- vector("list", n_starts)
    multistart_started <- proc.time()[["elapsed"]]
    evaluate_start <- function(i) {
        opt <- try(optim(
            par = c(starts$alpha[i], starts$sigma[i]),
            fn = .blp_objective_runner(context),
            method = "L-BFGS-B", lower = c(alpha_bounds[1], 0),
            upper = c(alpha_bounds[2], sigma_upper),
            control = optimizer_control
        ), silent = TRUE)
        if (inherits(opt, "try-error")) return(NULL)
        list(optim = opt, converged = isTRUE(opt$convergence == 0), value = opt$value)
    }
    fits[evaluated_indices] <- lapply(evaluated_indices, evaluate_start)
    valid <- vapply(fits, function(x) !is.null(x) && x$converged &&
                        is.finite(x$value) && x$value < 1e100, logical(1))

    pilot_convergence_count <- sum(valid[evaluated_indices])
    pilot_values <- vapply(fits[evaluated_indices], function(x) {
        if (is.null(x) || !is.finite(x$value) || x$value >= 1e100) NA_real_ else x$value
    }, numeric(1))
    pilot_values <- pilot_values[is.finite(pilot_values)]
    boundary_tolerance <- 1e-6
    artificial_boundary <- any(vapply(fits[evaluated_indices], function(x) {
        if (is.null(x) || !isTRUE(x$converged) || length(x$optim$par) != 2L) return(FALSE)
        alpha_contact <- any(abs(x$optim$par[1] - alpha_bounds) <= boundary_tolerance *
                                 pmax(1, abs(alpha_bounds)))
        sigma_contact <- abs(x$optim$par[2] - sigma_upper) <= boundary_tolerance *
            max(1, abs(sigma_upper))
        isTRUE(alpha_contact || sigma_contact)
    }, logical(1)))
    invalid_count <- context$metrics$invalid_evaluations_multistart
    if (is.null(invalid_count)) invalid_count <- 0L
    decision <- .blp_multistart_decision(
        convergence_count = pilot_convergence_count,
        objective_values = pilot_values,
        artificial_boundary = artificial_boundary,
        invalid_evaluations = invalid_count,
        strategy = strategy
    )
    fallback_reasons <- decision$reasons
    if (identical(strategy, "adaptive") && decision$fallback) {
        if (length(fallback_reasons)) {
            remaining <- setdiff(seq_len(n_starts), evaluated_indices)
            fits[remaining] <- lapply(remaining, evaluate_start)
            evaluated_indices <- sort(unique(c(evaluated_indices, remaining)))
            valid <- vapply(fits, function(x) !is.null(x) && x$converged &&
                                is.finite(x$value) && x$value < 1e100, logical(1))
        }
    }
    if (!any(valid)) stop("BLP calibration did not converge from any deterministic starting value.")
    multistart_elapsed <- proc.time()[["elapsed"]] - multistart_started
    values <- vapply(fits[valid], `[[`, numeric(1), "value")
    best_index <- which(valid)[which.min(values)]
    best <- fits[[best_index]]$optim
    context$phase <- "final"
    details <- .blp_objective(best$par, context, details = TRUE)
    model <- details$model
    ## Candidate models are not validity-checked inside the outer optimizer;
    ## validate the selected fitted model once, after calibration.
    validObject(model)
    model@mcPre <- calcMC(model, preMerger = TRUE)
    model@mcPost <- calcMC(model, preMerger = FALSE)
    model@pricePre <- as.numeric(prices)
    model@pricePost <- as.numeric(prices)

    starts_report <- data.frame(
        alphaMean = starts$alpha[evaluated_indices], sigma = starts$sigma[evaluated_indices],
        objective = vapply(fits[evaluated_indices], function(x) if (is.null(x)) NA_real_ else x$value, numeric(1)),
        convergence = vapply(fits[evaluated_indices], function(x) if (is.null(x)) NA_integer_ else x$optim$convergence, integer(1))
    )
    residuals <- details$residuals
    rmse <- sqrt(mean((residuals[context$moment_index])^2))
    wrong_sign <- if (best$par[2] == 0) 0 else if (output) {
        1 - pnorm((0 - best$par[1]) / best$par[2])
    } else {
        pnorm((0 - best$par[1]) / best$par[2])
    }
    profile <- .blp_profile_function(context, alpha_bounds)
    identification_started <- proc.time()[["elapsed"]]
    identification <- .blp_identification(
        context, best$par, context$moment_index
    )
    identification_elapsed <- proc.time()[["elapsed"]] - identification_started
    metric_snapshot <- as.list(context$metrics)
    final_convergence_count <- sum(valid[evaluated_indices])
    fallback <- identical(strategy, "adaptive") && length(fallback_reasons) > 0L
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
        ## Observed-margin moments are the pre-merger conduct equations
        ## evaluated at observed prices. Keep this explicit for diagnostics;
        ## the units are those of the calibration objective.
        preMergerFOCResidual = max(abs(residuals[context$moment_index])),
        focResidualUnits = "proportional margin",
        marginMoments = length(context$moment_index), weights = context$weights,
        s0 = s0, alphaMean = best$par[1], sigma = best$par[2],
        wrongSignProbability = wrong_sign,
        sigmaOnBoundary = isTRUE(all.equal(best$par[2], 0)),
        starts = starts_report,
        multistart = list(
            strategy = strategy,
            evaluatedStarts = length(evaluated_indices),
            totalStarts = n_starts,
            pilotStarts = length(pilot_indices),
            pilotConvergedStarts = pilot_convergence_count,
            convergedStarts = final_convergence_count,
            pilotObjectiveAgreement = decision$objective_agreement,
            pilotArtificialBoundaryContact = artificial_boundary,
            pilotInvalidEvaluationCount = invalid_count,
            fallback = fallback,
            fallbackReason = if (length(fallback_reasons)) fallback_reasons else NA_character_
        ),
        ## The default profile is intentionally on demand. These retained
        ## placeholders keep the historical diagnostic names without storing
        ## an unevaluated grid.
        profile_sigma_grid = NULL,
        profile_sigma_values = NULL,
        profile_sigma = profile,
        identification = identification,
        performance = list(
            multistartSeconds = multistart_elapsed,
            identificationSeconds = identification_elapsed,
            objectiveEvaluations = metric_snapshot$objective_evaluations,
            objectiveEvaluationsMultistart = metric_snapshot$objective_evaluations_multistart,
            objectiveEvaluationsIdentification = metric_snapshot$objective_evaluations_identification,
            contractionCalls = metric_snapshot$contraction_calls,
            contractionIterations = metric_snapshot$contraction_iterations,
            metrics = metric_snapshot
        ),
        optimizer = list(method = "L-BFGS-B", convergence = best$convergence,
                         message = best$message)
    )
    params <- model@slopes
    ## Calibration returns scalar structural parameters as plain numeric
    ## values, matching the established fit@parameters contract.  The model
    ## slots themselves retain the legacy names used by parameterized sim().
    params$alphaMean <- unname(params$alphaMean)
    params$sigma <- unname(params$sigma)
    if (spec$conduct == "bargaining") {
        params$bargpowerPre <- model@bargpowerPre
        params$bargpowerPost <- model@bargpowerPost
    }
    new(
        "AntitrustFit", spec = spec, model = model, parameters = params,
        observed = list(prices = prices, shares = shares, margins = margins,
                        ownerPre = ownerPre, s0 = s0),
        diagnostics = c(diagnostics,
                        if (identical(spec$conduct, "moncom"))
                            .moncom_diagnostics(model) else list())
    )
}


.specify_blp_conduct_fit <- function(spec, prices, parameters, ownerPre,
                                     shares, margins, insideSize, output,
                                     dots, specification_args) {
    if (is.null(output)) output <- TRUE
    .blp_parameters_validate(parameters, output)
    alpha <- if (!is.null(parameters$alphaMean)) parameters$alphaMean else
        if (!is.null(parameters$alpha)) parameters$alpha else parameters$alpha_mean
    sigma <- parameters$sigma
    integration_dots <- dots
    for (name in intersect(names(parameters), c(
        "draws", "consDraws", "drawWeights", "integrationWeights",
        "integration", "nNodes", "nDraws"
    ))) {
        if (is.null(integration_dots[[name]])) {
            integration_dots[[name]] <- parameters[[name]]
        }
    }
    integration <- .blp_integration(integration_dots)
    delta <- parameters$meanval
    contraction_messages <- character()
    price_outside <- if (is.null(dots[["priceOutside"]])) 0 else dots[["priceOutside"]]
    requested_s0 <- dots[["s0"]]
    ## With no supplied s0, preserve the historical outside-good default.
    ## An explicit s0 = 0 instead selects the no-outside-good normalization
    ## when reconstructing shares from supplied meanval.
    outside <- if (is.null(requested_s0)) {
        TRUE
    } else {
        is.numeric(requested_s0) && length(requested_s0) == 1L &&
            is.finite(requested_s0) && requested_s0 > 0
    }
    if (is.null(shares) && is.null(delta)) {
        stop("'shares' must be supplied for BLP parameter loading unless 'meanval' is supplied.")
    }
    if (is.null(shares)) {
        ## A supplied mean utility already identifies the observed inside
        ## shares, conditional on the supplied integration rule and outside
        ## good.  Construct that market state directly; forcing callers to
        ## repeat shares here would make specify() needlessly different from
        ## the parameterized simulation boundary.
        if (!is.numeric(delta) || length(delta) != length(prices) ||
            any(!is.finite(delta))) {
            stop("BLP 'meanval' must be a finite length-k vector when supplied.")
        }
        predicted <- .blp_stable_shares(
            delta, prices, alpha + sigma * integration$draws,
            integration$draws, integration$weights,
            priceOutside = price_outside, outside = outside
        )$aggregate
        shares <- as.numeric(predicted)
        names(shares) <- names(delta)
    }
    ## With observed shares, s0 retains the historical share-derived default;
    ## with supplied meanval, shares were just constructed from the same
    ## outside-good normalization and imply this identical value.
    s0 <- if (is.null(dots$s0)) 1 - sum(shares) else dots$s0
    .blp_validate_inputs(prices, shares, if (is.null(margins)) rep(1 / length(prices), length(prices)) else margins,
                         ownerPre, s0, output)
    if (is.null(delta)) {
        contraction_tol <- if (is.null(dots[["contractionTol"]])) {
            1e-10
        } else {
            dots[["contractionTol"]]
        }
        contraction_max_iter <- if (is.null(dots[["contractionMaxIter"]])) {
            1200L
        } else {
            dots[["contractionMaxIter"]]
        }
        contracted <- .blp_contract(
            prices, shares, alpha, sigma, integration$draws,
            integration$weights, s0,
            priceOutside = price_outside,
            tol = contraction_tol, maxIter = contraction_max_iter,
            dampFactor = .5
        )
        if (!contracted$converged) stop("BLP supplied-parameter contraction did not converge.")
        delta <- contracted$delta
        contraction_messages <- c(
            "Note: 'meanval' (delta) not provided for BLP. It will be recovered via BLP contraction from observed shares/prices.",
            paste0("Running BLP contraction (tol=", sprintf("%.0e", contraction_tol),
                   ", maxIter=", contraction_max_iter, ")..."),
            paste0("BLP contraction converged in ", contracted$iterations, " iterations")
        )
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
        dots = c(dots, list(insideSize = insideSize,
                            priceOutside = price_outside)),
        bargpowerPre = barg_pre, bargpowerPost = dots$bargpowerPost,
        weights = dots$weights, validate = TRUE
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
        diagnostics = c(
            list(status = "completed", source = "specified", route = "specify",
                 model_class = class(model)[[1]], specification_args = specification_args,
                 messages = contraction_messages,
                 integration = list(rule = integration$rule, nodes = integration$draws,
                                    weights = integration$weights),
                 s0 = s0,
                 wrongSignProbability = if (sigma == 0) 0 else if (output) {
                     1 - stats::pnorm(-alpha / sigma)
                 } else {
                     stats::pnorm(-alpha / sigma)
                 }),
            if (identical(spec$conduct, "moncom"))
                .moncom_diagnostics(model) else list()
        ))
}


#' Profile the BLP calibration objective over random-coefficient scale
#'
#' For a calibrated price-random-coefficient BLP fit, evaluate the minimum
#' margin-distance objective at each supplied value of \code{sigma},
#' re-optimizing \code{alphaMean}. This is an identification diagnostic, not
#' an econometric standard-error calculation.
#'
#' @param fit An \code{AntitrustFit} produced by BLP calibration.
#' @param sigma_grid A finite, non-negative numeric vector of values at which
#'   to profile the objective.
#' @return A data frame with \code{sigma}, \code{alphaMean}, \code{objective},
#'   and optimizer convergence columns.
#' @export
profileBLP <- function(fit, sigma_grid) {
    if (!methods::is(fit, "AntitrustFit") || fit@spec$demand != "blp") {
        stop("'fit' must be an AntitrustFit for BLP demand.")
    }
    if (is.null(fit@diagnostics$profile_sigma)) {
        stop("the BLP fit does not retain profile diagnostics.")
    }
    fit@diagnostics$profile_sigma(sigma_grid)
}
