## Model-aware synthetic markets and their neutral design are deliberately
## implemented in antitrust. The design layer does not dispatch demand or
## supply models; realization remains model-specific.

.antitrust_synthetic_or <- function(x, y) if (is.null(x)) y else x

.antitrust_synthetic_validate_scalar <- function(x, name, positive = FALSE) {
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
        (positive && x <= 0)) {
        stop("'", name, "' must be a single ",
             if (positive) "strictly positive " else "finite ", "number")
    }
    invisible(x)
}

.antitrust_synthetic_margin_input <- function(demand, conduct, markup,
                                              reference_price) {
    ## antitrust's standard output-market demand calibrators use proportional
    ## margins. Auction and bargaining calibrators use level margins.
    level_conduct <- conduct %in% c("auction2nd", "bargaining", "bargaining2nd")
    if (level_conduct) markup else markup / reference_price
}

.antitrust_synthetic_product_counts <- function(n_firms, n_products) {
    n_products_input <- as.numeric(n_products)
    if (length(n_products_input) != 1L &&
        length(n_products_input) != n_firms) {
        stop("'n_products' must be a positive integer scalar or a vector of length n_firms")
    }
    if (any(!is.finite(n_products_input)) ||
        any(n_products_input != as.integer(n_products_input)) ||
        any(n_products_input < 1)) {
        stop("'n_products' must contain positive integers")
    }
    products_per_firm <- if (length(n_products_input) == 1L) {
        rep(as.integer(n_products_input), n_firms)
    } else {
        as.integer(n_products_input)
    }
    if (sum(products_per_firm) > 2147483646) {
        stop("the number of inside products is too large")
    }
    list(
        input = if (length(n_products_input) == 1L) {
            as.integer(n_products_input)
        } else {
            products_per_firm
        },
        products_per_firm = products_per_firm,
        n_inside = as.integer(sum(products_per_firm)),
        n_total = as.integer(sum(products_per_firm) + 1L)
    )
}

.antitrust_synthetic_complete_parameters <- function(spec, parameters, shares,
                                                     prices, reference_product,
                                                     output = TRUE) {
    if (!is.list(parameters) || is.null(names(parameters)) ||
        any(!nzchar(names(parameters)))) {
        stop("'parameters' must be a named list in primitives mode")
    }
    result <- parameters
    if (spec$demand %in% c("logit", "logit_nests", "logit_cap")) {
        alpha <- .antitrust_synthetic_or(
            result$alpha,
            .antitrust_synthetic_or(result$alphaMean, result$alpha_mean)
        )
        alpha_valid <- isTRUE(output) && alpha < 0 ||
            !isTRUE(output) && alpha > 0
        if (is.null(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
            !alpha_valid) {
            stop("primitives mode requires a finite ",
                 if (isTRUE(output)) "negative" else "positive", " 'alpha' for ",
                 spec$demand)
        }
        result$alpha <- as.numeric(alpha)
    }

    if (spec$demand == "logit" && is.null(result$meanval)) {
        result$meanval <- log(shares / shares[reference_product]) -
            result$alpha * (prices - prices[reference_product])
        result$meanval[reference_product] <- 0
    }
    if (spec$demand == "ces" && is.null(result$meanval)) {
        gamma <- result$gamma
        if (is.null(gamma) || length(gamma) != 1L || !is.finite(gamma)) {
            stop("primitives mode requires 'gamma' to derive CES meanval")
        }
        result$meanval <- (shares * prices^gamma) /
            (shares[reference_product] * prices[reference_product]^gamma)
        result$meanval[reference_product] <- 1
    }
    if (spec$demand == "logit_nests" && is.null(result$meanval)) {
        stop("primitives mode requires 'meanval' for logit_nests; ",
             "nest-specific mean utilities are not inferred")
    }
    if (spec$demand == "ces_nests" && is.null(result$meanval)) {
        stop("primitives mode requires 'meanval' for ces_nests; ",
             "nest-specific mean values are not inferred")
    }
    if (spec$demand %in% c("logit", "logit_nests") &&
        !is.null(result$meanval)) {
        if (length(result$meanval) != length(shares) ||
            any(!is.finite(result$meanval))) {
            stop("'meanval' must be a finite vector with one element per product")
        }
        result$meanval <- result$meanval - result$meanval[reference_product]
    }
    if (spec$demand %in% c("ces", "ces_nests") &&
        !is.null(result$meanval)) {
        if (length(result$meanval) != length(shares) ||
            any(!is.finite(result$meanval)) ||
            result$meanval[reference_product] == 0) {
            stop("'meanval' must be finite and non-zero at the reference product")
        }
        result$meanval <- result$meanval / result$meanval[reference_product]
    }
    result
}

.antitrust_synthetic_oracle_sign <- function(model) {
    if (!"output" %in% methods::slotNames(model) ||
        length(model@output) != 1L || is.na(model@output)) return(NA_real_)
    if (isTRUE(model@output)) 1 else -1
}

.antitrust_synthetic_oracle_parameters <- function(model, truth) {
    slopes <- if ("slopes" %in% methods::slotNames(model) &&
                  is.list(model@slopes)) model@slopes else list()
    if (!is.list(truth)) truth <- list()
    generated <- if (is.list(truth$generated_parameters)) {
        truth$generated_parameters
    } else list()
    get_one <- function(...) {
        values <- list(...)
        for (value in values) {
            if (!is.null(value) && length(value) == 1L &&
                is.numeric(value) && is.finite(value)) return(as.numeric(value))
        }
        NULL
    }
    list(
        alpha = get_one(truth$alpha, truth$alphaMean, truth$alpha_mean,
                        generated$alpha, generated$alphaMean,
                        slopes$alpha, slopes$alphaMean),
        alphaMean = get_one(truth$alphaMean, truth$alpha, truth$alpha_mean,
                            generated$alphaMean, generated$alpha,
                            slopes$alphaMean, slopes$alpha),
        sigma = get_one(truth$sigma, generated$sigma, slopes$sigma),
        gamma = get_one(truth$gamma, generated$gamma, slopes$gamma),
        meanval = if (!is.null(truth$meanval)) truth$meanval else
            if (!is.null(generated$meanval)) generated$meanval else
                if (!is.null(slopes$meanval)) as.numeric(slopes$meanval) else NULL,
        draws = if (!is.null(truth$draws)) truth$draws else
            if (!is.null(truth$consDraws)) truth$consDraws else
                if (!is.null(generated$draws)) generated$draws else
                    if (!is.null(slopes$consDraws)) as.numeric(slopes$consDraws) else
                        if (!is.null(slopes$draws)) as.numeric(slopes$draws) else NULL,
        drawWeights = if (!is.null(truth$drawWeights)) truth$drawWeights else
            if (!is.null(truth$integrationWeights)) truth$integrationWeights else
                if (!is.null(generated$drawWeights)) generated$drawWeights else
                    if (!is.null(slopes$drawWeights)) as.numeric(slopes$drawWeights) else
                        if (!is.null(slopes$integrationWeights)) {
                            as.numeric(slopes$integrationWeights)
                        } else NULL
    )
}

.antitrust_synthetic_oracle_blp_draws <- function(prices, meanval,
                                                   alpha_mean, sigma, draws,
                                                   weights, price_outside,
                                                   outside_share) {
    n <- length(prices)
    if (length(meanval) != n || any(!is.finite(meanval)) ||
        length(draws) < 1L || any(!is.finite(draws)) ||
        length(weights) != length(draws) || any(!is.finite(weights)) ||
        any(weights < 0) || sum(weights) <= 0 ||
        length(alpha_mean) != 1L || !is.finite(alpha_mean) ||
        length(sigma) != 1L || !is.finite(sigma) || sigma < 0 ||
        length(price_outside) != 1L || !is.finite(price_outside) ||
        length(outside_share) != 1L || !is.finite(outside_share) ||
        outside_share < 0 || outside_share >= 1) return(NULL)
    weights <- as.numeric(weights / sum(weights))
    alpha <- as.numeric(alpha_mean + sigma * draws)
    if (any(!is.finite(alpha))) return(NULL)

    ## Keep the supplied points and weights exactly. This is deliberately a
    ## small, direct draw calculation rather than calcShares(), so the oracle
    ## does not share the production demand-output path it validates.
    utility <- outer(alpha, prices - price_outside, "*")
    utility <- sweep(utility, 2L, as.numeric(meanval), "+")
    max_utility <- apply(utility, 1L, max)
    outside <- outside_share > 1e-12
    if (outside) max_utility <- pmax(0, max_utility)
    exp_inside <- exp(utility - max_utility)
    denominator <- rowSums(exp_inside)
    if (outside) denominator <- denominator + exp(-max_utility)
    draw_shares <- sweep(exp_inside, 1L, denominator, "/")
    if (any(!is.finite(draw_shares))) return(NULL)
    list(alpha = alpha, weights = weights, draw_shares = t(draw_shares),
         aggregate = as.vector(crossprod(weights, draw_shares)))
}

.antitrust_synthetic_oracle_truth <- function(fit, market, reference_markup,
                                               mode, parameter_truth) {
    model <- fit@model
    spec <- fit@spec
    retained_truth <- if (is.list(market$truth)) market$truth else list()
    prices <- if (!is.null(retained_truth$generated_prices)) {
        as.numeric(retained_truth$generated_prices)
    } else as.numeric(market$prices)
    shares <- if (!is.null(retained_truth$generated_shares)) {
        as.numeric(retained_truth$generated_shares)
    } else as.numeric(market$shares)
    ownership <- if (!is.null(retained_truth$generated_ownership)) {
        as.matrix(retained_truth$generated_ownership)
    } else as.matrix(market$ownership)
    n <- length(shares)
    ref <- market$design$reference_product
    output_sign <- .antitrust_synthetic_oracle_sign(model)
    unavailable <- function(reason) list(
        status = "unavailable", reason = reason, costs = rep(NA_real_, n),
        markup = rep(NA_real_, n), quantity = rep(NA_real_, n),
        derivative = matrix(NA_real_, nrow = n, ncol = n),
        expected_markup = rep(NA_real_, n), parameter = NULL
    )
    if (length(output_sign) != 1L || !is.finite(output_sign)) {
        return(unavailable("model output/input orientation is unavailable"))
    }
    if (length(ref) != 1L || !is.finite(ref) || ref != as.integer(ref) ||
        ref < 1L || ref > n || length(prices) != n ||
        any(!is.finite(prices)) || any(prices <= 0) ||
        length(shares) != n || any(!is.finite(shares)) || any(shares <= 0)) {
        return(unavailable("synthetic market primitives are invalid"))
    }
    ref <- as.integer(ref)
    params <- .antitrust_synthetic_oracle_parameters(model, parameter_truth)
    derivative <- matrix(NA_real_, nrow = n, ncol = n)
    quantity <- rep(NA_real_, n)
    signed_markup <- rep(NA_real_, n)
    parameter <- NULL

    ## The retained standard realization uses the same full system G as the
    ## independent oracle below, but the diagnostic derivative and cost truth
    ## are rebuilt here from the primitives and do not call calcMargins().
    if (identical(spec$demand, "logit") &&
        identical(spec$conduct, "bertrand") &&
        identical(spec$variant, "standard")) {
        if (nrow(ownership) != n || ncol(ownership) != n ||
            any(!is.finite(ownership)) || any(!(ownership %in% c(0, 1)))) {
            return(unavailable("ownership is not a finite zero-one matrix"))
        }
        alpha <- params$alpha
        d_unit <- diag(shares) - tcrossprod(shares)
        G <- t(ownership * d_unit)
        z <- tryCatch(solve(G, shares), error = function(e) NULL)
        if (is.null(z) || any(!is.finite(z))) {
            return(unavailable("ownership-adjusted Logit FOC is singular"))
        }
        if (mode == "observed") {
            if (length(reference_markup) != 1L || !is.finite(reference_markup) ||
                reference_markup <= 0) {
                return(unavailable("observed reference markup is unavailable"))
            }
            alpha <- -z[ref] / (output_sign * reference_markup)
        }
        if (length(alpha) != 1L || !is.finite(alpha) || alpha == 0 ||
            (isTRUE(output_sign == 1) && alpha >= 0) ||
            (isTRUE(output_sign == -1) && alpha <= 0)) {
            return(unavailable("Logit price coefficient is invalid"))
        }
        quantity <- shares
        derivative <- alpha * d_unit
        signed_markup <- -z / alpha
        parameter <- list(alpha = alpha)
    } else if (identical(spec$demand, "logit") &&
               identical(spec$conduct, "moncom") &&
               identical(spec$variant, "standard")) {
        alpha <- params$alpha
        if (mode == "observed") {
            if (length(reference_markup) != 1L || !is.finite(reference_markup) ||
                reference_markup <= 0) {
                return(unavailable("observed reference markup is unavailable"))
            }
            alpha <- -1 / (output_sign * reference_markup)
        }
        if (length(alpha) != 1L || !is.finite(alpha) || alpha == 0 ||
            (isTRUE(output_sign == 1) && alpha >= 0) ||
            (isTRUE(output_sign == -1) && alpha <= 0)) {
            return(unavailable("MonCom Logit price coefficient is invalid"))
        }
        quantity <- shares
        derivative <- diag(alpha * quantity)
        signed_markup <- rep(-1 / alpha, n)
        parameter <- list(alpha = alpha)
    } else if (identical(spec$demand, "ces") &&
               identical(spec$conduct, "moncom") &&
               identical(spec$variant, "standard")) {
        gamma <- params$gamma
        if (mode == "observed") {
            if (length(reference_markup) != 1L || !is.finite(reference_markup) ||
                reference_markup <= 0) {
                return(unavailable("observed reference markup is unavailable"))
            }
            gamma <- prices[ref] / (output_sign * reference_markup)
        }
        if (length(gamma) != 1L || !is.finite(gamma) || gamma == 0 ||
            (isTRUE(output_sign == 1) && gamma <= 0) ||
            (isTRUE(output_sign == -1) && gamma >= 0)) {
            return(unavailable("MonCom CES curvature is invalid"))
        }
        ## CES synthetic shares are the revenue-share primitive.  Quantity
        ## is therefore proportional to revenue share / price; the common
        ## market-size factor cancels from the FOC but is retained explicitly.
        quantity <- shares / prices
        derivative <- diag(-gamma * quantity / prices)
        signed_markup <- prices / gamma
        parameter <- list(gamma = gamma)
    } else if (identical(spec$demand, "blp") &&
               identical(spec$conduct, "moncom") &&
               identical(spec$variant, "standard")) {
        alpha_mean <- params$alphaMean
        sigma <- params$sigma
        draws <- params$draws
        weights <- params$drawWeights
        if (is.null(weights) && !is.null(draws)) {
            weights <- rep(1 / length(draws), length(draws))
        }
        if (is.null(params$meanval) || is.null(draws) || is.null(weights) ||
            length(params$meanval) != n) {
            return(unavailable("MonCom BLP requires stored mean values, points, and weights"))
        }
        price_outside <- if ("priceOutside" %in% methods::slotNames(model)) {
            as.numeric(model@priceOutside)
        } else 0
        outside_share <- max(0, 1 - sum(shares))
        integrated <- .antitrust_synthetic_oracle_blp_draws(
            prices, params$meanval, alpha_mean, sigma, draws, weights,
            price_outside, outside_share
        )
        if (is.null(integrated)) {
            return(unavailable("MonCom BLP integration primitives are invalid"))
        }
        scale <- if ("insideSize" %in% methods::slotNames(model) &&
                     length(model@insideSize) == 1L &&
                     is.finite(model@insideSize) && model@insideSize > 0) {
            as.numeric(model@insideSize)
        } else 1
        quantity <- scale * integrated$aggregate
        direct <- as.vector(integrated$draw_shares %*%
                            (integrated$weights * integrated$alpha))
        derivative <- diag(scale * direct)
        if (any(!is.finite(quantity)) || any(!is.finite(direct)) ||
            any(abs(direct) < 1e-12)) {
            return(unavailable("MonCom BLP has a singular integrated own derivative"))
        }
        signed_markup <- -quantity / direct / scale
        parameter <- list(alphaMean = alpha_mean, sigma = sigma,
                          draws = draws, drawWeights = integrated$weights)
    } else {
        return(unavailable("no independent synthetic oracle is implemented for this model"))
    }

    if (any(!is.finite(quantity)) || any(!is.finite(derivative)) ||
        any(!is.finite(signed_markup))) {
        return(unavailable("synthetic oracle produced non-finite values"))
    }
    costs <- prices - signed_markup
    if (any(!is.finite(costs))) {
        return(unavailable("synthetic oracle produced non-finite costs"))
    }
    list(status = "available", reason = NULL, costs = unname(costs),
         markup = unname(output_sign * signed_markup),
         signed_markup = unname(signed_markup), quantity = unname(quantity),
         derivative = derivative, expected_markup = unname(output_sign * signed_markup),
         parameter = parameter)
}

.antitrust_synthetic_diagnostics <- function(fit, market, reference_markup,
                                             mode, parameter_truth = list()) {
    model <- fit@model
    retained_truth <- if (is.list(market$truth)) market$truth else list()
    shares <- if (!is.null(retained_truth$generated_shares)) {
        as.numeric(retained_truth$generated_shares)
    } else as.numeric(market$shares)
    candidate_prices <- as.numeric(market$prices)
    n <- length(shares)
    price_pre <- if ("pricePre" %in% methods::slotNames(model)) {
        as.numeric(model@pricePre)
    } else rep(NA_real_, n)
    price_residual <- if (length(price_pre) == n && length(candidate_prices) == n &&
                          all(is.finite(price_pre)) && all(is.finite(candidate_prices))) {
        max(abs(price_pre - candidate_prices))
    } else NA_real_
    oracle <- .antitrust_synthetic_oracle_truth(
        fit, market, reference_markup, mode, parameter_truth
    )
    unavailable <- function(reason = oracle$reason) {
        list(
            status = "unavailable", mode = mode, price_residual = price_residual,
            retained_cost = rep(NA_real_, n), actual_markup = rep(NA_real_, n),
            implied_markup = rep(NA_real_, n), reference_markup_error = NA_real_,
            ownership_residual = NA_real_,
            market_ownership_residual = NA_real_,
            quantity = rep(NA_real_, n), derivative = matrix(NA_real_, nrow = n, ncol = n),
            foc = rep(NA_real_, n), foc_residual = NA_real_, foc_tolerance = 1e-8,
            foc_status = "unavailable", equilibrium_check = NA,
            oracle_reason = reason
        )
    }
    if (!identical(oracle$status, "available")) return(unavailable())
    retained_cost <- market$products$cost
    if (is.null(retained_cost) || length(retained_cost) != n ||
        any(!is.finite(retained_cost))) {
        return(unavailable("retained structural costs are unavailable"))
    }
    orientation <- .antitrust_synthetic_oracle_sign(model)
    ## Demand derivatives are derivatives with respect to the model price.
    ## Keep p - c signed here (negative in the input convention); the
    ## orientation is applied only when reporting the positive level markup.
    actual_signed_markup <- price_pre - retained_cost
    foc <- rep(NA_real_, n)
    ownership_residual <- NA_real_
    market_ownership_residual <- NA_real_
    if (length(price_pre) == n && all(is.finite(actual_signed_markup))) {
        if (identical(fit@spec$conduct, "bertrand")) {
            owner <- if ("ownerPre" %in% methods::slotNames(model)) {
                as.matrix(model@ownerPre)
            } else matrix(NA_real_, nrow = n, ncol = n)
            reference_owner <- if (!is.null(retained_truth$generated_ownership)) {
                as.matrix(retained_truth$generated_ownership)
            } else as.matrix(market$ownership)
            market_owner <- as.matrix(market$ownership)
            if (nrow(owner) == n && ncol(owner) == n &&
                nrow(reference_owner) == n && ncol(reference_owner) == n &&
                all(is.finite(owner)) && all(is.finite(reference_owner))) {
                ownership_residual <- max(abs(owner - reference_owner))
            }
            if (nrow(market_owner) == n && ncol(market_owner) == n &&
                nrow(reference_owner) == n && ncol(reference_owner) == n &&
                all(is.finite(market_owner)) && all(is.finite(reference_owner))) {
                market_ownership_residual <- max(abs(market_owner - reference_owner))
            }
            if (nrow(owner) == n && ncol(owner) == n) {
                foc <- as.vector(oracle$quantity +
                                 t(owner * oracle$derivative) %*%
                                 actual_signed_markup)
            }
        } else if (identical(fit@spec$conduct, "moncom")) {
            foc <- as.vector(oracle$quantity +
                             diag(oracle$derivative) * actual_signed_markup)
        }
    }
    foc_residual <- if (all(is.finite(foc))) max(abs(foc)) else NA_real_
    expected_markup <- oracle$expected_markup
    actual_markup <- unname(orientation * actual_signed_markup)
    reference_error <- if (length(reference_markup) == 1L &&
                           is.finite(reference_markup)) {
        actual_markup[market$design$reference_product] - reference_markup
    } else NA_real_
    supported <- is.finite(foc_residual) && is.finite(price_residual) &&
        (!identical(fit@spec$conduct, "bertrand") ||
         (is.finite(ownership_residual) &&
          is.finite(market_ownership_residual)))
    list(
        status = if (supported) "completed" else "unavailable", mode = mode,
        price_residual = price_residual,
        retained_cost = unname(retained_cost),
        actual_markup = actual_markup,
        implied_markup = expected_markup,
        reference_markup_error = unname(reference_error),
        ownership_residual = ownership_residual,
        market_ownership_residual = market_ownership_residual,
        quantity = oracle$quantity, derivative = oracle$derivative,
        foc = unname(foc), foc_residual = foc_residual, foc_tolerance = 1e-8,
        foc_status = if (supported) "verified" else "unavailable",
        equilibrium_check = if (supported) {
            isTRUE(price_residual < 1e-8 && foc_residual < 1e-8 &&
                   (!identical(fit@spec$conduct, "bertrand") ||
                    (ownership_residual < 1e-8 &&
                     market_ownership_residual < 1e-8)))
        } else NA,
        oracle_reason = if (supported) NULL else
            "price or ownership state is unavailable for the independent FOC check"
    )
}

.antitrust_synthetic_parameter_error <- function(truth, recovered) {
    if (!length(truth) || !length(recovered)) return(list())
    common <- intersect(names(truth), names(recovered))
    common <- common[vapply(common, function(name) {
        is.numeric(truth[[name]]) && is.numeric(recovered[[name]]) &&
            length(truth[[name]]) == length(recovered[[name]])
    }, logical(1))]
    stats::setNames(lapply(common, function(name) {
        unname(recovered[[name]] - truth[[name]])
    }), common)
}

.antitrust_synthetic_attach <- function(fit, market, mode, reference_markup,
                                        parameter_truth = list(),
                                        foc_diagnostics = list()) {
    oracle_truth <- .antitrust_synthetic_oracle_truth(
        fit, market, reference_markup, mode, parameter_truth
    )
    if (identical(oracle_truth$status, "available")) {
        ## Retain the independently generated structural cost state on the
        ## synthetic market. Diagnostics always use this state instead of
        ## asking production calcMC()/calcMargins() to recreate its truth.
        market$products$cost <- unname(oracle_truth$costs)
        market$products$markup <- unname(oracle_truth$markup)
        market$products$margin <- unname(oracle_truth$markup / market$prices)
        market$costs <- unname(oracle_truth$costs)
        market$markups <- unname(oracle_truth$markup)
        market$truth$generated_cost <- unname(oracle_truth$costs)
        market$truth$generated_markup <- unname(oracle_truth$markup)
        market$truth$generated_parameters <- oracle_truth$parameter
        market$truth$generated_shares <- unname(market$shares)
        market$truth$generated_prices <- unname(market$prices)
        market$truth$generated_ownership <- market$ownership
        market$diagnostics$oracle_status <- "available"
    } else {
        market$diagnostics$oracle_status <- "unavailable"
        market$diagnostics$oracle_reason <- oracle_truth$reason
    }
    fit@observed$shares <- market$shares
    fit@observed$prices <- market$prices
    fit@observed$ownerPre <- market$products$firm_id
    fit@observed$synthetic_market <- market
    fit@observed$reference_product <- market$design$reference_product
    fit@observed$reference_share <- market$design$outside_share
    fit@observed$reference_price <- market$design$reference_price
    fit@observed$outside_margin <- reference_markup
    fit@observed$markup_units <- "level price difference"
    fit@diagnostics$synthetic_market <- market
    fit@diagnostics$synthetic_truth <- parameter_truth
    fit@diagnostics$synthetic_recovered <- fit@parameters
    fit@diagnostics$synthetic_parameter_error <-
        .antitrust_synthetic_parameter_error(parameter_truth, fit@parameters)
    synthetic_diagnostics <- .antitrust_synthetic_diagnostics(
        fit, market, reference_markup, mode, parameter_truth = parameter_truth
    )
    if (length(foc_diagnostics)) {
        synthetic_diagnostics$foc_rank <- foc_diagnostics$rank
        synthetic_diagnostics$foc_condition_number <-
            foc_diagnostics$condition_number
        synthetic_diagnostics$foc_condition_limit <-
            foc_diagnostics$condition_limit
    }
    fit@diagnostics$synthetic <- synthetic_diagnostics
    fit
}

#' Generate a model-consistent synthetic antitrust market
#'
#' `synthetic_market()` is the antitrust-native entry point for fake markets.
#' It draws only a product-share and ownership design. The active reference
#' product has a real positive price. In observed mode, the selected
#' demand/supply implementation uses the reference-product level markup and
#' all product shares to recover structural demand parameters. In primitives
#' mode, the supplied parameters are passed through [specify()] without hidden
#' recalibration.
#'
#' `n_firms` counts inside firms. The active reference product is an additional
#' one-product firm and is included in `ownerPre`. `n_products` may be a scalar
#' shared by all inside firms or a vector with one count per firm, and
#' `dirichlet_alpha` is product-level.
#'
#' Prices and margins are not independently drawn. `reference_price` is the
#' positive price normalization for the market; if `prices` is supplied, it
#' must contain the complete positive product-price vector and its last value
#' must equal `reference_price`. The selected economic implementation derives
#' the remaining markups and marginal costs.
#'
#' @param demand Demand-system name, with `"logit"` as the default.
#' @param supply Supply/conduct name, with `"bertrand"` as the default.
#' @param mode Either `"observed"` or `"primitives"`.
#' @param n_firms Number of inside firms; the reference firm is additional.
#' @param n_products Number of products per inside firm. A scalar is recycled
#' across firms; a vector must have length `n_firms`.
#' @param dirichlet_alpha Positive product-level Dirichlet parameters, one per
#' inside product. If omitted, all shapes equal one.
#' @param outside_beta Positive Beta shape parameters for the reference share.
#' @param reference_price A positive level price for the reference product.
#' @param prices Optional complete positive price vector. It is not randomly
#'   generated and the final element is the reference-product price.
#' @param outside_margin Optional level reference-product markup in observed
#'   mode. Otherwise the open numerical implementation of `U(0, 100)` is used.
#' @param parameters Named model-specific primitives in primitives mode.
#' @param seed Optional explicit integer seed.
#' @param ... Model-specific arguments forwarded to [calibrate()] or
#'   [specify()], such as `nests`, `capacities`, or solver controls.
#' @return An [AntitrustFit] directly accepted by [simulate()].
#' @export
synthetic_market <- function(
    demand = "logit", supply = "bertrand",
    mode = c("observed", "primitives"),
    n_firms = 3L, n_products = 1L,
    dirichlet_alpha = NULL,
    outside_beta = c(2, 8), reference_price = 100,
    prices = NULL, outside_margin = NULL, parameters = NULL,
    seed = NULL, ...) {
    mode <- match.arg(mode)
    .antitrust_synthetic_validate_scalar(reference_price,
                                         "reference_price", positive = TRUE)
    if (!is.numeric(n_firms) || length(n_firms) != 1L ||
        n_firms != as.integer(n_firms) || n_firms < 1L) {
        stop("'n_firms' must be a positive integer")
    }
    n_firms <- as.integer(n_firms)
    counts <- .antitrust_synthetic_product_counts(n_firms, n_products)
    n <- counts$n_total
    if (!is.null(dirichlet_alpha) &&
        (!is.numeric(dirichlet_alpha) || length(dirichlet_alpha) != n - 1L ||
         any(!is.finite(dirichlet_alpha)) || any(dirichlet_alpha <= 0))) {
        stop("'dirichlet_alpha' must be a finite, strictly positive vector of length ",
             n - 1L)
    }
    if (!is.numeric(outside_beta) || length(outside_beta) != 2L ||
        any(!is.finite(outside_beta)) || any(outside_beta <= 0)) {
        stop("'outside_beta' must be a finite, strictly positive vector of length 2")
    }
    if (!is.null(outside_margin) &&
        (length(outside_margin) != 1L || !is.numeric(outside_margin) ||
         !is.finite(outside_margin) || outside_margin <= 0 ||
         outside_margin >= 100)) {
        stop("'outside_margin' must lie strictly inside the level support (0, 100)")
    }
    if (mode == "primitives" && is.null(parameters)) {
        stop("primitives mode requires a named 'parameters' list")
    }
    if (!is.null(parameters) && !is.list(parameters)) {
        stop("'parameters' must be a list")
    }
    if (!is.null(prices) &&
        (!is.numeric(prices) || length(prices) != n ||
         any(!is.finite(prices)) || any(prices <= 0) ||
         !isTRUE(all.equal(unname(prices[n]), unname(reference_price))))) {
        stop("'prices' must be a finite, strictly positive all-product vector whose reference price equals 'reference_price'")
    }
    dots <- list(...)
    duplicate <- intersect(names(dots), c("prices", "shares", "margins",
                                           "ownerPre", "parameters", "demand",
                                           "supply", "conduct", "variant",
                                           "priceOutside", "labels"))
    if (length(duplicate)) {
        stop("argument(s) supplied more than once: ",
             paste(duplicate, collapse = ", "))
    }

    spec <- model_spec(demand, supply)
    design <- fake_market(
        mode = if (mode == "observed") "observed" else "primitives",
        n_firms = n_firms, n_products = n_products,
        dirichlet_alpha = dirichlet_alpha, outside_beta = outside_beta,
        prices = prices, price_level = reference_price,
        reference_price = reference_price,
        outside_margin = if (mode == "observed") outside_margin else NULL,
        parameters = if (mode == "primitives") parameters else list(),
        seed = seed
    )
    shares <- design$shares
    prices <- design$prices
    owner <- design$products$firm_id
    ref <- design$design$reference_product
    drawn_markup <- design$observed$outside_margin
    markup <- if (mode == "observed") drawn_markup else NA_real_
    if (mode == "observed" && markup >= reference_price) {
        stop("the reference markup implies a non-positive reference cost; use a larger 'reference_price' or supply a smaller 'outside_margin'")
    }
    foc_diagnostics <- list()
    if (identical(spec$demand, "logit") &&
        identical(spec$conduct, "bertrand") &&
        identical(spec$variant, "standard")) {
        foc_diagnostics <- .antitrust_synthetic_foc_solution(
            shares, design$ownership
        )
    }

    if (mode == "observed") {
        margin_input <- rep(NA_real_, n)
        margin_input[ref] <- .antitrust_synthetic_margin_input(
            spec$demand, spec$conduct, markup, reference_price
        )
        if (is.null(dots$control.slopes)) {
            dots$control.slopes <- list(reltol = 1e-12)
        }
        calibration_args <- c(
            list(demand = spec$demand, conduct = spec$conduct,
                 prices = prices, shares = shares, margins = margin_input,
                 ownerPre = owner), dots
        )
        fit <- do.call(calibrate, calibration_args)
        return(.antitrust_synthetic_attach(
            fit, design, mode, markup, parameter_truth = list(),
            foc_diagnostics = foc_diagnostics
        ))
    }

    parameters <- .antitrust_synthetic_complete_parameters(
        spec, parameters, shares, prices, ref,
        output = if (is.null(dots$output)) TRUE else dots$output
    )
    specification_args <- c(
        list(demand = spec$demand, conduct = spec$conduct,
             prices = prices, parameters = parameters, ownerPre = owner,
             priceOutside = reference_price,
             labels = paste0("Prod", seq_len(n))), dots
    )
    if (spec$demand == "blp") specification_args$shares <- shares
    fit <- do.call(specify, specification_args)
    .antitrust_synthetic_attach(
        fit, design, mode, NA_real_, parameter_truth = parameters,
        foc_diagnostics = foc_diagnostics
    )
}
