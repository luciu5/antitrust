# Local translations between the flat and nested Logit and CES demand systems.
#
# These helpers intentionally live above the existing S4 economic
# implementations.  They construct target parameters and then use specify()
# to recover the target supply state; they do not alter any legacy calibration
# or equilibrium equations.

.translation_slot <- function(object, name, default = NULL) {
    if (name %in% methods::slotNames(object)) methods::slot(object, name) else default
}

.translation_nest_values <- function(values, nests) {
    values <- as.numeric(values)
    nest_levels <- levels(factor(nests))
    if (length(values) == 1L) values <- rep(values, length(nest_levels))
    if (length(values) != length(nest_levels)) {
        stop("nesting parameters do not match the number of nests")
    }
    if (!is.null(names(values)) && all(nest_levels %in% names(values))) {
        values <- values[nest_levels]
    }
    names(values) <- nest_levels
    values
}

.translation_market_elasticity <- function(model) {
    value <- try(elast(model, preMerger = TRUE, market = TRUE), silent = TRUE)
    if (inherits(value, "try-error")) NA_real_ else as.numeric(value)
}

.translation_source_state <- function(fit) {
    model <- fit@model
    prices <- as.numeric(model@pricePre)
    quantity_shares <- as.numeric(calcShares(model, preMerger = TRUE,
                                             revenue = FALSE))
    revenue_shares <- as.numeric(calcShares(model, preMerger = TRUE,
                                            revenue = TRUE))
    quantities <- as.numeric(calcQuantities(model, preMerger = TRUE))
    revenues <- as.numeric(calcRevenues(model, preMerger = TRUE))
    ## Some calibrated legacy calls intentionally omit market-size inputs.
    ## Their downstream quantity methods then return NA even though the
    ## normalized baseline shares are available.  Preserve that historical
    ## normalized-share convention for translation rather than rejecting an
    ## otherwise valid fitted model.
    if (any(!is.finite(quantities))) {
        inside_size <- .translation_slot(model, "insideSize", NA_real_)
        if (length(inside_size) && is.finite(inside_size)) {
            quantities <- quantity_shares * inside_size
        } else {
            quantities <- quantity_shares
        }
    }
    if (any(!is.finite(revenues))) revenues <- quantities * prices
    if (any(!is.finite(prices)) || any(!is.finite(quantity_shares)) ||
        any(!is.finite(revenue_shares)) || any(!is.finite(quantities)) ||
        any(!is.finite(revenues))) {
        stop("respecify() requires finite baseline prices, shares, and quantities")
    }

    nests <- if ("nests" %in% methods::slotNames(model)) model@nests else NULL
    list(
        model = model,
        prices = prices,
        quantity_shares = quantity_shares,
        revenue_shares = revenue_shares,
        quantities = quantities,
        revenues = revenues,
        has_outside_quantity = sum(quantity_shares) < 1 - 1e-8,
        has_outside_revenue = sum(revenue_shares) < 1 - 1e-8,
        norm_index = .translation_slot(model, "normIndex", NA_integer_),
        price_outside = .translation_slot(model, "priceOutside", NA_real_),
        labels = .translation_slot(model, "labels", names(prices)),
        owner = .translation_slot(model, "ownerPre"),
        nests = nests,
        elasticity = as.matrix(elast(model, preMerger = TRUE)),
        market_elasticity = .translation_market_elasticity(model)
    )
}

.translation_reference <- function(state) {
    index <- as.integer(state$norm_index)[1]
    if (is.na(index) || index < 1L || index > length(state$prices)) 1L else index
}

.translation_logit_meanval <- function(shares, prices, alpha,
                                       has_outside, price_outside, reference) {
    if (any(shares <= 0) || any(shares >= 1)) {
        stop("local demand translation requires strictly positive baseline shares")
    }
    if (has_outside) {
        outside <- 1 - sum(shares)
        if (outside <= 0) stop("the target Logit model has no positive outside share")
        meanval <- log(shares / outside) - alpha * (prices - price_outside)
    } else {
        meanval <- log(shares / shares[reference]) -
            alpha * (prices - prices[reference])
        meanval[reference] <- 0
    }
    names(meanval) <- names(shares)
    meanval
}

.translation_ces_meanval <- function(shares, prices, gamma,
                                     has_outside, price_outside, reference) {
    if (any(shares <= 0) || any(shares >= 1)) {
        stop("local demand translation requires strictly positive baseline shares")
    }
    if (has_outside) {
        outside <- 1 - sum(shares)
        if (outside <= 0 || price_outside <= 0) {
            stop("the target CES model requires a positive outside-good price")
        }
        meanval <- (shares / outside) *
            (prices / price_outside)^(gamma - 1)
    } else {
        meanval <- (shares / shares[reference]) *
            (prices[reference] / prices)^(1 - gamma)
        meanval[reference] <- 1
    }
    names(meanval) <- names(shares)
    meanval
}

.translation_ces_nest_meanval <- function(shares, prices, gamma, sigma,
                                          nests, has_outside, price_outside,
                                          reference) {
    sigma <- .translation_nest_values(sigma, nests)
    nest_index <- as.character(nests)
    sigma_product <- unname(sigma[nest_index])
    if (any(abs(1 - sigma_product) < 1e-8)) {
        stop("local nested Logit/CES translation requires CES nesting parameters below 1")
    }
    nest_shares <- as.numeric(tapply(shares, nests, sum))
    names(nest_shares) <- levels(factor(nests))
    conditional <- shares / nest_shares[nest_index]
    exponent <- (1 - gamma) / (1 - sigma)

    if (has_outside) {
        outside <- 1 - sum(shares)
        if (outside <= 0 || price_outside <= 0) {
            stop("the target nested CES model requires a positive outside-good price")
        }
        nest_raw <- (nest_shares / outside)^(1 / exponent)
    } else {
        reference_nest <- nest_index[reference]
        reference_raw <- prices[reference]^(1 - sigma_product[reference]) /
            conditional[reference]
        reference_between <- reference_raw^exponent[reference_nest]
        scale <- reference_between / nest_shares[reference_nest]
        nest_raw <- (scale * nest_shares)^(1 / exponent)
    }

    meanval <- conditional * nest_raw[nest_index] /
        prices^(1 - sigma_product)
    if (!has_outside) meanval[reference] <- 1
    names(meanval) <- names(shares)
    meanval
}

.translation_logit_nest_meanval <- function(shares, prices, alpha, sigma,
                                            nests, has_outside, reference) {
    sigma <- .translation_nest_values(sigma, nests)
    nest_index <- as.character(nests)
    sigma_product <- unname(sigma[nest_index])
    if (any(sigma_product <= 0) || any(sigma_product > 1)) {
        stop("the target nested Logit nesting parameters must be in (0, 1]")
    }
    nest_shares <- as.numeric(tapply(shares, nests, sum))
    names(nest_shares) <- levels(factor(nests))
    conditional <- shares / nest_shares[nest_index]

    if (has_outside) {
        outside <- 1 - sum(shares)
        if (outside <= 0) stop("the target nested Logit model has no outside share")
        nest_index_value <- (nest_shares / outside)^(1 / sigma)
    } else {
        reference_nest <- nest_index[reference]
        reference_index <- exp(alpha * prices[reference] /
                                   sigma_product[reference]) /
            conditional[reference]
        reference_between <- reference_index^sigma[reference_nest]
        scale <- reference_between / nest_shares[reference_nest]
        nest_index_value <- (scale * nest_shares)^(1 / sigma)
    }

    meanval <- sigma_product * log(conditional *
                                       nest_index_value[nest_index]) -
        alpha * prices
    if (!has_outside) meanval[reference] <- 0
    names(meanval) <- names(shares)
    meanval
}

.translation_parameter_names <- function(demand) {
    if (demand %in% c("logit", "logit_nests")) "alpha" else "gamma"
}

.translation_target_args <- function(state, target, parameters,
                                     target_shares, target_inside_size,
                                     target_price_outside) {
    args <- list(
        demand = target$demand,
        conduct = target$conduct,
        variant = target$variant,
        prices = state$prices,
        parameters = parameters,
        ownerPre = state$owner,
        shares = target_shares,
        insideSize = target_inside_size,
        priceOutside = target_price_outside,
        labels = state$labels
    )
    if (target$demand %in% c("logit_nests", "ces_nests")) {
        args$nests <- state$nests
    }
    args
}

.translation_build_antitrust <- function(state, target, parameters,
                                         target_shares, target_inside_size,
                                         target_price_outside) {
    args <- .translation_target_args(state, target, parameters,
                                     target_shares, target_inside_size,
                                     target_price_outside)
    result <- do.call(specify, args)
    ## Several legacy nested methods use a strict `< 1` test for an outside
    ## option.  Normalize the target's no-outside metadata after the analytic
    ## share construction so round-off does not manufacture an outside good.
    if (sum(target_shares) >= 1 - 1e-8) {
        result@model@shares <- result@model@shares / sum(result@model@shares)
        if ("shareInside" %in% methods::slotNames(result@model)) {
            result@model@shareInside <- 1
        }
    }
    if ("mktSize" %in% methods::slotNames(result@model)) {
        result@model@mktSize <- target_inside_size / sum(target_shares)
    }
    result
}

.translation_distance <- function(source_elasticity, target_fit) {
    target_elasticity <- as.matrix(elast(target_fit@model, preMerger = TRUE))
    difference <- target_elasticity - source_elasticity
    finite <- is.finite(difference)
    if (!any(finite)) return(Inf)
    mean(difference[finite]^2)
}

.translation_optimise_gamma <- function(state, target, target_shares,
                                        target_inside_size, target_price_outside,
                                        sigma, has_outside, reference) {
    sigma <- .translation_nest_values(sigma, state$nests)
    lower <- max(1e-4, max(1 - sigma) + 1e-5)
    upper <- 1 - 1e-5
    if (lower >= upper) {
        stop("no admissible nested CES curvature region was found for this transition")
    }

    objective <- function(gamma) {
        target_sigma <- 1 + (gamma - 1) / sigma
        if (any(!is.finite(target_sigma)) || any(target_sigma <= 0) ||
            any(target_sigma >= 1)) return(1e12)
        meanval <- try(.translation_ces_nest_meanval(
            target_shares, state$prices, gamma, target_sigma, state$nests,
            has_outside, target_price_outside, reference
        ), silent = TRUE)
        if (inherits(meanval, "try-error")) return(1e12)
        parameters <- list(gamma = gamma, meanval = meanval,
                           ## Pass the full vector.  The legacy nested
                           ## constructor expands scalar sigma values before
                           ## constructing parmsStart, which is invalid for
                           ## its constrained validity path.
                           sigma = target_sigma,
                           shareInside = sum(target_shares))
        candidate <- try(suppressWarnings(.translation_build_antitrust(
            state, target, parameters, target_shares, target_inside_size,
            target_price_outside
        )), silent = TRUE)
        if (inherits(candidate, "try-error")) return(1e12)
        .translation_distance(state$elasticity, candidate)
    }

    result <- stats::optimize(objective, c(lower, upper))
    gamma <- result$minimum
    target_sigma <- 1 + (gamma - 1) / sigma
    list(
        gamma = gamma,
        sigma = target_sigma,
        objective = result$objective,
        convergence = TRUE,
        mapping = list(
            formula = "sigmaCES = 1 + (gamma - 1) / lambdaLogit",
            source_lambda = sigma,
            target_gamma = gamma
        )
    )
}

.translation_optimise_nested_logit <- function(state, target, target_shares,
                                               target_inside_size,
                                               target_price_outside, gamma,
                                               sigma, has_outside, reference) {
    sigma <- .translation_nest_values(sigma, state$nests)
    nest_levels <- levels(factor(state$nests))
    counts <- table(factor(state$nests, levels = nest_levels))
    free_nests <- nest_levels[as.numeric(counts) > 1]
    constrained <- isTRUE(state$model@constraint)
    n_free <- if (constrained) 1L else length(free_nests)
    mapped <- rep(.5, length(nest_levels))
    names(mapped) <- nest_levels
    mapped[as.numeric(counts) == 1] <- 1
    valid_mapping <- abs(sigma - 1) > 1e-8
    mapped[valid_mapping] <- (gamma - 1) / (sigma[valid_mapping] - 1)
    mapped[!is.finite(mapped) | mapped <= 1e-4 | mapped >= 1 - 1e-4] <- .5
    mapped[as.numeric(counts) == 1] <- 1
    initial_sigma <- if (constrained) mean(mapped[as.numeric(counts) > 1]) else mapped[free_nests]
    if (!length(initial_sigma) || any(!is.finite(initial_sigma))) initial_sigma <- .5

    weighted_price <- sum(state$prices * target_shares) / sum(target_shares)
    alpha_initial <- if (isTRUE(state$model@output)) -abs(gamma / weighted_price) else abs(gamma / weighted_price)
    alpha_initial <- min(max(alpha_initial, -100), 100)
    alpha_lower <- if (isTRUE(state$model@output)) -1e6 else 1e-5
    alpha_upper <- if (isTRUE(state$model@output)) -1e-5 else 1e6
    start <- c(alpha_initial, initial_sigma)
    lower <- c(alpha_lower, rep(1e-4, n_free))
    upper <- c(alpha_upper, rep(1 - 1e-4, n_free))

    unpack <- function(par) {
        lambda <- mapped
        if (constrained) {
            lambda[as.numeric(counts) > 1] <- par[2]
        } else {
            lambda[free_nests] <- par[-1]
        }
        lambda
    }
    objective <- function(par) {
        lambda <- unpack(par)
        meanval <- try(.translation_logit_nest_meanval(
            target_shares, state$prices, par[1], lambda, state$nests,
            has_outside, reference
        ), silent = TRUE)
        if (inherits(meanval, "try-error")) return(1e12)
        parameters <- list(alpha = par[1], meanval = meanval,
                           sigma = if (constrained) lambda[1] else lambda)
        candidate <- try(suppressWarnings(.translation_build_antitrust(
            state, target, parameters, target_shares, target_inside_size,
            target_price_outside
        )), silent = TRUE)
        if (inherits(candidate, "try-error")) return(1e12)
        .translation_distance(state$elasticity, candidate)
    }

    result <- stats::optim(start, objective, method = "L-BFGS-B",
                           lower = lower, upper = upper)
    lambda <- unpack(result$par)
    list(
        alpha = result$par[1],
        sigma = lambda,
        objective = result$value,
        convergence = result$convergence,
        mapping = list(
            formula = "lambdaLogit = (gamma - 1) / (sigmaCES - 1)",
            source_gamma = gamma,
            source_sigma = sigma,
            initial_lambda = mapped
        )
    )
}

.translate_antitrust_demand <- function(fit, target) {
    state <- .translation_source_state(fit)
    source_demand <- fit@spec$demand
    target_demand <- target$demand
    source_flat <- source_demand %in% c("logit", "ces")
    target_flat <- target_demand %in% c("logit", "ces")
    source_nested <- source_demand %in% c("logit_nests", "ces_nests")
    target_nested <- target_demand %in% c("logit_nests", "ces_nests")
    if (!((source_flat && target_flat) || (source_nested && target_nested))) {
        stop("local demand translation requires flat-to-flat or nested-to-nested demand systems")
    }
    if (source_nested && is.null(state$nests)) {
        stop("nested demand translation requires retained nest assignments")
    }

    to_ces <- target_demand %in% c("ces", "ces_nests")
    target_shares <- if (to_ces) state$revenue_shares else state$quantity_shares
    target_inside_size <- if (to_ces) sum(state$revenues) else sum(state$quantities)
    target_has_outside <- if (to_ces) state$has_outside_revenue else state$has_outside_quantity
    target_price_outside <- if (to_ces) 1 else 0
    reference <- .translation_reference(state)
    source_parameters <- fit@parameters
    source_sigma <- if (source_nested) {
        source_parameters$sigma
    } else NULL
    if (source_nested && is.null(source_sigma)) {
        source_sigma <- .translation_slot(state$model, "slopes")$sigma
    }

    if (source_flat) {
        if (to_ces) {
            objective_gamma <- function(gamma) {
                meanval <- try(.translation_ces_meanval(
                    target_shares, state$prices, gamma, target_has_outside,
                    target_price_outside, reference
                ), silent = TRUE)
                if (inherits(meanval, "try-error")) return(1e12)
                parameters <- list(gamma = gamma, meanval = meanval,
                                   shareInside = sum(target_shares))
                candidate <- try(suppressWarnings(.translation_build_antitrust(
                    state, target, parameters, target_shares, target_inside_size,
                    target_price_outside
                )), silent = TRUE)
                if (inherits(candidate, "try-error")) return(1e12)
                .translation_distance(state$elasticity, candidate)
            }
            bounds <- if (isTRUE(state$model@output)) c(1 + 1e-5, 100) else c(1e-5, 100)
            result <- stats::optimize(objective_gamma, bounds)
            gamma <- result$minimum
            meanval <- .translation_ces_meanval(
                target_shares, state$prices, gamma, target_has_outside,
                target_price_outside, reference
            )
            parameters <- list(gamma = gamma, meanval = meanval,
                               shareInside = sum(target_shares))
            mapping <- list(
                formula = "mu_j proportional to revenue_share_j * price_j^(gamma - 1)",
                optimizer = "stats::optimize over gamma"
            )
            optimizer <- list(convergence = TRUE, objective = result$objective)
        } else {
            source_gamma <- as.numeric(source_parameters$gamma)[1]
            weighted_price <- sum(state$prices * target_shares) / sum(target_shares)
            objective_alpha <- function(alpha) {
                meanval <- try(.translation_logit_meanval(
                    target_shares, state$prices, alpha, target_has_outside,
                    target_price_outside, reference
                ), silent = TRUE)
                if (inherits(meanval, "try-error")) return(1e12)
                parameters <- list(alpha = alpha, meanval = meanval)
                candidate <- try(suppressWarnings(.translation_build_antitrust(
                    state, target, parameters, target_shares, target_inside_size,
                    target_price_outside
                )), silent = TRUE)
                if (inherits(candidate, "try-error")) return(1e12)
                .translation_distance(state$elasticity, candidate)
            }
            sign <- if (isTRUE(state$model@output)) -1 else 1
            initial <- sign * max(abs(source_gamma) / weighted_price, 1e-3)
            bounds <- if (sign < 0) c(-1e6, -1e-5) else c(1e-5, 1e6)
            result <- stats::optimize(objective_alpha, bounds,
                                      tol = .Machine$double.eps^0.25)
            alpha <- result$minimum
            meanval <- .translation_logit_meanval(
                target_shares, state$prices, alpha, target_has_outside,
                target_price_outside, reference
            )
            parameters <- list(alpha = alpha, meanval = meanval)
            mapping <- list(
                formula = "meanval_j = log(quantity_share_j / reference_share) - alpha * price_difference",
                initial_alpha = initial,
                optimizer = "stats::optimize over alpha"
            )
            optimizer <- list(convergence = TRUE, objective = result$objective)
        }
    } else if (to_ces) {
        source_lambda <- .translation_nest_values(source_sigma, state$nests)
        result <- .translation_optimise_gamma(
            state, target, target_shares, target_inside_size,
            target_price_outside, source_lambda, target_has_outside, reference
        )
        meanval <- .translation_ces_nest_meanval(
            target_shares, state$prices, result$gamma, result$sigma,
            state$nests, target_has_outside, target_price_outside, reference
        )
        parameters <- list(
            gamma = result$gamma,
            meanval = meanval,
            sigma = result$sigma,
            shareInside = sum(target_shares)
        )
        mapping <- result$mapping
        optimizer <- list(convergence = result$convergence,
                          objective = result$objective)
    } else {
        source_gamma <- as.numeric(source_parameters$gamma)[1]
        result <- .translation_optimise_nested_logit(
            state, target, target_shares, target_inside_size,
            target_price_outside, source_gamma, source_sigma,
            target_has_outside, reference
        )
        meanval <- .translation_logit_nest_meanval(
            target_shares, state$prices, result$alpha, result$sigma,
            state$nests, target_has_outside, reference
        )
        parameters <- list(
            alpha = result$alpha,
            meanval = meanval,
            sigma = result$sigma
        )
        mapping <- result$mapping
        optimizer <- list(convergence = result$convergence,
                          objective = result$objective)
    }

    target_fit <- .translation_build_antitrust(
        state, target, parameters, target_shares, target_inside_size,
        target_price_outside
    )
    target_elasticity <- as.matrix(elast(target_fit@model, preMerger = TRUE))
    share_check <- as.numeric(calcShares(target_fit@model, TRUE,
                                         revenue = to_ces))
    quantity_check <- as.numeric(calcQuantities(target_fit@model, TRUE))
    share_difference <- share_check - target_shares
    quantity_difference <- quantity_check - state$quantities
    elasticity_difference <- target_elasticity - state$elasticity
    diagnostics <- list(
        source_demand = source_demand,
        target_demand = target_demand,
        transition_type = "local-demand-translation",
        baseline_share_discrepancy = max(abs(share_difference), na.rm = TRUE),
        baseline_quantity_discrepancy = max(abs(quantity_difference), na.rm = TRUE),
        elasticity_rmse = sqrt(mean(elasticity_difference^2, na.rm = TRUE)),
        maximum_absolute_elasticity_difference = max(abs(elasticity_difference), na.rm = TRUE),
        source_market_elasticity = state$market_elasticity,
        target_market_elasticity = .translation_market_elasticity(target_fit@model),
        parameter_mapping = mapping,
        optimizer = optimizer,
        target_parameters = parameters
    )
    list(fit = target_fit, state = state, parameters = parameters,
         shares = target_shares, inside_size = target_inside_size,
         diagnostics = diagnostics)
}
