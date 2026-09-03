# Deterministic demand-system translations for fitted Antitrust models.
#
# These helpers sit above the existing S4 economic implementations. They
# construct target primitives analytically (or take them explicitly from the
# caller), then use specify() to rebuild the target supply state. They never
# call a margin-based calibration routine and never optimize a target demand
# parameter.

.translation_slot <- function(object, name, default = NULL) {
    if (name %in% methods::slotNames(object)) methods::slot(object, name) else default
}

.translation_nest_values <- function(values, nests) {
    if (is.null(values)) return(NULL)
    supplied_names <- names(values)
    values <- as.numeric(values)
    nest_levels <- levels(factor(nests))
    if (length(values) == 1L) values <- rep(values, length(nest_levels))
    if (length(values) != length(nest_levels)) {
        stop("nesting parameters do not match the number of nests")
    }
    if (!is.null(supplied_names) && all(nest_levels %in% supplied_names)) {
        values <- values[match(nest_levels, supplied_names)]
    }
    names(values) <- nest_levels
    values
}

.translation_validate_nests <- function(nests, n, source = NULL) {
    if (is.null(nests)) {
        if (is.null(source)) stop("target 'nests' must be supplied for this respecify() transition")
        nests <- source
    }
    if (length(nests) != n) {
        stop("target 'nests' must have the same length as the fitted prices")
    }
    nests <- factor(nests, levels = unique(nests))
    if (nlevels(nests) < 2L) {
        stop("target 'nests' must contain at least two nests")
    }
    if (!is.null(source)) {
        source <- factor(source, levels = levels(factor(source)))
        if (!identical(as.character(nests), as.character(source))) {
            stop("nested-to-nested respecification must preserve the source nest assignments")
        }
    }
    nests
}

.translation_validate_sigma <- function(sigma, nests, target) {
    sigma <- .translation_nest_values(sigma, nests)
    if (any(!is.finite(sigma)) || any(sigma <= 0) || any(sigma > 1)) {
        stop("target nesting parameters must be numeric values in (0, 1]")
    }
    if (identical(target, "ces_nests") && any(sigma >= 1 - 1e-10)) {
        stop("target nested CES parameters must be strictly below 1")
    }
    sigma
}

.translation_market_elasticity <- function(model) {
    value <- try(elast(model, preMerger = TRUE, market = TRUE), silent = TRUE)
    if (inherits(value, "try-error")) NA_real_ else as.numeric(value)
}

.translation_source_state <- function(fit) {
    model <- fit@model
    prices <- as.numeric(model@pricePre)
    quantities <- as.numeric(calcQuantities(model, preMerger = TRUE))
    if (fit@spec$demand %in% c("linear", "loglin") &&
        is.numeric(fit@observed$quantities) &&
        length(fit@observed$quantities) == length(prices)) {
        quantities <- as.numeric(fit@observed$quantities)
    }
    quantity_shares <- as.numeric(calcShares(model, preMerger = TRUE,
                                             revenue = FALSE))
    revenue_shares <- as.numeric(calcShares(model, preMerger = TRUE,
                                            revenue = TRUE))
    if (any(!is.finite(quantities)) && is.numeric(fit@observed$quantities)) {
        quantities <- as.numeric(fit@observed$quantities)
    }
    if (any(!is.finite(quantities))) {
        quantities <- quantity_shares / sum(quantity_shares)
    }
    if (any(!is.finite(quantity_shares))) {
        quantity_shares <- quantities / sum(quantities)
    }
    revenues <- quantities * prices
    if (any(!is.finite(revenue_shares))) {
        revenue_shares <- revenues / sum(revenues)
    }
    if (any(!is.finite(prices)) || any(!is.finite(quantity_shares)) ||
        any(!is.finite(revenue_shares)) || any(!is.finite(quantities)) ||
        any(!is.finite(revenues))) {
        stop("respecify() requires finite baseline prices, shares, and quantities")
    }

    nests <- if ("nests" %in% methods::slotNames(model)) model@nests else NULL
    owner <- fit@observed$ownerPre
    if (is.null(owner)) owner <- .translation_slot(model, "ownerPre")
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
        price_outside = .translation_slot(model, "priceOutside", 0),
        labels = .translation_slot(model, "labels", names(prices)),
        owner = owner,
        nests = nests,
        elasticity = as.matrix(elast(model, preMerger = TRUE)),
        market_elasticity = .translation_market_elasticity(model)
    )
}

.translation_parameter <- function(fit, name, default = NULL) {
    parameters <- fit@parameters
    if (!is.null(parameters[[name]])) return(parameters[[name]])
    slopes <- .translation_slot(fit@model, "slopes", list())
    if (is.list(slopes) && !is.null(slopes[[name]])) slopes[[name]] else default
}

.translation_reference <- function(state) {
    index <- as.integer(state$norm_index)[1]
    if (is.na(index) || index < 1L || index > length(state$prices)) 1L else index
}

.translation_logit_meanval <- function(shares, prices, alpha,
                                       has_outside, price_outside, reference) {
    if (any(shares <= 0)) {
        stop("demand translation requires strictly positive baseline shares")
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
    if (any(shares <= 0)) {
        stop("demand translation requires strictly positive baseline shares")
    }
    if (has_outside) {
        outside <- 1 - sum(shares)
        if (outside <= 0 || price_outside <= 0) {
            stop("the target CES model requires a positive outside-good price")
        }
        meanval <- (shares / outside) *
            (price_outside / prices)^(1 - gamma)
    } else {
        meanval <- (shares / shares[reference]) *
            (prices[reference] / prices)^(1 - gamma)
        meanval[reference] <- 1
    }
    names(meanval) <- names(shares)
    meanval
}

.translation_logit_nest_meanval <- function(shares, prices, alpha, sigma,
                                            nests, has_outside, price_outside,
                                            reference) {
    sigma <- .translation_nest_values(sigma, nests)
    nest_index <- as.character(nests)
    nest_sigma <- unname(sigma[nest_index])
    if (any(!is.finite(sigma)) || any(sigma <= 0) || any(sigma > 1)) {
        stop("target nested Logit nesting parameters must be in (0, 1]")
    }
    nest_shares <- as.numeric(tapply(shares, nests, sum))
    names(nest_shares) <- levels(factor(nests))
    conditional <- shares / nest_shares[nest_index]
    raw <- if (has_outside) {
        outside <- 1 - sum(shares)
        if (outside <= 0) stop("the target nested Logit model has no outside share")
        (nest_shares / outside)^(1 / sigma)
    } else {
        reference_nest <- nest_index[reference]
        reference_raw <- exp(alpha * (prices[reference] - price_outside) /
                                 nest_sigma[reference]) /
            conditional[reference]
        scale <- reference_raw^sigma[reference_nest] /
            nest_shares[reference_nest]
        (scale * nest_shares)^(1 / sigma)
    }
    meanval <- nest_sigma * log(conditional * raw[nest_index]) -
        alpha * (prices - price_outside)
    if (!has_outside) meanval[reference] <- 0
    names(meanval) <- names(shares)
    meanval
}

.translation_ces_nest_meanval <- function(shares, prices, gamma, sigma,
                                          nests, has_outside, price_outside,
                                          reference) {
    sigma <- .translation_nest_values(sigma, nests)
    nest_index <- as.character(nests)
    nest_sigma <- unname(sigma[nest_index])
    if (any(!is.finite(sigma)) || any(sigma <= 0) || any(sigma >= 1)) {
        stop("target nested CES parameters must be strictly below 1")
    }
    nest_shares <- as.numeric(tapply(shares, nests, sum))
    names(nest_shares) <- levels(factor(nests))
    conditional <- shares / nest_shares[nest_index]
    exponent <- (1 - gamma) / (1 - sigma)
    if (any(!is.finite(exponent)) || any(abs(exponent) < 1e-10)) {
        stop("target nested CES parameters imply an invalid between-nest exponent")
    }
    raw <- if (has_outside) {
        outside <- 1 - sum(shares)
        if (outside <= 0 || price_outside <= 0) {
            stop("the target nested CES model requires a positive outside-good price")
        }
        (nest_shares / outside)^(1 / exponent)
    } else {
        reference_nest <- nest_index[reference]
        reference_raw <- prices[reference]^(1 - nest_sigma[reference]) /
            conditional[reference]
        scale <- reference_raw^exponent[reference_nest] /
            nest_shares[reference_nest]
        (scale * nest_shares)^(1 / exponent)
    }
    meanval <- conditional * raw[nest_index] /
        prices^(1 - nest_sigma)
    if (!has_outside) meanval[reference] <- 1
    names(meanval) <- names(shares)
    meanval
}

.translation_target_price_outside <- function(state, target, has_outside) {
    value <- as.numeric(state$price_outside)[1]
    if (target %in% c("ces", "ces_nests") && (!is.finite(value) || value <= 0)) {
        ## Logit's historical default is a zero-priced outside option. CES
        ## requires a positive outside price for its revenue normalization;
        ## one is the existing CES normalization in that case.
        return(1)
    }
    if (!has_outside) return(if (is.finite(value) && value >= 0) value else 1)
    if (!is.finite(value) || value < 0) 0 else value
}

.translation_target_args <- function(state, target, parameters,
                                     target_shares, target_inside_size,
                                     target_price_outside) {
    n <- length(state$prices)
    if (target$demand %in% c("linear", "loglin")) {
        placeholder <- matrix(0.5 / max(1, n - 1), nrow = n, ncol = n)
        diag(placeholder) <- -1
        return(list(
            demand = target$demand,
            conduct = target$conduct,
            variant = target$variant,
            prices = state$prices,
            parameters = parameters,
            ownerPre = state$owner,
            quantities = state$quantities,
            margins = rep(0, n),
            labels = state$labels,
            diversions = placeholder
        ))
    }
    args <- list(
        demand = target$demand,
        conduct = target$conduct,
        variant = target$variant,
        prices = state$prices,
        parameters = parameters,
        ownerPre = state$owner,
        shares = target_shares,
        margins = rep(0, n),
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
    ## The target's observed margins are not transition primitives. The
    ## target demand slopes and ownership are already sufficient for calcMC();
    ## leave the constructor's finite placeholder in the historical slot.
    methods::validObject(result@model)
    result
}

.translation_transition_parameters <- function(fit, state, target, transition,
                                               supplied) {
    source <- fit@spec$demand
    target_demand <- target$demand
    n <- length(state$prices)
    source_nested <- source %in% c("logit_nests", "ces_nests")
    target_nested <- target_demand %in% c("logit_nests", "ces_nests")

    missing <- setdiff(transition$required_arguments, names(supplied))
    if (length(missing)) {
        stop("respecify() transition from '", source, "' to '", target_demand,
             "' requires explicit target primitive(s): ", paste(missing, collapse = ", "))
    }

    nests <- NULL
    if (target_nested) {
        nests <- .translation_validate_nests(
            supplied$nests,
            n,
            source = if (source_nested) state$nests else NULL
        )
        state$nests <- nests
    }

    source_alpha <- as.numeric(.translation_parameter(fit, "alpha"))[1]
    source_gamma <- as.numeric(.translation_parameter(fit, "gamma"))[1]
    alpha <- if (target_demand %in% c("logit", "logit_nests")) {
        if (source %in% c("logit", "logit_nests")) source_alpha else supplied$alpha
    } else NULL
    gamma <- if (target_demand %in% c("ces", "ces_nests")) {
        if (source %in% c("ces", "ces_nests")) source_gamma else supplied$gamma
    } else NULL

    source_sigma <- if (source_nested) {
        .translation_nest_values(.translation_parameter(fit, "sigma"), state$nests)
    } else NULL
    sigma <- NULL
    mapping <- list()
    if (target_nested) {
        sigma <- if (target_demand == "logit_nests" && source == "ces_nests") {
            if (!is.finite(source_gamma) || is.null(source_sigma)) {
                stop("source nested CES fit does not contain gamma and sigma needed for translation")
            }
            (1 - source_gamma) / (1 - source_sigma)
        } else if (target_demand == "ces_nests" && source == "logit_nests") {
            if (!is.finite(gamma) || is.null(source_sigma)) {
                stop("source nested Logit fit does not contain lambda needed for translation")
            }
            1 + (gamma - 1) / source_sigma
        } else if (source_nested) {
            source_sigma
        } else {
            supplied$sigma
        }
        sigma <- .translation_validate_sigma(sigma, nests, target_demand)
        ## The legacy constructors use a scalar sigma to represent a common
        ## nest parameter.  Preserve that representation for flat-to-nested
        ## transitions when the caller supplied a scalar; this is important
        ## for markets containing singleton nests, whose historical
        ## `parmsStart` validation counts only non-singleton nests.
        sigma_argument <- if (!source_nested && length(supplied$sigma) == 1L) {
            as.numeric(supplied$sigma)
        } else {
            sigma
        }
        mapping <- list(nests = nests, sigma = sigma)
        if (source == "logit_nests" && target_demand == "ces_nests") {
            mapping$formula <- "sigmaCES = 1 + (gamma - 1) / lambdaLogit"
        } else if (source == "ces_nests" && target_demand == "logit_nests") {
            mapping$formula <- "lambdaLogit = (1 - gamma) / (1 - sigmaCES)"
        }
    }

    shares <- if (target_demand %in% c("ces", "ces_nests")) {
        state$revenue_shares
    } else if (target_demand %in% c("linear", "loglin")) {
        state$quantities / sum(state$quantities)
    } else {
        state$quantity_shares
    }
    target_has_outside <- if (target_demand %in% c("ces", "ces_nests")) {
        state$has_outside_revenue
    } else if (target_demand %in% c("linear", "loglin")) {
        FALSE
    } else {
        state$has_outside_quantity
    }
    ## The legacy nested-Logit constructor tests `shareInside < 1` directly,
    ## rather than using a tolerance.  Normalize no-outside shares before
    ## construction so round-off in a source model cannot accidentally create
    ## an outside option in the target model.
    if (!target_has_outside) {
        shares <- shares / sum(shares)
        shares[length(shares)] <- 1 - sum(shares[-length(shares)])
    }
    target_price_outside <- .translation_target_price_outside(
        state, target_demand, target_has_outside
    )
    reference <- .translation_reference(state)
    target_inside_size <- if (target_demand %in% c("ces", "ces_nests")) {
        sum(state$revenues)
    } else {
        sum(state$quantities)
    }

    if (target_demand == "logit") {
        meanval <- .translation_logit_meanval(
            shares, state$prices, alpha, target_has_outside,
            target_price_outside, reference
        )
        parameters <- list(alpha = alpha, meanval = meanval)
    } else if (target_demand == "ces") {
        meanval <- .translation_ces_meanval(
            shares, state$prices, gamma, target_has_outside,
            target_price_outside, reference
        )
        parameters <- list(gamma = gamma, meanval = meanval,
                           shareInside = sum(shares))
    } else if (target_demand == "logit_nests") {
        meanval <- .translation_logit_nest_meanval(
            shares, state$prices, alpha, sigma, nests,
            target_has_outside, target_price_outside, reference
        )
        parameters <- list(alpha = alpha, sigma = sigma_argument,
                           meanval = meanval)
    } else if (target_demand == "ces_nests") {
        meanval <- .translation_ces_nest_meanval(
            shares, state$prices, gamma, sigma, nests,
            target_has_outside, target_price_outside, reference
        )
        parameters <- list(gamma = gamma, sigma = sigma_argument, meanval = meanval,
                           shareInside = sum(shares))
    } else if (target_demand %in% c("linear", "loglin")) {
        quantities <- state$quantities
        jacobian <- state$elasticity * outer(quantities, 1 / state$prices)
        if (target_demand == "linear") {
            parameters <- list(
                slopes = jacobian,
                intercepts = quantities - as.vector(jacobian %*% state$prices)
            )
        } else {
            parameters <- list(
                slopes = state$elasticity,
                intercepts = as.vector(log(quantities) -
                    state$elasticity %*% log(state$prices))
            )
        }
        mapping <- list(
            source_jacobian = jacobian,
            source_elasticity = state$elasticity
        )
    } else {
        stop("respecify() does not have a deterministic target construction for '",
             target_demand, "'")
    }

    list(
        parameters = parameters,
        nests = nests,
        shares = shares,
        inside_size = target_inside_size,
        price_outside = target_price_outside,
        mapping = mapping,
        discarded = transition$discarded
    )
}

.translation_diagnostics <- function(state, target_fit, target, transition,
                                    translated) {
    target_model <- target_fit@model
    target_q <- as.numeric(calcQuantities(target_model, preMerger = TRUE))
    target_e <- as.matrix(elast(target_model, preMerger = TRUE))
    target_market_elasticity <- .translation_market_elasticity(target_model)
    target_share_kind <- target$demand %in% c("ces", "ces_nests")
    target_shares <- as.numeric(calcShares(
        target_model, preMerger = TRUE, revenue = target_share_kind
    ))
    share_difference <- target_shares - translated$shares
    quantity_difference <- target_q - state$quantities
    price_difference <- as.numeric(target_model@pricePre) - state$prices
    elasticity_difference <- target_e - state$elasticity
    source_jacobian <- state$elasticity * outer(state$quantities, 1 / state$prices)
    target_jacobian <- target_e * outer(target_q, 1 / state$prices)
    finite_e <- is.finite(elasticity_difference)
    finite_j <- is.finite(target_jacobian - source_jacobian)
    list(
        source_demand = state$source_demand,
        target_demand = target$demand,
        transition_kind = transition$kind,
        baseline_price_discrepancy = max(abs(price_difference), na.rm = TRUE),
        baseline_quantity_discrepancy = max(abs(quantity_difference), na.rm = TRUE),
        baseline_share_discrepancy = max(abs(share_difference), na.rm = TRUE),
        required_arguments = transition$required_arguments,
        derived_parameters = translated$parameters,
        discarded_parameters = translated$discarded,
        target_parameter_validity = isTRUE(methods::validObject(target_model, test = TRUE)),
        source_elasticity = state$elasticity,
        target_elasticity = target_e,
        source_market_elasticity = state$market_elasticity,
        target_market_elasticity = target_market_elasticity,
        elasticity_discrepancy = elasticity_difference,
        elasticity_rmse = if (any(finite_e)) sqrt(mean(elasticity_difference[finite_e]^2)) else Inf,
        maximum_absolute_elasticity_difference = if (any(finite_e)) {
            max(abs(elasticity_difference[finite_e]))
        } else Inf,
        source_jacobian = source_jacobian,
        target_jacobian = target_jacobian,
        jacobian_discrepancy = target_jacobian - source_jacobian,
        jacobian_rmse = if (any(finite_j)) sqrt(mean((target_jacobian - source_jacobian)[finite_j]^2)) else Inf,
        parameter_mapping = translated$mapping
    )
}

.translate_antitrust_demand <- function(fit, target, transition, supplied) {
    state <- .translation_source_state(fit)
    state$source_demand <- fit@spec$demand
    translated <- .translation_transition_parameters(
        fit, state, target, transition, supplied
    )
    if (!is.null(translated$nests)) state$nests <- translated$nests
    target_fit <- .translation_build_antitrust(
        state, target, translated$parameters, translated$shares,
        translated$inside_size, translated$price_outside
    )
    translated$fit <- target_fit
    translated$state <- state
    translated$diagnostics <- .translation_diagnostics(
        state, target_fit, target, transition, translated
    )
    translated
}
