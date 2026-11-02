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

.antitrust_synthetic_nested_meanval <- function(shares, prices, alpha,
                                                nests, sigma, ref) {
    nest_share <- as.numeric(tapply(shares, nests, sum)[as.character(nests)])
    nest_sigma <- as.numeric(sigma[as.character(nests)])
    utility <- nest_sigma * log(shares / nest_share) + log(nest_share)
    meanval <- utility - alpha * prices
    unname(meanval - meanval[ref])
}

.antitrust_synthetic_nested_ces_meanval <- function(shares, prices, gamma,
                                                    nests, sigma, ref) {
    nest_share <- tapply(shares, nests, sum)
    within <- shares / as.numeric(nest_share[nests])
    reference_nest <- nests[ref]
    reference_A <- prices[ref]^(1 - sigma[reference_nest]) / within[ref]
    log_reference_T <- (1 - gamma) / (1 - sigma[reference_nest]) *
        log(reference_A)
    log_A <- (log_reference_T + log(nest_share / nest_share[reference_nest])) *
        (1 - sigma[names(nest_share)]) / (1 - gamma)
    meanval <- within * exp(log_A[nests] -
                            (1 - sigma[nests]) * log(prices))
    unname(meanval / meanval[ref])
}

.antitrust_synthetic_missing_anchors <- function(spec) {
    if (spec$variant == "alm") {
        return(switch(spec$demand,
            logit_nests = "nest assignments and named nesting curvatures, plus a validated nested-ALM inversion using the passive outside share",
            logit_cap = "capacity limits, binding-status or shadow-value information, and the passive outside share",
            "a validated conduct-specific ALM inversion using the passive outside share"))
    }
    switch(spec$demand,
        aids = "an AIDS slope/diversion structure satisfying adding-up and symmetry restrictions",
        pcaids = "the documented proportional-diversion and market-elasticity conventions",
        pcaids_nests = "a PCAIDS diversion structure, nests, and nesting curvature",
        blp = "fixed heterogeneity dispersion and integration draws/weights for a validated joint inversion",
        logit_cap = "capacities and the binding active set or shadow-value information",
        auction2nd_cap = "capacity limits, seller-cost distribution support, and the binding regime",
        linear = "a relative inverse-demand shape, plant-level cost functions, and any leader state",
        loglin = "a relative elasticity shape, plant-level cost functions, and any leader state",
        "additional model-specific demand and conduct anchors")
}

.antitrust_synthetic_observed_solution <- function(spec, market,
                                                   bargpowerPre = NULL,
                                                   parameters = list(),
                                                   nests = NULL) {
    shares <- market$shares
    costs <- market$observed$costs
    margin <- market$observed$reference_margin
    ref <- market$design$reference_product
    handler <- .model_registry_entry(spec$demand, spec$conduct,
                                     spec$variant)$observed_synthetic
    if (identical(handler, "unsupported")) {
        stop("observed synthetic mode is underidentified or unsupported for ",
             spec$id, "; it requires ",
             .antitrust_synthetic_missing_anchors(spec),
             ". No validated observed inversion is currently implemented")
    }
    if (!is.numeric(costs) || length(costs) != length(shares) ||
        any(!is.finite(costs)) || any(costs <= 0)) {
        stop("observed synthetic mode requires finite positive supplied costs")
    }
    if (length(margin) != 1L || !is.finite(margin) ||
        margin <= 0 || margin >= 1) {
        stop("the reference margin must lie strictly between zero and one")
    }
    reference_price <- costs[ref] / (1 - margin)
    reference_markup <- reference_price - costs[ref]
    condition <- NA_real_
    rank <- NA_integer_
    firm_share <- as.vector(market$ownership %*% shares)
    if (grepl("bargaining", handler, fixed = TRUE)) {
        if (is.null(bargpowerPre)) {
            bargpowerPre <- rep(0.5, length(shares))
        }
        if (!is.numeric(bargpowerPre) || length(bargpowerPre) != length(shares) ||
            any(!is.finite(bargpowerPre)) ||
            any(bargpowerPre <= 0 | bargpowerPre >= 1)) {
            stop(spec$id, " requires an explicit all-product 'bargpowerPre' vector in (0, 1)")
        }
    }
    if (grepl("_alm$", handler)) {
        s0 <- market$observed$passive_outside_share
        if (!is.numeric(s0) || length(s0) != 1L ||
            !is.finite(s0) || s0 <= 0 || s0 >= 1) {
            stop("ALM observed mode requires a passive outside share in (0, 1)")
        }
        actual_shares <- (1 - s0) * shares
        actual_firm_share <- as.vector(market$ownership %*% actual_shares)
        if (handler == "logit_nests_bertrand_alm") {
            sigma <- parameters$sigma
            if (is.null(nests) || length(nests) != length(shares) ||
                anyNA(nests)) {
                stop("nested Logit ALM observed mode requires all-product 'nests'")
            }
            nests <- as.character(nests)
            if (!is.numeric(sigma) || is.null(names(sigma)) ||
                !setequal(names(sigma), unique(nests)) ||
                any(!is.finite(sigma)) || any(sigma <= 0 | sigma > 1)) {
                stop("nested Logit ALM requires named 'parameters$sigma' in (0,1] for every nest")
            }
            singleton <- table(nests) == 1L
            if (any(sigma[names(singleton)[singleton]] != 1)) {
                stop("singleton nests require sigma = 1")
            }
            within <- shares /
                as.numeric(tapply(shares, nests, sum)[nests])
            n <- length(shares)
            derivative <- diag(1 / sigma[nests], n) +
                outer(nests, nests, `==`) *
                matrix(rep((1 - 1 / sigma[nests]) * within,
                           each = n), n, n) -
                matrix(rep(actual_shares, each = n), n, n)
            derivative <- diag(actual_shares, n) %*% derivative
            G <- t(market$ownership * derivative)
            rank <- qr(G)$rank
            condition <- kappa(G, exact = TRUE)
            if (rank < n || !is.finite(condition) ||
                condition > 1 / sqrt(.Machine$double.eps)) {
                stop("nested Logit ALM ownership FOC is singular or ill-conditioned")
            }
            coefficient <- as.vector(solve(G, actual_shares))
            if (any(!is.finite(coefficient)) || any(coefficient <= 0)) {
                stop("nested Logit ALM implies invalid markups")
            }
            alpha <- -coefficient[ref] / reference_markup
            markups <- coefficient / (-alpha)
            prices <- costs + markups
            parameter <- list(alpha = unname(alpha), sigma = sigma,
                              passive_outside_share = s0)
        } else if (startsWith(handler, "logit_")) {
            if (handler == "logit_bertrand_alm") {
                system <- .antitrust_synthetic_foc_solution(
                    actual_shares, market$ownership)
                coefficient <- system$z
                condition <- system$condition_number
                rank <- system$rank
            } else if (handler == "logit_cournot_alm") {
                coefficient <- 1 + actual_firm_share / s0
            } else if (handler == "logit_auction2nd_alm") {
                coefficient <- -log1p(-actual_firm_share) /
                    actual_firm_share
            } else {
                relative_power <- bargpowerPre / (1 - bargpowerPre)
                M <- -market$ownership * actual_shares
                diag(M) <- diag(market$ownership) + diag(M)
                rank <- qr(M)$rank
                condition <- kappa(M, exact = TRUE)
                if (rank < length(shares) || !is.finite(condition) ||
                    condition > 1 / sqrt(.Machine$double.eps)) {
                    stop("ALM bargaining ownership system is singular or ill-conditioned")
                }
                direct <- -log1p(-actual_shares) /
                    (relative_power * actual_shares /
                         (1 - actual_shares) - log1p(-actual_shares))
                coefficient <- as.vector(solve(t(M), direct))
            }
            if (any(!is.finite(coefficient)) || any(coefficient <= 0)) {
                stop("ALM Logit conduct equations imply invalid markups")
            }
            alpha <- -coefficient[ref] / reference_markup
            markups <- coefficient / (-alpha)
            prices <- costs + markups
            parameter <- list(alpha = unname(alpha),
                              passive_outside_share = s0)
        } else {
            if (handler == "ces_cournot_alm") {
                alpha_ces <- s0 / (1 - s0)
                cournot_margin <- function(gamma) {
                    1 / gamma + (gamma - 1) * (1 + alpha_ces) /
                        (gamma * (1 + gamma * alpha_ces)) *
                        actual_firm_share
                }
                upper <- 2
                while (cournot_margin(upper)[ref] > margin && upper < 1e8) {
                    upper <- upper * 2
                }
                if (cournot_margin(upper)[ref] >= margin) {
                    stop("ALM CES Cournot reference margin does not identify finite admissible curvature")
                }
                gamma <- stats::uniroot(function(g) cournot_margin(g)[ref] -
                    margin, c(1 + 1e-8, upper), tol = 1e-12)$root
                all_margin <- cournot_margin(gamma)
            } else if (handler == "ces_auction2nd_alm") {
                gamma <- 1 + log1p(-actual_firm_share[ref]) /
                    log1p(-margin)
                all_margin <- 1 - (1 - actual_firm_share)^(1 / (gamma - 1))
            } else {
                adjusted <- if (handler == "ces_bargaining_alm") {
                    margin / (1 - bargpowerPre[ref])
                } else margin
                gamma <- (1 / adjusted - actual_firm_share[ref]) /
                    (1 - actual_firm_share[ref])
                all_margin <- 1 / (gamma * (1 - actual_firm_share) +
                                  actual_firm_share)
                if (handler == "ces_bargaining_alm") {
                    all_margin <- all_margin * (1 - bargpowerPre)
                }
            }
            if (!is.finite(gamma) || gamma <= 1 ||
                any(!is.finite(all_margin)) ||
                any(all_margin <= 0 | all_margin >= 1)) {
                stop("ALM CES outside share and reference margin imply inadmissible curvature or margins")
            }
            prices <- costs / (1 - all_margin)
            markups <- prices - costs
            parameter <- list(gamma = unname(gamma),
                              passive_outside_share = s0)
        }
    } else if (handler == "pcaids_bertrand") {
        known_elasticity <- -1 / margin
        probe <- calibrate("pcaids", "bertrand", shares = shares,
            prices = rep(reference_price, length(shares)),
            knownElast = known_elasticity, knownElastIndex = ref,
            mktElast = -1, ownerPre = market$products$firm_id)
        all_margin <- as.numeric(calcMargins(probe@model, TRUE))
        if (any(!is.finite(all_margin)) ||
            any(all_margin <= 0 | all_margin >= 1) ||
            abs(all_margin[ref] - margin) > 1e-7) {
            stop("PCAIDS reference margin and proportional-diversion calibration disagree or imply inadmissible product margins")
        }
        prices <- costs / (1 - all_margin)
        markups <- prices - costs
        E <- as.matrix(elast(probe@model, TRUE))
        G <- t(E) * market$ownership
        rank <- qr(G)$rank
        condition <- kappa(G, exact = TRUE)
        if (rank < length(shares) || !is.finite(condition) ||
            condition > 1 / sqrt(.Machine$double.eps)) {
            stop("PCAIDS ownership-adjusted elasticity system is singular or ill-conditioned")
        }
        parameter <- list(knownElast = known_elasticity,
                          mktElast = -1, slopes = probe@model@slopes)
    } else if (handler == "loglin_bertrand") {
        shape <- parameters$slopes
        n <- length(shares)
        if (!is.matrix(shape) || !is.numeric(shape) ||
            !identical(dim(shape), c(n, n)) || any(!is.finite(shape)) ||
            any(diag(shape) >= 0) || any(shape[row(shape) != col(shape)] < 0)) {
            stop("log-linear Bertrand observed mode requires 'parameters$slopes' as a finite all-product relative elasticity matrix with negative diagonal and nonnegative cross effects")
        }
        scale <- -1 / (margin * shape[ref, ref])
        slopes <- scale * shape
        G <- diag(shares) +
            t(market$ownership * slopes) %*% diag(shares)
        rank <- qr(G)$rank
        condition <- kappa(G, exact = TRUE)
        if (rank < n || !is.finite(condition) ||
            condition > 1 / sqrt(.Machine$double.eps)) {
            stop("log-linear Bertrand ownership-adjusted elasticity system is singular or ill-conditioned")
        }
        markups <- -as.vector(solve(G, shares * costs))
        prices <- costs + markups
        if (any(!is.finite(markups)) || any(markups <= 0) ||
            any(!is.finite(prices)) || any(prices <= 0) ||
            abs(markups[ref] - reference_markup) > 1e-8) {
            stop("log-linear elasticity shape implies inadmissible markups")
        }
        intercepts <- log(shares) - as.vector(slopes %*% log(prices))
        parameter <- list(slopes = slopes, intercepts = intercepts)
    } else if (handler == "linear_bertrand") {
        shape <- parameters$slopes
        n <- length(shares)
        if (!is.matrix(shape) || !is.numeric(shape) ||
            !identical(dim(shape), c(n, n)) || any(!is.finite(shape)) ||
            any(diag(shape) >= 0) || any(shape[row(shape) != col(shape)] < 0)) {
            stop("linear Bertrand observed mode requires 'parameters$slopes' as a finite all-product relative slope matrix with negative diagonal and nonnegative cross effects")
        }
        G <- t(market$ownership * shape)
        rank <- qr(G)$rank
        condition <- kappa(G, exact = TRUE)
        if (rank < n || !is.finite(condition) ||
            condition > 1 / sqrt(.Machine$double.eps)) {
            stop("linear Bertrand ownership-adjusted slope matrix is singular or ill-conditioned")
        }
        coefficient <- -as.vector(solve(G, shares))
        if (any(!is.finite(coefficient)) || any(coefficient <= 0)) {
            stop("linear slope shape implies nonpositive markups")
        }
        scale <- coefficient[ref] / reference_markup
        slopes <- scale * shape
        markups <- coefficient / scale
        prices <- costs + markups
        intercepts <- shares - as.vector(slopes %*% prices)
        if (any(!is.finite(intercepts)) || any(intercepts < 0)) {
            stop("linear observed inputs imply inadmissible negative demand intercepts")
        }
        parameter <- list(slopes = slopes, intercepts = intercepts)
    } else if (handler == "logit_nests_bertrand") {
        sigma <- parameters$sigma
        if (is.null(nests) || length(nests) != length(shares) || anyNA(nests)) {
            stop("nested Logit observed mode requires an all-product 'nests' vector")
        }
        nests <- as.character(nests)
        if (!is.numeric(sigma) || is.null(names(sigma)) ||
            !setequal(names(sigma), unique(nests)) ||
            any(!is.finite(sigma)) || any(sigma <= 0 | sigma > 1)) {
            stop("nested Logit observed mode requires named 'parameters$sigma' in (0,1] for every nest")
        }
        singleton <- table(nests) == 1L
        if (any(sigma[names(singleton)[singleton]] != 1)) {
            stop("singleton nests require sigma = 1")
        }
        probe_prices <- rep(reference_price, length(shares))
        probe_meanval <- .antitrust_synthetic_nested_meanval(
            shares, probe_prices, -1, nests, sigma, ref)
        probe <- specify("logit_nests", "bertrand",
            prices = probe_prices,
            parameters = list(alpha = -1, sigma = sigma,
                              meanval = probe_meanval),
            ownerPre = market$products$firm_id, nests = nests,
            priceOutside = reference_price)
        D <- numDeriv::jacobian(function(candidate_prices) {
            changed <- probe@model
            changed@pricePre <- candidate_prices
            as.numeric(calcShares(changed, TRUE))
        }, probe_prices)
        G <- t(market$ownership * D)
        rank <- qr(G)$rank
        condition <- kappa(G, exact = TRUE)
        if (rank < length(shares) || !is.finite(condition) ||
            condition > 1 / sqrt(.Machine$double.eps)) {
            stop("nested Logit ownership FOC system is singular or ill-conditioned")
        }
        coefficient <- -as.vector(solve(G, shares))
        native <- as.numeric(calcMargins(probe@model, TRUE, level = TRUE))
        if (any(!is.finite(coefficient)) || any(coefficient <= 0) ||
            max(abs(coefficient - native)) > 1e-6) {
            stop("nested Logit derivative and native markup equations disagree")
        }
        alpha <- -coefficient[ref] / reference_markup
        markups <- coefficient / (-alpha)
        prices <- costs + markups
        parameter <- list(alpha = unname(alpha), sigma = sigma)
    } else if (handler == "ces_nests_bertrand") {
        sigma <- parameters$sigma
        if (is.null(nests) || length(nests) != length(shares) || anyNA(nests)) {
            stop("nested CES observed mode requires an all-product 'nests' vector")
        }
        nests <- as.character(nests)
        if (!is.numeric(sigma) || is.null(names(sigma)) ||
            !setequal(names(sigma), unique(nests)) ||
            any(!is.finite(sigma)) || any(sigma <= 1)) {
            stop("nested CES observed mode requires named 'parameters$sigma' above one for every nest")
        }
        if (any(table(nests) == 1L)) {
            stop("nested CES observed mode requires at least two products per nest; singleton nesting curvature is not identified")
        }
        nest_share <- tapply(shares, nests, sum)
        within <- shares / as.numeric(nest_share[nests])
        nested_margin <- function(gamma, details = FALSE) {
            n <- length(shares)
            E <- matrix(rep((gamma - 1) * shares, each = n), n, n)
            E <- E + outer(nests, nests, `==`) *
                matrix(rep((sigma[nests] - gamma) * within,
                           each = n), n, n)
            diag(E) <- diag(E) - sigma[nests]
            B <- t(E) * market$ownership
            if (qr(B)$rank < n ||
                !is.finite(kappa(B, exact = TRUE)) ||
                kappa(B, exact = TRUE) > 1 / sqrt(.Machine$double.eps)) {
                return(if (details) NULL else NA_real_)
            }
            m <- -as.vector(solve(B, shares)) / shares
            if (details) list(margin = m, matrix = B) else m[ref]
        }
        upper <- min(sigma) - 1e-7
        lower <- 1 + 1e-7
        if (upper <= lower) stop("nested CES requires sigma_g > gamma > 1")
        grid <- seq(lower, upper, length.out = 61L)
        residual <- vapply(grid, function(g) nested_margin(g) - margin,
                           numeric(1))
        exact <- which(is.finite(residual) & abs(residual) < 1e-10)
        crossings <- which(is.finite(residual[-length(grid)]) &
                           is.finite(residual[-1L]) &
                           residual[-length(grid)] * residual[-1L] < 0)
        if (length(exact) == 1L &&
            all(crossings %in% c(exact - 1L, exact))) {
            gamma <- grid[exact]
        } else if (length(exact) == 0L && length(crossings) == 1L) {
            gamma <- stats::uniroot(function(g) nested_margin(g) - margin,
                                    grid[crossings + 0:1], tol = 1e-12)$root
        } else {
            stop("nested CES reference margin does not uniquely identify gamma in (1, min(sigma)); supply admissible shares, nesting curvatures, and margin")
        }
        details <- nested_margin(gamma, details = TRUE)
        all_margin <- details$margin
        condition <- kappa(details$matrix, exact = TRUE)
        rank <- qr(details$matrix)$rank
        if (any(!is.finite(all_margin)) ||
            any(all_margin <= 0 | all_margin >= 1)) {
            stop("nested CES shares and curvature imply inadmissible product margins")
        }
        prices <- costs / (1 - all_margin)
        markups <- prices - costs
        parameter <- list(gamma = unname(gamma), sigma = sigma)
    } else if (handler == "logit_bertrand") {
        system <- .antitrust_synthetic_foc_solution(shares, market$ownership)
        z <- system$z
        alpha <- -z[ref] / reference_markup
        markups <- -z / alpha
        prices <- costs + markups
        parameter <- list(alpha = unname(alpha))
        condition <- system$condition_number
        rank <- system$rank
    } else if (handler == "logit_moncom") {
        alpha <- -1 / reference_markup
        markups <- rep(reference_markup, length(shares))
        prices <- costs + markups
        parameter <- list(alpha = unname(alpha))
    } else if (handler == "logit_cournot") {
        coefficient <- 1 + firm_share / shares[ref]
        alpha <- -coefficient[ref] / reference_markup
        markups <- coefficient / (-alpha)
        prices <- costs + markups
        parameter <- list(alpha = unname(alpha))
    } else if (handler %in% c("logit_auction2nd", "logit_bargaining",
                               "logit_bargaining2nd")) {
        auction_coefficient <- -log1p(-firm_share) / firm_share
        if (handler == "logit_auction2nd") {
            coefficient <- auction_coefficient
        } else if (handler == "logit_bargaining2nd") {
            coefficient <- (1 - bargpowerPre) * auction_coefficient
        } else {
            relative_power <- bargpowerPre / (1 - bargpowerPre)
            M <- -market$ownership * shares
            diag(M) <- diag(market$ownership) + diag(M)
            rank <- qr(M)$rank
            condition <- kappa(M, exact = TRUE)
            if (rank < length(shares) || !is.finite(condition) ||
                condition > 1 / sqrt(.Machine$double.eps)) {
                stop("bargaining Logit ownership system is singular or ill-conditioned")
            }
            direct <- -log1p(-shares) /
                (relative_power * shares / (1 - shares) - log1p(-shares))
            coefficient <- as.vector(solve(t(M), direct))
        }
        if (any(!is.finite(coefficient)) || any(coefficient <= 0)) {
            stop("the Logit conduct equations imply invalid markups")
        }
        alpha <- -coefficient[ref] / reference_markup
        markups <- coefficient / (-alpha)
        prices <- costs + markups
        parameter <- list(alpha = unname(alpha))
    } else {
        if (handler == "ces_bertrand") {
            gamma <- (1 / margin - firm_share[ref]) / (1 - firm_share[ref])
            all_margin <- 1 / (gamma * (1 - firm_share) + firm_share)
        } else if (handler == "ces_moncom") {
            gamma <- 1 / margin
            all_margin <- rep(margin, length(shares))
        } else if (handler == "ces_cournot") {
            if (margin <= firm_share[ref]) {
                stop("CES Cournot requires reference margin above the reference firm's share for gamma > 1")
            }
            gamma <- (1 - firm_share[ref]) / (margin - firm_share[ref])
            all_margin <- firm_share + (1 - firm_share) / gamma
        } else if (handler == "ces_auction2nd" ||
                   handler == "ces_bargaining2nd") {
            adjusted <- if (handler == "ces_bargaining2nd") {
                margin / (1 - bargpowerPre[ref])
            } else margin
            if (adjusted <= 0 || adjusted >= 1) {
                stop("reference margin exceeds the share available after bargaining")
            }
            gamma <- 1 + log1p(-firm_share[ref]) / log1p(-adjusted)
            all_margin <- 1 - (1 - firm_share)^(1 / (gamma - 1))
            if (handler == "ces_bargaining2nd") {
                all_margin <- all_margin * (1 - bargpowerPre)
            }
        } else {
            adjusted <- margin / (1 - bargpowerPre[ref])
            gamma <- (1 / adjusted - firm_share[ref]) /
                (1 - firm_share[ref])
            all_margin <- (1 - bargpowerPre) /
                (gamma * (1 - firm_share) + firm_share)
        }
        if (!is.finite(gamma) || gamma <= 1 ||
            any(!is.finite(all_margin)) || any(all_margin <= 0 | all_margin >= 1)) {
            stop("the reference margin and shares imply an inadmissible CES curvature or product margin")
        }
        prices <- costs / (1 - all_margin)
        markups <- prices - costs
        parameter <- list(gamma = unname(gamma))
    }
    if (any(!is.finite(prices)) || any(prices <= 0) ||
        any(!is.finite(markups)) || any(markups <= 0)) {
        stop("observed inputs imply non-finite or inadmissible equilibrium prices")
    }
    list(prices = unname(prices), markups = unname(markups),
         parameters = parameter, reference_price = reference_price,
         reference_markup = reference_markup, handler = handler,
         foc_rank = rank, foc_condition_number = condition)
}

.antitrust_synthetic_attach_observed <- function(fit, market, solution) {
    n <- length(market$shares)
    prices <- solution$prices
    costs <- market$observed$costs
    target_shares <- if (fit@spec$variant == "alm") {
        market$observed$unconditional_shares
    } else market$shares
    native_shares <- as.numeric(calcShares(fit@model, TRUE,
                                           revenue = fit@spec$demand %in%
                                               c("ces", "ces_nests", "pcaids")))
    native_markups <- as.numeric(calcMargins(fit@model, TRUE, level = TRUE))
    native_costs <- as.numeric(fit@model@mcPre)
    share_residual <- max(abs(native_shares - target_shares))
    foc_residual <- max(abs(native_markups - (prices - costs)))
    cost_residual <- max(abs(native_costs - costs))
    reference_residual <- (prices[n] - costs[n]) / prices[n] -
        market$observed$reference_margin
    if (any(!is.finite(c(share_residual, foc_residual, cost_residual,
                         reference_residual))) ||
        max(share_residual, abs(reference_residual)) > 1e-8 ||
        max(foc_residual, cost_residual) > 1e-6) {
        stop("observed synthetic realization failed share, cost, reference-margin, or native FOC validation: ",
             paste(signif(c(share_residual, cost_residual, foc_residual,
                            reference_residual), 4), collapse = ", "))
    }
    market$prices <- prices
    market$reference_price <- prices[n]
    market$costs <- costs
    market$markups <- solution$markups
    market$products$price <- prices
    market$products$cost <- costs
    market$products$markup <- solution$markups
    market$products$margin <- solution$markups / prices
    market$design$reference_price <- prices[n]
    market$observed$prices <- prices
    market$observed$reference_price <- prices[n]
    market$diagnostics$equilibrium_status <- "verified"
    market$diagnostics$share_residual <- share_residual
    market$diagnostics$foc_residual <- foc_residual
    fit@observed$synthetic_market <- market
    fit@observed$shares <- market$shares
    fit@observed$unconditional_shares <- target_shares
    fit@observed$passive_outside_share <-
        market$observed$passive_outside_share
    fit@observed$prices <- prices
    fit@observed$costs <- costs
    fit@observed$reference_margin <- market$observed$reference_margin
    fit@observed$reference_price <- prices[n]
    fit@diagnostics$synthetic_market <- market
    fit@diagnostics$synthetic_recovered <- fit@parameters
    fit@diagnostics$synthetic <- list(
        status = "completed", mode = "observed", handler = solution$handler,
        target_shares = target_shares,
        conditional_shares = market$shares,
        passive_outside_share = market$observed$passive_outside_share,
        supplied_costs = costs,
        solved_prices = prices, implied_markup = solution$markups,
        implied_margin = solution$markups / prices,
        recovered_parameters = fit@parameters,
        reference_margin_residual = reference_residual,
        share_residual = share_residual, foc_residual = foc_residual,
        cost_residual = cost_residual, foc_rank = solution$foc_rank,
        foc_condition_number = solution$foc_condition_number,
        equilibrium_status = "verified", equilibrium_check = TRUE)
    fit
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
#' Structural demand parameters can be difficult to choose directly. Observed
#' mode lets users reason in terms of shares, ownership, marginal costs, and a
#' reference proportional margin, then recovers the structural parameter and
#' equilibrium prices implied by the selected model. This is an intuition and
#' calibration experiment, not an empirical data-generating process. In primitives
#' mode, the supplied parameters are passed through [specify()] without hidden
#' recalibration.
#'
#' `n_firms` counts inside firms. The active reference product is an additional
#' one-product firm and is included in `ownerPre`. `n_products` may be a scalar
#' shared by all inside firms or a vector with one count per firm, and
#' `dirichlet_alpha` is product-level.
#'
#' In observed mode, `reference_margin = (p_r-c_r)/p_r` pins down the
#' reference price `p_r=c_r/(1-reference_margin)`. Costs remain supplied
#' design inputs. Standard Logit and CES support Bertrand, Cournot,
#' monopolistic competition, second-score auctions, bargaining, and
#' second-score bargaining. Nested Logit and nested CES Bertrand additionally
#' require `nests` and named `parameters$sigma`; nested CES requires at least
#' two products per nest. Linear and log-linear Bertrand require a relative
#' `parameters$slopes` matrix. Bargaining power defaults to 0.5 per product;
#' `bargpowerPre` in `...` overrides it. Other registered combinations raise
#' an explicit unsupported or underidentification error. Primitives mode
#' retains `prices` and `reference_price`.
#' PCAIDS Bertrand uses its legacy proportional-to-share diversion rule and
#' aggregate market elasticity of -1; the reference margin identifies the
#' reference own-price elasticity.
#' Logit and CES ALM Bertrand, Cournot, second-score auction, and bargaining
#' routes use a separate passive outside-option share.
#' Nested Logit ALM Bertrand also requires `nests` and named
#' `parameters$sigma` with singleton-nest sigma fixed at one.
#' Product shares are conditional on choosing an active product, including
#' the reference product. If `passive_outside_share` is omitted, it is drawn
#' uniformly between 0.1 and 0.5; unconditional product shares are recorded.
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
#' @param shares Optional complete all-product share vector summing to one.
#' @param costs Optional complete positive cost vector in observed mode.
#' @param cost_rule `"common"` or `"uniform"` observed cost design.
#' @param cost_level Positive common cost, default 80.
#' @param cost_range Positive endpoints for uniform heterogeneous costs.
#' @param reference_margin Proportional reference-product margin in `(0,1)`;
#'   if omitted, drawn from the design's `margin_range`.
#' @param passive_outside_share Optional ALM passive outside-option share in
#'   `(0,1)`. If omitted for ALM, draw uniformly from
#'   `passive_outside_range`.
#' @param passive_outside_range ALM passive outside-share draw endpoints,
#'   default 0.1 and 0.5.
#' @param reference_price Positive reference price in primitives mode only.
#' @param prices Optional complete positive price vector in primitives mode.
#' @param outside_margin Retired price-first observed argument; use
#'   `reference_margin` with known costs.
#' @param parameters Named model-specific primitives in primitives mode.
#'   Observed nested Logit and nested CES need named `sigma` by nest;
#'   observed Linear and log-linear Bertrand need a relative all-product
#'   `slopes` matrix.
#' @param seed Optional explicit integer seed.
#' @param ... Model-specific arguments forwarded to [calibrate()] or
#'   [specify()], such as `nests`, `capacities`, or solver controls.
#' @return An [AntitrustFit] directly accepted by [simulate()].
#' @export
synthetic_market <- function(
    demand = "logit", supply = "bertrand",
    mode = c("observed", "primitives"),
    n_firms = 3L, n_products = 1L,
    dirichlet_alpha = NULL, outside_beta = c(2, 8), shares = NULL,
    costs = NULL, cost_rule = c("common", "uniform"), cost_level = 80,
    cost_range = c(50, 100), reference_margin = NULL,
    passive_outside_share = NULL,
    passive_outside_range = c(0.1, 0.5),
    reference_price = 100, prices = NULL, outside_margin = NULL,
    parameters = NULL,
    seed = NULL, ...) {
    mode <- match.arg(mode)
    if (mode == "observed") {
        if (!missing(reference_price) || !is.null(prices) ||
            !is.null(outside_margin)) {
            stop("observed mode uses 'costs' and proportional 'reference_margin'; 'prices', 'reference_price', and 'outside_margin' are primitives-mode/retired price-first inputs")
        }
    } else {
        .antitrust_synthetic_validate_scalar(reference_price,
                                             "reference_price", positive = TRUE)
        if (!is.null(passive_outside_share) ||
            !missing(passive_outside_range)) {
            stop("passive outside-share arguments are only valid in observed mode")
        }
    }
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
    if (mode == "primitives" && is.null(parameters)) {
        stop("primitives mode requires a named 'parameters' list")
    }
    if (!is.null(parameters) && !is.list(parameters)) {
        stop("'parameters' must be a list")
    }
    if (mode == "primitives" && !is.null(prices) &&
        (!is.numeric(prices) || length(prices) != n ||
         any(!is.finite(prices)) || any(prices <= 0) ||
         !isTRUE(all.equal(unname(prices[n]), unname(reference_price))))) {
        stop("'prices' must be a finite, strictly positive all-product vector whose reference price equals 'reference_price'")
    }
    dots <- list(...)
    duplicate <- intersect(names(dots), c("prices", "shares", "costs", "margins",
                                           "ownerPre", "parameters", "demand",
                                           "supply", "conduct", "variant",
                                           "priceOutside", "labels"))
    if (length(duplicate)) {
        stop("argument(s) supplied more than once: ",
             paste(duplicate, collapse = ", "))
    }

    spec <- model_spec(demand, supply)
    if (mode == "observed") {
        if (spec$variant != "alm" && !is.null(passive_outside_share) &&
            !(is.numeric(passive_outside_share) &&
              length(passive_outside_share) == 1L &&
              is.finite(passive_outside_share) &&
              passive_outside_share == 0)) {
            stop("'passive_outside_share' is only supported for ALM observed models")
        }
        if (spec$conduct %in% c("bargaining", "bargaining2nd") &&
            is.null(dots$bargpowerPre)) dots$bargpowerPre <- rep(0.5, n)
        design <- fake_market(
            mode = "observed", n_firms = n_firms, n_products = n_products,
            dirichlet_alpha = dirichlet_alpha, outside_beta = outside_beta,
            shares = shares, costs = costs, cost_rule = cost_rule,
            cost_level = cost_level, cost_range = cost_range,
            reference_margin = reference_margin,
            passive_outside_share = if (spec$variant == "alm") {
                passive_outside_share
            } else 0,
            passive_outside_range = passive_outside_range,
            seed = seed)
        solution <- .antitrust_synthetic_observed_solution(
            spec, design, bargpowerPre = dots$bargpowerPre,
            parameters = parameters, nests = dots$nests)
        if (spec$demand %in% c("linear", "loglin", "pcaids")) {
            meanval <- NULL
        } else if (spec$demand == "logit_nests") {
            meanval <- .antitrust_synthetic_nested_meanval(
                design$shares, solution$prices, solution$parameters$alpha,
                dots$nests, solution$parameters$sigma, n)
        } else if (spec$demand == "ces_nests") {
            meanval <- .antitrust_synthetic_nested_ces_meanval(
                design$shares, solution$prices, solution$parameters$gamma,
                as.character(dots$nests), solution$parameters$sigma, n)
        } else if (spec$demand == "logit") {
            meanval <- log(design$shares / design$shares[n])
            if (!spec$conduct %in% c("auction2nd", "bargaining2nd")) {
                meanval <- meanval - solution$parameters$alpha *
                    (solution$prices - solution$prices[n])
            }
            meanval[n] <- 0
        } else {
            meanval <- (design$shares / design$shares[n]) /
                (solution$prices / solution$prices[n])^
                    (1 - solution$parameters$gamma)
            meanval[n] <- 1
        }
        params <- c(solution$parameters,
                    if (is.null(meanval)) list() else list(meanval = meanval))
        if (spec$demand %in% c("ces", "ces_nests")) params$alpha <- 0
        if (spec$variant == "alm" && spec$demand == "logit_nests") {
            s0 <- design$observed$passive_outside_share
            sigma <- solution$parameters$sigma
            nest_order <- unique(as.character(dots$nests))
            non_singleton <- nest_order[
                table(dots$nests)[nest_order] > 1L]
            fit <- withCallingHandlers(calibrate("logit_nests", "bertrand", variant = "alm",
                shares = design$shares, prices = solution$prices,
                margins = solution$markups / solution$prices,
                ownerPre = design$products$firm_id,
                nests = dots$nests, constraint = FALSE,
                parmsStart = c(solution$parameters$alpha, s0,
                               sigma[non_singleton]),
                priceOutside = solution$prices[n]),
                warning = function(w) {
                    if (startsWith(conditionMessage(w),
                                   "Some nests contain only one product")) {
                        invokeRestart("muffleWarning")
                    }
                })
            if (abs((1 - fit@model@shareInside) - s0) > 1e-6 ||
                max(abs(fit@model@slopes$sigma[names(sigma)] - sigma)) >
                    1e-6) {
                stop("nested Logit ALM calibration did not retain passive outside share or nesting curvature")
            }
        } else if (spec$variant == "alm") {
            s0 <- design$observed$passive_outside_share
            mkt_elast <- if (spec$demand == "logit") {
                solution$parameters$alpha *
                    sum(design$shares * solution$prices) * s0
            } else (1 - solution$parameters$gamma) * s0 - 1
            calibration_args <- list(demand = spec$demand,
                conduct = spec$conduct, variant = "alm",
                shares = design$shares, prices = solution$prices,
                margins = if (spec$conduct == "auction2nd" &&
                              spec$demand == "logit") solution$markups
                    else solution$markups / solution$prices,
                mktElast = mkt_elast,
                parmsStart = c(if (spec$demand == "logit") {
                    solution$parameters$alpha
                } else solution$parameters$gamma, s0),
                priceOutside = solution$prices[n],
                ownerPre = design$products$firm_id)
            if (spec$conduct == "bargaining") {
                calibration_args$bargpowerPre <- dots$bargpowerPre
            }
            if (spec$conduct == "auction2nd" && spec$demand == "logit") {
                calibration_args$priceOutside <- NULL
            }
            fit <- do.call(calibrate, calibration_args)
            if (abs((1 - fit@model@shareInside) - s0) > 1e-6) {
                stop("ALM calibration did not retain the supplied passive outside share")
            }
        } else if (spec$demand == "pcaids") {
            fit <- calibrate("pcaids", "bertrand", shares = design$shares,
                prices = solution$prices,
                knownElast = solution$parameters$knownElast,
                knownElastIndex = n, mktElast = -1,
                ownerPre = design$products$firm_id)
        } else {
            fit <- do.call(specify, c(list(
                demand = spec$demand, conduct = spec$conduct,
                prices = solution$prices, parameters = params,
                ownerPre = design$products$firm_id,
                quantities = if (spec$demand %in% c("linear", "loglin"))
                    design$shares else NULL,
                priceOutside = solution$prices[n],
                labels = paste0("Prod", seq_len(n))), dots))
        }
        recovered <- if (spec$demand %in% c("linear", "loglin", "pcaids")) {
            fit@model@slopes
        } else if (spec$demand %in% c("logit", "logit_nests")) {
            fit@model@slopes$alpha
        } else fit@model@slopes$gamma
        target <- if (spec$demand %in% c("linear", "loglin", "pcaids")) {
            solution$parameters$slopes
        } else if (spec$demand %in% c("logit", "logit_nests")) {
            solution$parameters$alpha
        } else solution$parameters$gamma
        if (!is.numeric(recovered) || length(recovered) != length(target) ||
            any(!is.finite(recovered)) ||
            max(abs(recovered - target)) > 1e-6 * max(1, abs(target))) {
            stop("native calibration disagrees with the observed synthetic FOC inversion")
        }
        return(.antitrust_synthetic_attach_observed(fit, design, solution))
    }
    design <- fake_market(
        mode = "primitives",
        n_firms = n_firms, n_products = n_products,
        dirichlet_alpha = dirichlet_alpha, outside_beta = outside_beta,
        shares = shares,
        prices = prices, price_level = reference_price,
        reference_price = reference_price,
        parameters = if (mode == "primitives") parameters else list(),
        seed = seed
    )
    shares <- design$shares
    prices <- design$prices
    owner <- design$products$firm_id
    ref <- design$design$reference_product
    drawn_markup <- design$observed$outside_margin
    markup <- if (mode == "observed") drawn_markup else NA_real_
    foc_diagnostics <- list()
    if (identical(spec$demand, "logit") &&
        identical(spec$conduct, "bertrand") &&
        identical(spec$variant, "standard")) {
        foc_diagnostics <- .antitrust_synthetic_foc_solution(
            shares, design$ownership
        )
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
