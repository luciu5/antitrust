## Model-aware synthetic markets are deliberately implemented in antitrust.
## iopolicy supplies only the neutral share and ownership design; it does not
## dispatch demand or supply models.

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
    level_conduct <- conduct %in% c("auction2nd", "bargaining", "bargaining2nd",
                                    "vertical_bargaining")
    if (level_conduct) markup else markup / reference_price
}

.antitrust_synthetic_complete_parameters <- function(spec, parameters, shares,
                                                     prices, reference_product) {
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
        if (is.null(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
            alpha >= 0) {
            stop("primitives mode requires a finite negative 'alpha' for ",
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

.antitrust_synthetic_diagnostics <- function(fit, shares, prices,
                                             reference_product,
                                             reference_markup,
                                             mode) {
    model <- fit@model
    n <- length(shares)
    price_pre <- if ("pricePre" %in% methods::slotNames(model)) {
        as.numeric(model@pricePre)
    } else {
        rep(NA_real_, n)
    }
    price_residual <- if (all(is.finite(price_pre))) {
        max(abs(price_pre - prices))
    } else {
        NA_real_
    }
    implied_markup <- tryCatch(
        as.numeric(calcMargins(model, preMerger = TRUE, level = TRUE)),
        error = function(e) rep(NA_real_, n)
    )

    foc <- rep(NA_real_, n)
    foc_residual <- NA_real_
    ## This is the complete product-level Bertrand Logit FOC, with the
    ## reference product included in the same strategic system. Other model
    ## families retain their package-specific solver diagnostics below.
    if (methods::is(model, "Logit") && !methods::is(model, "LogitCournot") &&
        !methods::is(model, "LogitCap") &&
        !methods::is(model, "LogitNests")) {
        alpha <- model@slopes$alpha
        owner <- model@ownerPre
        mu <- implied_markup
        if (length(alpha) == 1L && is.finite(alpha) &&
            length(mu) == n && nrow(owner) == n) {
            derivative <- alpha * (diag(shares) - tcrossprod(shares))
            foc <- shares + t(owner * derivative) %*% mu
            foc_residual <- max(abs(foc))
        }
    }
    list(
        status = "completed",
        mode = mode,
        price_residual = price_residual,
        implied_markup = unname(implied_markup),
        reference_markup_error = unname(implied_markup[reference_product] -
                                            reference_markup),
        foc = unname(foc),
        foc_residual = foc_residual,
        foc_tolerance = 1e-8,
        equilibrium_check = if (is.finite(price_residual)) {
            price_residual < 1e-7
        } else {
            NA
        }
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
                                        parameter_truth = list()) {
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
    fit@diagnostics$synthetic <- .antitrust_synthetic_diagnostics(
        fit, market$shares, market$prices,
        market$design$reference_product, reference_markup, mode
    )
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
#' one-product firm and is included in `ownerPre`. All inside firms have
#' `n_products` products, and `dirichlet_alpha` is product-level with length
#' `n_firms * n_products`.
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
#' @param n_products Number of products per inside firm.
#' @param dirichlet_alpha Positive product-level Dirichlet parameters.
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
    dirichlet_alpha = rep(1, n_firms * n_products),
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
    if (!is.numeric(n_products) || length(n_products) != 1L ||
        n_products != as.integer(n_products) || n_products < 1L) {
        stop("'n_products' must be a positive integer")
    }
    n_firms <- as.integer(n_firms)
    n_products <- as.integer(n_products)
    n <- n_firms * n_products + 1L
    if (!is.numeric(dirichlet_alpha) || length(dirichlet_alpha) != n - 1L ||
        any(!is.finite(dirichlet_alpha)) || any(dirichlet_alpha <= 0)) {
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
    design <- iopolicy::fake_market(
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
            fit, design, mode, markup, parameter_truth = list()
        ))
    }

    parameters <- .antitrust_synthetic_complete_parameters(
        spec, parameters, shares, prices, ref
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
        fit, design, mode, NA_real_, parameter_truth = parameters
    )
}
