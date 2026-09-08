# Neutral synthetic-market design infrastructure owned by antitrust.
#
# This layer generates shares, prices, ownership, and observed reference
# markups. Demand equations, conduct FOCs, parameter inversion, and cost
# recovery remain in model-specific antitrust realization code.

.antitrust_market_max_seed <- 2147483646L

.antitrust_market_assert_scalar <- function(x, name, integer = FALSE) {
    if (length(x) != 1L || is.na(x) || !is.finite(x) ||
        (integer && x != as.integer(x))) {
        stop("'", name, "' must be a single finite ",
             if (integer) "integer" else "number")
    }
    invisible(x)
}

.antitrust_market_validate_positive_vector <- function(x, name, length = NULL) {
    if (!is.numeric(x) || (!is.null(length) && base::length(x) != length) ||
        any(!is.finite(x)) || any(x <= 0)) {
        requirement <- if (is.null(length)) "finite and strictly positive" else
            paste0("a finite, strictly positive vector of length ", length)
        stop("'", name, "' must be ", requirement)
    }
    invisible(x)
}

.antitrust_market_restore_rng <- function(had_seed, old_seed) {
    if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
    }
}

.antitrust_market_begin_rng <- function(seed) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    if (is.null(seed)) {
        actual_seed <- sample.int(.antitrust_market_max_seed, 1L)
    } else {
        .antitrust_market_assert_scalar(seed, "seed", integer = TRUE)
        actual_seed <- as.integer(seed)
        if (actual_seed < 1L || actual_seed > .antitrust_market_max_seed) {
            stop("'seed' must be between 1 and ", .antitrust_market_max_seed)
        }
    }
    set.seed(actual_seed)
    list(had_seed = had_seed, old_seed = old_seed, seed = actual_seed)
}

.antitrust_market_price_vector <- function(n_total_products, n_inside, prices,
                                           price_rule, price_level, price_range,
                                           reference_price) {
    if (!is.null(prices)) {
        if (!is.numeric(prices) || !all(is.finite(prices)) || any(prices <= 0) ||
            !length(prices) %in% c(n_inside, n_total_products)) {
            stop("'prices' must be a finite, strictly positive vector of inside or all-product length")
        }
        if (length(prices) == n_total_products) {
            if (!is.null(reference_price) &&
                !isTRUE(all.equal(prices[n_total_products], reference_price))) {
                stop("'reference_price' conflicts with the supplied all-product 'prices'")
            }
            return(list(values = unname(prices), rule = "user-supplied"))
        }
        ref <- if (is.null(reference_price)) price_level else reference_price
        if (length(ref) != 1L || !is.finite(ref) || ref <= 0) {
            stop("'reference_price' must be a finite, strictly positive number")
        }
        return(list(values = c(unname(prices), ref), rule = "user-supplied-inside"))
    }

    price_rule <- match.arg(price_rule, c("common", "uniform"))
    if (price_rule == "common") {
        .antitrust_market_assert_scalar(price_level, "price_level")
        if (price_level <= 0) stop("'price_level' must be strictly positive")
        return(list(values = rep(price_level, n_total_products), rule = "common"))
    }
    if (!is.numeric(price_range) || length(price_range) != 2L ||
        any(!is.finite(price_range)) || any(price_range <= 0) ||
        price_range[1] >= price_range[2]) {
        stop("'price_range' must contain two finite positive values in increasing order")
    }
    list(values = stats::runif(n_total_products, min = price_range[1],
                               max = price_range[2]), rule = "uniform")
}

.antitrust_market_open_uniform <- function(n, support) {
    if (!is.numeric(support) || length(support) != 2L ||
        any(!is.finite(support)) || support[1] < 0 || support[1] >= support[2]) {
        stop("'markup_range' must contain two finite values with 0 <= lower < upper")
    }
    width <- diff(support)
    eps <- max(.Machine$double.eps * max(1, width), 1e-12)
    if (width <= 2 * eps) stop("'markup_range' is too narrow for an open-boundary draw")
    stats::runif(n, min = support[1] + eps, max = support[2] - eps)
}

#' Generate a reproducible synthetic antitrust market design
#'
#' The design layer draws shares, prices, ownership, and optionally one
#' observed reference-product markup. It does not infer demand or conduct
#' primitives. The active reference product is an additional one-product firm
#' and participates in the ownership matrix; its mean utility is normalized by
#' the economic realization adapter.
#'
#' `n_products` may be a scalar, in which case all inside firms have the same
#' number of products, or a vector of length `n_firms`, in which case it gives
#' each firm's product count. Shares are drawn at the product level, so
#' `dirichlet_alpha` has one element for each inside product.
#'
#' @param mode Either `"observed"` or `"primitives"`.
#' @param n_firms Number of inside firms; the reference firm is additional.
#' @param n_products Number of products owned by each inside firm. A scalar is
#'   recycled across firms; a vector must have length `n_firms`.
#' @param dirichlet_alpha A positive finite vector with one element for each
#'   inside product. If omitted, all product-level Dirichlet shapes equal one.
#' @param outside_beta A positive finite length-two vector of Beta shapes.
#' @param prices Optional positive prices for inside products or all products.
#' @param reference_price Price appended when only inside prices are supplied.
#' @param price_rule Either `"common"` or `"uniform"` when prices are not
#'   supplied.
#' @param price_level Common positive price and default reference price.
#' @param price_range Two positive endpoints for the uniform price rule.
#' @param outside_margin Optional observed reference-product level margin.
#' @param observed_markup Backward-compatible alias for `outside_margin`.
#' @param reference_markup Backward-compatible alias for `outside_margin`.
#' @param markup_range Two endpoints for the observed markup draw.
#' @param parameters A list of model-specific known primitives in primitives
#'   mode, such as `list(alpha = -1)`.
#' @param alpha Optional shorthand for `parameters$alpha`.
#' @param seed An optional explicit integer seed.
#' @return A `SyntheticMarket` object.
#' @export
fake_market <- function(
    mode = "observed", n_firms = 3L, n_products = 1L,
    dirichlet_alpha = NULL, outside_beta = c(2, 8), prices = NULL,
    reference_price = NULL, price_rule = c("common", "uniform"),
    price_level = 100, price_range = c(50, 150), outside_margin = NULL,
    observed_markup = NULL, reference_markup = NULL, markup_range = c(0, 100),
    parameters = list(), alpha = NULL, seed = NULL) {
    mode <- match.arg(mode, c("observed", "primitives",
                              "observed_information", "known_primitives"))
    mode <- switch(mode, observed_information = "observed",
                   known_primitives = "primitives", mode)
    n_firms <- as.numeric(n_firms)
    .antitrust_market_assert_scalar(n_firms, "n_firms", integer = TRUE)
    n_firms <- as.integer(n_firms)
    if (n_firms < 1L) stop("'n_firms' must be at least one")
    n_products_input <- as.numeric(n_products)
    if (length(n_products_input) != 1L && length(n_products_input) != n_firms) {
        stop("'n_products' must be a positive integer scalar or a vector of length n_firms")
    }
    if (any(!is.finite(n_products_input)) ||
        any(n_products_input != as.integer(n_products_input)) ||
        any(n_products_input < 1)) {
        stop("'n_products' must contain positive integers")
    }
    products_per_firm <- if (length(n_products_input) == 1L) {
        rep(as.integer(n_products_input), n_firms)
    } else as.integer(n_products_input)
    if (sum(products_per_firm) > .antitrust_market_max_seed) {
        stop("the number of inside products is too large")
    }
    n_inside <- as.integer(sum(products_per_firm))
    n_total_products <- n_inside + 1L
    n_products_design <- if (length(n_products_input) == 1L) {
        as.integer(n_products_input)
    } else products_per_firm
    if (is.null(dirichlet_alpha)) dirichlet_alpha <- rep(1, n_inside)
    .antitrust_market_validate_positive_vector(dirichlet_alpha,
                                               "dirichlet_alpha", n_inside)
    .antitrust_market_validate_positive_vector(outside_beta, "outside_beta", 2L)
    if (!is.numeric(markup_range) || length(markup_range) != 2L ||
        any(!is.finite(markup_range)) || markup_range[1] < 0 ||
        markup_range[1] >= markup_range[2]) {
        stop("'markup_range' must contain two finite values with 0 <= lower < upper")
    }
    if (!is.list(parameters)) stop("'parameters' must be a list")
    if (!is.null(alpha)) {
        if (length(alpha) != 1L || !is.finite(alpha)) stop("'alpha' must be a finite scalar")
        if (!is.null(parameters$alpha) && !isTRUE(all.equal(parameters$alpha, alpha))) {
            stop("'alpha' conflicts with 'parameters$alpha'")
        }
        parameters$alpha <- alpha
    }
    if (mode == "primitives" && !length(parameters)) {
        stop("'parameters' must contain the known structural primitives in primitives mode")
    }
    if (mode == "primitives" && (!is.null(outside_margin) ||
        !is.null(observed_markup) || !is.null(reference_markup))) {
        stop("outside margin arguments are only valid in observed mode")
    }
    if (!is.null(outside_margin) && !is.null(observed_markup) &&
        !isTRUE(all.equal(outside_margin, observed_markup))) {
        stop("'outside_margin' conflicts with 'observed_markup'")
    }
    if (!is.null(outside_margin) && !is.null(reference_markup) &&
        !isTRUE(all.equal(outside_margin, reference_markup))) {
        stop("'outside_margin' conflicts with 'reference_markup'")
    }
    if (!is.null(reference_markup)) {
        if (!is.null(observed_markup) &&
            !isTRUE(all.equal(observed_markup, reference_markup))) {
            stop("'reference_markup' conflicts with 'observed_markup'")
        }
        observed_markup <- reference_markup
    }
    if (!is.null(outside_margin)) observed_markup <- outside_margin
    if (!is.null(observed_markup) &&
        (length(observed_markup) != 1L || !is.finite(observed_markup) ||
         observed_markup <= markup_range[1] || observed_markup >= markup_range[2])) {
        stop("'outside_margin' must lie strictly inside 'markup_range'")
    }

    rng <- .antitrust_market_begin_rng(seed)
    on.exit(.antitrust_market_restore_rng(rng$had_seed, rng$old_seed), add = TRUE)
    relative_product_shares <- stats::rgamma(n_inside, shape = dirichlet_alpha, rate = 1)
    relative_product_shares <- relative_product_shares / sum(relative_product_shares)
    outside_share <- stats::rbeta(1L, shape1 = outside_beta[1], shape2 = outside_beta[2])
    product_shares_inside <- (1 - outside_share) * relative_product_shares
    inside_firm_id <- rep(seq_len(n_firms), times = products_per_firm)
    firm_shares <- vapply(seq_len(n_firms), function(firm) {
        sum(product_shares_inside[inside_firm_id == firm])
    }, numeric(1))
    product_shares <- c(product_shares_inside, outside_share)
    price_info <- .antitrust_market_price_vector(
        n_total_products, n_inside, prices, price_rule, price_level, price_range,
        reference_price
    )
    price_values <- price_info$values
    if (mode == "observed") {
        if (is.null(observed_markup)) {
            observed_markup <- .antitrust_market_open_uniform(1L, markup_range)
            markup_rule <- "uniform-open-U(0,100)"
        } else markup_rule <- "user-supplied"
    } else {
        observed_markup <- NULL
        markup_rule <- "not-drawn"
    }
    observed_markup_product <- if (is.null(observed_markup)) NA_real_ else unname(observed_markup)
    firm_id <- c(inside_firm_id, n_firms + 1L)
    product_id <- seq_len(n_total_products)
    reference_product <- n_total_products
    firm_share_by_product <- c(rep(firm_shares, times = products_per_firm), outside_share)
    reference_firm <- n_firms + 1L
    ownership <- outer(firm_id, firm_id, FUN = "==") * 1
    dimnames(ownership) <- list(product_id, product_id)
    products <- data.frame(
        product_id = product_id, firm_id = firm_id,
        firm_share = unname(firm_share_by_product),
        product_share = unname(product_shares), price = unname(price_values),
        cost = rep(NA_real_, n_total_products), markup = rep(NA_real_, n_total_products),
        observed_markup = c(rep(NA_real_, n_inside), observed_markup_product),
        outside_margin = c(rep(NA_real_, n_inside), observed_markup_product),
        reference_product = product_id == reference_product,
        stringsAsFactors = FALSE
    )
    firms <- data.frame(
        firm_id = seq_len(reference_firm),
        reference_firm = seq_len(reference_firm) == reference_firm,
        firm_share = c(firm_shares, outside_share),
        n_products = c(products_per_firm, 1L), stringsAsFactors = FALSE
    )
    design <- list(
        mode = mode, seed = rng$seed, n_firms = n_firms,
        n_inside_products = n_inside, n_products = n_products_design,
        products_per_firm = unname(products_per_firm), n_total_products = n_total_products,
        reference_product = reference_product, reference_firm = reference_firm,
        dirichlet_alpha = unname(dirichlet_alpha), outside_beta = unname(outside_beta),
        outside_share = unname(outside_share), relative_product_shares = unname(relative_product_shares),
        firm_shares = unname(firm_shares), ownership_map = data.frame(product_id = product_id, firm_id = firm_id),
        price_rule = price_info$rule, price_level = price_level,
        price_range = unname(price_range), reference_price = unname(price_values[reference_product]),
        markup_rule = markup_rule, markup_range = unname(markup_range),
        outside_margin = observed_markup, observed_reference_markup = observed_markup,
        parameters = parameters
    )
    observed <- list(
        shares = unname(product_shares), prices = unname(price_values), ownership = ownership,
        reference_product = reference_product, reference_share = unname(outside_share),
        reference_price = unname(price_values[reference_product]), outside_margin = observed_markup,
        reference_markup = observed_markup
    )
    truth <- if (mode == "primitives") parameters else list()
    diagnostics <- list(equilibrium_status = "design-only", foc_residual = NA_real_,
                        rejection_reason = NULL, costs_status = "not-realized")
    metadata <- list(class = "SyntheticMarket", shares_sum = sum(product_shares),
                     units = list(price = "price level", markup = "price level"),
                     reference_normalization = "mean utility only; reference price is real")
    SyntheticMarket(design = design, firms = firms, products = products,
                   ownership = ownership, shares = product_shares, prices = price_values,
                   observed = observed, truth = truth, diagnostics = diagnostics,
                   metadata = metadata)
}

.antitrust_market_derived_seed <- function(seed, index) {
    as.integer(((as.double(seed) + as.double(index) - 1) %%
                .antitrust_market_max_seed) + 1)
}

#' Generate deterministic synthetic-market replications
#'
#' @param n Number of requested replications.
#' @param generator A function such as [fake_market].
#' @param seed Optional base integer seed.
#' @param max_attempts Maximum attempts per requested replication.
#' @param realizer Optional model-specific function applied to each market.
#' @param ... Arguments passed to `generator`.
#' @return A `SyntheticMarketBatch` list with markets, seeds, and diagnostics.
#' @export
simulate_markets <- function(n, generator = fake_market, seed = NULL,
                             max_attempts = 1L, realizer = NULL, ...) {
    .antitrust_market_assert_scalar(n, "n", integer = TRUE)
    .antitrust_market_assert_scalar(max_attempts, "max_attempts", integer = TRUE)
    n <- as.integer(n); max_attempts <- as.integer(max_attempts)
    if (n < 1L) stop("'n' must be at least one")
    if (max_attempts < 1L) stop("'max_attempts' must be at least one")
    if (!is.function(generator)) stop("'generator' must be a function")
    if (!is.null(realizer) && !is.function(realizer)) {
        stop("'realizer' must be NULL or a function")
    }
    rng <- .antitrust_market_begin_rng(seed)
    on.exit(.antitrust_market_restore_rng(rng$had_seed, rng$old_seed), add = TRUE)
    base_seed <- rng$seed
    markets <- vector("list", n); used_seeds <- vector("list", n)
    rejection_reasons <- vector("list", n)
    for (i in seq_len(n)) {
        reasons <- character()
        for (attempt in seq_len(max_attempts)) {
            attempt_index <- (i - 1L) * max_attempts + attempt
            draw_seed <- .antitrust_market_derived_seed(base_seed, attempt_index)
            used_seeds[[i]] <- draw_seed
            result <- tryCatch({
                generated <- do.call(generator, c(list(seed = draw_seed), list(...)))
                if (is.null(realizer)) generated else realizer(generated)
            }, error = function(e) e)
            if (!inherits(result, "error")) {
                markets[[i]] <- result
                break
            }
            reasons <- c(reasons, conditionMessage(result))
        }
        rejection_reasons[[i]] <- reasons
    }
    rejected <- vapply(markets, is.null, logical(1))
    structure(list(
        markets = markets, seeds = unlist(used_seeds, use.names = FALSE),
        diagnostics = list(n_requested = n, n_success = sum(!rejected),
                           n_rejected = sum(rejected), rejection_rate = mean(rejected),
                           rejection_reasons = rejection_reasons, base_seed = base_seed,
                           max_attempts = max_attempts, realized = !is.null(realizer))
    ), class = c("SyntheticMarketBatch", "list"))
}

#' @export
print.SyntheticMarketBatch <- function(x, ...) {
    d <- x$diagnostics
    cat("SyntheticMarketBatch\n")
    cat("  requested: ", d$n_requested, "\n", sep = "")
    cat("  successful: ", d$n_success, "\n", sep = "")
    cat("  rejected: ", d$n_rejected, " (", format(d$rejection_rate), ")\n", sep = "")
    invisible(x)
}
