#' A model-independent synthetic market design
#'
#' `SyntheticMarket` is a lightweight S3 representation of a generated market
#' design. It stores observables and experimental-design metadata without
#' imposing a demand system or conduct model.
#'
#' @param design A list describing the random design and reproducibility
#'   inputs.
#' @param firms A data frame of firm-level structure.
#' @param products A data frame of product-level structure.
#' @param ownership A product-level ownership matrix.
#' @param shares A vector of all product shares, including the reference
#'   product.
#' @param prices A vector of all product prices.
#' @param observed A list of observed inputs retained for calibration QA.
#' @param truth A list of supplied structural primitives, when known.
#' @param diagnostics A list of realization and numerical diagnostics.
#' @param metadata A list of non-economic object metadata.
#' @return An object with class `SyntheticMarket`.
#' @export
SyntheticMarket <- function(design, firms, products, ownership, shares,
                             prices, observed = list(), truth = list(),
                             diagnostics = list(), metadata = list()) {
    if (!is.list(design) || !is.data.frame(firms) || !is.data.frame(products) ||
        !is.matrix(ownership) || !is.numeric(shares) || !is.numeric(prices)) {
        stop("design, firms, products, ownership, shares, and prices have invalid types")
    }
    n <- nrow(products)
    if (n < 1L || nrow(ownership) != n || ncol(ownership) != n ||
        length(shares) != n || length(prices) != n) {
        stop("products, ownership, shares, and prices must describe the same products")
    }
    if (any(!is.finite(shares)) || any(shares <= 0) ||
        any(!is.finite(prices)) || any(prices <= 0)) {
        stop("shares and prices must be finite and strictly positive")
    }
    if (abs(sum(shares) - 1) > 100 * .Machine$double.eps) {
        stop("product shares, including the reference product, must sum to one")
    }
    if (any(!is.finite(ownership)) || any(!(ownership %in% c(0, 1))) ||
        !isTRUE(all.equal(ownership, t(ownership)))) {
        stop("ownership must be a finite symmetric zero-one matrix")
    }
    if (!all(c("firm_id", "product_id", "product_share", "price") %in%
             names(products))) {
        stop("products must contain firm_id, product_id, product_share, and price")
    }
    structure(list(
        design = design, firms = firms, products = products,
        ownership = ownership, owner = ownership, shares = unname(shares),
        prices = unname(prices), reference_product = design$reference_product,
        reference_share = if (!is.null(design$reference_product)) {
            unname(shares[design$reference_product])
        } else NA_real_,
        reference_price = if (!is.null(design$reference_product)) {
            unname(prices[design$reference_product])
        } else NA_real_,
        observed_markup = if (is.list(observed)) observed$reference_markup else NULL,
        costs = if ("cost" %in% names(products)) products$cost else NULL,
        markups = if ("markup" %in% names(products)) products$markup else NULL,
        observed = observed, truth = truth, diagnostics = diagnostics,
        metadata = metadata
    ), class = c("SyntheticMarket", "list"))
}

#' @rdname SyntheticMarket
#' @param x A `SyntheticMarket` object.
#' @export
is_SyntheticMarket <- function(x) inherits(x, "SyntheticMarket")

#' Print a synthetic market
#' @param x A `SyntheticMarket` object.
#' @param ... Ignored.
#' @export
print.SyntheticMarket <- function(x, ...) {
    cat("SyntheticMarket\n")
    cat("  mode: ", x$design$mode, "\n", sep = "")
    cat("  inside firms: ", x$design$n_firms, "\n", sep = "")
    cat("  products: ", nrow(x$products), "\n", sep = "")
    cat("  reference product: ", x$design$reference_product, "\n", sep = "")
    cat("  seed: ", x$design$seed, "\n", sep = "")
    if (!is.null(x$diagnostics$equilibrium_status)) {
        cat("  realization: ", x$diagnostics$equilibrium_status, "\n", sep = "")
    }
    invisible(x)
}

#' Return synthetic market design metadata
#' @param x A `SyntheticMarket` object.
#' @return The market's design metadata.
#' @export
market_design <- function(x) {
    if (!is_SyntheticMarket(x)) stop("'x' must be a SyntheticMarket")
    x$design
}
