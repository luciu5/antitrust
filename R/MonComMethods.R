# Monopolistic-competition methods for the existing flat Logit and CES
# demand classes.  The atomistic/competitive-fringe convention holds the
# relevant aggregate demand index fixed and uses only the own-product
# derivative.  For flat Logit this is D^MC_j = alpha * s_j, so an output
# market has markup -1 / alpha, rather than Bertrand's
# -1 / (alpha * (1 - s_j)).  CES uses the own derivative already encoded by
# the mature CES elasticity method.  Demand shares and counterfactual
# promotion remain inherited from those implementations.

.moncom_output_sign <- function(object) {
  if (isTRUE(object@output)) -1 else 1
}

.moncom_foc_residual <- function(object, preMerger = TRUE) {
  prices <- if (preMerger) object@pricePre else object@pricePost
  mc <- if (preMerger) object@mcPre else object@mcPost
  active <- if (preMerger) rep(TRUE, length(prices)) else object@subset
  if (!length(prices) || !any(active)) return(NA_real_)
  actual <- if (isTRUE(object@output)) prices - mc else mc - prices
  implied <- calcMargins(object, preMerger = preMerger, level = TRUE)
  max(abs(actual[active] - implied[active]), na.rm = TRUE)
}

.moncom_diagnostics <- function(object) {
  list(
    conduct = "moncom",
    foc_residual_pre = .moncom_foc_residual(object, preMerger = TRUE),
    foc_residual_post = .moncom_foc_residual(object, preMerger = FALSE),
    ownership_irrelevant = TRUE
  )
}

#' @rdname Params-Methods
#' @export
setMethod(
  f = "calcSlopes",
  signature = "MonComLogit",
  definition = function(object) {
    prices <- as.numeric(object@prices)
    margins <- as.numeric(object@margins)
    shares <- as.numeric(object@shares)
    observed <- which(is.finite(prices) & prices > 0 &
      is.finite(margins) & margins > 0)
    if (!length(observed)) {
      stop("MonCom Logit calibration requires at least one positive observed margin")
    }

    ## The own-product FOC is markup_j = -1 / alpha for output markets and
    ## markup_j = 1 / alpha for input markets.  A weighted least-squares
    ## average accommodates the rounding of reported margins without
    ## replacing the conduct equation with the Bertrand ownership system.
    markup <- margins[observed] * prices[observed]
    weights <- .product_weights(object)[observed]
    common_markup <- stats::weighted.mean(markup, w = weights)
    if (!is.finite(common_markup) || common_markup <= 0) {
      stop("MonCom Logit calibration requires positive finite observed markups")
    }
    alpha <- .moncom_output_sign(object) / common_markup

    idx <- object@normIndex
    share_inside <- object@shareInside
    if (is.na(idx)) {
      idx_share <- 1 - share_inside
      idx_price <- object@priceOutside
      if (!is.finite(idx_share) || idx_share <= 0) {
        stop("MonCom Logit calibration requires a positive outside share when normIndex is NA")
      }
      meanval <- log(shares / idx_share) - alpha * (prices - idx_price)
    } else {
      idx <- as.integer(idx)[1]
      if (idx < 1L || idx > length(shares) || shares[idx] <= 0) {
        stop("MonCom Logit calibration requires a valid positive-share normalization product")
      }
      idx_share <- shares[idx]
      idx_price <- prices[idx]
      meanval <- log(shares / idx_share) - alpha * (prices - idx_price)
      meanval[idx] <- 0
    }

    names(meanval) <- object@labels
    names(alpha) <- "alpha"
    object@slopes <- list(alpha = alpha, meanval = meanval)
    object@priceOutside <- idx_price
    if (is.finite(object@insideSize) && share_inside > 0) {
      object@mktSize <- object@insideSize / share_inside
    }
    object
  }
)

#' @rdname Params-Methods
#' @export
setMethod(
  f = "calcSlopes",
  signature = "MonComCES",
  definition = function(object) {
    prices <- as.numeric(object@prices)
    margins <- as.numeric(object@margins)
    shares <- as.numeric(object@shares)
    observed <- which(is.finite(prices) & prices > 0 &
      is.finite(margins) & margins > 0 & shares < 1)
    if (!length(observed)) {
      stop("MonCom CES calibration requires at least one positive observed margin with share below one")
    }

    ## CES own-price elasticity is e_jj = -gamma + (gamma - 1) s_j.
    ## MonCom uses only this own derivative, so gamma is identified from
    ## own margins and shares rather than the ownership-adjusted Bertrand
    ## elasticity matrix.
    own_elast <- .moncom_output_sign(object) / margins[observed]
    x <- 1 - shares[observed]
    weights <- .product_weights(object)[observed]
    gamma <- -sum(weights * x * (own_elast + shares[observed])) /
      sum(weights * x^2)
    if (!is.finite(gamma)) stop("MonCom CES calibration could not identify gamma")
    if (isTRUE(object@output) && gamma <= 1) {
      stop("MonCom output-market CES requires gamma > 1")
    }
    if (!isTRUE(object@output) &&
      any(gamma >= -shares / pmax(1 - shares, .Machine$double.eps))) {
      stop("MonCom input-market CES requires gamma below every own-elasticity singularity")
    }

    share_inside <- object@shareInside
    alpha <- if (is.finite(share_inside) && share_inside > 0 && share_inside < 1) {
      1 / share_inside - 1
    } else {
      NULL
    }
    idx <- object@normIndex
    if (is.na(idx)) {
      idx_share <- 1 - share_inside
      idx_price <- object@priceOutside
      if (!is.finite(idx_share) || idx_share <= 0 || idx_price <= 0) {
        stop("MonCom CES calibration requires a positive outside share and price when normIndex is NA")
      }
      meanval <- shares / (prices / idx_price)^(1 - gamma) / idx_share
    } else {
      idx <- as.integer(idx)[1]
      if (idx < 1L || idx > length(shares) || shares[idx] <= 0 || prices[idx] <= 0) {
        stop("MonCom CES calibration requires a valid positive-share normalization product")
      }
      idx_share <- shares[idx]
      idx_price <- prices[idx]
      meanval <- shares / (prices / idx_price)^(1 - gamma)
      meanval <- meanval / meanval[idx]
      meanval[idx] <- 1
    }

    names(meanval) <- object@labels
    names(gamma) <- "gamma"
    object@slopes <- list(alpha = alpha, gamma = gamma, meanval = meanval)
    object@priceOutside <- idx_price
    if (is.finite(object@insideSize)) {
      object@mktSize <- object@insideSize * (1 + if (is.null(alpha)) 0 else alpha)
    }
    object
  }
)

## calcMargins() is the only conduct-specific demand derivative used by the
## inherited Logit price solver.  Because it ignores ownerPre/ownerPost, the
## same stable solver now solves independent own-product problems.
#' @rdname Margins-Methods
#' @export
setMethod(
  f = "calcMargins",
  signature = "MonComLogit",
  definition = function(object, preMerger = TRUE, level = FALSE) {
    prices <- if (preMerger) object@pricePre else object@pricePost
    active <- if (preMerger) rep(TRUE, length(prices)) else object@subset
    alpha <- object@slopes$alpha
    if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha == 0) {
      stop("MonCom Logit requires a finite non-zero alpha")
    }
    result <- rep(NA_real_, length(prices))
    result[active] <- .moncom_output_sign(object) /
      (alpha * prices[active])
    if (level) result[active] <- result[active] * prices[active]
    names(result) <- object@labels
    result
  }
)

#' @rdname Margins-Methods
#' @export
setMethod(
  f = "calcMargins",
  signature = "MonComCES",
  definition = function(object, preMerger = TRUE, level = FALSE) {
    prices <- if (preMerger) object@pricePre else object@pricePost
    active <- if (preMerger) rep(TRUE, length(prices)) else object@subset
    own_elast <- diag(elast(object, preMerger = preMerger))
    result <- rep(NA_real_, length(prices))
    if (any(!is.finite(own_elast[active]) | abs(own_elast[active]) < 1e-12)) {
      stop("MonCom CES has a singular own-product elasticity")
    }
    result[active] <- .moncom_output_sign(object) / own_elast[active]
    if (level) result[active] <- result[active] * prices[active]
    names(result) <- object@labels
    result
  }
)

# The mature Logit solver evaluates calcMargins() inside its FOCs.  Dispatch
# therefore gives MonCom the same robust nleqslv/BBsolve behavior while the
# conduct-specific derivative above removes strategic cross-product effects.
#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "MonComLogit",
  definition = function(object, preMerger = TRUE, ...) {
    callNextMethod(object, preMerger = preMerger, ...)
  }
)

#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "MonComCES",
  definition = function(object, preMerger = TRUE, ...) {
    callNextMethod(object, preMerger = preMerger, ...)
  }
)

#' @rdname Diagnostics-Methods
#' @export
setMethod(
  f = "calcDiagnostics",
  signature = "MonComLogit",
  definition = function(object, labels = object@labels) {
    result <- callNextMethod(object, labels = labels)
    attr(result, "moncom") <- .moncom_diagnostics(object)
    result
  }
)

#' @rdname Diagnostics-Methods
#' @export
setMethod(
  f = "calcDiagnostics",
  signature = "MonComCES",
  definition = function(object, labels = object@labels) {
    result <- callNextMethod(object, labels = labels)
    attr(result, "moncom") <- .moncom_diagnostics(object)
    result
  }
)
