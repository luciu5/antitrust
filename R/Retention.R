#' Product revenue retention
#'
#' Revenue retention is the share of each product's gross revenue received by
#' its seller. A value below one represents an ad valorem charge and a value
#' above one represents a subsidy. Economic ownership is stored separately.
#' `setRetention()` changes metadata only: it does not recover costs,
#' recalibrate, or solve prices. Prefer `revenueRetentionPre` in `specify()` or
#' `calibrate()` and `revenueRetentionPost` in `simulate()` for fitted models.
#' When changing retention on an already fitted legacy model, refresh its
#' stored `mcPre` and `mcPost` from `calcMC()` before calling `calcPrices()`;
#' price solvers use those stored effective costs.
#'
#' @param object An antitrust model or `AntitrustFit`.
#' @param retentionPre,retentionPost Positive finite product-level vectors.
#'   Scalars are recycled. Named vectors are matched to product labels.
#'   Omitted values retain their existing state; on a new model the post state
#'   defaults to the pre state.
#' @return `setRetention()` returns the updated object; `getRetention()`
#'   returns the product-level vector, with an all-one default.
#' @export
setRetention <- function(object, retentionPre = NULL, retentionPost = NULL) {
  if (methods::is(object, "AntitrustFit")) {
    object@model <- setRetention(object@model, retentionPre, retentionPost)
    return(object)
  }
  if (!methods::is(object, "Antitrust")) {
    stop("'object' must be an Antitrust model or AntitrustFit.")
  }
  n <- length(object@prices)
  labels <- object@labels
  normalize <- function(x, name) {
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x <= 0)) {
      stop("'", name, "' must contain positive finite numbers.")
    }
    if (length(x) == 1L) x <- rep(x, n)
    if (length(x) != n) stop("'", name, "' must have one value per product.")
    if (!is.null(names(x))) {
      if (anyDuplicated(names(x)) || !setequal(names(x), labels)) {
        stop("'", name, "' names must match product labels.")
      }
      x <- x[labels]
    }
    as.numeric(x)
  }
  existing <- attr(object, "antitrust_revenue_retention", exact = TRUE)
  pre <- if (is.null(existing)) rep(1, n) else existing$pre
  post <- if (is.null(existing)) pre else existing$post
  if (!is.null(retentionPre)) {
    pre <- normalize(retentionPre, "retentionPre")
    if (is.null(existing) && is.null(retentionPost)) post <- pre
  }
  if (!is.null(retentionPost)) post <- normalize(retentionPost, "retentionPost")
  attr(object, "antitrust_revenue_retention") <- list(pre = pre, post = post)
  object
}

#' @rdname setRetention
#' @param preMerger Return the pre-merger vector if `TRUE`.
#' @export
getRetention <- function(object, preMerger = TRUE) {
  if (methods::is(object, "AntitrustFit")) object <- object@model
  if (!methods::is(object, "Antitrust")) {
    stop("'object' must be an Antitrust model or AntitrustFit.")
  }
  value <- attr(object, "antitrust_revenue_retention", exact = TRUE)
  if (is.null(value)) return(rep(1, length(object@prices)))
  if (preMerger) value$pre else value$post
}

.retention_owner_bertrand <- function(owner, retention) {
  if (is.null(dim(owner))) {
    owner <- matrix(owner, nrow = length(retention),
                    ncol = length(retention))
  }
  sweep(sweep(owner, 2L, retention, "*"), 1L, retention, "/")
}

.retention_owner_bargaining_logit <- function(owner, retention) {
  if (is.null(dim(owner))) {
    owner <- matrix(owner, nrow = length(retention),
                    ncol = length(retention))
  }
  sweep(sweep(owner, 1L, retention, "*"), 2L, retention, "/")
}

.require_uniform_auction_retention <- function(object, preMerger, subset) {
  if (.mixed_firm_retention(object, preMerger, subset)) {
    stop("mixed revenue retention within an auction firm's active portfolio requires a product-level bidding game")
  }
  invisible(TRUE)
}

.mixed_firm_retention <- function(object, preMerger, subset) {
  retention <- getRetention(object, preMerger)[subset]
  owner <- if (preMerger) object@ownerPre else object@ownerPost
  owner <- owner[subset, subset, drop = FALSE]
  any((owner > 0) &
      (abs(outer(log(retention), log(retention), "-")) > 1e-10))
}

.auction_effective_cost_delta <- function(object) {
  pre_r <- getRetention(object, TRUE)
  post_r <- getRetention(object, FALSE)
  if (all(pre_r == 1) && all(post_r == 1)) return(object@mcDelta)
  state <- attr(object, "antitrust_cost_state", exact = TRUE)
  if (is.list(state) && !is.null(state$base)) {
    return(.persistent_mc(object, FALSE) - .persistent_mc(object, TRUE))
  }
  pre <- object@mcPre
  if (length(pre) != length(pre_r) || any(!is.finite(pre))) {
    pre <- calcMC(object, TRUE)
  }
  if (methods::is(object, "Auction2ndLogit") &&
      !methods::is(object, "Auction2ndCES")) {
    (pre * pre_r + object@mcDelta) / post_r - pre
  } else {
    pre * (1 + object@mcDelta) * pre_r / post_r - pre
  }
}
