#' @title Monopolistic-competition calibration
#' @name MonCom-Functions
#' @aliases moncom.logit moncom.ces
#' @description Calibrate flat Logit or CES demand under differentiated-product
#' monopolistic competition.  Each product uses its own perceived demand
#' derivative and does not internalize cross-product effects.  For CES, the
#' aggregate CES index is held fixed, so the perceived own elasticity is
#' `-gamma`, rather than the full share-adjusted CES elasticity.  The BLP
#' MonCom model is exposed through the general `calibrate()` and `specify()`
#' lifecycle using the validated integration engine.
#' @param prices A length-k vector of observed product prices.
#' @param shares A length-k vector of quantity shares for Logit or revenue
#'   shares for CES.
#' @param margins A length-k vector of observed positive margins. At least one
#'   margin is required; additional margins are used as precision weights.
#' @param ownerPre Pre-merger ownership vector or matrix. It is retained as
#'   model metadata but does not enter the MonCom FOC.
#' @param ownerPost Post-merger ownership vector or matrix.
#' @param output `TRUE` for an output market and `FALSE` for an input market.
#' @param weights A length-k vector of non-negative product weights used in
#'   calibration.
#' @param mcDelta A vector of proportional marginal-cost changes.
#' @param subset A logical vector identifying products active in the post-
#'   counterfactual market.
#' @param insideSize Total inside quantity (Logit) or revenue (CES).
#' @param priceOutside Outside-good price when the market has an outside good.
#' @param normIndex Normalization product when inside shares exhaust the market.
#' @param priceStart Starting values for the post-counterfactual price solver.
#' @param isMax If `TRUE`, request the legacy local-profit maximum diagnostic.
#' @param control.slopes Optional legacy calibration controls.
#' @param control.equ Optional legacy equilibrium-solver controls.
#' @param labels Product labels.
#' @param ... Additional legacy model options.
#' @return A `MonComLogit` or `MonComCES` object.
#' @export
moncom.logit <- function(
    prices, shares, margins, ownerPre, ownerPost,
    output = TRUE, weights = rep(1, length(shares)),
    normIndex = ifelse(isTRUE(all.equal(sum(shares), 1,
                                        check.names = FALSE)), 1, NA),
    mcDelta = rep(0, length(prices)), subset = rep(TRUE, length(prices)),
    insideSize = NA_real_, priceOutside = 0, priceStart = prices,
    isMax = FALSE, control.slopes, control.equ,
    labels = paste("Prod", 1:length(prices), sep = ""), ...) {
  if (missing(ownerPost)) ownerPost <- ownerPre
  diversions <- matrix(NA_real_, nrow = length(prices), ncol = length(prices))
  result <- new(
    "MonComLogit", prices = prices, shares = shares, margins = margins,
    diversion = diversions, normIndex = normIndex, ownerPre = ownerPre,
    ownerPost = ownerPost, insideSize = insideSize, output = output,
    mcDelta = mcDelta, subset = subset, weights = weights,
    priceOutside = priceOutside, priceStart = priceStart,
    shareInside = sum(shares), labels = labels
  )
  if (!missing(control.slopes)) result@control.slopes <- control.slopes
  if (!missing(control.equ)) result@control.equ <- control.equ
  result@ownerPre <- ownerToMatrix(result, TRUE)
  result@ownerPost <- ownerToMatrix(result, FALSE)
  result <- calcSlopes(result)
  result@pricePre <- prices
  result@mcPre <- calcMC(result, TRUE)
  result@mcPost <- calcMC(result, FALSE)
  result@pricePost <- calcPrices(result, FALSE, isMax = isMax, ...)
  result
}

#' @rdname MonCom-Functions
#' @export
moncom.ces <- function(
    prices, shares, margins, ownerPre, ownerPost,
    output = TRUE, weights = rep(1, length(shares)),
    normIndex = ifelse(isTRUE(all.equal(sum(shares), 1,
                                        check.names = FALSE)), 1, NA),
    mcDelta = rep(0, length(prices)), subset = rep(TRUE, length(prices)),
    insideSize = NA_real_, priceOutside = 1, priceStart = prices,
    isMax = FALSE, control.slopes, control.equ,
    labels = paste("Prod", 1:length(prices), sep = ""), ...) {
  if (missing(ownerPost)) ownerPost <- ownerPre
  diversions <- matrix(NA_real_, nrow = length(prices), ncol = length(prices))
  result <- new(
    "MonComCES", prices = prices, shares = shares, margins = margins,
    diversion = diversions, normIndex = normIndex, ownerPre = ownerPre,
    ownerPost = ownerPost, insideSize = insideSize, output = output,
    mcDelta = mcDelta, subset = subset, weights = weights,
    priceOutside = priceOutside, priceStart = priceStart,
    shareInside = sum(shares), labels = labels
  )
  if (!missing(control.slopes)) result@control.slopes <- control.slopes
  if (!missing(control.equ)) result@control.equ <- control.equ
  result@ownerPre <- ownerToMatrix(result, TRUE)
  result@ownerPost <- ownerToMatrix(result, FALSE)
  result <- calcSlopes(result)
  result@pricePre <- prices
  result@mcPre <- calcMC(result, TRUE)
  result@mcPost <- calcMC(result, FALSE)
  result@pricePost <- calcPrices(result, FALSE, isMax = isMax, ...)
  result
}
