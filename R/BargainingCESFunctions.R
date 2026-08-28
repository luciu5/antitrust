#' @title Nash Bargaining Model with CES Demand
#' @name BargainingCES-Functions
#' @aliases bargaining.ces bargaining2nd.ces bargaining.ces.alm
#' @description Calibrates consumer demand using Constant Elasticity of Substitution (CES) and then
#' simulates the price effect of a merger between two firms
#' under the assumption that firms and customers in the market are playing a
#' differentiated products Nash Bargaining game.
#' @description Let k denote the number of products produced by all firms playing the
#' Nash Bargaining game below.
#' @param prices A length k vector of product prices.
#' @param shares A length k vector of product revenue shares. Values must be
#' between 0 and 1.
#' @param margins A length k vector of product margins (in levels, not percents), some of which may
#' equal NA.
#' @param normIndex An integer equalling the index (position) of the
#' inside product whose price will be normalized to 1. Default
#' is 1, unless \sQuote{shares} sum to less than 1, in which case the default is
#' NA and an outside good is assumed to exist.
#' @param ownerPre EITHER a vector of length k whose values
#' indicate which firm produced a product pre-merger OR
#' a k x k matrix of pre-merger ownership shares.
#' @param ownerPost EITHER a vector of length k whose values
#' indicate which firm produced a product after the merger OR
#' a k x k matrix of post-merger ownership shares.
#' @param output a length 1 logical vector equal to TRUE if merger simulation performed
#' on output market, FALSE if simulation performed on input market. Default TRUE.
#' @param bargpowerPre A length k vector of pre-merger bargaining power parameters. Values
#' must be between 0 (sellers have the power) and 1 (buyers the power). NA values are allowed,
#' though must be calibrated from additional margin and share data. Default is 0.5.
#' @param bargpowerPost A length k vector of post-merger bargaining power parameters. Values
#' must be between 0 (sellers have the power) and 1 (buyers the power). NA values are allowed,
#' though must be calibrated from additional margin and share data. Default is \sQuote{bargpowerPre}.
#' @param insideSize An integer equal to total pre-merger revenue spent inside market.
#' If shares sum to one, this also equals the size of the market.
#' @param mcDelta A vector of length k where each element equals the
#' (level) change in a product's marginal costs due to
#' the merger. Default is 0, which assumes that the merger does not
#' affect any products' marginal cost.
#' @param mcDeltaOutside A length 1 vector indicating the change in the marginal cost of the
#' outside good. Default is 0.
#' @param subset A vector of length k where each element equals TRUE if
#' the product indexed by that element should be included in the
#' post-merger simulation and FALSE if it should be excluded. Default is a
#' length k vector of TRUE.
#' @param priceOutside A positive real number equal to the price of the outside good.
#' Default is 1.
#' @param priceStart A length k vector of starting values for price calculation.
#' @param parmsStart A length 2 vector of starting values for ALM demand parameter estimation.
#' @param weights A length k vector of product weights. Default is rep(1, length(shares)).
#' @param solver A length-1 character vector specifying the solver algorithm used to calculate pre- and post-merger price equilibria. Options are \code{"nleqslv"} (default) or \code{"ag"} (Aggregative Games).
#' @param control.slopes A list of control parameters for non-linear optimization routine.
#' @param control.equ A list of control parameters for non-linear equation solver routine.
#' @param labels A length-k vector of product labels.
#'
#' @details
#' \code{bargaining.ces} uses observed product prices, revenue shares, margins, and pre-/post-merger
#' ownership matrices to calibrate a Constant Elasticity of Substitution (CES) demand system and
#' predict post-merger equilibrium outcomes under Nash Bargaining.
#'
#' \code{bargaining2nd.ces} calculates post-merger price effects under a 2nd-score auction Nash bargaining model.
#'
#' \code{bargaining.ces.alm} estimates the unobserved outside share alongside CES parameters.
#'
#' @return \code{bargaining.ces} returns an instance of \code{\linkS4class{BargainingCES}},
#' \code{bargaining2nd.ces} returns an instance of \code{\linkS4class{Bargaining2ndCES}},
#' \code{bargaining.ces.alm} returns an instance of \code{\linkS4class{BargainingCESALM}}.
#'
#' @author Charles Taragin
#' @seealso \code{\link{ces}} for Bertrand pricing under CES demand, and \code{\link{bargaining.logit}} for Nash Bargaining under Logit demand.
#' @export
bargaining.ces <- function(prices, shares, margins,
                           ownerPre, ownerPost,
                           bargpowerPre = rep(0.5, length(prices)),
                           bargpowerPost = bargpowerPre,
                           output = TRUE,
                           weights = rep(1, length(shares)),
                           normIndex = ifelse(isTRUE(all.equal(sum(shares), 1, check.names = FALSE)), 1, NA),
                           mcDelta = rep(0, length(prices)),
                           subset = rep(TRUE, length(prices)),
                           priceStart = prices,
                           insideSize = NA_real_,
                           priceOutside = 1,
                           control.slopes,
                           control.equ,
                           solver = c("nleqslv", "ag"),
                           labels = paste("Prod", 1:length(prices), sep = "")) {
  solver <- match.arg(solver)

  ## Create BargainingCES container to store relevant data
  result <- new("BargainingCES",
    prices = prices, shares = shares,
    margins = margins,
    normIndex = normIndex,
    ownerPre = ownerPre,
    ownerPost = ownerPost,
    bargpowerPre = bargpowerPre,
    bargpowerPost = bargpowerPost,
    output = output,
    insideSize = insideSize,
    mcDelta = mcDelta,
    subset = subset,
    weights = weights,
    priceOutside = priceOutside,
    shareInside = ifelse(isTRUE(all.equal(sum(shares), 1, check.names = FALSE)), 1, sum(shares)),
    priceStart = priceStart,
    labels = labels,
    cls = "BargainingCES"
  )

  if (!missing(control.slopes)) {
    result@control.slopes <- control.slopes
  }

  if (!missing(control.equ)) {
    result@control.equ <- control.equ
  }

  ## Convert ownership vectors to ownership matrices
  result@ownerPre <- ownerToMatrix(result, TRUE)
  result@ownerPost <- ownerToMatrix(result, FALSE)

  ## Calculate Demand Slope Coefficients
  result <- calcSlopes(result)

  ## Calculate marginal cost
  result@mcPre <- calcMC(result, TRUE)
  result@mcPost <- calcMC(result, FALSE)

  ## Solve Non-Linear System for Price Changes
  if (solver == "ag") {
    result@pricePre  <- calcPricesAG(result, preMerger = TRUE)
    result@pricePost <- calcPricesAG(result, preMerger = FALSE)
  } else {
    result@pricePre  <- calcPrices(result, preMerger = TRUE)
    result@pricePost <- calcPrices(result, preMerger = FALSE)
  }

  return(result)
}


#' @rdname BargainingCES-Functions
#' @export
bargaining2nd.ces <- function(prices, shares, margins,
                             ownerPre, ownerPost,
                             bargpowerPre = rep(0.5, length(prices)),
                             bargpowerPost = bargpowerPre,
                             output = TRUE,
                             weights = rep(1, length(shares)),
                             normIndex = ifelse(isTRUE(all.equal(sum(shares), 1, check.names = FALSE)), 1, NA),
                             mcDelta = rep(0, length(prices)),
                             subset = rep(TRUE, length(prices)),
                             insideSize = NA_real_,
                             mcDeltaOutside = 0,
                             control.slopes,
                             labels = paste("Prod", 1:length(prices), sep = "")) {
  if (missing(prices)) {
    prices <- rep(NA_integer_, length(shares))
  }
  ## Create Bargaining2ndCES container to store relevant data
  result <- new("Bargaining2ndCES",
    prices = prices, shares = shares,
    margins = margins,
    normIndex = normIndex,
    ownerPre = ownerPre,
    ownerPost = ownerPost,
    bargpowerPre = bargpowerPre,
    bargpowerPost = bargpowerPost,
    output = output,
    insideSize = insideSize,
    mcDelta = mcDelta,
    subset = subset,
    weights = weights,
    priceOutside = 1,
    priceStart = prices,
    shareInside = ifelse(isTRUE(all.equal(sum(shares), 1, check.names = FALSE)), 1, sum(shares)),
    labels = labels,
    cls = "Bargaining2ndCES"
  )

  if (!missing(control.slopes)) {
    result@control.slopes <- control.slopes
  }

  ## Convert ownership vectors to ownership matrices
  result@ownerPre <- ownerToMatrix(result, TRUE)
  result@ownerPost <- ownerToMatrix(result, FALSE)

  ## Calculate Demand Slope Coefficients
  result <- calcSlopes(result)

  ## Calculate marginal cost
  result@mcPre <- calcMC(result, TRUE)
  result@mcPost <- calcMC(result, FALSE)

  ## Solve Non-Linear System for Price Changes
  result@pricePre <- calcPrices(result, preMerger = TRUE)
  result@pricePost <- calcPrices(result, preMerger = FALSE)

  return(result)
}


#' @rdname BargainingCES-Functions
#' @export
bargaining.ces.alm <- function(prices, shares, margins,
                               ownerPre, ownerPost,
                               bargpowerPre = rep(0.5, length(prices)),
                               bargpowerPost = bargpowerPre,
                               output = TRUE,
                               weights = rep(1, length(shares)),
                               parmsStart = c(0.1, 0.5),
                               mcDelta = rep(0, length(prices)),
                               subset = rep(TRUE, length(prices)),
                               priceStart = prices,
                               insideSize = NA_real_,
                               priceOutside = 1,
                               control.slopes,
                               control.equ,
                               labels = paste("Prod", 1:length(prices), sep = "")) {
  ## Create BargainingCESALM container to store relevant data
  result <- new("BargainingCESALM",
    prices = prices, shares = shares,
    margins = margins,
    normIndex = NA_integer_,
    ownerPre = ownerPre,
    ownerPost = ownerPost,
    bargpowerPre = bargpowerPre,
    bargpowerPost = bargpowerPost,
    output = output,
    insideSize = insideSize,
    mcDelta = mcDelta,
    subset = subset,
    weights = weights,
    priceOutside = priceOutside,
    parmsStart = parmsStart,
    shareInside = 1,
    priceStart = priceStart,
    labels = labels,
    cls = "BargainingCESALM"
  )

  if (!missing(control.slopes)) {
    result@control.slopes <- control.slopes
  }

  if (!missing(control.equ)) {
    result@control.equ <- control.equ
  }

  ## Convert ownership vectors to ownership matrices
  result@ownerPre <- ownerToMatrix(result, TRUE)
  result@ownerPost <- ownerToMatrix(result, FALSE)

  ## Calculate Demand Slope Coefficients
  result <- calcSlopes(result)

  ## Calculate marginal cost
  result@mcPre <- calcMC(result, TRUE)
  result@mcPost <- calcMC(result, FALSE)

  ## Solve Non-Linear System for Price Changes
  result@pricePre <- calcPrices(result, preMerger = TRUE)
  result@pricePost <- calcPrices(result, preMerger = FALSE)

  return(result)
}
