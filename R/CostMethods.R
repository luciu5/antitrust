#' @title Methods for Calculating marginal and Variable Costs
#' @name Cost-Methods
#' @docType methods
#'
#' @aliases calcMC calcMC,ANY-method calcMC,Bertrand-method calcMC,Auction2ndLogit-method calcMC,Cournot-method calcMC,Auction2ndCap-method calcdMC calcdMC,ANY-method calcdMC,Stackelberg-method calcVC calcVC,ANY-method calcVC,Cournot-method
#'
#' @description
#' For Auction2ndCap, calcMC calculates (constant) marginal cost for each
#' product. For those classes that do not require prices, returns a
#' length-k vector of NAs when prices are not supplied.
#'
#' For Bertrand, calcMC computes either pre- or post-merger marginal costs.
#' Marginal costs are assumed to be constant. Post-merger marginal costs are
#' equal to pre-merger marginal costs multiplied by 1+\sQuote{mcDelta}, a
#' length-k vector of marginal cost changes. \sQuote{mcDelta} will typically
#' be between 0 and 1. The second-score Logit auction retains its legacy
#' additive cost-level interpretation of \sQuote{mcDelta}.
#'
#' For Auction2ndLogit, calcMC computes constant marginal costs impied by the model.
#'
#' For Cournot, calcMC calculates marginal cost for each product.
#'
#' calcdMC computes the derivative of either pre- or post-merger marginal costs. The derivative of Marginal costs
#' is assumed to be constant. Post-merger marginal costs are equal to
#' pre-merger marginal costs multiplied by 1+\sQuote{mcDelta}, a length-k
#' vector of marginal cost changes. \sQuote{mcDelta} will typically be between 0 and 1.
#'
#' calcVC computes either pre- or post-merger variable costs. Variable costs
#' are assumed to be quadratic by default. Post-merger variable costs are equal to
#' pre-merger variable costs multiplied by 1+\sQuote{mcDelta}, a length-k
#' vector of marginal cost changes. \sQuote{mcDelta} will typically be between 0
#' and 1.
#'
#' @param object An instance of the respective class (see description for the classes)
#' @param  preMerger If TRUE, the pre-merger ownership structure is used. If FALSE, the post-merger ownership structure is used.
#' Default is TRUE.
#' @param exAnte If \sQuote{exAnte} equals TRUE then the
#' \emph{ex ante} expected result for each firm is produced, while FALSE produces the
#' expected result conditional on each firm winning the auction. Default is FALSE.
#' @param t The capacity profile of each supplier. Default is \sQuote{preMerger} capacities.
#'
#' @include PricesMethods.R
#' @keywords methods
NULL

setGeneric(
  name = "calcMC",
  def = function(object, ...) {
    standardGeneric("calcMC")
  }
)

setGeneric(
  name = "calcdMC",
  def = function(object, ...) {
    standardGeneric("calcdMC")
  }
)

setGeneric(
  name = "calcVC",
  def = function(object, ...) {
    standardGeneric("calcVC")
  }
)

## A calibrated model has two different kinds of cost information.  The
## legacy `mcPre` slot is the value implied by the calibration observations,
## while a counterfactual needs the structural cost primitive that remains
## fixed when ownership, quality, or the active product set changes.  Keep the
## latter as an ordinary attribute so the public S4 classes and their mature
## constructors do not change.  The attribute is copied by S4 value
## semantics, including when a result is promoted into the next path step.
.cost_state_attribute <- "antitrust_cost_state"

.cost_state <- function(object) {
  attr(object, .cost_state_attribute, exact = TRUE)
}

.cost_state_mode <- function(object) {
  ## The second-score Logit family historically treats `mcDelta` as an
  ## additive cost level.  All other constant-cost demand classes use the
  ## documented proportional change.  Keep that one-step convention intact.
  if (methods::is(object, "Auction2ndLogit") &&
      !methods::is(object, "Auction2ndCES")) {
    "additive"
  } else {
    "multiplicative"
  }
}

.set_cost_state <- function(object, state) {
  attr(object, .cost_state_attribute) <- state
  object
}

## Capture the calibrated cost level after a legacy constructor has completed.
## For Cournot/Stackelberg the cost functions themselves are the structural
## primitive; their realized value remains quantity-dependent and is therefore
## deliberately not frozen here.
#' Initialize persistent structural cost state
#'
#' @param object A fitted structural model.
#' @return The model with its structural cost state initialized.
#' @export
initialize_cost_state <- function(object) {
  ## Cournot and Stackelberg marginal costs are functions of equilibrium
  ## quantities.  Their function closures and derivative closures already
  ## live in the S4 slots and are promoted with the quantity state; freezing
  ## the calibration-time evaluation here would change their economics.
  if (methods::is(object, "Cournot")) return(object)
  slots <- methods::slotNames(object)
  if (!"mcPre" %in% slots) return(object)
  base <- methods::slot(object, "mcPre")
  if (!is.numeric(base) || !length(base)) return(object)
  state <- list(
    base = as.numeric(base),
    mode = .cost_state_mode(object)
  )
  .set_cost_state(object, state)
}

.initialize_cost_state <- function(object) initialize_cost_state(object)

#' Compound structural marginal-cost shocks
#'
#' Combines a new marginal-cost change with a cumulative change using the
#' cost convention of `object`. Most structural models use proportional
#' changes; second-score Logit auctions retain their historical additive
#' convention.
#'
#' @param object A fitted structural model whose cost convention is used.
#' @param current The current cumulative cost change, or `NULL`.
#' @param change The next cost change, or `NULL`.
#' @return The compounded cost-change vector, or `current` when `change` is
#'   `NULL`.
#' @export
compound_cost_shocks <- function(object, current = NULL, change = NULL) {
  state <- .cost_state(object)
  mode <- if (!is.null(state) && !is.null(state$mode)) {
    state$mode
  } else {
    .cost_state_mode(object)
  }
  .compound_costs(current, change, mode = mode)
}

## Return a persistent cost level when the model carries one.  A NULL return
## means that the object is a direct legacy S4 object, for which calcMC keeps
## its original inference behavior.
.persistent_mc <- function(object, preMerger) {
  state <- .cost_state(object)
  if (is.null(state) || is.null(state$base)) return(NULL)
  base <- state$base
  if (!is.numeric(base) || !length(base)) return(NULL)
  if (!preMerger) {
    delta <- if ("mcDelta" %in% methods::slotNames(object)) object@mcDelta else NULL
    if (is.null(delta) || length(delta) != length(base)) return(NULL)
    if (identical(state$mode, "additive")) {
      base <- base + delta
    } else {
      base <- base * (1 + delta)
    }
  }
  names(base) <- if ("labels" %in% methods::slotNames(object) &&
                     is.character(object@labels)) object@labels else names(base)
  base
}

## Create a method to recover marginal cost using
## demand parameters and supplied prices
#' @rdname Cost-Methods
#' @export
setMethod(
  f = "calcMC",
  signature = "Bertrand",
  definition = function(object, preMerger = TRUE) {
    persistent <- .persistent_mc(object, preMerger)
    if (!is.null(persistent)) {
      isNegMC <- persistent < 0
      if (preMerger && any(isNegMC, na.rm = TRUE)) {
        warning(paste("Negative marginal costs were calibrated for the following firms:",
                      paste(object@labels[isNegMC], collapse = ",")))
      }
      return(persistent)
    }
    output <- object@output

    object@pricePre <- object@prices


    marginPre <- calcMargins(object, preMerger = TRUE, level = FALSE)

    if (output) {
      mc <- (1 - marginPre) * object@prices
    } else {
      mc <- (1 + marginPre) * object@prices
    }

    if (!preMerger) {
      mc <- mc * (1 + object@mcDelta)
    }

    # mc <- as.vector(mc)

    names(mc) <- object@labels


    isNegMC <- mc < 0

    if (preMerger && any(isNegMC, na.rm = TRUE)) {
      warning(paste("Negative marginal costs were calibrated for the following firms:", paste(object@labels[isNegMC], collapse = ",")))
    }

    return(mc)
  }
)

#' @rdname Cost-Methods
#' @export
setMethod(
  f = "calcMC",
  signature = "Auction2ndCap",
  definition = function(object, t, preMerger = TRUE, exAnte = TRUE) {
    cdfF <- match.fun(object@sellerCostCDF)
    pdfF <- object@sellerCostPDF
    sellerCostBounds <- object@sellerCostBounds


    if (preMerger) {
      capacities <- object@capacities
      r <- object@reservePre
    } else {
      capacities <- tapply(object@capacities * (1 + object@mcDelta), object@ownerPost, sum)
      r <- object@reservePost
    }

    totCap <- sum(capacities)

    if (missing(t)) {
      t <- capacities
    }


    ## The expected production cost
    ecIntegrand <- function(c, t) {
      sellerCostParms <- c(list(c), as.list(object@sellerCostParms))

      fc <- do.call(pdfF, sellerCostParms)

      sellerCostParms <- c(sellerCostParms,
        lower.tail = as.list(object@sellerCostCDFLowerTail)
      )
      Fc <- do.call(cdfF, sellerCostParms)

      retval <- t * c * fc * (1 - Fc)^(totCap - 1)
      retval <- ifelse(is.finite(retval), retval, 0)

      return(retval)
    }

    result <- sapply(
      t,
      function(t.i) {
        if (r < sellerCostBounds[2]) {
          retval <- integrate(ecIntegrand, lower = sellerCostBounds[1], upper = r, stop.on.error = FALSE, t = t.i)$value
        } else {
          retval <- integrate(ecIntegrand, lower = sellerCostBounds[1], upper = sellerCostBounds[2], stop.on.error = FALSE, t = t.i)$value
        }

        return(retval)
      }
    )


    if (!preMerger && length(t) > 1) {
      temp <- rep(NA, length(object@ownerPre))
      temp[object@ownerPre == object@ownerPost] <- result
      result <- temp
    }

    if (!exAnte) {
      result <- result / calcShares(object, preMerger = preMerger, exAnte = TRUE)
    }

    return(result)
  }
)


#' @rdname Cost-Methods
#' @export
setMethod(
  f = "calcMC",
  signature = "Cournot",
  definition = function(object, preMerger = TRUE) {
    if (preMerger) {
      quantity <- object@quantityPre
      mcfun <- object@mcfunPre
      cap <- object@capacitiesPre
    } else {
      quantity <- object@quantityPost
      mcfun <- object@mcfunPost
      cap <- object@capacitiesPost
    }

    plantQuant <- rowSums(quantity, na.rm = TRUE)


    nplants <- nrow(quantity)

    mc <- rep(NA, nplants)

    for (f in 1:nplants) {
      mc[f] <- mcfun[[f]](quantity[f, ])
    }

    if (!preMerger) {
      mc <- mc * (1 + object@mcDelta)
    }

    mc <- mc + 1 / (100 * (pmax(cap - plantQuant, 1e-16))) + 1 / (100 * (pmax(1e-16, plantQuant)))
    # mc <- ifelse(plantQuant <= cap & plantQuant >= 0 , mc, max(mc,na.rm=TRUE) * 1e3)


    names(mc) <- object@labels[[1]]

    return(mc)
  }
)


## Create a method to recover marginal cost using
## demand parameters and supplied prices
#' @rdname Cost-Methods
#' @export
setMethod(
  f = "calcMC",
  signature = "Auction2ndLogit",
  definition = function(object, preMerger = TRUE, exAnte = FALSE) {
    persistent <- .persistent_mc(object, preMerger)
    if (!is.null(persistent)) {
      ## Auction2ndLogit's public legacy method reports the conditional cost
      ## level by default and multiplies by the product share for exAnte.
      mc <- persistent
      if (exAnte) mc <- mc * calcShares(object, preMerger = preMerger)
      return(as.vector(mc))
    }
    prices <- object@prices
    output <- object@output

    marginPre <- calcMargins(object, preMerger = TRUE, level = TRUE)

    if (output) {
      mc <- prices - marginPre
    } else {
      mc <- marginPre + prices
    }

    if (!preMerger) {
      mc <- mc + object@mcDelta
    }

    if (exAnte) {
      mc <- mc * calcShares(object, preMerger = preMerger)
    }

    names(mc) <- object@labels

    mc <- as.vector(mc)

    isNegMC <- mc < 0

    if (preMerger && any(isNegMC, na.rm = TRUE)) {
      warning(paste("Negative marginal costs were calibrated for the following firms:", paste(object@labels[isNegMC], collapse = ",")))
    }

    return(mc)
  }
)

#' @rdname Cost-Methods
#' @export
setMethod(
  f = "calcMC",
  signature = "Auction2ndCES",
  definition = function(object, preMerger = TRUE) {
    persistent <- .persistent_mc(object, preMerger)
    if (!is.null(persistent)) return(as.vector(persistent))
    prices <- object@prices
    output <- object@output

    object@pricePre <- prices

    marginPre <- calcMargins(object, preMerger = TRUE, level = FALSE)

    if (output) {
      mc <- prices * (1 - marginPre)
    } else {
      mc <- prices * (1 + marginPre)
    }

    if (!preMerger) {
      mc <- mc * (1 + object@mcDelta)
    }

    names(mc) <- object@labels

    mc <- as.vector(mc)

    isNegMC <- mc < 0

    if (preMerger && any(isNegMC, na.rm = TRUE)) {
      warning(paste("Negative marginal costs were calibrated for the following firms:", paste(object@labels[isNegMC], collapse = ",")))
    }

    return(mc)
  }
)

#' @rdname Cost-Methods
#' @export
setMethod(
  f = "calcdMC",
  signature = "Stackelberg",
  definition = function(object, preMerger = TRUE) {
    if (preMerger) {
      quantity <- object@quantityPre
      dmcfun <- object@dmcfunPre
    } else {
      quantity <- object@quantityPost
      dmcfun <- object@dmcfunPost
    }


    nplants <- nrow(quantity)

    dmc <- rep(NA, nplants)

    for (f in 1:nplants) {
      dmc[f] <- dmcfun[[f]](quantity[f, ])
    }

    if (!preMerger) {
      dmc <- dmc * (1 + object@mcDelta)
    }

    names(dmc) <- object@labels[[1]]

    return(dmc)
  }
)


#' @rdname Cost-Methods
#' @export
setMethod(
  f = "calcVC",
  signature = "Cournot",
  definition = function(object, preMerger = TRUE) {
    if (preMerger) {
      quantity <- object@quantityPre
      vcfun <- object@vcfunPre
    } else {
      quantity <- object@quantityPost
      vcfun <- object@vcfunPost
    }


    nplants <- nrow(quantity)

    vc <- rep(NA, nplants)

    for (f in 1:nplants) {
      vc[f] <- vcfun[[f]](quantity[f, ])
    }

    if (!preMerger) {
      vc <- vc * (1 + object@mcDelta)
    }

    names(vc) <- object@labels[[1]]

    return(vc)
  }
)
