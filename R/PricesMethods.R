#' @title \dQuote{Calculating Prices} Methods
#' @name Prices-Methods
#' @docType methods
#' @aliases calcPrices,ANY-method calcPrices,Auction2ndCap-method calcPrices,Cournot-method calcPrices,Linear-method calcPrices,Logit-method calcPrices,LogitBLP-method calcPrices,LogLin-method calcPrices,AIDS-method calcPrices,LogitCap-method calcPrices,Auction2ndLogit-method calcPrices,BargainingBLP-method calcPrices
#' @param object An instance of the respective class (see description for the classes)
#' @param  preMerger If TRUE, the pre-merger ownership structure is used. If FALSE, the post-merger ownership structure is used.
#' Default is TRUE.
#' @param exAnte If \sQuote{exAnte} equals TRUE then the
#' \emph{ex ante} expected result for each firm is produced, while FALSE produces the
#' expected result conditional on each firm winning the auction. Default is FALSE.
#' @param subset A vector of length k where each element equals TRUE if
#' the product indexed by that element should be included in the
#' post-merger simulation and FALSE if it should be excluded. Default is a
#' length k vector of TRUE.
#' @param isMax If TRUE, a check is run to determine if the calculated equilibrium price vector locally maximizes profits.
#' Default is FALSE.
#'
#' @param ... For Logit, additional values that may be used to change the
#' default values of \code{\link[BB]{BBsolve}}, the non-linear equation solver.
#'
#' For others, additional values that may be used to change the default values of \code{\link[stats]{constrOptim}}, the non-linear
#' equation solver used to enforce non-negative equilibrium quantities.
#' @description For Auction2ndCap, the calcPrices method computes the expected price that the buyer pays,
#' conditional on the buyer purchasing from a particular firm.
#' @description For Logit, the calcPrices method computes either pre-merger or post-merger equilibrium prices under the assumptions
#' that consumer demand is Logit and firms play a differentiated product Bertrand Nash pricing game.
#' @description For LogitCap, the calcPrices method computes either pre-merger or post-merger equilibrium shares under the assumptions that
#' consumer demand is Logit and firms play a differentiated product Bertrand Nash pricing game with capacity constraints.
#' @description For Logit, the calcPrices method computes either pre-merger or post-merger equilibrium prices under the assumptions
#' that consumer demand is Logit and firms play a differentiated product Bertrand Nash pricing game.
#' @description For LogLin, the calcPrices method computes either pre-merger or post-merger equilibrium prices under the assumptions
#' that consumer demand is Log-Linear and firms play a differentiated product Bertrand Nash pricing game.
#' @description For AIDS, the calcPrices method computes either pre-merger or post-merger equilibrium prices under the assumptions
#' that consumer demand is AIDS and firms play a differentiated product Bertrand Nash pricing game.
#' It returns a length-k vector of NAs if the user did not supply prices.
#'
#' @include AuctionCapMethods.R
#' @include BargainingClasses.R
#' @keywords methods
NULL
setGeneric(
  name = "calcPrices",
  def = function(object, ...) {
    standardGeneric("calcPrices")
  }
)
#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "Cournot",
  definition = function(object, preMerger = TRUE) {
    if (preMerger) {
      quantities <- object@quantityPre
    } else {
      quantities <- object@quantityPost
    }
    intercepts <- object@intercepts
    slopes <- object@slopes
    mktQuant <- colSums(quantities, na.rm = TRUE)
    prices <- ifelse(object@demand == "linear",
      intercepts + slopes * mktQuant,
      exp(intercepts) * mktQuant^slopes
    )
    names(prices) <- object@labels[[2]]
    return(prices)
  }
)
#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "Auction2ndCap",
  definition = function(object, preMerger = TRUE, exAnte = TRUE) {
    val <- calcProducerSurplus(object, preMerger = preMerger, exAnte = exAnte) + calcMC(object, preMerger = preMerger, exAnte = exAnte)
    names(val) <- object@labels
    return(val)
  }
)
#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "Logit",
  definition = function(object, preMerger = TRUE, isMax = FALSE, subset, ...) {
    output <- object@output


    priceStart <- object@priceStart
    outSign <- ifelse(output, -1, 1)
    nprods <- length(object@shares)

    if (preMerger) {
      owner <- object@ownerPre
      mc <- object@mcPre
    } else {
      owner <- object@ownerPost
      mc <- object@mcPost
    }

    if (missing(subset)) {
      subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    }
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }

    ## A homogeneous CES market with no outside good has constant total
    ## revenue.  If one firm owns every product, its Bertrand profit is
    ## unbounded above and there is no finite price equilibrium.  Returning a
    ## very large stalled price vector is more dangerous than a clear error.
    owner_check <- if (preMerger) object@ownerPre else object@ownerPost
    plain_ces <- class(object)[1] %in% c("CES", "CESNests")
    if (plain_ces && object@output && !is.na(object@normIndex) &&
      nrow(owner_check) > 1 && all(owner_check == 1)) {
      stop("CES demand without an outside good has no finite Bertrand equilibrium when one firm owns every product.")
    }


    priceStart <- priceStart[subset]


    priceEst <- rep(NA, nprods)

    FOC <- function(priceCand) {
      thisPrice <- priceEst
      thisPrice[subset] <- priceCand
      if (preMerger) {
        object@pricePre <- thisPrice
        mc <- object@mcPre[subset]
      } else {
        object@pricePost <- thisPrice
        mc <- object@mcPost[subset]
      }

      if (output) {
        margins <- priceCand - mc
      } else {
        margins <- mc - priceCand
      }

      predMargin <- calcMargins(object, preMerger, level = TRUE)[subset]

      thisFOC <- margins - predMargin
      return(thisFOC)
    }


    ## Find price changes that set FOCs equal to 0
    ## Try nleqslv first (faster, more reliable for smooth FOCs)
    nleqslv_maxit <- as.integer(object@control.equ$maxit)
    if (length(nleqslv_maxit) == 0 || is.na(nleqslv_maxit[1]) || nleqslv_maxit[1] < 1) nleqslv_maxit <- 150L
    minResult <- nleqslv::nleqslv(priceStart, FOC,
      method = "Broyden",
      control = list(
        ftol = object@control.equ$tol,
        maxit = nleqslv_maxit
      )
    )

    ## Fallback to BBsolve if nleqslv fails
    if (minResult$termcd > 2) {
      minResult <- BBsolve(priceStart, FOC, quiet = TRUE, control = object@control.equ, ...)
      priceEst_solution <- minResult$par
      if (minResult$convergence != 0) {
        warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBsolve' reports: '", minResult$message, "'")
      }
    } else {
      priceEst_solution <- minResult$x
      if (minResult$termcd > 1) {
        warning("'calcPrices' may not have fully converged. 'nleqslv' termcd: ", minResult$termcd)
      }
    }
    if (isMax) {
      hess <- genD(FOC, priceEst_solution) # compute the numerical approximation of the FOC hessian at optimium
      hess <- hess$D[, 1:hess$p]
      hess <- hess * (owner > 0) # 0 terms not under the control of a common owner
      state <- ifelse(preMerger, "Pre-merger", "Post-merger")
      if (any(eigen(hess)$values > 0)) {
        warning("Hessian of first-order conditions is not positive definite. ", state, " price vector may not maximize profits. Consider rerunning 'calcPrices' using different starting values")
      }
    }
    priceEst[subset] <- priceEst_solution

    names(priceEst) <- object@labels
    return(priceEst)
  }
)


#' CES Cournot solves a positive-price interior equilibrium. If high
#' product-level retention differences require a boundary quantity or
#' endogenous product exit, the method reports that no valid interior root
#' was found; callers can specify an explicit product subset.
#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "CESCournot",
  definition = function(object, preMerger = TRUE, isMax = FALSE, subset, ...) {
    n <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, n) else object@subset
    if (!is.logical(subset) || length(subset) != n || !any(subset)) {
      stop("'subset' must select at least one product.")
    }
    mc <- if (preMerger) object@mcPre else object@mcPost
    if (any(!is.finite(mc[subset]))) {
      stop("CES Cournot has non-finite effective costs.")
    }
    residual <- function(prices) {
      if (length(prices) != sum(subset) || any(!is.finite(prices)) ||
          any(prices <= 0)) return(rep(1e6, sum(subset)))
      candidate <- rep(NA_real_, n)
      candidate[subset] <- prices
      trial <- object
      if (preMerger) trial@pricePre <- candidate else trial@pricePost <- candidate
      margin <- try(calcMargins(trial, preMerger, level = TRUE)[subset],
                    silent = TRUE)
      if (inherits(margin, "try-error") || any(!is.finite(margin))) {
        return(rep(1e6, sum(subset)))
      }
      (prices - mc[subset] - margin) / pmax(prices, 1)
    }
    if (.mixed_firm_retention(object, preMerger, subset)) {
      start <- if (!preMerger && all(is.finite(object@pricePre[subset]))) {
        object@pricePre[subset]
      } else object@priceStart[subset]
      start <- pmax(start, .Machine$double.eps^0.25)
      foc <- function(z) {
        if (any(!is.finite(z)) || any(z > 100)) {
          return(rep(1e6, length(z)))
        }
        residual(exp(z))
      }
      maxit <- as.integer(object@control.equ$maxit)
      if (length(maxit) != 1L || !is.finite(maxit) || maxit < 1L) maxit <- 300L
      tol <- object@control.equ$tol
      if (length(tol) != 1L || !is.finite(tol) || tol <= 0) tol <- 1e-10
      solution <- try(nleqslv::nleqslv(log(start), foc,
          control = list(ftol = tol, maxit = maxit)), silent = TRUE)
      z <- if (inherits(solution, "try-error")) log(start) else solution$x
      if (any(!is.finite(foc(z))) || max(abs(foc(z))) > 1e-8) {
        alternative <- try(BB::BBsolve(z, foc,
            control = list(tol = tol, maxit = maxit), quiet = TRUE),
            silent = TRUE)
        if (!inherits(alternative, "try-error") &&
            all(is.finite(alternative$par)) &&
            max(abs(foc(alternative$par))) < max(abs(foc(z)))) {
          z <- alternative$par
        }
      }
      prices <- rep(NA_real_, n)
      prices[subset] <- exp(z)
    } else {
      prices <- callNextMethod(object, preMerger = preMerger,
                              isMax = isMax, subset = subset, ...)
    }
    if (any(!is.finite(prices[subset])) || any(prices[subset] <= 0) ||
        max(abs(residual(prices[subset]))) > 1e-7) {
      stop("CES Cournot failed to find a valid positive interior equilibrium; a boundary or product-exit solution may be required.")
    }
    names(prices) <- object@labels
    prices
  }
)

#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "Auction2ndLogit",
  definition = function(object, preMerger = TRUE, exAnte = FALSE) {
    nprods <- length(object@shares)
    output <- object@output

    if (preMerger) {
      owner <- object@ownerPre
      mc <- object@mcPre
    } else {
      owner <- object@ownerPost
      mc <- object@mcPost
    }
    margins <- calcMargins(object, preMerger, exAnte = FALSE)

    if (output) {
      prices <- margins + mc
    } else {
      prices <- mc - margins
    }

    if (exAnte) {
      prices <- prices * calcShares(object, preMerger = preMerger, revenue = FALSE)
    }
    names(prices) <- object@labels
    return(prices)
  }
)
#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "Auction2ndCES",
  definition = function(object, preMerger = TRUE, exAnte = FALSE) {
    nprods <- length(object@shares)
    output <- object@output

    if (preMerger) {
      owner <- object@ownerPre
      mc <- object@mcPre
      subset <- rep(TRUE, nprods)
    } else {
      owner <- object@ownerPost
      mc <- object@mcPost
      subset <- object@subset
    }

    priceStart <- object@priceStart
    if (!preMerger) {
      cand <- object@pricePre
      if (length(cand) == nprods && all(is.finite(cand))) priceStart <- cand
    }
    priceStart <- priceStart[subset]

    FOC <- function(priceCand) {
      if (preMerger) {
        object@pricePre[subset] <- priceCand
      } else {
        object@pricePost[subset] <- priceCand
      }
      margins_prop <- calcMargins(object, preMerger, exAnte = FALSE, level = FALSE)[subset]

      if (output) {
        return(priceCand - mc[subset] / (1 - margins_prop))
      } else {
        return(priceCand - mc[subset] / (1 + margins_prop))
      }
    }

    nleqslv_maxit <- as.integer(object@control.equ$maxit)
    if (length(nleqslv_maxit) == 0 || is.na(nleqslv_maxit[1]) || nleqslv_maxit[1] < 1) nleqslv_maxit <- 150L
    minResult <- nleqslv::nleqslv(priceStart, FOC,
      method = "Broyden",
      control = list(
        ftol = object@control.equ$tol,
        maxit = nleqslv_maxit
      )
    )

    if (minResult$termcd > 2) {
      minResult <- BBsolve(priceStart, FOC, quiet = TRUE, control = object@control.equ)
      priceEst_solution <- minResult$par
      if (minResult$convergence != 0) {
        warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBsolve' reports: '", minResult$message, "'")
      }
    } else {
      priceEst_solution <- minResult$x
      if (minResult$termcd > 1) {
        warning("'calcPrices' may not have fully converged. 'nleqslv' termcd: ", minResult$termcd)
      }
    }

    prices <- rep(NA, nprods)
    prices[subset] <- priceEst_solution

    if (exAnte) {
      prices <- prices * calcShares(object, preMerger = preMerger, revenue = FALSE)
    }
    names(prices) <- object@labels
    return(prices)
  }
)
#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "LogitCap",
  definition = function(object, preMerger = TRUE, isMax = FALSE, subset, ...) {
    output <- object@output
    # We'll pick a robust starting vector after we know nprods; default to priceStart
    priceStart <- object@priceStart
    # For post-merger, provide a stronger starting point from pre-merger prices
    priceStart <- if (preMerger) object@priceStart else object@pricePre
    if (preMerger) {
      owner <- object@ownerPre
      mc <- object@mcPre
      capacities <- object@capacitiesPre
    } else {
      owner <- object@ownerPost
      mc <- object@mcPost
      capacities <- object@capacitiesPost
    }
    nprods <- length(object@shares)
    # For post-merger, prefer pre-merger prices if available and finite; otherwise fall back
    if (!preMerger) {
      cand <- object@pricePre
      if (length(cand) == nprods && all(is.finite(cand))) {
        priceStart <- cand
      }
    }
    if (missing(subset)) {
      subset <- rep(TRUE, nprods)
    }
    if (!is.logical(subset) || length(subset) != nprods) {
      stop("'subset' must be a logical vector the same length as 'shares'")
    }
    if (any(!subset)) {
      owner <- owner[subset, subset]
      mc <- mc[subset]
      priceStart <- priceStart[subset]
      capacities <- capacities[subset]
    }
    owner <- .retention_owner_bertrand(owner,
                                      getRetention(object, preMerger)[subset])
    priceEst <- rep(NA, nprods)
    ## Define system of FOC as a function of prices
    FOC <- function(priceCand) {
      if (preMerger) {
        object@pricePre[subset] <- priceCand
      } else {
        object@pricePost[subset] <- priceCand
      }

      quantities <- calcQuantities(object, preMerger = preMerger)
      quantities <- quantities[subset]
      if (output) {
        margins <- 1 - mc / priceCand
        revenues <- calcShares(object, preMerger = preMerger,
                               revenue = TRUE)[subset]
        elasticities <- elast(object, preMerger)[subset, subset]
        thisFOC <- revenues * diag(owner) +
          as.vector((t(elasticities) * owner) %*% (margins * revenues))
      } else {
        shares <- calcShares(object, preMerger = preMerger,
                             revenue = FALSE)[subset]
        demand_jacobian <- object@mktSize * object@slopes$alpha *
          (diag(shares, length(shares)) - tcrossprod(shares))
        thisFOC <- quantities -
          as.vector((t(demand_jacobian) * owner) %*% (mc - priceCand))
      }
      constraint <- ifelse(is.finite(capacities), (quantities - capacities) / object@insideSize, 0)
      ## Fischer-Burmeister complementarity residual for finite capacities.
      ## A positive infinite capacity is the unconstrained-product sentinel,
      ## so those products must retain the ordinary FOC rather than a
      ## complementarity residual with a zero placeholder constraint.
      finite_capacity <- is.finite(capacities)
      complementarity <- thisFOC + constraint +
        sqrt(thisFOC^2 + constraint^2)
      measure <- ifelse(finite_capacity, complementarity, thisFOC)
      return(measure)
    }
    ## Find price changes that set FOCs equal to 0
    ## Try nleqslv first (faster for smooth FOCs); fallback to BBsolve for non-smooth constraints
    nleqslv_maxit <- as.integer(object@control.equ$maxit)
    if (length(nleqslv_maxit) == 0 || is.na(nleqslv_maxit[1]) || nleqslv_maxit[1] < 1) nleqslv_maxit <- 150L
    minResult <- nleqslv::nleqslv(priceStart, FOC,
      method = "Broyden",
      control = list(
        ftol = object@control.equ$tol,
        maxit = nleqslv_maxit
      )
    )

    if (minResult$termcd > 2) {
      bb_control <- object@control.equ
      bb_control$price_domain <- NULL
      minResult <- BBsolve(priceStart, FOC, quiet = TRUE, control = bb_control, ...)
      priceEst_solution <- minResult$par
      if (minResult$convergence != 0) {
        warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBsolve' reports: '", minResult$message, "'")
      }
    } else {
      priceEst_solution <- minResult$x
      if (minResult$termcd > 1) {
        warning("'calcPrices' may not have fully converged. 'nleqslv' termcd: ", minResult$termcd)
      }
    }
    if (any(!is.finite(priceEst_solution))) {
      stop("'calcPrices' returned non-finite LogitCap prices; no valid equilibrium was found.")
    }
    signed_input <- !isTRUE(object@output) &&
      identical(object@control.equ$price_domain, "real")
    if (!signed_input && any(priceEst_solution <= 0)) {
      if (!isTRUE(object@output)) {
        i <- which.min(priceEst_solution)
        stop(structure(list(
          message = "Input LogitCap equilibrium has a nonpositive rate; use price_domain = 'real' to retain signed rates.",
          call = NULL, category = "positive_domain_violation",
          price_domain = "positive", minimum_rate = priceEst_solution[[i]],
          minimum_product = as.character(object@labels[which(subset)[[i]]]),
          preMerger = preMerger),
          class = c("antitrust_price_domain_error", "error", "condition")))
      }
      stop("'calcPrices' returned non-positive or non-finite LogitCap prices; no valid equilibrium was found.")
    }
    priceEst[subset] <- priceEst_solution
    names(priceEst) <- object@labels
    if (isMax) {
      hess <- genD(FOC, priceEst_solution) # compute the numerical approximation of the FOC hessian at optimium
      hess <- hess$D[, 1:hess$p]
      hess <- hess * (owner > 0) # 0 terms not under the control of a common owner
      state <- ifelse(preMerger, "Pre-merger", "Post-merger")
      if (any(eigen(hess)$values > 0)) {
        warning("Hessian of first-order conditions is not positive definite. ", state, " price vector may not maximize profits. Consider rerunning 'calcPrices' using different starting values")
      }
    }
    return(priceEst)
  }
)
#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "Linear",
  definition = function(object, preMerger = TRUE, subset, ...) {
    slopes <- object@slopes
    intercept <- object@intercepts
    if (preMerger) {
      owner <- object@ownerPre
      mc <- object@mcPre
    } else {
      owner <- object@ownerPost
      mc <- object@mcPost
    }
    nprods <- length(object@quantities)
    if (missing(subset)) {
      subset <- rep(TRUE, nprods)
    }
    if (!is.logical(subset) || length(subset) != nprods) {
      stop("'subset' must be a logical vector the same length as 'quantities'")
    }
    owner <- .retention_owner_bertrand(owner,
                                      getRetention(object, preMerger))

    ## First try the closed-form Bertrand FOC solution for linear demand.
    analytic <- try(
      solve((slopes * diag(owner)) + (t(slopes) * owner)) %*%
        ((t(slopes) * owner) %*% mc - (intercept * diag(owner))),
      silent = TRUE
    )

    if (!any(class(analytic) == "try-error")) {
      prices <- as.vector(analytic)
      quantities <- as.vector(intercept + slopes %*% prices)

      if (all(subset) && all(is.finite(prices)) && all(quantities >= -1e-8, na.rm = TRUE)) {
        names(prices) <- object@labels
        return(prices)
      }
    }

    FOC <- function(priceCand) {
      if (preMerger) {
        object@pricePre <- priceCand
      } else {
        object@pricePost <- priceCand
      }
      margins <- priceCand - mc
      quantities <- calcQuantities(object, preMerger)
      thisFOC <- quantities * diag(owner) + (t(slopes) * owner) %*% margins
      thisFOC[!subset] <- quantities[!subset] # set quantity equal to 0 for firms not in subset
      return(as.vector(crossprod(thisFOC)))
    }

    isInterior <- function(priceCand) {
      all(is.finite(priceCand)) &&
        all(as.vector(intercept + slopes %*% priceCand) > 1e-8, na.rm = TRUE)
    }

    candidates <- list(object@priceStart)
    if (!any(class(analytic) == "try-error")) candidates <- c(candidates, list(as.vector(analytic)))
    if (length(object@prices) == nprods) candidates <- c(candidates, list(object@prices))

    targetQuantities <- pmax(as.vector(intercept + slopes %*% object@priceStart), 1e-4)
    feasibleFromQuantities <- try(as.vector(solve(slopes, targetQuantities - intercept)), silent = TRUE)
    if (any(class(feasibleFromQuantities) == "try-error")) {
      feasibleFromQuantities <- as.vector(MASS::ginv(slopes) %*% (targetQuantities - intercept))
    }
    candidates <- c(candidates, list(feasibleFromQuantities))

    startIndex <- which(vapply(candidates, isInterior, logical(1)))[1]
    if (is.na(startIndex)) {
      stop("Unable to find a feasible starting price vector for constrained linear equilibrium solve")
    }

    minResult <- constrOptim(candidates[[startIndex]], FOC, grad = NULL, ui = slopes, ci = -intercept, ...)
    if (!isTRUE(all.equal(minResult$convergence, 0, check.names = FALSE))) {
      warning("'calcPrices' solver may not have successfully converged.'constrOptim' reports: '", minResult$message, "'")
    }
    prices <- minResult$par
    names(prices) <- object@labels
    return(prices)
  }
)
#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "LogLin",
  definition = function(object, preMerger = TRUE, subset, ...) {
    # Use better starting values for post-merger; fallback to priceStart if missing/invalid
    priceStart <- object@priceStart
    if (preMerger) {
      owner <- object@ownerPre
      mc <- object@mcPre
    } else {
      owner <- object@ownerPost
      mc <- object@mcPost
    }
    nprods <- length(object@quantities)
    if (missing(subset)) {
      subset <- rep(TRUE, nprods)
    }
    if (!is.logical(subset) || length(subset) != nprods) {
      stop("'subset' must be a logical vector the same length as 'quantities'")
    }
    owner <- .retention_owner_bertrand(owner,
                                      getRetention(object, preMerger))
    if (!preMerger) {
      cand <- object@pricePre
      if (length(cand) == nprods && all(is.finite(cand))) priceStart <- cand
    }
    if (object@output && any(diag(object@slopes) >= -1)) {
      stop("LogLin output demand requires every own-price slope to be strictly below -1 for a finite interior price equilibrium.")
    }
    if (object@output && any(!is.finite(mc) | mc <= 0)) {
      stop("LogLin output demand requires strictly positive finite marginal costs for a finite interior price equilibrium.")
    }
    FOC <- function(priceCand) {
      if (preMerger) {
        object@pricePre <- priceCand
      } else {
        object@pricePost <- priceCand
      }
      margins <- 1 - mc / priceCand
      quantities <- calcQuantities(object, preMerger, revenue = TRUE)
      revenues <- priceCand * quantities
      elasticities <- t(elast(object, preMerger))
      thisFOC <- revenues * diag(owner) + as.vector((elasticities * owner) %*% (margins * revenues))
      thisFOC[!subset] <- revenues[!subset] # set quantity equal to 0 for firms not in subset
      return(thisFOC)
    }
    ## Try nleqslv first, fallback to BBsolve
    nleqslv_maxit <- as.integer(object@control.equ$maxit)
    if (length(nleqslv_maxit) == 0 || is.na(nleqslv_maxit[1]) || nleqslv_maxit[1] < 1) nleqslv_maxit <- 150L
    minResult <- nleqslv::nleqslv(priceStart, FOC,
      method = "Broyden",
      control = list(
        ftol = object@control.equ$tol,
        maxit = nleqslv_maxit
      )
    )

    if (minResult$termcd > 2) {
      minResult <- BBsolve(priceStart, FOC, quiet = TRUE, control = object@control.equ, ...)
      priceEst <- minResult$par
      if (minResult$convergence != 0) {
        warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBSolve' reports: '", minResult$message, "'")
      }
    } else {
      priceEst <- minResult$x
      if (minResult$termcd > 1) {
        warning("'calcPrices' may not have fully converged. 'nleqslv' termcd: ", minResult$termcd)
      }
    }
    names(priceEst) <- object@labels
    return(priceEst)
  }
)
#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "AIDS",
  definition = function(object, preMerger = TRUE, ...) {
    ## if(any(is.na(object@prices)){warning("'prices' contains missing values. AIDS can only predict price changes, not price levels")}
    if (preMerger) {
      prices <- object@prices
    } else {
      prices <- object@prices * (1 + object@priceDelta)
    }
    names(prices) <- object@labels
    return(prices)
  }
)
#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "BargainingLogit",
  definition = function(object, preMerger = TRUE, isMax = FALSE, subset, ...) {
    # Better starting values improve convergence for post-merger; guard against missing/NA
    priceStart <- object@priceStart

    alpha <- object@slopes$alpha

    if (preMerger) {
      owner <- object@ownerPre
      mc <- object@mcPre
      barg <- object@bargpowerPre
    } else {
      owner <- object@ownerPost
      mc <- object@mcPost
      barg <- object@bargpowerPost
    }

    barg <- barg / (1 - barg)

    nprods <- length(object@shares)
    if (missing(subset)) {
      subset <- rep(TRUE, nprods)
    }

    if (!is.logical(subset) || length(subset) != nprods) {
      stop("'subset' must be a logical vector the same length as 'shares'")
    }

    # If solving post-merger, try to initialize from pre-merger prices when valid
    if (!preMerger) {
      cand <- object@pricePre
      if (length(cand) == nprods && all(is.finite(cand))) {
        priceStart <- cand
      }
    }

    if (any(!subset)) {
      owner <- owner[subset, subset]
      mc <- mc[subset]
      priceStart <- priceStart[subset]
      barg <- barg[subset]
    }
    owner <- .retention_owner_bargaining_logit(owner,
        getRetention(object, preMerger)[subset])

    priceEst <- rep(NA, nprods)


    ## Define system of FOC as a function of prices
    FOC <- function(priceCand) {
      if (preMerger) {
        object@pricePre[subset] <- priceCand
      } else {
        object@pricePost[subset] <- priceCand
      }


      shares <- calcShares(object, preMerger = preMerger, revenue = FALSE)[subset]
      levelMargin <- if (object@output) {
        priceCand - mc
      } else {
        mc - priceCand
      }
      outSign <- ifelse(object@output, -1, 1)

      elastInv <- owner
      # diag(elastInv) <- -1*diag(elastInv)
      elastInv <- -elastInv * shares
      diag(elastInv) <- diag(owner) + diag(elastInv)

      tmp <- try(solve(t(elastInv)), silent = TRUE)
      if (any(class(tmp) == "try-error")) {
        elastInv <- MASS::ginv(t(elastInv))
      } else {
        elastInv <- tmp
      }

      thisFOC <- levelMargin - elastInv %*% ((log(1 - shares) * diag(owner)) / (-1 * outSign * alpha * (barg * shares / (1 - shares) -
        log(1 - shares))))

      return(as.vector(thisFOC))
    }

    ## Find price changes that set FOCs equal to 0
    ## Try nleqslv first, fallback to BBsolve
    nleqslv_maxit <- as.integer(object@control.equ$maxit)
    if (length(nleqslv_maxit) == 0 || is.na(nleqslv_maxit[1]) || nleqslv_maxit[1] < 1) nleqslv_maxit <- 150L
    minResult <- nleqslv::nleqslv(priceStart, FOC,
      method = "Broyden",
      control = list(
        ftol = object@control.equ$tol,
        maxit = nleqslv_maxit
      )
    )

    if (minResult$termcd > 2) {
      minResult <- BBsolve(priceStart, FOC, quiet = TRUE, control = object@control.equ, ...)
      priceEst_solution <- minResult$par
      if (minResult$convergence != 0) {
        warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBsolve' reports: '", minResult$message, "'")
      }
    } else {
      priceEst_solution <- minResult$x
      if (minResult$termcd > 1) {
        warning("'calcPrices' may not have fully converged. 'nleqslv' termcd: ", minResult$termcd)
      }
    }


    if (isMax) {
      hess <- genD(FOC, priceEst_solution) # compute the numerical approximation of the FOC hessian at optimium
      hess <- hess$D[, 1:hess$p]
      hess <- hess * (owner > 0) # 0 terms not under the control of a common owner

      state <- ifelse(preMerger, "Pre-merger", "Post-merger")

      if (any(eigen(hess)$values > 0)) {
        warning("Hessian of first-order conditions is not positive definite. ", state, " price vector may not maximize profits. Consider rerunning 'calcPrices' using different starting values")
      }
    }


    priceEst[subset] <- priceEst_solution
    names(priceEst) <- object@labels

    return(priceEst)
  }
)

#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "LogitBLP",
  definition = function(object, preMerger = TRUE, isMax = FALSE, ...) {
    # 1 Set prices and ownership
    output <- object@output
    priceStart <- object@priceStart
    nprods <- length(object@shares)

    if (preMerger) {
      subset <- rep(TRUE, nprods)
      owner <- object@ownerPre
      mc <- object@mcPre
    } else {
      subset <- object@subset
      owner <- object@ownerPost
      mc <- object@mcPost
    }


    priceStart <- priceStart[subset]

    # 2 Define FOCs function (Unified for Root-Finding and Fixed-Point)
    FOC <- function(priceCand, as_fp = FALSE) {
      thisobj <- object
      thisPrice <- rep(0, nprods)
      thisPrice[subset] <- priceCand

      if (preMerger) {
        thisobj@pricePre <- thisPrice
      } else {
        thisobj@pricePost <- thisPrice
      }

      # Predicted Margins (Omega^-1 * S)
      predMargin <- calcMargins(thisobj, preMerger, level = TRUE)[subset]

      if (as_fp) {
        # Fixed Point Update: P = mc +/- Omega^-1 * S
        if (output) {
          return(mc[subset] + predMargin)
        } else {
          return(mc[subset] - predMargin)
        }
      } else {
        # Residual: (Actual Margin) - (Predicted Margin)
        margins <- if (output) priceCand - mc[subset] else mc[subset] - priceCand[subset]
        return(margins - predMargin)
      }
    }


    # 4 Solve FOCs

    # Strategy 1: SQUAREM (Contraction Mapping)
    sq_control <- list(tol = object@control.equ$tol, maxiter = 1500)
    if (!is.null(object@control.equ$maxit) && !is.na(object@control.equ$maxit)) {
      sq_control$maxiter <- object@control.equ$maxit
    }

    minResult <- tryCatch(
      {
        SQUAREM::squarem(par = priceStart, fixptfn = FOC, as_fp = TRUE, control = sq_control)
      },
      error = function(e) NULL
    )

    success <- FALSE
    if (!is.null(minResult) && minResult$convergence) {
      priceEst_solution <- minResult$par
      success <- TRUE
    }


    # Strategy 2: BBsolve (Fallback)
    if (!success) {
      # Fallback to original FOC (margin based) with BBsolve
      minResult <- BB::BBsolve(
        par = priceStart, fn = FOC,
        control = list(tol = object@control.equ$tol, maxit = 1500),
        quiet = TRUE
      )
      if (minResult$convergence == 0) {
        priceEst_solution <- minResult$par
        success <- TRUE
      }
    }

    # Strategy 3: nleqslv (Newton with Analytic Jacobian)
    if (!success) {
      nleqslv_maxit <- as.integer(object@control.equ$maxit)
      if (nprods <= 30 && (length(nleqslv_maxit) == 0 || is.na(nleqslv_maxit[1]) || nleqslv_maxit[1] < 1)) {
        nleqslv_maxit <- 150L
      }

      minResult <- nleqslv::nleqslv(
        x = priceStart, fn = FOC,
        method = "Newton",
        control = list(ftol = object@control.equ$tol, maxit = nleqslv_maxit)
      )

      priceEst_solution <- minResult$x
      if (minResult$termcd == 1) success <- TRUE
    }

    if (!success) {
      warning("'calcPrices' solver hierarchy (SQUAREM -> BBsolve -> nleqslv) may not have fully converged.")
    }

    # 5 Optional Hessian check for maximization
    if (isMax) {
      hess <- numDeriv::jacobian(FOC, priceEst_solution)
      hess <- hess * (owner[subset, subset] > 0)
      if (any(eigen(hess)$values > 0)) {
        warning("Hessian not positive definite at solution.")
      }
    }

    # 6 Return final prices
    result <- rep(NA, nprods)
    result[subset] <- priceEst_solution


    names(result) <- object@labels
    return(result)
  }
)
