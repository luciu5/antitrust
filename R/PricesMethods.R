#' @title \dQuote{Calculating Prices} Methods
#' @name Prices-Methods
#' @docType methods
#' @aliases calcPrices,ANY-method
#' calcPrices,Auction2ndCap-method
#' calcPrices,Cournot-method
#' calcPrices,Linear-method
#' calcPrices,Logit-method
#' calcPrices,LogitBLP-method
#' calcPrices,LogLin-method
#' calcPrices,AIDS-method
#' calcPrices,LogitCap-method
#' calcPrices,Auction2ndLogit-method
#' calcPrices,VertBargBertLogit-method
#' calcPrices,VertBarg2ndLogit-method
#' calcPrices
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
#' For others, additional values that may be used to change the default values of \code{\link{constrOptim}}, the non-linear
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
#' @include AuctionCapMethods.R VerticalClasses.R
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
    # Start from pre-merger prices when available; fall back to priceStart
    priceStart <- object@priceStart
    outSign <- ifelse(output, -1, 1)
    alpha <- object@slopes$alpha

    if (preMerger) {
      owner <- object@ownerPre
      mc <- object@mcPre
    } else {
      owner <- object@ownerPost
      mc <- object@mcPost
    }
    nprods <- length(object@shares)
    if (!preMerger) {
      cand <- object@pricePre
      if (length(cand) == nprods && all(is.finite(cand))) priceStart <- cand
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
    }
    priceEst <- rep(NA, nprods)
    ## Define system of FOC as a function of prices
    FOC <- function(priceCand) {
      if (preMerger) {
        object@pricePre[subset] <- priceCand
      } else {
        object@pricePost[subset] <- priceCand
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
#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "LogitBLP",
  definition = function(object, preMerger = TRUE, isMax = FALSE, subset, ...) {
    # 1 Set prices and ownership
    output <- object@output
    priceStart <- object@priceStart

    if (preMerger) {
      owner <- object@ownerPre
      mc <- object@mcPre
    } else {
      owner <- object@ownerPost
      mc <- object@mcPost
    }

    nprods <- length(object@shares)
    if (!preMerger) {
      cand <- object@pricePre
      if (length(cand) == nprods && all(is.finite(cand))) priceStart <- cand
    }
    if (missing(subset)) subset <- rep(TRUE, nprods)
    if (!is.logical(subset) || length(subset) != nprods) {
      stop("'subset' must be a logical vector the same length as 'shares'")
    }

    if (any(!subset)) {
      owner <- owner[subset, subset]
      mc <- mc[subset]
      priceStart <- priceStart[subset]
    }

    # 2 Define FOCs function (Unified for Root-Finding and Fixed-Point)
    FOC <- function(priceCand, as_fp = FALSE) {
      if (preMerger) {
        object@pricePre[subset] <- priceCand
      } else {
        object@pricePost[subset] <- priceCand
      }

      # Predicted Margins (Omega^-1 * S)
      predMargin <- calcMargins(object, preMerger, level = TRUE)[subset]

      if (as_fp) {
        # Fixed Point Update: P = mc +/- Omega^-1 * S
        if (output) {
          return(mc + predMargin)
        } else {
          return(mc - predMargin)
        }
      } else {
        # Residual: (Actual Margin) - (Predicted Margin)
        margins <- if (output) priceCand - mc else mc - priceCand
        return(margins - predMargin)
      }
    }

    # 3 Define analytic Jacobian using individual shares
    # Helper to compute Share Jacobian J_s
    calc_Js <- function(priceCand) {
      tmpObject <- object
      if (preMerger) {
        tmpObject@pricePre[subset] <- priceCand
      } else {
        tmpObject@pricePost[subset] <- priceCand
      }

      shares_draw <- calcShares(tmpObject, preMerger = preMerger, aggregate = FALSE)[subset, , drop = FALSE]
      shares_draw <- pmax(pmin(shares_draw, 1 - 1e-12), 1e-12)

      alphas <- tmpObject@slopes$alphas
      sigmaNest <- tmpObject@slopes$sigmaNest
      if (is.null(sigmaNest)) sigmaNest <- 1
      R <- length(alphas)

      S_g <- colSums(shares_draw)
      S_jg <- sweep(shares_draw, 2, S_g, "/")

      A <- sweep(shares_draw, 2, alphas, "*")
      invSigma <- 1 / sigmaNest
      term_nest <- 1 - invSigma

      part1 <- tcrossprod(A, S_jg) * (term_nest / R)
      part2 <- tcrossprod(A, shares_draw) * (-1 / R)
      diag_boost <- rowMeans(A) * invSigma

      J_s <- part1 + part2
      diag(J_s) <- diag(J_s) + diag_boost
      return(J_s)
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

    # Strategy 2: nleqslv (Newton with Analytic Jacobian)
    if (!success) {
      # Define Share-based FOC and Jacobian for Newton
      FOC_Newton <- function(priceCand) {
        if (preMerger) {
          object@pricePre[subset] <- priceCand
        } else {
          object@pricePost[subset] <- priceCand
        }

        S <- calcShares(object, preMerger, aggregate = TRUE)[subset]
        J_s <- calc_Js(priceCand)

        margins <- if (output) priceCand - mc else mc - priceCand

        # FOC: S + (t(J_s) * owner) %*% margins = 0
        return(S + as.vector((t(J_s) * owner) %*% margins))
      }

      JAC_Newton <- function(priceCand) {
        J_s <- calc_Js(priceCand)
        # Approx Jacobian: J_s + (t(J_s) * owner)
        return(J_s + (t(J_s) * owner))
      }

      nleqslv_maxit <- as.integer(object@control.equ$maxit)
      if (nprods <= 30 && (length(nleqslv_maxit) == 0 || is.na(nleqslv_maxit[1]) || nleqslv_maxit[1] < 1)) {
        nleqslv_maxit <- 150L
      }

      minResult <- nleqslv::nleqslv(
        x = priceStart, fn = FOC_Newton, jac = JAC_Newton,
        method = "Newton",
        control = list(ftol = object@control.equ$tol, maxit = nleqslv_maxit)
      )

      priceEst_solution <- minResult$x
      if (minResult$termcd == 1) success <- TRUE
    }

    # Strategy 3: BBsolve (Fallback)
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

    if (!success) {
      warning("'calcPrices' solver hierarchy (SQUAREM -> BBsolve -> nleqslv) may not have fully converged.")
    }

    # 5 Optional Hessian check for maximization
    if (isMax) {
      hess <- numDeriv::jacobian(FOC, priceEst_solution)
      hess <- hess * (owner > 0)
      if (any(eigen(hess)$values > 0)) {
        warning("Hessian not positive definite at solution.")
      }
    }

    # 6 Return final prices
    priceEst <- rep(NA, nprods)
    priceEst[subset] <- priceEst_solution
    names(priceEst) <- object@labels
    return(priceEst)
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
    priceEst <- rep(NA, nprods)
    ## Define system of FOC as a function of prices
    FOC <- function(priceCand) {
      if (preMerger) {
        object@pricePre[subset] <- priceCand
      } else {
        object@pricePost[subset] <- priceCand
      }

      if (output) {
        margins <- 1 - mc / priceCand
      } else {
        margins <- mc / priceCand - 1
      }

      revenues <- calcShares(object, preMerger = preMerger, revenue = TRUE)
      quantities <- calcQuantities(object, preMerger = preMerger)
      revenues <- revenues[subset]
      elasticities <- elast(object, preMerger)[subset, subset]
      thisFOC <- revenues * diag(owner) + as.vector(t(elasticities * owner) %*% (margins * revenues))
      constraint <- ifelse(is.finite(capacities), (quantities - capacities) / object@insideSize, 0)
      measure <- ifelse(constraint != 0, thisFOC + constraint + sqrt(thisFOC^2 + constraint^2), thisFOC)
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
    ## first try the analytic solution
    #      prices <-
    #          solve((slopes*diag(owner)) + (t(slopes)*owner)) %*%
    #              ((t(slopes)*owner) %*% mc - (intercept*diag(owner)))
    #
    #       prices <- as.vector(prices)
    #       quantities <- as.vector(intercept + t(slopes) %*% prices)
    #
    #      ##use the numeric solution if analytic solution yields negative quantities
    #      if(any(subset) || any(quantities<0)){
    #
    #          if(any(quantities<0)){
    #           warning("Equilibrium prices yield negative equilibrium quantities. Recomputing equilibrium prices  under the restriction that equilbrium quantities must be non-negative")
    #                               }
    #          else if(any(subset)){
    #            warning("Elements of 'subset' are  FALSE. Computing equilbrium under the restriction that these products have 0 sales")
    #          }
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
    ## Find starting value that always meets boundary conditions
    ## startParm <- as.vector(solve(slopes) %*% (-intercept + 1))
    minResult <- constrOptim(object@priceStart, FOC, grad = NULL, ui = slopes, ci = -intercept, ...)
    if (!isTRUE(all.equal(minResult$convergence, 0, check.names = FALSE))) {
      warning("'calcPrices' solver may not have successfully converged.'constrOptim' reports: '", minResult$message, "'")
    }
    prices <- minResult$par
    # }
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
    if (!preMerger) {
      cand <- object@pricePre
      if (length(cand) == nprods && all(is.finite(cand))) priceStart <- cand
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
    }

    priceEst <- rep(NA, nprods)


    ## Define system of FOC as a function of prices
    FOC <- function(priceCand) {
      if (preMerger) {
        object@pricePre[subset] <- priceCand
      } else {
        object@pricePost[subset] <- priceCand
      }


      shares <- calcShares(object, preMerger = preMerger, revenue = FALSE)[subset]

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

      thisFOC <- (priceCand - mc) - elastInv %*% ((log(1 - shares) * diag(owner)) / (alpha * (barg * shares / (1 - shares) -
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
  signature = "VertBargBertLogit",
  definition = function(object, preMerger = TRUE, ...) {
    up <- object@up
    down <- object@down

    chain_level <- object@chain_level

    upPricePre <- up@prices
    downPricePre <- down@prices

    # isHorizontal <- object@isHorizontal


    subsetDown <- down@subset
    alpha <- down@slopes$alpha
    meanval <- down@slopes$meanval[subsetDown]

    priceStartUp <- up@priceStart
    priceStartDown <- down@priceStart

    if (chain_level == "full") {
      priceStart <- c(priceStartUp[subsetDown], priceStartDown[subsetDown])
    } else if (chain_level == "retailer") {
      priceStart <- priceStartDown[subsetDown]
    } else if (chain_level == "wholesaler") {
      priceStart <- priceStartUp[subsetDown]
    }


    nprods <- length(meanval)

    ownerUp <- ownerToMatrix(up, preMerger = preMerger)[subsetDown, subsetDown]
    ownerDown <- ownerToMatrix(down, preMerger = preMerger)[subsetDown, subsetDown]

    if (preMerger) {
      # ownerUp <- up@ownerPre[subsetDown,subsetDown]
      # ownerDown <- down@ownerPre[subsetDown,subsetDown]
      up@ownerPre <- ownerUp
      down@ownerPre <- ownerDown
      bargparm <- up@bargpowerPre[subsetDown]
      mcUp <- up@mcPre[subsetDown]
      mcDown <- down@mcPre[subsetDown]
      # ownerVDown <- object@ownerDownPre[subsetDown,subsetDown]
      # ownerDownLambda <- object@ownerDownLambdaPre[subsetDown,subsetDown]
      # ownerUpLambda <- object@ownerUpLambdaPre[subsetDown,subsetDown]
    } else {
      # ownerUp <- up@ownerPost[subsetDown,subsetDown]
      # ownerDown <- down@ownerPost[subsetDown,subsetDown]
      up@ownerPost <- ownerUp
      down@ownerPost <- ownerDown
      bargparm <- up@bargpowerPost[subsetDown]
      mcUp <- up@mcPost[subsetDown]
      mcDown <- down@mcPost[subsetDown]
      # ownerVDown <- object@ownerDownPost[subsetDown,subsetDown]
      # ownerDownLambda <- object@ownerDownLambdaPost[subsetDown,subsetDown]
      # ownerUpLambda <- object@ownerUpLambdaPost[subsetDown,subsetDown]
    }


    FOC <- function(priceCand) {
      thisobj <- object


      if (chain_level == "full") {
        priceCandUp <- rep(NA, 1, nprods)[subsetDown]
        priceCandUp <- priceCand[1:length(priceCandUp)]
        priceCandDown <- priceCand[-(1:length(priceCandUp))]
      } else if (chain_level == "retailer") {
        priceCandUp <- up@prices[subsetDown]
        priceCandDown <- priceCand
      } else if (chain_level == "wholesaler") {
        priceCandDown <- down@prices[subsetDown]
        priceCandUp <- priceCand
      }


      if (preMerger) {
        up@pricePre <- priceCandUp
        down@pricePre <- priceCandDown
      } else {
        up@pricePost <- priceCandUp
        down@pricePost <- priceCandDown
      }

      thisobj@up <- up
      thisobj@down <- down
      # shareCandDown  <- calcShares(down, preMerger=preMerger,revenue=FALSE)
      #
      #
      # elast <-  -alpha*tcrossprod(shareCandDown)
      # diag(elast) <- alpha*shareCandDown + diag(elast)
      # elast.inv <- try(solve(ownerDown * elast),silent=TRUE)
      # if(any(class(elast.inv) == "try-error")){elast.inv <- MASS::ginv(ownerDown * elast)}
      # marginsDownCand <- -as.vector(elast.inv %*% shareCandDown)
      #
      #
      # div <- tcrossprod(1/(1-shareCandDown),shareCandDown)*shareCandDown
      # diag(div) <- -shareCandDown
      # div <- as.vector(div)
      #
      # #if(!isHorizontal){
      #   marginsDownCand <-  marginsDownCand - elast.inv %*% ( (ownerVDown * elast) %*% (priceCandUp-mcUp) )
      #
      # #}
      # marginsUpCand <- as.vector(solve(ownerUpLambda * div) %*% (ownerDownLambda * div) %*% marginsDownCand)
      #

      theseMargins <- calcMargins(thisobj, preMerger = preMerger, level = TRUE)

      upFOC <- priceCandUp - mcUp - theseMargins$up
      downFOC <- priceCandDown - priceCandUp - mcDown - theseMargins$down


      thisFOC <- c(upFOC, downFOC)

      if (chain_level == "full") {
        thisFOC <- c(upFOC, downFOC)
      } else if (chain_level == "retailer") {
        thisFOC <- downFOC
      } else if (chain_level == "wholesaler") {
        thisFOC <- upFOC
      }

      return(thisFOC)
    }


    ## Try nleqslv first, fallback to BBsolve
    nleqslv_maxit <- as.integer(down@control.equ$maxit)
    if (length(nleqslv_maxit) == 0 || is.na(nleqslv_maxit[1]) || nleqslv_maxit[1] < 1) nleqslv_maxit <- 150L
    minResult <- nleqslv::nleqslv(priceStart, FOC,
      method = "Broyden",
      control = list(
        ftol = down@control.equ$tol,
        maxit = nleqslv_maxit
      )
    )

    if (minResult$termcd > 2) {
      minResult <- BBsolve(priceStart, FOC, quiet = TRUE, control = down@control.equ, ...)
      priceResult <- minResult$par
      if (minResult$convergence != 0) {
        warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBsolve' reports: '", minResult$message, "'")
      }
    } else {
      priceResult <- minResult$x
      if (minResult$termcd > 1) {
        warning("'calcPrices' may not have fully converged. 'nleqslv' termcd: ", minResult$termcd)
      }
    }

    minResultUp <- rep(NA, length(priceStartUp))
    minResultDown <- rep(NA, length(priceStartDown))


    if (chain_level == "full") {
      minResultUp[subsetDown] <- priceResult[1:length(priceStartUp[subsetDown])]
      minResultDown[subsetDown] <- priceResult[-(1:length(priceStartUp[subsetDown]))]
    } else if (chain_level == "retailer") {
      minResultUp[subsetDown] <- upPricePre[subsetDown]
      minResultDown[subsetDown] <- priceResult
    } else if (chain_level == "wholesaler") {
      minResultUp[subsetDown] <- priceResult
      minResultDown[subsetDown] <- downPricePre[subsetDown]
    }


    return(list(up = minResultUp, down = minResultDown))
  }
)
#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "VertBarg2ndLogit",
  definition = function(object, preMerger = TRUE, ...) {
    up <- object@up
    down <- object@down

    alpha <- down@slopes$alpha
    meanval <- down@slopes$meanval


    subsetDown <- down@subset


    priceStartUp <- up@priceStart
    priceStartDown <- down@priceStart
    priceStart <- c(priceStartUp[subsetDown], priceStartDown[subsetDown])


    nprods <- length(priceStartDown)


    ownerDown <- ownerToMatrix(down, preMerger = preMerger)[subsetDown, subsetDown]

    priceOutside <- down@priceOutside

    if (preMerger) {
      bargparm <- up@bargpowerPre[subsetDown]
      down@ownerPre <- ownerDown
      mcUp <- up@mcPre[subsetDown]
      mcDown <- down@mcPre[subsetDown]
      ownerDownLambda <- object@ownerDownLambdaPre[subsetDown, subsetDown]
      ownerUpLambda <- object@ownerUpLambdaPre[subsetDown, subsetDown]
    } else {
      down@ownerPost <- ownerDown
      bargparm <- up@bargpowerPost[subsetDown]
      mcUp <- up@mcPost[subsetDown]
      mcDown <- down@mcPost[subsetDown]
      ownerDownLambda <- object@ownerDownLambdaPost[subsetDown, subsetDown]
      ownerUpLambda <- object@ownerUpLambdaPost[subsetDown, subsetDown]
    }


    FOC <- function(priceCand) {
      priceCandUp <- rep(NA, 1, nprods)[subsetDown]
      priceCandUp <- priceCand[1:length(priceCandUp)]
      priceCandDown <- priceCand[-(1:length(priceCandUp))]


      down@slopes$meanval <- meanval + alpha * (priceCandUp - priceOutside)

      shareCandDown <- calcShares(down, preMerger = preMerger, revenue = FALSE)
      marginsDownCand <- calcMargins(down, preMerger = preMerger, level = TRUE)


      div <- tcrossprod(1 / (1 - shareCandDown), shareCandDown) * shareCandDown
      diag(div) <- -shareCandDown
      div <- as.vector(div)

      # marginsUpCand <- as.vector(solve(ownerUpLambda * div) %*% ((ownerDownLambda * div) %*% marginsDownCand))


      upFOC <- (ownerUpLambda * div) %*% (priceCandUp - mcUp) - ((ownerDownLambda * div) %*% marginsDownCand)

      downFOC <- priceCandDown - priceCandUp - mcDown - marginsDownCand


      thisFOC <- c(
        downFOC,
        as.vector(upFOC)
      )

      return(thisFOC)
    }


    ## Try nleqslv first, fallback to BBsolve
    nleqslv_maxit <- as.integer(down@control.equ$maxit)
    if (length(nleqslv_maxit) == 0 || is.na(nleqslv_maxit[1]) || nleqslv_maxit[1] < 1) nleqslv_maxit <- 150L
    minResult <- nleqslv::nleqslv(priceStart, FOC,
      method = "Broyden",
      control = list(
        ftol = down@control.equ$tol,
        maxit = nleqslv_maxit
      )
    )

    if (minResult$termcd > 2) {
      minResult <- BB::BBsolve(priceStart, FOC, quiet = TRUE, control = down@control.equ, ...)
      priceResult <- minResult$par
      if (minResult$convergence != 0) {
        warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBsolve' reports: '", minResult$message, "'")
      }
    } else {
      priceResult <- minResult$x
      if (minResult$termcd > 1) {
        warning("'calcPrices' may not have fully converged. 'nleqslv' termcd: ", minResult$termcd)
      }
    }

    minResultUp <- rep(NA, length(priceStartUp))
    minResultDown <- rep(NA, length(priceStartDown))


    minResultUp[subsetDown] <- priceResult[1:length(priceStartUp[subsetDown])]
    minResultDown[subsetDown] <- priceResult[-(1:length(priceStartUp[subsetDown]))]

    down@slopes$meanval <- meanval + alpha * (minResultUp - priceOutside)

    # marginsDown <- calcMargins(down, preMerger=preMerger,level=TRUE)

    # minResultDown <- marginsDown  + minResultUp + mcDown

    return(list(up = minResultUp, down = minResultDown))
  }
)
