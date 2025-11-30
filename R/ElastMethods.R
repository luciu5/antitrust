#' @title Methods For Calculating Own and Cross-Price Elasticities
#' @name Elast-Methods
#' @docType methods
#'
#' @aliases elast-methods
#' elast
#' elast,ANY-method
#' elast,AIDS-method
#' elast,CES-method
#' elast,CESNests-method
#' elast,Linear-method
#' elast,LogLin-method
#' elast,Logit-method
#' elast,LogitNests-method
#' elast,Cournot-method
#' elast,VertBargBertLogit-method
#' elast,LogitBLP-method
#'
#' @description Calculate the own and cross-price elasticity between any two products in the market.
#' @param object An instance of one of the classes listed above.
#' @param preMerger If TRUE, calculates pre-merger price elasticities. If
#' FALSE, calculates post-merger price elasticities. Default is TRUE.
#' @param market If TRUE, calculates the market (aggregate) elasticity. If
#' FALSE, calculates matrix of own- and cross-price elasticities. Default is FALSE.
#'
#' @details When \sQuote{market} is FALSE, this method computes the matrix
#' of own and cross-price elasticities. Element i,j of this matrix is
#' the percentage change in the demand for good i from a small change in
#' the price of good j. When \sQuote{market} is TRUE, this method computes the
#' market (aggregate) elasticities using share-weighted prices.
#'
#' When \sQuote{preMerger} is TRUE, elasticities are
#' calculated at pre-merger equilibrium prices and shares, and when \sQuote{preMerger} is FALSE, they
#' are calculated at post-merger equilibrium prices and shares.
#'
#' @return returns a k x k matrix of own- and cross-price elasticities,
#' where k is the number of products in the market.
#'
#' @include DiversionMethods.R
#' @keywords methods
NULL

setGeneric(
  name = "elast",
  def = function(object, ...) {
    standardGeneric("elast")
  }
)

#' @rdname Elast-Methods
#' @export
setMethod(
  f = "elast",
  signature = "Cournot",
  definition = function(object, preMerger = TRUE, market = FALSE) {
    isLinearD <- object@demand == "linear"
    slopes <- object@slopes
    intercepts <- object@intercepts


    if (preMerger) {
      quantities <- object@quantityPre
      owner <- object@ownerPre
    } else {
      quantities <- object@quantityPost
      owner <- object@ownerPost
    }

    quantities[is.na(quantities)] <- 0

    quantOwner <- owner %*% quantities

    prices <- calcPrices(object, preMerger = preMerger)

    mktQuant <- colSums(quantities, na.rm = TRUE)

    ## dPdQ
    partial <- ifelse(isLinearD,
      slopes,
      exp(intercepts) * slopes * mktQuant^(slopes - 1)
    )

    ## dQdP
    partial <- 1 / partial


    elast <- partial * prices / mktQuant


    if (!market) {
      sharesOwner <- t(quantOwner) / mktQuant


      elast <- t(elast / sharesOwner)


      dimnames(elast) <- object@labels
    }

    return(elast)
  }
)


#' @rdname Elast-Methods
#' @export
setMethod(
  f = "elast",
  signature = "Linear",
  definition = function(object, preMerger = TRUE, market = FALSE) {
    if (preMerger) {
      prices <- object@pricePre
    } else {
      prices <- object@pricePost
    }

    slopes <- object@slopes


    quantities <- calcQuantities(object, preMerger)

    if (market) {
      elast <- sum(slopes) / sum(quantities) * sum(quantities * prices / sum(quantities))
    } else {
      elast <- slopes * tcrossprod(1 / quantities, prices)
      dimnames(elast) <- list(object@labels, object@labels)
    }

    return(elast)
  }
)

#' @rdname Elast-Methods
#' @export
setMethod(
  f = "elast",
  signature = "Logit",
  definition = function(object, preMerger = TRUE, market = FALSE) {
    if (preMerger) {
      prices <- object@pricePre
    } else {
      prices <- object@pricePost
    }


    labels <- object@labels

    alpha <- object@slopes$alpha

    shares <- calcShares(object, preMerger = preMerger)

    if (market) {
      elast <- alpha * sum(shares / sum(shares) * prices, na.rm = TRUE) * (1 - sum(shares, na.rm = TRUE))

      names(elast) <- NULL
    } else {
      nprods <- length(shares)

      elast <- -alpha * matrix(prices * shares, ncol = nprods, nrow = nprods, byrow = TRUE)
      diag(elast) <- alpha * prices + diag(elast)

      dimnames(elast) <- list(labels, labels)
    }

    return(elast)
  }
)

#' @rdname Elast-Methods
#' @export
setMethod(
  f = "elast",
  signature = "LogitBLP",
  definition = function(object, preMerger = TRUE, market = FALSE) {
    if (preMerger) {
      prices <- object@pricePre
    } else {
      prices <- object@pricePost
    }

    # Retrieve random coefficients on price created by calcSlopes(LogitBLP)
    alphas <- object@slopes$alphas
    nDraws <- length(alphas)

    # Outside-good nesting parameter: sigmaNest in (0,1]
    # sigmaNest = 1 is flat logit
    sigmaNest <- object@slopes$sigmaNest
    if (is.null(sigmaNest)) sigmaNest <- 1

    shares <- calcShares(object, preMerger = preMerger)

    if (market) {
      # Calculate market elasticity by aggregating the individual product elasticities
      # This ensures consistency with the random coefficients and nesting structure

      # Get the full elasticity matrix (recursive call with market=FALSE)
      elast_mat <- elast(object, preMerger = preMerger, market = FALSE)

      # Re-normalize to inside shares for weighting
      inside_shares <- shares / sum(shares, na.rm = TRUE)

      # Market Elasticity = sum( w_j * sum_k(epsilon_jk) )
      mkt_elast <- sum(inside_shares * rowSums(elast_mat), na.rm = TRUE)

      names(mkt_elast) <- NULL
      return(mkt_elast)
    } else {
      nprods <- length(shares)

      # For individual elasticities, we need to account for the distribution of random coefficients
      # Calculate the matrix of share derivatives for each consumer type and then average
      # Elasticity = (P/S) * dS/dP
      # dS/dP = mean( d s_r / d p )
      
      derivDraws <- array(0, dim = c(nprods, nprods, nDraws))
      sharesDraws <- matrix(0, nrow = nprods, ncol = nDraws)

      delta <- object@slopes$meanval

      # Vectorized calculation over draws
      for (r in 1:nDraws) {
        # Utility: V_j = delta_j + alpha_r * (p_j - p_0)
        util <- matrix(rep(delta, 1), ncol = 1) + outer(prices - object@priceOutside, alphas[r])

        # Inclusive Value Calculation
        # D_g = sum(exp(V_k / sigma))
        expUtil_sigma <- exp(util / sigmaNest)
        D_g <- sum(expUtil_sigma)

        # Shares
        # s_j|g = exp(V_j / sigma) / D_g
        s_jg <- expUtil_sigma / D_g

        # s_g = D_g^sigma / (1 + D_g^sigma)
        D_g_sigma <- D_g^sigmaNest
        s_g <- D_g_sigma / (1 + D_g_sigma)

        # s_j = s_j|g * s_g
        sInd <- s_jg * s_g
        sharesDraws[, r] <- sInd

        # Elasticities (Standard Nested Logit)
        # Own-price: alpha * p_j * [ 1/sigma + s_jg * (1 - 1/sigma - s_g) ]
        # Cross-price: alpha * p_k * s_k * [ 1/s_g - 1 - 1/(sigma * s_g) ]

        # Pre-calculate common terms
        inv_sigma <- 1 / sigmaNest
        term2 <- 1 - inv_sigma - s_g
        term3 <- term2 / s_g

        # Vectorized Own-Price
        diag_elast <- alphas[r] * prices * (inv_sigma + s_jg * term2)

        # Vectorized Cross-Price
        cross_elast_row <- alphas[r] * prices * sInd * term3
        elast_matrix <- matrix(cross_elast_row, nrow = nprods, ncol = nprods, byrow = TRUE)
        diag(elast_matrix) <- diag_elast

        # Convert to derivatives: d s_r / d p = E_r * s_r / p
        # Element (j,k): E_jkr * s_jr / p_k
        deriv_matrix <- elast_matrix * (sInd %*% t(1/prices))
        
        derivDraws[, , r] <- deriv_matrix
      }

      # Average derivatives and shares across consumer types
      avg_deriv <- apply(derivDraws, c(1, 2), mean)
      avg_shares <- rowMeans(sharesDraws)
      
      # Compute Aggregate Elasticity
      # E_jk = (p_k / S_j) * (dS_j / dP_k)
      elast <- avg_deriv * ( (1/avg_shares) %*% t(prices) )
      
      dimnames(elast) <- list(object@labels, object@labels)

      return(elast)
    }
  }
)

#' @rdname Elast-Methods
#' @export
setMethod(
  f = "elast",
  signature = "LogLin",
  definition = function(object, preMerger = TRUE, market = FALSE) {
    if (market) {
      quantities <- calcQuantities(object, preMerger)
      prices <- calcPrices(object, preMerger)
      elast <- sum(t(t(object@slopes * quantities) * 1 / prices)) / sum(quantities) * sum(quantities * prices / sum(quantities))
    } else {
      elast <- object@slopes
      dimnames(elast) <- list(object@labels, object@labels)
    }

    return(elast)
  }
)

#' @rdname Elast-Methods
#' @export
setMethod(
  f = "elast",
  signature = "AIDS",
  definition = function(object, preMerger = TRUE, market = FALSE) {
    if (market) {
      return(object@mktElast)
    } else {
      shares <- calcShares(object, preMerger)

      elast <- t(object@slopes / shares) + shares * (object@mktElast + 1) # Caution: returns TRANSPOSED elasticity matrix
      diag(elast) <- diag(elast) - 1
      dimnames(elast) <- list(object@labels, object@labels)

      return(t(elast))
    }
  }
)

#' @rdname Elast-Methods
#' @export
setMethod(
  f = "elast",
  signature = "LogitNests",
  definition = function(object, preMerger = TRUE, market = FALSE) {
    if (preMerger) {
      prices <- object@pricePre
    } else {
      prices <- object@pricePost
    }


    nests <- object@nests
    alpha <- object@slopes$alpha
    sigma <- object@slopes$sigma
    meanval <- object@slopes$meanval

    shares <- calcShares(object, preMerger, revenue = FALSE)

    if (market) {
      elast <- alpha * sum(shares * prices) * (1 - sum(shares))
      names(elast) <- NULL
    } else {
      sharesNests <- shares / tapply(shares, nests, sum, na.rm = TRUE)[nests]


      nprods <- length(shares)

      elast <- diag((1 / sigma - 1) * alpha)
      elast <- elast[nests, nests]
      elast <- elast * matrix(sharesNests * prices, ncol = nprods, nrow = nprods, byrow = TRUE)
      elast <- -1 * (elast + alpha * matrix(shares * prices, ncol = nprods, nrow = nprods, byrow = TRUE))
      diag(elast) <- diag(elast) + (1 / sigma[nests]) * alpha * prices

      dimnames(elast) <- list(object@labels, object@labels)
    }
    return(elast)
  }
)

#' @rdname Elast-Methods
#' @export
setMethod(
  f = "elast",
  signature = "CES",
  definition = function(object, preMerger = TRUE, market = FALSE) {
    gamma <- object@slopes$gamma

    shares_r <- calcShares(object, preMerger, revenue = TRUE)
    shares_q <- calcShares(object, preMerger, revenue = FALSE)

    if (market) {
      if (preMerger) {
        prices <- object@pricePre
      } else {
        prices <- object@pricePost
      }

      # avgPrice <- sum(prices*shares_q)/sum(shares_q)

      elast <- (1 - gamma) * (1 - sum(shares_r)) - 1

      # elast <- -1 * elast
      # elast <- -gamma * ( 1 - sum(shares_r))*avgPrice

      names(elast) <- NULL
    } else {
      nprods <- length(shares_r)
      elast <- (gamma - 1) * matrix(shares_r, ncol = nprods, nrow = nprods, byrow = TRUE)
      diag(elast) <- -gamma + diag(elast)

      dimnames(elast) <- list(object@labels, object@labels)
    }
    return(elast)
  }
)

#' @rdname Elast-Methods
#' @export
setMethod(
  f = "elast",
  signature = "CESNests",
  definition = function(object, preMerger = TRUE, market = FALSE) {
    nests <- object@nests
    gamma <- object@slopes$gamma
    sigma <- object@slopes$sigma
    meanval <- object@slopes$meanval

    shares <- calcShares(object, preMerger, revenue = TRUE)
    sharesNests <- shares / tapply(shares, nests, sum, na.rm = TRUE)[nests]

    if (market) {
      alpha <- object@slopes$alpha
      if (is.null(alpha)) {
        stop("'shareInside' must be between 0 and 1 to  calculate Market Elasticity")
      }
      elast <- (1 + alpha) * (1 - gamma) * sum(shares) * (1 - sum(shares))
      names(elast) <- NULL
    } else {
      nprods <- length(shares)

      elast <- diag(sigma - gamma)
      elast <- elast[nests, nests]
      elast <- elast * matrix(sharesNests, ncol = nprods, nrow = nprods, byrow = TRUE)
      elast <- elast + (gamma - 1) * matrix(shares, ncol = nprods, nrow = nprods, byrow = TRUE)
      diag(elast) <- diag(elast) - sigma[nests]

      dimnames(elast) <- list(object@labels, object@labels)
    }
    return(elast)
  }
)


#' @rdname Elast-Methods
#' @export
setMethod(
  f = "elast",
  signature = "VertBargBertLogit",
  definition = function(object, preMerger = TRUE, market = FALSE) {
    result <- elast(object@down, preMerger = preMerger, market = market)

    return(result)
  }
)
