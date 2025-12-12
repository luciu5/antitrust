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
#' @param partial if TRUE returns the partial derivative rather than the elasticity. Default is FALSE
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
  definition = function(object, preMerger = TRUE, market = FALSE, partial = FALSE) {
   
    outSign <- ifelse(object@output, 1, -1)  # 1 for output, -1 for input 
     # Get subsetting info and prices
    if (preMerger) {
      subset <- rep(TRUE, length(object@labels))
      prices <- object@pricePre
    } else {
      subset <- object@subset
      prices <- object@pricePost
    }
    
    # Get parameters
    alphas <- object@slopes$alphas
    sigmaNest <- object@slopes$sigmaNest
    if (is.null(sigmaNest)) sigmaNest <- 1
    
    # Get individual shares by draw (full size with NAs for excluded products)
    shares_draw <- calcShares(object, preMerger = preMerger, aggregate = FALSE)
    nDraws <- ncol(shares_draw)
    nprods <- nrow(shares_draw)
    
    # Replace NAs with zeros for calculation purposes
    shares_draw[is.na(shares_draw)] <- 0
    
    # Get aggregate shares
    shares <- calcShares(object, preMerger = preMerger)
    shares[is.na(shares)] <- 0
    
    # Group shares for each draw (only include active products)
    s_g <- colSums(shares_draw * matrix(subset, nrow=nprods, ncol=nDraws))
    
    # Conditional shares within group for each draw
    s_jg <- sweep(shares_draw, 2, s_g, "/")
    s_jg[is.na(s_jg) | is.infinite(s_jg)] <- 0  # Handle division by zero
    
    # Common terms for derivatives
    inv_sigma <- 1/sigmaNest
    term2 <- 1 - inv_sigma - s_g  # Length nDraws
    term3 <- 1/s_g - 1 - inv_sigma/s_g  # Length nDraws
    
    # Replace infinities with zeros
    term3[is.infinite(term3)] <- 0
    
    # Initialize partial derivatives matrix
    partial_deriv <- matrix(0, nrow = nprods, ncol = nprods)
    
    # Calculate average partial derivatives across draws
    for (r in 1:nDraws) {
      # Only calculate for active products
      active_r <- subset & (!is.na(shares_draw[,r]))
      
      if (sum(active_r) > 0) {  # Skip if no active products in this draw
        # Outer product of shares for cross-derivatives
        share_outer <- outer(shares_draw[,r], shares_draw[,r])
        
        # Cross-derivatives for this draw
        cross_deriv <- -1 * outSign * alphas[r] * share_outer * term3[r]
        
        # Own-derivatives for this draw
        own_deriv <-  -1 * outSign * shares_draw[,r] * (inv_sigma + s_jg[,r] * term2[r])
        
        # Combine into full matrix for this draw
        deriv_r <- cross_deriv
        diag(deriv_r) <- own_deriv
        
        # Mask out inactive products
        deriv_r <- deriv_r * outer(active_r, active_r)
        
        # Add to average
        partial_deriv <- partial_deriv + deriv_r / nDraws
      }
    }
    
    # If market elasticity is requested
    if (market) {
      # Calculate elasticities from partials
      shares_for_elast <- pmax(shares, 1e-10)  # Avoid division by zero
      elast_mat <- partial_deriv * (outer(1/shares_for_elast, prices))
      
      # Mask out inactive products
      elast_mat <- elast_mat * outer(subset, subset)
      
      # Calculate market elasticity as weighted average of row sums
      active_shares <- shares[subset]
      inside_shares <- active_shares / sum(active_shares)
      active_elast <- elast_mat[subset, subset]
      mkt_elast <- sum(inside_shares * rowSums(active_elast))
      
      names(mkt_elast) <- NULL
      return(mkt_elast)
    }
    
    # Return partials if requested
    if (partial) {
      # Zero out entries for inactive products
      partial_deriv <- partial_deriv * outer(subset, subset)
      dimnames(partial_deriv) <- list(object@labels, object@labels)
      return(partial_deriv)
    }
    
    # Calculate elasticities
    shares_for_elast <- pmax(shares, 1e-10)  # Avoid division by zero
    elast <- partial_deriv * (outer(1/shares_for_elast, prices))
    
    # Zero out entries for inactive products
    elast <- elast * outer(subset, subset)
    
    dimnames(elast) <- list(object@labels, object@labels)
    return(elast)
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
