#' @title Methods For Calculating Own and Cross-Price Elasticities
#' @name Elast-Methods
#' @docType methods
#'
#' @aliases elast-methods elast elast,ANY-method elast,AIDS-method elast,CES-method elast,CESNests-method elast,Linear-method elast,LogLin-method elast,Logit-method elast,LogitNests-method elast,Cournot-method elast,LogitBLP-method elast,Auction2ndBLP-method elast,BargainingBLP-method
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

    # Get aggregate shares
    shares <- calcShares(object, preMerger = preMerger)
    shares[is.na(shares)] <- 0

    # Initialize partial derivatives matrix
    draw_weights <- .blp_draw_weights(object, nDraws)

    if (sigmaNest == 1) {
      ## For the flat model, aggregate the draw-level Logit Jacobians
      ## directly.  Keep the active mask separate from the zero-filled
      ## kernel matrix so excluded products and missing draw shares retain
      ## the loop's zero contribution.
      active_draw <- !is.na(shares_draw) &
        matrix(subset, nrow = nprods, ncol = nDraws)
      S <- shares_draw
      S[!active_draw] <- 0
      a <- draw_weights * alphas
      partial_deriv <- diag(drop(S %*% a)) -
        sweep(S, 2L, a, "*") %*% t(S)
    } else {
      ## Nested (one-nest) draw Jacobians, aggregated over draws by matrix
      ## products. Per draw r, with inside share s_g and within-nest shares
      ## s_j|g: d s_j / d p_k = alpha_r * s_j * s_k * term3 off the diagonal,
      ## where term3 = 1/s_g - 1 - (1/sigma)/s_g, and
      ## alpha_r * s_j * (1/sigma + s_j|g * term2) on it, term2 = 1 - 1/sigma - s_g.
      ## Since s_j s_k term3 = s_j|g s_k|g s_g term2, the cross term needs no
      ## division by s_g, so a draw whose inside shares underflow to zero
      ## contributes zero rather than NaN.
      S <- shares_draw
      S[is.na(S)] <- 0
      S[!subset, ] <- 0
      s_g <- colSums(S)
      S_jg <- sweep(S, 2L, s_g, "/")
      S_jg[!is.finite(S_jg)] <- 0
      inv_sigma <- 1 / sigmaNest
      term2 <- 1 - inv_sigma - s_g
      a <- draw_weights * alphas
      partial_deriv <- tcrossprod(sweep(S_jg, 2L, a * s_g * term2, "*"), S_jg)
      diag(partial_deriv) <- drop(
        (S * (inv_sigma + S_jg * matrix(term2, nrow = nprods, ncol = nDraws, byrow = TRUE))) %*% a
      )
      partial_deriv <- partial_deriv * outer(subset, subset)
    }

    # If market elasticity is requested
    if (market) {
      # Calculate elasticities from partials
      shares_for_elast <- pmax(shares, 1e-10) # Avoid division by zero
      elast_mat <- partial_deriv * (outer(1 / shares_for_elast, prices))

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
    shares_for_elast <- pmax(shares, 1e-10) # Avoid division by zero
    elast <- partial_deriv * (outer(1 / shares_for_elast, prices))

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
    output <- object@output

    shares_r <- calcShares(object, preMerger, revenue = TRUE)

    if (market) {
      outSign <- ifelse(output, 1, -1)
      elast <- outSign * ((1 - gamma) * (1 - sum(shares_r)) - 1)
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
    outSign <- ifelse(object@output, 1, -1)

    shares <- calcShares(object, preMerger, revenue = TRUE)
    sharesNests <- shares / tapply(shares, nests, sum, na.rm = TRUE)[nests]

    if (market) {
      alpha <- object@slopes$alpha
      if (is.null(alpha)) {
        stop("'shareInside' must be between 0 and 1 to  calculate Market Elasticity")
      }
      elast <- outSign * ((1 + alpha) * (1 - gamma) * sum(shares) * (1 - sum(shares)) - 1)
      names(elast) <- NULL
    } else {
      nprods <- length(shares)

      elast <- diag(sigma - gamma)
      elast <- elast[nests, nests]
      elast <- elast * matrix(sharesNests, ncol = nprods, nrow = nprods, byrow = TRUE)
      elast <- elast + (gamma - 1) * matrix(shares, ncol = nprods, nrow = nprods, byrow = TRUE)
      diag(elast) <- diag(elast) - sigma[nests]
      elast <- outSign * elast

      dimnames(elast) <- list(object@labels, object@labels)
    }
    return(elast)
  }
)
