## BLP mean-utility recovery.
##
## This generic + LogitBLP method implement the BLP contraction mapping that
## recovers mean utilities (delta) rationalizing observed shares given
## random-coefficient parameters. It is exported so that downstream packages
## building coordinated-effects models on top of LogitBLP (e.g. price
## leadership with BLP demand) can reuse the exact same contraction rather
## than reaching into antitrust internals.

#' @title Calculate BLP Mean Utilities via Contraction Mapping
#' @name calcMeanval
#' @rdname calcMeanval
#' @description \code{calcMeanval} runs the BLP contraction mapping to recover
#'   mean utilities (delta/meanval) that rationalize observed shares given
#'   random-coefficient parameters. It is defined for the \code{LogitBLP} class
#'   and inherited by its subclasses.
#' @param object A \code{LogitBLP} object (or a subclass) with BLP parameters
#'   stored in its \code{slopes} slot.
#' @param ... Additional arguments (currently unused).
#' @return The object with updated \code{slopes$meanval} (the delta vector) and
#'   related BLP integration bookkeeping.
#' @export
setGeneric(
  name = "calcMeanval",
  def = function(object, ...) {standardGeneric("calcMeanval")}
)

#' @rdname calcMeanval
#' @details Implements the BLP contraction mapping
#' delta^(t+1) = delta^(t) + log(s_observed) - log(s_predicted(delta^(t))),
#' using \code{BB::dfsane} for robust convergence with a fallback to a dampened
#' fixed point. Handles nesting (sigmaNest) and demographic interactions
#' (piDemog).
#' @export
setMethod(
  f = "calcMeanval",
  signature = "LogitBLP",
  definition = function(object){

    shares       <-  object@shares
    prices       <-  object@prices
    idx          <-  object@normIndex
    shareInside  <-  object@shareInside

    # Get parameters from object (accept aliases for alpha)
    alphaMean    <-  if(!is.null(object@slopes$alpha)) object@slopes$alpha else object@slopes$alphaMean
    if(is.null(alphaMean)){
      stop("calcBLPDelta: missing 'alpha' (or 'alphaMean') in slopes.")
    }
    sigma        <-  object@slopes$sigma
    if(is.null(sigma)){
      stop("calcBLPDelta: missing 'sigma' in slopes.")
    }
    piDemog      <-  object@slopes$piDemog
    if(is.null(piDemog)) piDemog <- numeric(0)
    nDemog       <-  object@slopes$nDemog
    if(is.null(nDemog)) nDemog <- length(piDemog)
    if(length(nDemog) != 1L || !is.finite(nDemog) || nDemog < 0 ||
       nDemog != as.integer(nDemog)) {
      stop("'nDemog' must be a non-negative integer.")
    }
    nDemog <- as.integer(nDemog)
    sigmaNest    <-  object@slopes$sigmaNest

    # Check if meanval (delta) is already provided
    deltaProvided <- "meanval" %in% names(object@slopes) && !is.null(object@slopes$meanval)

    ## Reuse the canonical integration and materialization paths.  In
    ## particular, a tensor quadrature matrix represents rows of consumer
    ## draws, rather than length(integration$draws) scalar draws.
    integration <- .blp_object_integration(object, legacy_default = "monte-carlo")
    prodChar <- object@slopes$prodChar
    beta <- object@slopes$beta
    sigmaChar <- object@slopes$sigmaChar
    pi <- object@slopes$pi
    hasChar <- is.matrix(prodChar) && nrow(prodChar) == length(shares)
    materialized <- .blp_materialize_draws(
      integration = integration, alphaMean = alphaMean, sigma = sigma,
      nDemog = nDemog, piDemog = piDemog,
      demogMean = object@slopes$demogMean,
      demogCov = object@slopes$demogCov,
      prodChar = prodChar, sigmaChar = sigmaChar, pi = pi,
      output = object@output,
      storedDemogDraws = object@slopes$demogDraws,
      storedCharDraws = object@slopes$charDraws
    )
    consDraws <- materialized$consDraws
    drawWeights <- materialized$weights
    nDraws <- length(drawWeights)
    demogDraws <- materialized$demogDraws
    alphas <- materialized$alphas
    charDraws <- materialized$charDraws
    char_random <- materialized$char_random

    # Use output slot to verify sign consistency
    output <- object@output
    expectedSign <- ifelse(output, -1, 1)

    # Check if alphaMean has the wrong sign for this market type
    if(sign(alphaMean) != expectedSign && alphaMean != 0){
      warning("Price coefficient sign inconsistent with market type (output=", output,
              "). Expected ", ifelse(output, "negative", "positive"), " alpha, got ",
              round(alphaMean, 4), ". Flipping sign to maintain consistency.")
      alphas <- -alphas
      alphaMean <- -alphaMean
    }

    # Ensure all individual alphas have the correct sign
    wrongSigns <- if(expectedSign > 0) sum(alphas <= 0) else sum(alphas >= 0)
    if(wrongSigns > 0){
      warning(wrongSigns, " out of ", length(alphas), " individual price coefficients have wrong sign. ",
              "Consider reducing sigma or adjusting alphaMean to keep all draws on correct side of zero.")
    }

    nprods <- length(shares)

    if(is.na(idx)){
      idxShare <- 1 - shareInside
      idxPrice <- object@priceOutside
    }
    else{
      idxShare <- shares[idx]
      idxPrice <- prices[idx]
    }

    # Set default sigmaNest if not provided: sigmaNest=1 is flat logit (no nesting)
    if(is.null(sigmaNest)) sigmaNest <- 1

    # When sigmaNest != 1, nesting requires an outside good
    if(!is.na(idx) && sigmaNest != 1){
      idx <- NA
      idxShare <- 1 - sum(shares)
      idxPrice <- object@priceOutside
    }

    if (deltaProvided) {
      delta <- object@slopes$meanval
      message("Using provided meanval (delta) for BLP - skipping contraction mapping")
    } else if (sigmaNest == 1 && nDemog == 0L && !hasChar) {
      ## Keep PriceLeadershipBLP on the same flat-Logit contraction path as
      ## the observed-data BLP adapter.  The branch below remains for nested
      ## and characteristic-rich legacy objects.
      contracted <- .blp_contract(
        prices = prices, shares = shares, alphaMean = alphaMean,
        sigma = sigma, draws = consDraws, weights = drawWeights,
        s0 = if (is.na(idx)) 1 - shareInside else 0,
        priceOutside = idxPrice,
        tol = ifelse(!is.null(object@slopes$contractionTol),
                     object@slopes$contractionTol, 1e-12),
        maxIter = ifelse(!is.null(object@slopes$contractionMaxIter),
                         object@slopes$contractionMaxIter, 1000L)
      )
      if (!contracted$converged) {
        warning("BLP contraction mapping did not converge within the configured iteration limit.")
      }
      delta <- contracted$delta
    } else {
      # Pre-compute constant price term for efficiency
      price_diff <- prices - object@priceOutside

      # Define the fixed point function
      fpFunction <- function(delta) {
        # Compute utilities: nDraws x nProducts matrix
        utilities <- outer(alphas, price_diff, "*")
        utilities <- sweep(utilities, 2, delta, "+")
        if (hasChar) {
          utilities <- utilities + char_random
        }

        # Prevent numeric overflow
        maxUtil <- 700
        utilities <- pmin(pmax(utilities/sigmaNest, -maxUtil), maxUtil)
        expUtil <- exp(utilities)

        sumExpUtil <- rowSums(expUtil)
        insideIV <- sumExpUtil^sigmaNest
        if (is.na(idx)) {
          denom <- 1 + insideIV
        } else {
          denom <- insideIV
        }
        predShares <- as.vector(crossprod(drawWeights, expUtil / denom))

        return(delta + log(shares) - log(predShares))
      }

      # Initial guess
      delta <- log(shares)

      # Use BB::dfsane for fixed point iteration
      tol <- ifelse(!is.null(object@slopes$contractionTol),
                    object@slopes$contractionTol, 1e-12)
      maxIter <- ifelse(!is.null(object@slopes$contractionMaxIter),
                        object@slopes$contractionMaxIter, 1000)

      message("Running BLP contraction with BB::dfsane (tol=", sprintf("%.0e", tol),
              ", maxIter=", maxIter, ")...")

      result <- try(BB::dfsane(par = delta, fn = function(x) fpFunction(x) - x,
                               control = list(tol = tol, maxit = maxIter, trace = FALSE)))

      if (class(result)[1] == "try-error" || result$convergence != 0) {
        # Fall back to dampened fixed point iteration if BB fails
        message("BB::dfsane failed to converge, falling back to dampened fixed point iteration")
        dampFactor <- 0.5
        delta <- log(shares) + log(object@insideSize)

        for (iter in 1:maxIter) {
          deltaNew <- fpFunction(delta)
          absDiff <- max(abs(deltaNew - delta))
          relDiff <- max(abs((deltaNew - delta)/(abs(delta) + 1e-8)))

          if (relDiff < tol || absDiff < tol*1e-2) {
            delta <- deltaNew
            message("BLP contraction converged in ", iter, " iterations")
            break
          }
          delta <- delta + dampFactor * (deltaNew - delta)
        }
        if (iter == maxIter) {
          warning("BLP contraction mapping did not converge within ", maxIter, " iterations")
        }
      } else {
        delta <- result$par
        message("BB::dfsane converged in ", result$iter, " iterations")
      }
    }

    # Store results
    names(delta) <- object@labels
    names(alphaMean) <- "alphaMean"
    names(sigma) <- "sigma"
    names(sigmaNest) <- "sigmaNest"
    if(length(piDemog) > 0){
      names(piDemog) <- paste0("pi_", 1:length(piDemog))
    }

    # Store both alpha and alphaMean for downstream compatibility
    object@slopes$alpha <- as.numeric(alphaMean)
    object@slopes$alphaMean <- alphaMean
    object@slopes$meanval <- delta
    object@slopes$sigma <- sigma
    object@slopes$sigmaNest <- sigmaNest
    object@slopes$piDemog <- piDemog
    object@slopes$nDemog <- nDemog
    object@slopes$alphas <- as.numeric(alphas)
    object@slopes$consDraws <- consDraws
    object@slopes$priceDraws <- materialized$priceDraws
    object@slopes$demogDraws <- demogDraws
    object@slopes$drawWeights <- drawWeights
    object@slopes$integrationWeights <- drawWeights
    object@slopes$integrationWeightsNormalized <- TRUE
    object@slopes$integration <- integration$rule
    object@slopes$integrationPoints <- integration$integrationPoints
    object@slopes$factorOrder <- integration$factorOrder
    object@slopes$nodesPerAxis <- integration$nodesPerAxis
    object@slopes$nNodes <- if(identical(integration$rule, "gauss-hermite")) {
      integration$nodesPerAxis
    } else NULL
    if(hasChar) {
      object@slopes$prodChar <- prodChar
      object@slopes$beta <- beta
      object@slopes$char_random <- char_random
      if(!is.null(sigmaChar)) object@slopes$sigmaChar <- sigmaChar
      if(!is.null(charDraws)) object@slopes$charDraws <- charDraws
      if(!is.null(pi)) object@slopes$pi <- pi
    }
    object@nDraws <- as.numeric(nDraws)

    return(object)
  }
)
