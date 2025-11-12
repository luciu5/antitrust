 #' @include BertrandRUMClasses.R
 #' @title Price Leadership Model (PLE) Calibration and Simulation
 #' @name PriceLeadership-Functions
 #' @aliases ple
 #' @description Calibrates consumer demand using Logit and simulates coordinated effects
 #' using the price leadership model of Mansley, Miller, Sheu & Weinberg (2023).
 #' In this framework, a leader firm announces a supermarkup above Bertrand prices,
 #' coalition firms follow the leader, and fringe firms best-respond. The model
 #' tests incentive compatibility constraints and calibrates the timing parameter
 #' from binding constraints.
 #'   between 0 and 1.
 #' @param margins A length k vector of product margins. At least one coalition margin
 #'   and one fringe margin must be supplied for calibration.
 #' @param diversions A k x k matrix of diversion ratios with diagonal
 #'   elements equal to -1. Default is missing.
 #' @param coalitionPre A numeric vector of product indices that are part of the
 #'   pre-merger coordinating coalition. Products not in this vector are treated as fringe.
 #'   Internally converted to a control matrix.
 #' @param coalitionPost Optional. A numeric vector of product indices in the post-merger
 #'   coalition. If not specified, automatically determined: if a merger brings a coalition
 #'   member together with a non-coalition member under common ownership, both join the
 #'   post-merger coalition. Internally converted to a control matrix.
 #' @param timingParam Optional. If provided, firm-specific timing/discount parameters
 #'   (delta) are fixed at these values rather than calibrated. Should be a named vector
 #'   where names match firm identifiers and values are between 0 and 1. If any IC constraint
 #'   binds, all firm-specific timing parameters are estimated and saved.
 #' @param insideSize An integer equal to total pre-merger units sold.
 #'   If shares sum to one, this also equals the size of the market.
 #' @param normIndex An integer equalling the index (position) of the
 #'   inside product whose mean valuation will be normalized to 1. Default
 #'   is 1, unless \sQuote{shares} sum to less than 1, in which case the default is
 #'   NA and an outside good is assumed to exist.
 #' @param ownerPre EITHER a vector of length k whose values
 #'   indicate which firm produced a product pre-merger OR
 #'   a k x k matrix of pre-merger ownership shares.
 #' @param ownerPost EITHER a vector of length k whose values
 #'   indicate which firm produced a product after the merger OR
 #'   a k x k matrix of post-merger ownership shares.
 #' @param output a length 1 logical vector equal to TRUE if simulation performed
 #'   on output market, FALSE if simulation performed on input market. Default TRUE.
 #' @param mcDelta A vector of length k where each element equals the
 #'   proportional change in a product's marginal costs due to
 #'   the merger. Default is 0, which assumes that the merger does not
 #'   affect any products' marginal cost.
 #' @param subset A vector of length k where each element equals TRUE if
 #'   the product indexed by that element should be included in the
 #'   post-merger simulation and FALSE if it should be excluded. Default is a
 #'   length k vector of TRUE.
 #' @param priceOutside A length 1 vector indicating the price of the
 #'   outside good. Default is 0.
 #' @param priceStart A length k vector of starting values used to solve for
 #'   equilibrium price. Default is the \sQuote{prices} vector.
 #' @param isMax If TRUE, checks to see whether computed price equilibrium
 #'   locally maximizes firm profits and returns a warning if not. Default is FALSE.
 #' @param control.slopes A list of  \code{\link{optim}}
 #'   control parameters passed to the calibration routine.
 #' @param control.equ A list of  \code{\link[BB]{BBsolve}} control parameters
#'   passed to the non-linear equation solver.
#' @param labels A k-length vector of labels. Default is "Prod#", where
#'   \sQuote{#} is a number between 1 and the length of \sQuote{prices}.
#' @param ... Additional options to feed to the optimizer.
#'
#' @details Using observed product prices, shares, and margins from both
#' coalition and fringe firms, \code{ple} calibrates a Logit
#' demand system and identifies:
#' \itemize{
#'   \item The price coefficient (alpha) from fringe firm margins
#'   \item The supermarkup (m) from coalition firm margins
#'   \item Firm-specific timing parameters (delta_f) from binding incentive compatibility constraints
#' }
#'
#' The price leadership equilibrium is characterized by:
#' \itemize{
#'   \item Coalition products priced at p = p_Bertrand + m
#'   \item Fringe products best-responding to coalition prices
#'   \item Incentive compatibility: No coalition firm wants to deviate
#' }
#'
#' If any coalition IC constraint binds, firm-specific timing parameters are estimated
#' for all coalition firms. If all constraints are slack (unconstrained equilibrium),
#' timing parameters cannot be identified from the data.
#'
#' @return An instance of class \code{\linkS4class{PriceLeadership}} with
#' calibrated parameters including supermarkup, firm-specific timing parameters 
#' (if any IC constrained), and slack function values for each coalition firm.
#'
#' @references Mansley, Miller, Sheu & Weinberg (2023), "A price leadership model
#' for coordinated effects in merger analysis"
#'
#' @seealso \code{\link{logit}} for standard Bertrand pricing
#' @seealso \code{\link{calcProducerSurplusGrimTrigger}} for full cartel Grim Trigger analysis
#' @examples
#' \dontrun{
#' # Beer industry example (MSW 2023)
#' prices <- c(0.93, 0.88, 1.10, 1.02)
#' shares <- c(0.35, 0.25, 0.25, 0.15)
#' margins <- c(0.40, NA, 0.25, NA)  # Coalition and fringe margins
#' 
#' result <- ple(
#'   prices = prices,
#'   shares = shares,
#'   margins = margins,
#'   ownerPre = c("AB", "AB", "MC", "Fringe"),
#'   ownerPost = c("AB", "AB", "MC", "Fringe"),
#'   coalitionPre = c(1, 2, 3),  # Products 1-3 in coalition, 4 is fringe
#' )
#' 
#' summary(result)
#' result@supermarkupPre   # Pre-merger equilibrium supermarkup
#' result@supermarkupPost  # Post-merger equilibrium supermarkup
#' result@timingParam      # Calibrated timing parameter
#' result@bindingFirm      # Which firm has binding IC
#' }
#' @include BertrandRUMClasses.R
#' @export
ple <- function(
  prices,
  shares,
  margins,
  diversions,
  ownerPre,
  ownerPost,
  coalitionPre,
  coalitionPost,
  timingParam = NA_real_,
  normIndex = ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE,tolerance=1e-3)),1,NA),
  insideSize = NA_real_,
  output = TRUE,
  mcDelta = rep(0,length(prices)),
  subset = rep(TRUE,length(prices)),
  priceOutside = 0,
  priceStart = prices,
  isMax = FALSE,
  control.slopes,
  control.equ,
  labels = paste("Prod",1:length(prices),sep=""),
  ...
){
  
  ## Check required arguments (before creating object)
  if(missing(coalitionPre)){
    stop("'coalitionPre' must be specified (vector of product indices in pre-merger coordinating coalition)")
  }
  
  nprods <- length(prices)
  
  # Check margin availability
  coalitionMargins <- margins[coalitionPre]
  fringe <- setdiff(1:nprods, coalitionPre)
  fringeMargins <- margins[fringe]
  
  if(all(is.na(coalitionMargins))){
    stop("At least one coalition margin must be supplied for calibration")
  }
  
  if(all(is.na(fringeMargins))){
    stop("At least one fringe margin must be supplied for calibration")
  }
  
  ## Determine post-merger coalition
  ## Logic: If merger brings a coalition member together with a non-coalition member,
  ## assume both join the coalition post-merger (common ownership enables coordination)
  
  ## Check if user specified post-merger coalition
  if(!missing(coalitionPost)){
    # User provided explicit post-merger coalition
    if(length(coalitionPost) < 2){
      stop("'coalitionPost' must contain at least 2 products")
    }
    # Use as-is
  } else {
    # Automatically determine post-merger coalition based on ownership changes
    coalitionPost <- coalitionPre  # Start with same coalition
    
    ## Detect merger changes by comparing ownership
    ## Convert ownership vectors to matrices for comparison (vectorized for performance)
    if(!is.matrix(ownerPre)){
      ownerPreMat <- outer(ownerPre, ownerPre, "==") * 1
    } else {
      ownerPreMat <- ownerPre
    }
    
    if(!is.matrix(ownerPost)){
      ownerPostMat <- outer(ownerPost, ownerPost, "==") * 1
    } else {
      ownerPostMat <- ownerPost
    }
    
    ## Find products that become commonly owned post-merger but weren't pre-merger
    newCommonOwnership <- (ownerPostMat > 0) & (ownerPreMat == 0)
    
    ## For each coalition product, check if it's now commonly owned with any non-coalition product
    for(coalProd in coalitionPre){
      # Find products that are now commonly owned with this coalition product
      newPartners <- which(newCommonOwnership[coalProd, ])
      
      # Add any non-coalition partners to the post-merger coalition
      for(partner in newPartners){
        if(!(partner %in% coalitionPost)){
          coalitionPost <- c(coalitionPost, partner)
          message("Product ", partner, " added to post-merger coalition (merged with coalition product ", coalProd, ")")
        }
      }
    }
    
    ## Sort coalition for consistency
    coalitionPost <- sort(coalitionPost)
  }
  
  ## Create object - validity method will check coalition, timingParam, etc.
  result <- new("PriceLeadership",
                prices = prices,
                shares = shares,
                margins = margins,
                normIndex = normIndex,
                ownerPre = ownerPre,
                ownerPost = ownerPost,
                output = output,
                insideSize = insideSize,
                mcDelta = mcDelta,
                subset = subset,
                priceOutside = priceOutside,
                priceStart = priceStart,
                shareInside = sum(shares),
                labels = labels,
                coalitionPre = coalitionPre,
                coalitionPost = coalitionPost,
                timingParam = timingParam
  )
  
  if(!missing(diversions)){result@diversions <- diversions}
  
  if(!missing(control.slopes)){
    result@control.slopes <- control.slopes
  }
  if(!missing(control.equ)){
    result@control.equ <- control.equ
  }
  
  ## Convert ownership vectors to ownership matrices
  result@ownerPre  <- ownerToMatrix(result, TRUE)
  result@ownerPost <- ownerToMatrix(result, FALSE)
  
  ## Convert coalition vectors to coalition matrices
  result@coalitionPre  <- ownerToMatrix(result, TRUE, control = TRUE)
  result@coalitionPost <- ownerToMatrix(result, FALSE, control = TRUE)
  
  ## Calculate Demand Slope Coefficients, marginal costs, and supermarkups
  ## (calcSlopes now handles both pre and post-merger calibration)
  result <- calcSlopes(result)
  
  ## Solve Non-Linear System for Price Changes
  result@pricePre  <- calcPrices(result, preMerger = TRUE, isMax = isMax, ...)
  result@pricePost <- calcPrices(result, preMerger = FALSE, isMax = isMax, subset = subset, ...)
  
  return(result)
}


#'@rdname PriceLeadership-Functions
#'@export
ple.blp <- function(
  prices,
  shares,
  ownerPre,
  ownerPost,
  coalitionPre,
  coalitionPost,
  timingParam = NA_real_,
  normIndex = ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE,tolerance=1e-3)),1,NA),
  insideSize = NA_real_,
  output = TRUE,
  mcDelta = rep(0,length(prices)),
  subset = rep(TRUE,length(prices)),
  priceOutside = 0,
  priceStart = prices,
  isMax = FALSE,
  nDraws = 500,
  slopes,
  control.slopes,
  control.equ,
  labels = paste("Prod",1:length(prices),sep=""),
  ...
){
  
  nprods <- length(prices)
  
  ## Determine post-merger coalition (if not specified)
  if(missing(coalitionPost)){
    coalitionPost <- .determineCoalitionPost(coalitionPre, ownerPre, ownerPost)
  }
  
  ## Create object - validity method will check coalition, slopes, set defaults, etc.
  ## Note: slopes contains pre-calibrated BLP demand parameters
  result <- new("PriceLeadershipBLP",
                prices = prices,
                shares = shares,
                margins = rep(1/nprods, nprods),  # Not used for BLP
                normIndex = normIndex,
                ownerPre = ownerPre,
                ownerPost = ownerPost,
                coalitionPre = coalitionPre,
                coalitionPost = coalitionPost,
                timingParam = timingParam,
                insideSize = insideSize,
                output = output,
                mcDelta = mcDelta,
                subset = subset,
                priceOutside = priceOutside,
                priceStart = priceStart,
                shareInside = ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE,tolerance=1e-3)),1,sum(shares)),
                slopes = slopes,
                nDraws = nDraws,
                labels = labels)
  
  if(!missing(control.slopes)){
    result@control.slopes <- control.slopes
  }
  if(!missing(control.equ)){
    result@control.equ <- control.equ
  }
  
  ## Convert ownership vectors to ownership matrices
  result@ownerPre  <- ownerToMatrix(result, TRUE)
  result@ownerPost <- ownerToMatrix(result, FALSE)
  
  ## Convert coalition vectors to coalition matrices
  result@coalitionPre  <- ownerToMatrix(result, TRUE, control = TRUE)
  result@coalitionPost <- ownerToMatrix(result, FALSE, control = TRUE)
  
  ## Calculate Demand Slope Coefficients, marginal costs, and supermarkups
  ## (calcSlopes now handles both pre and post-merger calibration)
  result <- calcSlopes(result)
  
  ## Solve Non-Linear System for Price Changes
  result@pricePre  <- calcPrices(result, preMerger = TRUE, isMax = isMax, ...)
  result@pricePost <- calcPrices(result, preMerger = FALSE, isMax = isMax, subset = subset, ...)
  
  return(result)
}


## Helper function: Extract owner vector from matrix or vector
## Used to avoid repeated apply() operations
.extractOwnerVec <- function(owner){
  if(is.vector(owner) || is.factor(owner)){
    return(owner)
  }
  # Extract from matrix - use row indices as firm IDs
  apply(owner, 1, function(row) which(row > 0)[1])
}


## Helper function: Determine post-merger coalition
## Automatically expands coalition if merger brings coalition member together with non-coalition member
.determineCoalitionPost <- function(coalitionPre, ownerPre, ownerPost){
  
  coalitionPost <- coalitionPre  # Start with same coalition
  
  ## Convert ownership vectors to matrices for comparison (vectorized for performance)
  if(!is.matrix(ownerPre)){
    ownerPreMat <- outer(ownerPre, ownerPre, "==") * 1
  } else {
    ownerPreMat <- ownerPre
  }
  
  if(!is.matrix(ownerPost)){
    ownerPostMat <- outer(ownerPost, ownerPost, "==") * 1
  } else {
    ownerPostMat <- ownerPost
  }
  
  ## Find products that become commonly owned post-merger but weren't pre-merger
  newCommonOwnership <- (ownerPostMat > 0) & (ownerPreMat == 0)
  
  ## For each coalition product, check if it's now commonly owned with any non-coalition product
  for(coalProd in coalitionPre){
    # Find products that are now commonly owned with this coalition product
    newPartners <- which(newCommonOwnership[coalProd, ])
    
    # Add any non-coalition partners to the post-merger coalition
    for(partner in newPartners){
      if(!(partner %in% coalitionPost)){
        coalitionPost <- c(coalitionPost, partner)
        message("Product ", partner, " added to post-merger coalition (merged with coalition product ", coalProd, ")")
      }
    }
  }
  
  ## Sort coalition for consistency
  return(sort(coalitionPost))
}


## Helper function: Extract owner vector from matrix or vector
## Used to avoid repeated apply() operations
.extractOwnerVec <- function(owner){
  if(is.vector(owner) || is.factor(owner)){
    return(owner)
  }
  # Extract from matrix - use row indices as firm IDs
  apply(owner, 1, function(row) which(row > 0)[1])
}


## Generic and method for calibrating price leadership parameters
## Shared by PriceLeadership and PriceLeadershipBLP classes

setGeneric(
  name = "calcPriceLeadershipParams",
  def = function(object, ...) {standardGeneric("calcPriceLeadershipParams")}
)

#' Calculate Price Leadership Parameters (Supermarkup and Timing)
#' 
#' @description Internal method that calibrates price leadership equilibrium parameters:
#' supermarkup and firm-specific timing parameters. This method implements the core
#' IC constraint logic shared by both PriceLeadership and PriceLeadershipBLP.
#' 
#' @param object A PriceLeadership or PriceLeadershipBLP object with calibrated demand
#' @param preMerger Logical, if TRUE calibrate pre-merger params, else post-merger
#' @param pricesColl Coordinated/collusive prices for the coalition
#' @param coalition Vector of product indices in the coalition
#' 
#' @return List with elements:
#'   \itemize{
#'     \item supermarkup - Equilibrium supermarkup above coordinated prices
#'     \item timingParam - Named vector of firm-specific timing parameters
#'     \item bindingFirm - Firm with binding IC constraint (or NA)
#'     \item slackValues - IC slack values by firm
#'   }
#' 
#' @details Strategy:
#' 1. Calculate profits at collusive prices (π^PL)
#' 2. Calculate Bertrand punishment profits (π^B)
#' 3. Calculate optimal deviation profits for each firm (π^D)
#' 4. Identify supermarkup and check IC constraints
#' 5. Estimate firm-specific timing parameters if any IC binds
#' 
#' @keywords internal
setMethod(
  f = "calcPriceLeadershipParams",
  signature = "PriceLeadership",
  definition = function(object, preMerger = TRUE, pricesColl, coalition){
    
    prices <- object@prices
    nprods <- length(prices)
    tol <- 1e-6
    
    ownerPre <- object@ownerPre
    
    # Supermarkup that best explains observed prices
    supermarkup_observed <- mean(prices[coalition] - pricesColl[coalition])
    
    ## Check IC constraints at full collusion
    ## IC: π^PL - π^D + δ/(1-δ) * (π^PL - π^B) ≥ 0
    
    # Calculate profits at full collusion prices
    if(preMerger){
      object@pricePre <- pricesColl
    } else {
      object@pricePost <- pricesColl
    }
    profitsPL <- calcProducerSurplus(object, preMerger = preMerger)
    
    # Calculate Bertrand and deviation profits
    # Temporarily swap owner matrix to break coordination
    ownerOriginal <- if(preMerger) object@ownerPre else object@ownerPost
    if(preMerger){
      object@ownerPre <- ownerToMatrix(object, preMerger = TRUE, control = FALSE)
    } else {
      object@ownerPost <- ownerToMatrix(object, preMerger = FALSE, control = FALSE)
    }
    
    # Bertrand profits
    pBertrand <- calcPrices(object, preMerger = preMerger, isMax = FALSE)
    if(preMerger){
      object@pricePre <- pBertrand
    } else {
      object@pricePost <- pBertrand
    }
    profitsBertrand <- calcProducerSurplus(object, preMerger = preMerger)
    
    # Restore original owner matrix
    if(preMerger){
      object@ownerPre <- ownerOriginal
    } else {
      object@ownerPost <- ownerOriginal
    }
    
    # Calculate optimal deviation profits for each coalition firm
    ownerVec_temp <- .extractOwnerVec(ownerOriginal)
    coalitionFirms_temp <- unique(ownerVec_temp[coalition])
    
    profitsDeviation_byFirm_temp <- sapply(coalitionFirms_temp, function(firm){
      firmProducts <- which(ownerVec_temp == firm)
      deviatingCoalitionProducts <- intersect(firmProducts, coalition)
      
      if(length(deviatingCoalitionProducts) == 0) return(0)
      
      subset <- rep(FALSE, nprods)
      subset[deviatingCoalitionProducts] <- TRUE
      
      optimalPrices <- calcPrices(object, preMerger = preMerger, subset = subset, regime = "deviation")
      
      if(preMerger){
        object@pricePre <- optimalPrices
      } else {
        object@pricePost <- optimalPrices
      }
      
      deviationProfit <- calcProducerSurplus(object, preMerger = preMerger)
      sum(deviationProfit[firmProducts])
    })
    
    # Get firm identifiers for aggregation
    ownerVec <- .extractOwnerVec(if(preMerger) ownerPre else object@ownerPost)
    coalitionFirms <- unique(ownerVec[coalition])
    nCoalitionFirms <- length(coalitionFirms)
    
    # Aggregate profits by firm
    profitsPL_byFirm <- sapply(coalitionFirms, function(firm){
      sum(profitsPL[ownerVec == firm])
    })
    
    profitsDeviation_byFirm <- profitsDeviation_byFirm_temp
    
    profitsBertrand_byFirm <- sapply(coalitionFirms, function(firm){
      sum(profitsBertrand[ownerVec == firm])
    })
    
    names(profitsPL_byFirm) <- coalitionFirms
    names(profitsDeviation_byFirm) <- as.character(coalitionFirms)
    names(profitsBertrand_byFirm) <- coalitionFirms
    
    # Check IC slack for each firm
    IC_slack_byFirm <- profitsPL_byFirm - profitsDeviation_byFirm
    
    ## Identify supermarkup and timing parameters based on IC constraints
    if(supermarkup_observed <= tol){
      # No supermarkup observed - just coordinated Bertrand
      supermarkup <- 0
      timingParam <- numeric()
      bindingFirm <- NA_integer_
      
    } else if(all(IC_slack_byFirm >= -tol)){
      # All ICs satisfied at δ=1 (static constraint holds)
      supermarkup <- supermarkup_observed
      timingParam <- rep(1, nCoalitionFirms)
      names(timingParam) <- coalitionFirms
      bindingFirm <- NA_integer_
      
    } else {
      # At least one IC binds - estimate firm-specific timing parameters
      # δ_f = (π^D_f - π^PL_f) / (π^D_f - π^B_f)
      
      timingParam <- sapply(coalitionFirms, function(firm){
        delta_f <- (profitsDeviation_byFirm[as.character(firm)] - profitsPL_byFirm[as.character(firm)]) / 
                   (profitsDeviation_byFirm[as.character(firm)] - profitsBertrand_byFirm[as.character(firm)])
        
        # Constrain to valid range
        if(is.na(delta_f) || delta_f < 0){
          return(0)
        } else if(delta_f >= 1){
          return(0.999)
        } else {
          return(delta_f)
        }
      })
      
      names(timingParam) <- coalitionFirms
      bindingFirm <- coalitionFirms[which.min(timingParam)]
      supermarkup <- supermarkup_observed
    }
    
    # Store slack values for diagnostics
    slackValues <- IC_slack_byFirm
    names(slackValues) <- coalitionFirms
    
    return(list(
      supermarkup = supermarkup,
      timingParam = timingParam,
      bindingFirm = bindingFirm,
      slackValues = slackValues
    ))
  }
)


setMethod(
  f = "calcPriceLeadershipParams",
  signature = "PriceLeadershipBLP",
  definition = function(object, preMerger = TRUE, pricesColl, coalition) {
    selectMethod("calcPriceLeadershipParams", "PriceLeadership")(object, preMerger, pricesColl, coalition)
  }
)

## Generic and method for BLP contraction mapping
## Shared by LogitBLP and PriceLeadershipBLP classes

setGeneric(
  name = "calcMeanval",
  def = function(object, ...) {standardGeneric("calcMeanval")}
)

#' Calculate Mean Utilities via Contraction Mapping
#' 
#' @description Internal method that runs the BLP contraction mapping to recover
#' mean utilities (delta/meanval) that rationalize observed shares given random
#' coefficient parameters. This method is shared by LogitBLP and PriceLeadershipBLP.
#' 
#' @param object A LogitBLP or PriceLeadershipBLP object with BLP parameters
#' 
#' @return The object with updated slopes$meanval (delta vector)
#' 
#' @details Implements the BLP contraction mapping:
#' delta^(t+1) = delta^(t) + log(s_observed) - log(s_predicted(delta^(t)))
#' 
#' Uses BB::dfsane for robust convergence with fallback to dampened fixed point.
#' Handles nesting parameters (sigmaNest) and demographic interactions (piDemog).
#' 
#' @keywords internal
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
    nDraws       <-  object@nDraws
    piDemog      <-  object@slopes$piDemog
    if(is.null(piDemog)) piDemog <- numeric(0)
    nDemog       <-  object@slopes$nDemog
    if(is.null(nDemog)) nDemog <- 0
    sigmaNest    <-  object@slopes$sigmaNest
    
    # Check if meanval (delta) is already provided
    deltaProvided <- "meanval" %in% names(object@slopes) && !is.null(object@slopes$meanval)
    
    # Check if draws already exist (to ensure consistency across calls)
    drawsExist <- "consDraws" %in% names(object@slopes) && !is.null(object@slopes$consDraws)
    
    if(drawsExist) {
      # Reuse existing draws
      consDraws <- object@slopes$consDraws
      demogDraws <- object@slopes$demogDraws
      if(!is.null(demogDraws) && ncol(demogDraws) > 0) {
        demogEffect <- demogDraws %*% piDemog
      } else {
        demogEffect <- 0
      }
    } else {
      # Generate new consumer heterogeneity draws (unobserved)
      consDraws <- rnorm(nDraws)
      
      # Generate demographic draws (observed heterogeneity)
      if(!is.null(nDemog) && nDemog > 0){
        demogDraws <- matrix(rnorm(nDraws * nDemog), nrow=nDraws, ncol=nDemog)
        demogEffect <- demogDraws %*% piDemog
      } else {
        demogDraws <- matrix(nrow=nDraws, ncol=0)
        demogEffect <- 0
        piDemog <- numeric(0)
        nDemog <- 0
      }
    }
    
    # Compute individual-specific price coefficients
    # alpha_i = alphaMean + sigma * nu_i + pi * d_i
    alphas <- alphaMean + sigma * consDraws + demogEffect
    
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
    }
    else {
      # Pre-compute constant price term for efficiency
      price_diff <- prices - object@priceOutside
      
      # Define the fixed point function
      fpFunction <- function(delta) {
        # Compute utilities: nDraws x nProducts matrix
        utilities <- outer(alphas, price_diff, "*")
        utilities <- sweep(utilities, 2, delta, "+")
        
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
        predShares <- colMeans(expUtil/denom)
        
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
    object@slopes$demogDraws <- demogDraws
    
    return(object)
  }
)


## ownerToMatrix method for PriceLeadership class
## Returns ownership matrix accounting for coalition coordination when control=TRUE

#' @rdname Ownership-methods
#' @param control For PriceLeadership class only. If TRUE, returns the effective 
#'   control matrix where coalition products are treated as coordinating (full 
#'   cross-ownership). If FALSE (default), returns traditional ownership matrix.
#' @export
setMethod(
  f = "ownerToMatrix",
  signature = "PriceLeadership",
  definition = function(object, preMerger = TRUE, control = FALSE){
    
    # Get traditional ownership matrix first
    if(preMerger){
      thisOwner <- object@ownerPre
      coalition <- object@coalitionPre
    } else {
      thisOwner <- object@ownerPost
      coalition <- object@coalitionPost
    }
    
    # Convert to matrix if needed (using parent method logic)
    if(is.vector(thisOwner) || is.factor(thisOwner)){
      nprod <- length(object@labels)
      owners <- as.numeric(factor(thisOwner, levels = unique(thisOwner)))
      ownerMat <- matrix(0, ncol = nprod, nrow = nprod)
      
      for(o in unique(owners)){
        ownerMat[owners == o, owners == o] <- 1
      }
    } else {
      ownerMat <- thisOwner
    }
    
    # If control=TRUE, modify to reflect coalition coordination
    if(control && length(coalition) > 0){
      # Coalition products coordinate pricing (act as single entity)
      # Set full cross-ownership for coalition products
      ownerMat[coalition, coalition] <- 1
    }
    
    return(ownerMat)
  }
)


#' @rdname Ownership-methods
#' @export
setMethod(
  f = "ownerToMatrix",
  signature = "PriceLeadershipBLP",
  definition = function(object, preMerger = TRUE, control = FALSE){
    # Reuse PriceLeadership implementation - logic is identical
    # Both classes have coalitionPre/Post and ownerPre/Post slots
    selectMethod("ownerToMatrix", "PriceLeadership")(object, preMerger, control)
  }
)


## calcSlopes method for PriceLeadership class
## Calibrates demand parameters (alpha, meanval), supermarkup, and timing parameter

#' @rdname Params-Methods
#' @export
setMethod(
  f = "calcSlopes",
  signature = "PriceLeadership",
  definition = function(object){
    
    ## Extract data
    prices <- object@prices
    shares <- object@shares
    margins <- object@margins
    
    # Extract coalition indices (handle both vector and matrix forms)
    if(is.matrix(object@coalitionPre)){
      coalition <- which(rowSums(object@coalitionPre) > 1)
    } else {
      coalition <- object@coalitionPre
    }
    
    fringe <- setdiff(1:length(prices), coalition)
    ownerPre <- object@ownerPre
    idx <- object@normIndex
    shareInside <- object@shareInside
    output <- object@output
    outSign <- ifelse(output, -1, 1)
    tol <- 1e-6  # Tolerance for IC constraint checking
    
    nprods <- length(prices)
    
    if(is.na(idx)){
      idxShare <- 1 - shareInside
      idxPrice <- object@priceOutside
    } else {
      idxShare <- shares[idx]
      idxPrice <- prices[idx]
    }
    
    ## Step 1: Calibrate alpha (price coefficient) using fringe margins
    ## Use minimum distance estimation with all observed fringe margins
    
    # Find all fringe products with observed margins
    fringeWithMargins <- fringe[!is.na(margins[fringe])]
    
    if(length(fringeWithMargins) == 0){
      stop("At least one fringe margin required to identify price coefficient")
    }
    
    ## Minimum distance objective function
    ## For fringe firms playing Bertrand: margin_f = -1 / (alpha * p_f * (1 - share_f))
    ## In output markets: alpha < 0, so margin = -1 / (alpha * p * (1-s)) = 1 / (|alpha| * p * (1-s))
    ## Find alpha that minimizes sum of squared differences between observed and predicted margins
    
    objectiveFn <- function(alpha_cand){
      # Predicted margins for fringe products
      # margin = -1 / (alpha * p * (1-s)) where alpha has sign from output convention
      predMargins <- -1 / (alpha_cand * prices[fringeWithMargins] * (1 - shares[fringeWithMargins]))
      
      # Observed margins
      obsMargins <- margins[fringeWithMargins]
      
      # Sum of squared errors
      sse <- sum((obsMargins - predMargins)^2)
      
      return(sse)
    }
    
    # Initial guess from first fringe margin
    # Solve: margin = -1 / (alpha * p * (1-s)) for alpha
    # alpha = -1 / (margin * p * (1-s))
    alpha_init <- -1 / (margins[fringeWithMargins[1]] * prices[fringeWithMargins[1]] * (1 - shares[fringeWithMargins[1]]))
    
    # Optimize
    # For output markets (output=TRUE): alpha should be negative
    # For input markets (output=FALSE): alpha should be positive
    if(output){
      lowerBound <- -100
      upperBound <- -0.001
    } else {
      lowerBound <- 0.001
      upperBound <- 100
    }
    
    optResult <- optim(
      par = alpha_init,
      fn = objectiveFn,
      method = "Brent",
      lower = lowerBound,
      upper = upperBound,
      control = object@control.slopes
    )
    
    if(optResult$convergence != 0){
      warning("Alpha calibration may not have converged. Using result anyway.")
    }
    
    alpha <- optResult$par
    
    ## Step 2: Recover mean valuations from shares
    ## shares_j = exp(meanval_j + alpha*(p_j - p_0)) / (1 + sum(exp(...)))
    
    meanval <- log(shares) - log(idxShare) - alpha * (prices - idxPrice)
    
    ## Store demand parameters (needed for MC and price calculations)
    object@slopes <- list(
      alpha = alpha,
      meanval = meanval
    )
    object@mktSize <- object@insideSize / shareInside
    
    ## Calculate marginal costs using observed margins
    ## This will properly handle multi-product firms and missing margins
    object@mcPre <- calcMC(object, preMerger = TRUE)
    
    ## Step 3: Identify supermarkup and timing parameter using shared method
    ## Strategy: 
    ## 1. Solve for coordinated Bertrand prices (coalition coordinates, m=0)
    ## 2. Calculate IC-constrained equilibrium using calcPriceLeadershipParams
    
    # Solve for coordinated Bertrand prices (coalition coordinates with m=0)
    pricesColl <- calcPrices(object, preMerger = TRUE, isMax = FALSE, regime = "coordination")
    
    # Calculate price leadership parameters (supermarkup, timing, IC slack)
    pleParams <- calcPriceLeadershipParams(object, preMerger = TRUE, 
                                           pricesColl = pricesColl, 
                                           coalition = coalition)
    
    object@supermarkupPre <- pleParams$supermarkup
    
    ## Store timing parameters (user can override by providing timingParam in constructor)
    if(length(object@timingParam) == 0){
      object@timingParam <- pleParams$timingParam
    }
    
    object@bindingFirm <- pleParams$bindingFirm
    object@slackValues <- pleParams$slackValues
    
    ## Step 5: Calculate post-merger marginal costs
    object@mcPost <- calcMC(object, preMerger = FALSE)
    
    ## Step 6: Calibrate post-merger supermarkup
    # Get post-merger coalition
    if(is.matrix(object@coalitionPost)){
      coalitionPost <- which(rowSums(object@coalitionPost) > 1)
    } else {
      coalitionPost <- object@coalitionPost
    }
    
    if(length(coalitionPost) == 0){
      # No post-merger coalition - no supermarkup
      object@supermarkupPost <- 0
    } else {
      # Calculate post-merger coordinated/collusive prices
      # With coalitionPost (control=TRUE), calcPrices with regime="coordination"
      # already gives monopoly prices (coalition treated as single firm)
      pricesCollPost <- calcPrices(object, preMerger = FALSE, isMax = FALSE, regime = "coordination")
      
      # Check if we identified timing parameters from pre-merger IC binding
      timingNonNA <- object@timingParam[!is.na(object@timingParam)]
      if(length(timingNonNA) > 0 && isTRUE(any(timingNonNA < 1))){
        # Timing parameters were identified - use them to constrain post-merger supermarkup
        # Solve for maximum sustainable supermarkup given firm-specific deltas
        object@supermarkupPost <- calcSupermarkup(
          object, 
          preMerger = FALSE, 
          constrained = TRUE
        )
        
      } else {
        # Timing parameters not identified (pre-merger ICs were slack or all delta=1)
        # No additional supermarkup beyond collusive prices
        object@supermarkupPost <- calcSupermarkup(
          object,
          preMerger = FALSE,
          constrained = FALSE
        )
      }
    }
    
    return(object)
  }
)


## calcProducerSurplus method for PriceLeadershipBLP class
## Reuses PriceLeadership implementation - logic is identical

#' @rdname PS-methods
#' @export
setMethod(
  f = "calcProducerSurplus",
  signature = "PriceLeadershipBLP",
  definition = function(object, preMerger = TRUE, regime = c("coordination", "bertrand"), ...){
    # Delegate to PriceLeadership implementation
    selectMethod("calcProducerSurplus", "PriceLeadership")(object, preMerger, regime, ...)
  }
)


## calcMargins method for PriceLeadershipBLP class
## Reuses PriceLeadership implementation - logic is identical

#' @rdname Margins-Methods
#' @export
setMethod(
  f = "calcMargins",
  signature = "PriceLeadershipBLP",
  definition = function(object, preMerger = TRUE, level = FALSE, 
                        regime = c("coordination", "bertrand", "deviation", "constrained")){
    # Delegate to PriceLeadership implementation
    selectMethod("calcMargins", "PriceLeadership")(object, preMerger, level, regime)
  }
)


## calcPrices method for PriceLeadershipBLP class
## Reuses PriceLeadership implementation - logic is identical

#' @rdname Prices-Methods
#' @export
setMethod(
  f = "calcPrices",
  signature = "PriceLeadershipBLP",
  definition = function(object, preMerger = TRUE, isMax = FALSE, subset, 
                        regime = c("coordination", "bertrand", "deviation", "constrained"), ...){
    # Delegate to PriceLeadership implementation
    selectMethod("calcPrices", "PriceLeadership")(object, preMerger, isMax, subset, regime, ...)
  }
)


## calcSlopes method for PriceLeadershipBLP class
## Takes BLP demand parameters as given (from slopes), calibrates only price leadership parameters

#' @rdname Params-Methods
#' @export
setMethod(
  f = "calcSlopes",
  signature = "PriceLeadershipBLP",
  definition = function(object){
    
    ## For LogitBLP with price leadership:
    ## 1. BLP demand parameters are PRE-CALIBRATED (provided in slopes list)
    ## 2. Use BLP contraction mapping to ensure meanval is consistent with shares
    ## 3. Calibrate only price leadership parameters: supermarkup, timing
    ##
    ## This differs from PriceLeadership which calibrates both demand AND leadership params
    
    nprods <- length(object@prices)
    
    # Extract coalition indices (handle both vector and matrix forms)
    if(is.matrix(object@coalitionPre)){
      coalition <- which(rowSums(object@coalitionPre) > 1)
    } else {
      coalition <- object@coalitionPre
    }
    
    ## Validate BLP parameters (deferred from validity method)
    if(is.null(object@slopes) || !is.list(object@slopes)){
      stop("BLP demand parameters must be provided in 'slopes' list")
    }
    
    if(!("sigma" %in% names(object@slopes))){
      stop("'slopes' must contain 'sigma' (std. dev. of random price coefficient)")
    }
    
    # Alpha validation - accept either alpha or alphaMean
    if(!("alpha" %in% names(object@slopes)) && !("alphaMean" %in% names(object@slopes))){
      stop("'slopes' must contain either 'alpha' or 'alphaMean' (mean price coefficient)")
    }
    
  ## Run BLP contraction mapping to get meanval (delta) via shared method
  ## This ensures delta is consistent with observed shares given random coefficients
  object <- calcMeanval(object)
    
    ## Set market size if not already set
    if(is.na(object@mktSize) || length(object@mktSize) == 0){
      shareInside <- object@shareInside
      object@mktSize <- object@insideSize / shareInside
    }
    
    ## Calculate marginal costs using BLP-calibrated demand
    object@mcPre <- calcMC(object, preMerger = TRUE)
    
    ## Identify supermarkup and timing parameter using shared method
    # Solve for coordinated Bertrand prices (coalition coordinates with m=0)
    pricesColl <- calcPrices(object, preMerger = TRUE, isMax = FALSE, regime = "coordination")
    
    # Calculate price leadership parameters (supermarkup, timing, IC slack)
    pleParams <- calcPriceLeadershipParams(object, preMerger = TRUE, 
                                           pricesColl = pricesColl, 
                                           coalition = coalition)
    
    object@supermarkupPre <- pleParams$supermarkup
    
    ## Store timing parameters (user can override by providing timingParam in constructor)
    if(length(object@timingParam) == 0){
      object@timingParam <- pleParams$timingParam
    }
    
    object@bindingFirm <- pleParams$bindingFirm
    object@slackValues <- pleParams$slackValues
    
    ## Calculate post-merger marginal costs
    object@mcPost <- calcMC(object, preMerger = FALSE)
    
    ## Post-merger supermarkup
    if(is.matrix(object@coalitionPost)){
      coalitionPost <- which(rowSums(object@coalitionPost) > 1)
    } else {
      coalitionPost <- object@coalitionPost
    }
    
    if(length(coalitionPost) == 0){
      object@supermarkupPost <- 0
    } else {
      pricesCollPost <- calcPrices(object, preMerger = FALSE, isMax = FALSE, regime = "coordination")
      
      timingNonNA <- object@timingParam[!is.na(object@timingParam)]
      if(length(timingNonNA) > 0 && isTRUE(any(timingNonNA < 1))){
        object@supermarkupPost <- calcSupermarkup(object, preMerger = FALSE, constrained = TRUE)
      } else {
        object@supermarkupPost <- calcSupermarkup(object, preMerger = FALSE, constrained = FALSE)
      }
    }
    
    return(object)
  }
)


## calcProducerSurplus method for PriceLeadership class
## Profit calculation supporting multiple pricing regimes

#' @rdname PS-methods
#' @param regime Character string specifying pricing regime for profit evaluation.
#'   Options:
#'   \itemize{
#'     \item "bertrand" - Solve for Bertrand prices then evaluate
#'     \item "coordination" - Solve for coordination prices then evaluate (default)
#'   }
#'   Note: For deviation profits, use \code{\link{calcDeviationProfit}} directly.
#'   To evaluate at current prices without solving, use the parent class method directly.
#' @export
setMethod(
  f = "calcProducerSurplus",
  signature = "PriceLeadership",
  definition = function(object, preMerger = TRUE,
                        regime = c("coordination", "bertrand"), ...){
    
    regime <- match.arg(regime)
    
    # Store original prices to restore later
    pricesOriginal <- if(preMerger) object@pricePre else object@pricePost
    
    # Solve for equilibrium prices based on regime
    newPrices <- calcPrices(object, preMerger = preMerger, 
                             regime = regime, isMax = FALSE, ...)
    
    # Set prices temporarily
    if(preMerger){
      object@pricePre <- newPrices
    } else {
      object@pricePost <- newPrices
    }
    
    # Calculate profits at solved prices using parent method
    # Returns product-level variable profits (producer surplus absent fixed costs)
    profits <- callNextMethod(object, preMerger = preMerger)
    
    # Restore original prices
    if(preMerger){
      object@pricePre <- pricesOriginal
    } else {
      object@pricePost <- pricesOriginal
    }
    
    return(profits)
  }
)


## Generic for calcSupermarkup
setGeneric(
  name = "calcSupermarkup",
  def = function(object, ...) {standardGeneric("calcSupermarkup")}
)

#' @rdname calcSupermarkup-methods
#' Calculate Equilibrium Supermarkup
#' 
#' @description Calculates the equilibrium supermarkup for the price leadership model.
#' Can compute either unconstrained (full collusion) or IC-constrained supermarkup.
#' 
#' @param object A PriceLeadership object
#' @param preMerger Logical, if TRUE use pre-merger structure, FALSE for post-merger
#' @param constrained Logical, if TRUE find maximum IC-constrained supermarkup,
#'   if FALSE return unconstrained (full collusion) supermarkup
#' @param ... Additional arguments (not currently used)
#' 
#' @return Scalar supermarkup value
#' 
#' @details Extracts coalition from object@@coalitionPre or object@@coalitionPost 
#' based on preMerger flag.
#' 
#' When constrained=FALSE, computes unconstrained monopoly supermarkup by:
#' 1. Solving for monopoly prices via calcPrices(regime="coordination")
#' 2. Comparing to observed prices (pre-merger) or coordinated prices (post-merger)
#' 3. Returning the average difference for coalition products
#' 
#' When constrained=TRUE, uses binary search to find maximum supermarkup m such that
#' all coalition firms' IC constraints are satisfied using firm-specific timing 
#' parameters from object@@timingParam:
#' IC_f: π^PL_f(m) - π^D_f(m) + δ_f/(1-δ_f) * [π^PL_f(m) - π^B_f] ≥ 0
#' 
#' @export
setMethod(
  f = "calcSupermarkup",
  signature = "PriceLeadership",
  definition = function(object, preMerger = TRUE, constrained = FALSE, ...){
    # Extract coalition indices from appropriate slot
    if(preMerger){
      if(is.matrix(object@coalitionPre)){
        coalition <- which(rowSums(object@coalitionPre) > 1)
      } else {
        coalition <- object@coalitionPre
    }
  } else {
    if(is.matrix(object@coalitionPost)){
      coalition <- which(rowSums(object@coalitionPost) > 1)
    } else {
      coalition <- object@coalitionPost
    }
  }
  
  if(!constrained){
    # Unconstrained case: calculate monopoly/collusive supermarkup
    # This is the difference between monopoly prices and coordinated Bertrand prices
    
    # Get monopoly/collusive prices (coalition acts as single firm)
    pricesMonopoly <- calcPrices(object, preMerger = preMerger, regime = "coordination")
    
    if(preMerger){
      # Pre-merger: compare monopoly to observed prices
      # Supermarkup that would explain observed behavior as monopoly pricing
      observedPrices <- object@prices
      supermarkup <- mean(observedPrices[coalition] - pricesMonopoly[coalition])
    } else {
      # Post-merger: monopoly prices ARE the coordinated prices
      # So unconstrained supermarkup is 0 (already at monopoly)
      supermarkup <- 0
    }
    
    return(supermarkup)
  }
  
  # Constrained case: find max IC-constrained supermarkup
  
  # Get firm-specific timing parameters from object slot
  timingParam <- object@timingParam
  
  if(length(timingParam) == 0){
    warning("Constrained supermarkup requested but no timing parameters identified. Returning 0.")
    return(0)
  }
  
  # Validate timing parameters are non-negative
  if(isTRUE(any(timingParam < 0, na.rm = TRUE))){
    stop("Timing parameters must be non-negative (>= 0)")
  }
  
  # Get Bertrand coordinated prices
  pricesBertrand <- calcPrices(object, preMerger = preMerger, regime = "coordination")
  
  # Get ownership structure to identify firms
  ownerPre <- object@ownerPre
  ownerPost <- object@ownerPost
  owner <- if(preMerger) ownerPre else ownerPost
  
  # Extract owner vector (cached for efficiency)
  ownerVec <- .extractOwnerVec(owner)
  
  # Get coalition firms
  coalitionFirms <- unique(ownerVec[coalition])
  
  # Calculate Bertrand profits ONCE (independent of supermarkup m)
  # This avoids recalculating in every binary search iteration
  priceBertrand_fixed <- calcPrices(object, preMerger = preMerger, regime = "bertrand")
  if(preMerger){
    object@pricePre <- priceBertrand_fixed
  } else {
    object@pricePost <- priceBertrand_fixed
  }
  profitsBertrand_fixed <- calcProducerSurplus(object, preMerger = preMerger)
  
  # Aggregate Bertrand profits by firm (reused in every iteration)
  profitsBertrand_byFirm <- sapply(coalitionFirms, function(firm){
    sum(profitsBertrand_fixed[ownerVec == firm])
  })
  names(profitsBertrand_byFirm) <- coalitionFirms
  
  # Objective: Find maximum m such that IC_f holds for all coalition firms f
  # IC_f: π^PL_f(m) - π^D_f(m) + δ_f/(1-δ_f) * (π^PL_f(m) - π^B_f) ≥ 0
  
  # Binary search for maximum sustainable supermarkup
  m_low <- 0
  m_high <- 5 * mean(pricesBertrand[coalition])  # Upper bound
  tol <- 1e-4
  
  while(m_high - m_low > tol){
    m_mid <- (m_low + m_high) / 2
    
    # Calculate PL profits at this supermarkup
    if(preMerger){
      object@pricePre <- pricesBertrand
      object@pricePre[coalition] <- pricesBertrand[coalition] + m_mid
      profitsPL <- calcProducerSurplus(object, preMerger = TRUE)
    } else {
      object@pricePost <- pricesBertrand
      object@pricePost[coalition] <- pricesBertrand[coalition] + m_mid
      profitsPL <- calcProducerSurplus(object, preMerger = FALSE)
    }
    
    # Calculate proper deviation profits for each coalition firm
    # Each firm optimizes while others maintain PL prices (pricesBertrand + m_mid)
    plPricesWithMarkup <- pricesBertrand
    plPricesWithMarkup[coalition] <- pricesBertrand[coalition] + m_mid
    
    # Store original prices
    pricesOriginal <- if(preMerger) object@pricePre else object@pricePost
    
    # Calculate optimal deviation for each coalition firm
    profitsDeviation_byFirm <- sapply(coalitionFirms, function(firm){
      # Set all products to PL prices
      if(preMerger){
        object@pricePre <- plPricesWithMarkup
      } else {
        object@pricePost <- plPricesWithMarkup
      }
      
      # Find which products this firm owns
      firmProducts <- which(ownerVec == firm)
      
      # Filter to only coalition products owned by deviating firm
      deviatingCoalitionProducts <- intersect(firmProducts, coalition)
      
      if(length(deviatingCoalitionProducts) == 0) return(0)
      
      # Create subset vector: TRUE only for deviating firm's coalition products
      subset <- rep(FALSE, length(ownerVec))
      subset[deviatingCoalitionProducts] <- TRUE
      
      # Solve for optimal deviation prices (deviation regime)
      optimalPrices <- calcPrices(object, preMerger = preMerger, subset = subset, regime = "deviation")
      
      # Set prices to optimal deviation
      if(preMerger){
        object@pricePre <- optimalPrices
      } else {
        object@pricePost <- optimalPrices
      }
      
      # Calculate firm's profit at optimal deviation
      deviationProfit <- calcProducerSurplus(object, preMerger = preMerger)
      
      return(sum(deviationProfit[firmProducts]))
    })
    
    names(profitsDeviation_byFirm) <- coalitionFirms
    
    # Restore original prices
    if(preMerger){
      object@pricePre <- pricesOriginal
    } else {
      object@pricePost <- pricesOriginal
    }
    
    # Aggregate PL profits by firm
    profitsPL_byFirm <- sapply(coalitionFirms, function(firm){
      sum(profitsPL[ownerVec == firm])
    })
    
    # Note: profitsBertrand_byFirm already calculated before loop
    
    # Check IC for each coalition firm using firm-specific timing parameter
    IC_slack_byFirm <- sapply(coalitionFirms, function(firm){
      firmName <- as.character(firm)
      delta_f <- timingParam[firmName]
      
      # If firm doesn't have timing param in vector, use mean
      if(is.na(delta_f)){
        delta_f <- mean(timingParam, na.rm = TRUE)
      }
      
      # IC_f: π^PL_f - π^D_f + δ_f/(1-δ_f) * (π^PL_f - π^B_f) ≥ 0
      slack_f <- (profitsPL_byFirm[firmName] - profitsDeviation_byFirm[firmName]) + 
        (delta_f / (1 - delta_f)) * (profitsPL_byFirm[firmName] - profitsBertrand_byFirm[firmName])
      
      return(slack_f)
    })
    
    # Check if all ICs hold
    if(all(IC_slack_byFirm >= 0)){
      # All ICs hold, try higher supermarkup
      m_low <- m_mid
    } else {
      # At least one IC violated, reduce supermarkup
      m_high <- m_mid
    }
  }
  
  return(m_low)
}
)
## calcSupermarkup method for PriceLeadershipBLP class
## Reuses PriceLeadership implementation - logic is identical
#' @rdname calcSupermarkup
#' @export
setMethod(
  f = "calcSupermarkup",
  signature = "PriceLeadershipBLP",
  definition = function(object, preMerger = TRUE, constrained = FALSE, ...){
    # Delegate to PriceLeadership implementation
    selectMethod("calcSupermarkup", "PriceLeadership")(object, preMerger, constrained, ...)
  }
)

## calcMargins method for PriceLeadership class
## Flexible markup calculation supporting 4 cases via 'regime' parameter:
## 1) "bertrand" - standard Bertrand (no coordination)
## 2) "deviation" - deviator best-responds, others at PLE prices  
## 3) "coordination" - coalition coordinates (unconstrained PLE)
## 4) "constrained" - coalition coordinates with IC binding (includes supermarkup)

#' @rdname Margins-Methods
#' @param regime Character string specifying pricing regime. Options:
#'   \itemize{
#'     \item "bertrand" - Standard Bertrand competition (no coordination)
#'     \item "deviation" - Deviation equilibrium (deviator optimizes, others fixed)
#'     \item "coordination" - Coordinated equilibrium (coalition coordinates, IC slack)
#'     \item "constrained" - Constrained coordination (IC binds, includes supermarkup)
#'   }
#'   Default is "coordination".
#' @param deviator For regime="deviation", numeric vector of coalition product indices
#'   that are deviating. Other coalition products' prices are held fixed.
#' @export
setMethod(
  f = "calcMargins",
  signature = "PriceLeadership",
  definition = function(object, preMerger = TRUE, level = FALSE, 
                        regime = c("coordination", "bertrand", "deviation", "constrained")){
    
    regime <- match.arg(regime)
    
    # Get coalition indices for supermarkup application
    if(preMerger){
      if(is.matrix(object@coalitionPre)){
        coalition <- which(rowSums(object@coalitionPre) > 1)
      } else {
        coalition <- object@coalitionPre
      }
    } else {
      if(is.matrix(object@coalitionPost)){
        coalition <- which(rowSums(object@coalitionPost) > 1)
      } else {
        coalition <- object@coalitionPost
      }
    }
    
    # Select ownership matrix based on regime
    # Note: ownerPre/Post and coalitionPre/Post are already populated in constructor
    if(preMerger){
      ownerMatrix <- switch(regime,
        # Case 1 & 2: Standard Bertrand or Deviation - use traditional ownership
        "bertrand" = object@ownerPre,
        "deviation" = object@ownerPre,
        
        # Case 3 & 4: Coordination (IC slack or binding) - use coalition control matrix
        "coordination" = object@coalitionPre,
        "constrained" = object@coalitionPre
      )
    } else {
      ownerMatrix <- switch(regime,
        # Case 1 & 2: Standard Bertrand or Deviation - use traditional ownership
        "bertrand" = object@ownerPost,
        "deviation" = object@ownerPost,
        
        # Case 3 & 4: Coordination (IC slack or binding) - use coalition control matrix
        "coordination" = object@coalitionPost,
        "constrained" = object@coalitionPost
      )
    }
    
    # Temporarily swap ownership matrix and call parent method
    if(preMerger){
      ownerOriginal <- object@ownerPre
      object@ownerPre <- ownerMatrix
      result <- callNextMethod(object, preMerger = preMerger, level = level)
      object@ownerPre <- ownerOriginal
    } else {
      ownerOriginal <- object@ownerPost
      object@ownerPost <- ownerMatrix
      result <- callNextMethod(object, preMerger = preMerger, level = level)
      object@ownerPost <- ownerOriginal
    }
    
    # For constrained regime: Add supermarkup to coalition products
    # This represents the binding IC constraint replacing the FOC
    if(regime == "constrained" && length(coalition) > 0){
      supermarkup <- if(preMerger) object@supermarkupPre else object@supermarkupPost
      
      if(level){
        # Level margins: add supermarkup directly
        result[coalition] <- result[coalition] + supermarkup
      } else {
        # Percentage margins: adjust proportionally
        prices <- if(preMerger) object@pricePre else object@pricePost
        result[coalition] <- result[coalition] + supermarkup / prices[coalition]
      }
    }
    
    return(result)
  }
)


## calcPrices method for PriceLeadership class
## Generic solver that finds prices satisfying: margin = calcMargins(object, regime)
## The regime parameter controls what equilibrium concept is used (via calcMargins)

#' @rdname Prices-Methods
#' @param regime Character string specifying pricing regime. Passed to calcMargins.
#'   Options: "bertrand", "coordination" (default), "deviation", "constrained"
#' @export
setMethod(
  f = "calcPrices",
  signature = "PriceLeadership",
  definition = function(object, preMerger = TRUE, isMax = FALSE, subset, 
                        regime = c("coordination", "bertrand", "deviation", "constrained"), ...){
    
    regime <- match.arg(regime)
    output <- object@output
    
    # Get marginal costs
    mc <- if(preMerger) object@mcPre else object@mcPost
    
    nprods <- length(object@shares)
    
    if(missing(subset)){
      subset <- rep(TRUE, nprods)
    }
    
    if(!is.logical(subset) || length(subset) != nprods){
      stop("'subset' must be a logical vector the same length as 'shares'")
    }
    
    mc_sub <- mc[subset]
    priceStart_sub <- object@priceStart[subset]
    
    # Get ownership matrix for Hessian check
    if(preMerger){
      owner <- if(regime == "bertrand") object@ownerPre else object@coalitionPre
    } else {
      owner <- if(regime == "bertrand") object@ownerPost else object@coalitionPost
    }
    owner_sub <- owner[subset, subset]
    
    priceEst <- rep(NA, nprods)
    
    ## Generic FOC: margin = calcMargins(object, regime)
    ## All regime-specific logic is in calcMargins
    FOC <- function(priceCand){
      if(preMerger){ 
        object@pricePre[subset] <- priceCand
      } else {
        object@pricePost[subset] <- priceCand
      }
      
      if(output){
        margins <- priceCand - mc_sub
      } else {
        margins <- mc_sub - priceCand
      }
      
      # Get predicted margins for this regime
      predMargin <- calcMargins(object, preMerger, level = TRUE, regime = regime)[subset]
      
      thisFOC <- margins - predMargin
      
      return(thisFOC)
    }
    
    ## Solve for equilibrium prices
    nleqslv_maxit <- as.integer(object@control.equ$maxit)
    if(length(nleqslv_maxit) == 0 || is.na(nleqslv_maxit[1]) || nleqslv_maxit[1] < 1){
      nleqslv_maxit <- 150L
    }
    
    minResult <- nleqslv::nleqslv(
      priceStart_sub, 
      FOC, 
      method = "Broyden",
      control = list(
        ftol = object@control.equ$tol,
        maxit = nleqslv_maxit
      )
    )
    
    ## Fallback to BBsolve if needed
    if(minResult$termcd > 2){
      minResult <- BBsolve(priceStart_sub, FOC, quiet = TRUE, 
                          control = object@control.equ, ...)
      priceEst_solution <- minResult$par
      if(minResult$convergence != 0){
        warning("'calcPrices' (PriceLeadership) nonlinear solver may not have converged. BBsolve reports: '",
                minResult$message, "'")
      }
    } else {
      priceEst_solution <- minResult$x
      if(minResult$termcd > 1){
        warning("'calcPrices' (PriceLeadership) may not have fully converged. nleqslv termcd: ",
                minResult$termcd)
      }
    }
    
    if(isMax && length(fringe_sub) > 0){
      hess <- genD(FOC, priceEst_solution)$D[, 1:length(priceEst_solution)]
      hess <- hess * (owner_sub > 0)
      
      state <- ifelse(preMerger, "Pre-merger", "Post-merger")
      
      if(any(eigen(hess)$values > 0)){
        warning("Hessian not positive definite. ", state, 
                " price vector may not maximize profits.")
      }
    }
    
    priceEst[subset] <- priceEst_solution
    names(priceEst) <- object@labels
    
    return(priceEst)
  }
)


#' Calculate Incentive Compatibility Slack Functions
#' 
#' @description Calculates the slack function g_f(m) for each coalition firm,
#' which measures the incentive compatibility constraint. A firm's IC constraint
#' is satisfied when g_f(m) >= 0.
#' 
#' @param object A PriceLeadership object
#' @param preMerger Logical, if TRUE use pre-merger structure
#' @param ... Additional arguments passed to solver
#' 
#' @return Named vector of slack values for each coalition firm
#' 
#' @details The slack function is:
#' g_f(m) = [π^PL_f - π^D_f] + δ/(1-δ) * [π^PL_f - π^B_f]
#' 
#' where:
#' \itemize{
#'   \item π^PL_f: Profit under price leadership
#'   \item π^D_f: Profit from optimal deviation
#'   \item π^B_f: Profit under Bertrand punishment
#'   \item δ: Timing/discount parameter
#' }
#' 
#' @export
calcSlack <- function(object, preMerger = TRUE, ...){
  
  supermarkup <- if(preMerger) object@supermarkupPre else object@supermarkupPost
  timingParam <- object@timingParam
  nprods <- length(object@prices)
  
  if(preMerger){
    # Extract coalition indices (handle both vector and matrix forms)
    if(is.matrix(object@coalitionPre)){
      coalition <- which(rowSums(object@coalitionPre) > 1)
    } else {
      coalition <- object@coalitionPre
    }
    owner <- object@ownerPre
  } else {
    # Extract coalition indices (handle both vector and matrix forms)
    if(is.matrix(object@coalitionPost)){
      coalition <- which(rowSums(object@coalitionPost) > 1)
    } else {
      coalition <- object@coalitionPost
    }
    owner <- object@ownerPost
  }
  
  # Get firm identifiers for coalition products
  if(is.matrix(owner)){
    ownerVec <- rownames(owner)
    if(is.null(ownerVec)){
      ownerVec <- 1:nrow(owner)
    }
  } else {
    ownerVec <- owner
  }
  
  # Unique firms in coalition
  coalitionOwners <- unique(ownerVec[coalition])
  
  # Storage for slack values
  slackValues <- rep(NA_real_, length(coalitionOwners))
  names(slackValues) <- coalitionOwners
  
  ## Step 1: Calculate price leadership profits (π^PL) at current prices
  plProfits <- calcProducerSurplus(object, preMerger = preMerger)
  
  # Aggregate by firm
  plProfitsByFirm <- sapply(coalitionOwners, function(firm){
    if(is.character(firm)){
      sum(plProfits[ownerVec == firm])
    } else {
      sum(plProfits[owner[firm, ] > 0])
    }
  })
  
  ## Step 2: Calculate Bertrand profits (π^B) - punishment
  # Save original owner matrix and current prices
  ownerOriginal <- if(preMerger) object@ownerPre else object@ownerPost
  currentPrices <- if(preMerger) object@pricePre else object@pricePost
  
  # Calculate Bertrand prices (no coordination) - swap owner matrix
  if(preMerger){
    object@ownerPre <- ownerToMatrix(object, preMerger = TRUE, control = FALSE)
  } else {
    object@ownerPost <- ownerToMatrix(object, preMerger = FALSE, control = FALSE)
  }
  bertrandPrices <- calcPrices(object, preMerger = preMerger, isMax = FALSE, ...)
  
  # Set prices to Bertrand and calculate profits
  if(preMerger){
    object@pricePre <- bertrandPrices
  } else {
    object@pricePost <- bertrandPrices
  }
  
  bertrandProfits <- calcProducerSurplus(object, preMerger = preMerger)
  
  bertrandProfitsByFirm <- sapply(coalitionOwners, function(firm){
    if(is.character(firm)){
      sum(bertrandProfits[ownerVec == firm])
    } else {
      sum(bertrandProfits[owner[firm, ] > 0])
    }
  })
  
  ## Step 3: Calculate deviation profits (π^D) for each coalition firm
  ## Owner matrix is already set to no-coordination from Step 2
  
  # Reset to current (PL) prices and calculate what each firm would earn if they deviated
  if(preMerger){
    object@pricePre <- currentPrices
  } else {
    object@pricePost <- currentPrices
  }
  
  deviationProfits <- calcProducerSurplus(object, preMerger = preMerger)
  
  deviationProfitsByFirm <- sapply(coalitionOwners, function(firm){
    if(is.character(firm)){
      sum(deviationProfits[ownerVec == firm])
    } else {
      sum(deviationProfits[owner[firm, ] > 0])
    }
  })
  
  # Restore original owner matrix
  if(preMerger){
    object@ownerPre <- ownerOriginal
  } else {
    object@ownerPost <- ownerOriginal
  }
  
  ## Step 4: Calculate slack function
  ## g_f(m) = [π^PL_f - π^D_f] + δ/(1-δ) * [π^PL_f - π^B_f]
  
  if(is.na(timingParam)){
    # Cannot calculate slack without timing parameter
    # Return immediate payoff difference only
    slackValues <- plProfitsByFirm - deviationProfitsByFirm
    warning("Timing parameter not specified. Slack values show only immediate payoff difference [π^PL - π^D].")
  } else {
    # Full slack calculation
    slackValues <- (plProfitsByFirm - deviationProfitsByFirm) + 
                   (timingParam / (1 - timingParam)) * (plProfitsByFirm - bertrandProfitsByFirm)
  }
  
  return(slackValues)
}



