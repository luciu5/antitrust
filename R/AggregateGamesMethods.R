#' @title Methods for Calculating Prices and Margins Using Aggregative Games
#' @name MarginsAG-Methods
#' @docType methods
#'
#' @aliases calcMarginsAG calcPricesAG
#' calcMarginsAG,ANY-method
#' calcPricesAG,ANY-method
#' calcMarginsAG,Logit-method
#' calcMarginsAG,CES-method
#' calcMarginsAG,BargainingLogit-method
#' calcMarginsAG,BargainingCES-method
#' calcMarginsAG,Cournot-method
#' calcPricesAG,Logit-method
#' calcPricesAG,CES-method
#' calcPricesAG,BargainingLogit-method
#' calcPricesAG,BargainingCES-method
#' calcMarginsAG,Auction2ndLogit-method
#' calcPricesAG,Auction2ndLogit-method
#' calcMarginsAG,Auction2ndCES-method
#' calcPricesAG,Auction2ndCES-method
#' calcMarginsAG,Bargaining2ndLogit-method
#' calcPricesAG,Bargaining2ndLogit-method
#' calcMarginsAG,Bargaining2ndCES-method
#' calcPricesAG,Bargaining2ndCES-method
#' calcPricesAG,Cournot-method
#'
#' @description Computes equilibrium product margins and prices using the aggregative games technique described in
#' Nocke and Schutz (2018). Assumes that firms are playing a
#' Nash-Bertrand pricing game with either Logit or CES demand; bargaining and
#' second-score wrappers reuse the corresponding model-specific first-order
#' conditions after obtaining an AG-compatible starting point.
#' @param object An instance of one of the classes listed above.
#' @param preMerger If TRUE, returns pre-merger outcome. If
#' FALSE, returns post-merger outcome.  Default is TRUE.
#' @param level IF TRUE, return margins in dollars. If FALSE, returns 
#' margins in proportions. Default for most classes is FALSE.
#' @param subset A vector whose length equals the number  of products where each element equals TRUE if
#' the product indexed by that element should be included in the
#' post-merger simulation and FALSE if it should be excluded. Default is a
#' length k vector of TRUE.
#' @param isMax If TRUE, a check is run to determine if the calculated equilibrium price vector locally maximizes profits.
#' Default is FALSE.
#' @include MarginsMethods.R
#' @keywords methods
#' @references 
#' Nocke, V. and Schutz, N. (2018), Multiproduct-Firm Oligopoly: An Aggregative Games Approach. Econometrica, 86: 523-557.\doi{10.3982/ECTA14720}/
NULL

setGeneric (
  name= "calcMarginsAG",
  def=function(object,...){standardGeneric("calcMarginsAG")}
)

setGeneric (
  name= "calcPricesAG",
  def=function(object,...){standardGeneric("calcPricesAG")}
)

.checkAGMax <- function(object, prices, preMerger, subset) {
  if (length(prices) == 0 || any(!is.finite(prices))) return(invisible(NULL))

  owner <- if (preMerger) object@ownerPre else object@ownerPost
  mc <- if (preMerger) object@mcPre else object@mcPost
  owner <- owner[subset, subset, drop = FALSE]
  mc <- mc[subset]
  outSign <- ifelse(object@output, 1, -1)

  foc <- function(priceCand) {
    old <- if (preMerger) object@pricePre else object@pricePost
    old[subset] <- priceCand
    if (preMerger) object@pricePre <- old else object@pricePost <- old
    margin <- calcMargins(object, preMerger = preMerger, level = TRUE)[subset]
    (priceCand - mc) - outSign * margin
  }

  hess <- try(numDeriv::genD(foc, prices[subset]), silent = TRUE)
  if (inherits(hess, "try-error")) return(invisible(NULL))
  hess <- hess$D[, 1:hess$p]
  hess <- hess * (owner > 0)
  if (any(eigen(hess, symmetric = TRUE, only.values = TRUE)$values > 0)) {
    state <- ifelse(preMerger, "Pre-merger", "Post-merger")
    warning("Hessian of first-order conditions is not positive definite. ",
      state, " price vector may not maximize profits. Consider rerunning the solver using different starting values")
  }
  invisible(NULL)
}

## compute margins using aggregate games method
## 
#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcMarginsAG",
  signature= "Logit",
  definition=function(object,preMerger=TRUE, level=FALSE, subset){
    
    output <- ifelse(object@output,-1,1)
    meanval <- object@slopes$meanval
    alpha <- object@slopes$alpha
    idx <- object@normIndex
    
    outPrice <- object@priceOutside
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }
    
    if( preMerger) {
      mc <- object@mcPre
      owner  <- object@ownerPre
    }
    
    else{
      mc <- object@mcPost
      owner  <- object@ownerPost
    }

    meanval <- meanval[subset]
    mc <- mc[subset]
    owner <- owner[subset, subset, drop = FALSE]
    
    ## Use log-sum-exp rescaling to prevent exp() overflow
    ## when alpha > 0 (input markets) and mc is large
    log_type <- meanval + alpha*mc
    log_H0 <- outPrice*alpha
    
    if (is.na(idx)) {
      max_log <- max(c(log_type, log_H0))
    } else {
      max_log <- max(log_type)
    }
    
    type <- exp(log_type - max_log)
    H0 <- is.na(idx) * exp(log_H0 - max_log)
    
    typeFirm <- as.numeric(owner %*% type)
    
    solve_mu_logit_vec <- function(typeFirm, H, max_iter = 8) {
      A <- typeFirm / H
      m <- pmax(1.0, log(pmax(A, 1.001)))
      for (i in 1:max_iter) {
        exp_neg_m <- exp(-m)
        f <- m * (1 - A * exp_neg_m) - 1
        f_prime <- 1 + (m - 1) * A * exp_neg_m
        m <- pmax(1.0, m - f / f_prime)
      }
      return(m)
    }
    
    Hfun <- function(h){
      mu <- solve_mu_logit_vec(typeFirm, h)
      price <- mc - mu/alpha
      val <- exp(meanval + alpha*price - max_log)
      omega <- sum(H0,val)/h
      return((omega - 1)^2)
    }
    
    HStart <- H0+sum(type)/exp(1)
    
    if (!is.finite(HStart) || HStart <= 0) {
      stop("calcMarginsAG (Logit): non-finite starting value for aggregative games optimizer. ",
           "This typically occurs in input markets where alpha > 0 causes exp(alpha*mc) to overflow.")
    }
    
    HBest <- optim(HStart,Hfun,method="Brent",lower=0,upper=1e6*HStart)
    
    margins <- solve_mu_logit_vec(typeFirm, HBest$par)
    
    margins <- output*margins/alpha
    
    if(!level) {
      price <- mc - output * margins
      margins <- margins / price }
    
    result <- rep(NA_real_, nprods)
    result[subset] <- as.numeric(margins)
    names(result) <- object@labels
    return(result)
  }
  
)

#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcMarginsAG",
  signature= "CES",
  definition=function(object,preMerger=TRUE, level=FALSE, subset){
    
    output <- ifelse(object@output,1,-1)
    meanval <- object@slopes$meanval
    gamma <- object@slopes$gamma
    idx <- object@normIndex
    
    outPrice <- object@priceOutside
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }

    ## The closed-form CES aggregation is numerically well behaved for output
    ## markets (gamma > 1).  For input markets, use the package's elasticity
    ## implementation, which handles the reversed price/demand sign directly.
    ## This keeps the AG API accurate while avoiding fractional-power branches
    ## that can produce spurious markups when gamma < 1.
    if (!object@output) {
      margins <- calcMargins(object, preMerger = preMerger, level = level)
      margins[!subset] <- NA_real_
      return(margins)
    }
    
    H0 <- is.na(idx)*(outPrice^(1-gamma))
    
    if( preMerger) {
      mc <- object@mcPre
      owner  <- object@ownerPre
    }
    
    else{
      mc <- object@mcPost
      owner  <- object@ownerPost
    }

    meanval <- meanval[subset]
    mc <- mc[subset]
    owner <- owner[subset, subset, drop = FALSE]
    
    
    type <-  meanval * mc^(1-gamma)
    
    typeFirm <- as.numeric(owner %*% type)
    
    mufun <- function(m,H){
      
      return(m*(1 - (gamma - 1)/gamma*(typeFirm/H)*abs(1 - m/gamma)^(gamma - 1)) - 1)
    }
    
    Hfun <- function(h){
      
      muStart <- pmax(1,gamma*(1 - (typeFirm/h)^(1/(gamma - 1))))
      
      mu <- BB::BBsolve(muStart,mufun,H=h,quiet=TRUE)
      
      price <-  mc/(1 - mu$par/gamma)
      
      val <- meanval * price^(1-gamma) 
      
      omega <- sum(H0,val)/h
      
      return((omega - 1)^2)
    }
    
    ## Use abs() to handle gamma < 1 where (1-1/gamma) is negative
    ## and raising to fractional power (gamma-1) would produce NaN
    HStart <- H0+sum(type)*abs(1 - 1/gamma)^(gamma - 1)
    
    if (!is.finite(HStart) || HStart <= 0) {
      stop("calcMarginsAG (CES): non-finite starting value for aggregative games optimizer. ",
           "This typically occurs in input markets where gamma < 1 causes (1-1/gamma)^(gamma-1) to overflow.")
    }
    
    HBest <- optim(HStart,Hfun,method="Brent",lower=0,upper=1e6*HStart)
    
    
    muStart <- pmax(1,gamma*(1 - (typeFirm/HBest$par)^(1/(gamma - 1)))) 
    margins <- BB::BBsolve(muStart,mufun,H=HBest$par,quiet=TRUE)$par
    
    margins <- abs(margins/gamma)
    
    if(level) {
      price <- mc/(1 - output*margins)
      margins <- margins * price }
    
    result <- rep(NA_real_, nprods)
    result[subset] <- as.numeric(margins)
    names(result) <- object@labels
    return(result)
  }
  
)


## compute prices using aggregate games method
## 
#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcPricesAG",
  signature= "Logit",
  definition=function(object,preMerger=TRUE,isMax=FALSE,subset){
    
    
    priceStart <- object@priceStart
    output    <-  object@output
    outSign <- ifelse(output,1,-1)
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }
    
    if(preMerger){
      owner <- object@ownerPre
      mc    <- object@mcPre
    }
    else{
      owner <- object@ownerPost
      mc    <- object@mcPost
    }
    
  
      margins <- calcMarginsAG(object,preMerger=preMerger,level=TRUE,subset=subset)
      priceEst <- rep(NA_real_, nprods)
      priceEst[subset] <- outSign*margins[subset] + mc[subset]
      names(priceEst) <- object@labels
      if (isMax) .checkAGMax(object, priceEst, preMerger, subset)
      return(priceEst)
    }
)


#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcPricesAG",
  signature= "CES",
  definition=function(object,preMerger=TRUE,isMax=FALSE,subset){
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }
    ## Use the standard CES FOC solver for input markets.  The aggregation
    ## formula above is an output-market construction; the standard solver is
    ## the sign-safe implementation for gamma < 1.
    if (!object@output) {
      return(calcPrices(object, preMerger = preMerger, isMax = isMax, subset = subset))
    }
    output <- ifelse(object@output, 1, -1)
    mc <- if (preMerger) object@mcPre else object@mcPost
    margins <- calcMarginsAG(object, preMerger = preMerger, level = TRUE, subset = subset)
    prices <- rep(NA_real_, nprods)
    prices[subset] <- mc[subset] + output * margins[subset]
    names(prices) <- object@labels
    if (isMax) .checkAGMax(object, prices, preMerger, subset)
    return(prices)
  }
)


#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcMarginsAG",
  signature= "BargainingLogit",
  definition=function(object,preMerger=TRUE, level=FALSE, subset){
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    barg <- if (preMerger) object@bargpowerPre else object@bargpowerPost
    margins <- methods::selectMethod("calcMarginsAG", "Logit")(
      object, preMerger = preMerger, level = level, subset = subset
    )
    margins[subset] <- (1 - barg[subset]) * margins[subset]
    names(margins) <- object@labels
    margins
  }
)

#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcPricesAG",
  signature= "BargainingLogit",
  definition=function(object,preMerger=TRUE,isMax=FALSE,subset){
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }
    ## Use the AG solution as a starting point, then enforce the bargaining
    ## first-order conditions with the established bargaining price solver.
    ## Bargaining changes the firm objective, so the plain Logit AG markup
    ## equation is not itself the bargaining equilibrium condition.
    agStart <- try(suppressWarnings(
      methods::selectMethod("calcPricesAG", "Logit")(
        object, preMerger = preMerger, isMax = FALSE, subset = subset
      )
    ), silent = TRUE)
    if (!inherits(agStart, "try-error") && length(agStart) == nprods &&
      all(is.finite(agStart[subset]))) {
      object@priceStart[subset] <- agStart[subset]
    }
    methods::selectMethod("calcPrices", "BargainingLogit")(
      object, preMerger = preMerger, isMax = isMax, subset = subset
    )
  }
)


#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcMarginsAG",
  signature= "BargainingCES",
  definition=function(object,preMerger=TRUE, level=FALSE, subset){
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    barg <- if (preMerger) object@bargpowerPre else object@bargpowerPost
    margins <- methods::selectMethod("calcMarginsAG", "CES")(
      object, preMerger = preMerger, level = level, subset = subset
    )
    margins[subset] <- (1 - barg[subset]) * margins[subset]
    names(margins) <- object@labels
    margins
  }
)

#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcPricesAG",
  signature= "BargainingCES",
  definition=function(object,preMerger=TRUE,isMax=FALSE,subset){
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }
    ## Use the CES AG solution to initialize the full bargaining FOCs.  The
    ## bargaining-specific objective is not the plain CES AG markup equation.
    agStart <- try(suppressWarnings(
      methods::selectMethod("calcPricesAG", "CES")(
        object, preMerger = preMerger, isMax = FALSE, subset = subset
      )
    ), silent = TRUE)
    if (!inherits(agStart, "try-error") && length(agStart) == nprods &&
      all(is.finite(agStart[subset]))) {
      object@priceStart[subset] <- agStart[subset]
    }
    methods::selectMethod("calcPrices", "Logit")(
      object, preMerger = preMerger, isMax = isMax, subset = subset
    )
  }
)

#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcMarginsAG",
  signature= "Auction2ndLogit",
  definition=function(object,preMerger=TRUE, level=FALSE, subset){
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }
    margins <- calcMargins(object, preMerger = preMerger, level = level)
    margins[!subset] <- NA_real_
    margins
  }
)

#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcPricesAG",
  signature= "Auction2ndLogit",
  definition=function(object,preMerger=TRUE,isMax=FALSE,subset){
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }
    if (!preMerger) object@subset <- subset
    prices <- calcPrices(object, preMerger = preMerger)
    prices[!subset] <- NA_real_
    prices
  }
)

#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcMarginsAG",
  signature= "Auction2ndCES",
  definition=function(object,preMerger=TRUE, level=FALSE, subset){
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }
    margins <- calcMargins(object, preMerger = preMerger, level = level)
    margins[!subset] <- NA_real_
    margins
  }
)

#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcPricesAG",
  signature= "Auction2ndCES",
  definition=function(object,preMerger=TRUE,isMax=FALSE,subset){
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }
    if (!preMerger) object@subset <- subset
    prices <- calcPrices(object, preMerger = preMerger)
    prices[!subset] <- NA_real_
    prices
  }
)

#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcMarginsAG",
  signature= "Bargaining2ndLogit",
  definition=function(object,preMerger=TRUE, level=FALSE, subset){
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }
    margins <- calcMargins(object, preMerger = preMerger, level = level)
    margins[!subset] <- NA_real_
    margins
  }
)

#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcPricesAG",
  signature= "Bargaining2ndLogit",
  definition=function(object,preMerger=TRUE,isMax=FALSE,subset){
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }
    if (!preMerger) object@subset <- subset
    prices <- calcPrices(object, preMerger = preMerger)
    prices[!subset] <- NA_real_
    prices
  }
)

#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcMarginsAG",
  signature= "Bargaining2ndCES",
  definition=function(object,preMerger=TRUE, level=FALSE, subset){
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }
    margins <- calcMargins(object, preMerger = preMerger, level = level)
    margins[!subset] <- NA_real_
    margins
  }
)

#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcPricesAG",
  signature= "Bargaining2ndCES",
  definition=function(object,preMerger=TRUE,isMax=FALSE,subset){
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }
    if (!preMerger) object@subset <- subset
    prices <- calcPrices(object, preMerger = preMerger)
    prices[!subset] <- NA_real_
    prices
  }
)


#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcMarginsAG",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE, level=FALSE, subset){
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }
    margins <- calcMargins(object, preMerger = preMerger, level = level)
    margins[!subset] <- NA_real_
    margins
  }
)

#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcPricesAG",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE,isMax=FALSE,subset){
    nprods <- length(object@shares)
    if (missing(subset)) subset <- if (preMerger) rep(TRUE, nprods) else object@subset
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
      stop("'subset' must be a logical vector the same length as 'shares' with at least one TRUE value")
    }
    prices <- calcPrices(object, preMerger = preMerger)
    prices[!subset] <- NA_real_
    prices
  }
)
