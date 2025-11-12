#'@title \dQuote{Bertrand RUM} Classes
#'@name BertrandRUM-Classes
#'@aliases Logit-class PriceLeadership-class

#'@description Each class contains all the information needed to calibrate a specific type of demand system and
#'perform a merger simulation analysis under the assumption that firms are playing a differentiated products Bertrand pricing game. Also included is \dQuote{LogitCournot}, for modelling firms playing a differentiated products Cournot quantity game.
#'
#'@description The \dQuote{Logit} class has the information for a Logit demand system.
#'@description The \dQuote{PriceLeadership} class extends \dQuote{Logit} to model coordinated effects
#'using the price leadership framework of Mansley, Miller, Sheu \& Weinberg (2023). In this model,
#'a leader firm announces a supermarkup above Bertrand prices, coalition firms follow the leader,
#'and fringe firms best-respond. The model calibrates the supermarkup and timing parameter from
#'observed margins and tests incentive compatibility constraints.
#'@description The \dQuote{LogitCap} class has the information for a Logit demand system and assumes that
#'firms are playing a differentiated products Bertrand pricing game with capacity constraints.
#'\dQuote{LogitCapALM} extends \dQuote{LogitCap} to allow for an unobserved outside share.
#'@description The \dQuote{LogitNests} class has the information for a nested Logit
#'@description The \dQuote{LogitNestsALM} class has the information for a nested Logit
#'demand system under the assumption that the share of the outside product is not known.
#'Once the model parameters have been calibrated, methods exist that perform a merger simulation analysis under the assumption that
#'firms are playing a differentiated products Bertrand pricing game.
#'@description The \dQuote{LogitALM} class has the information for a Logit demand system
#'assuming that firms are playing a differentiated products Bertrand pricing game with unknown market elasticity.
#'@description The \dQuote{CES} class has the information for a CES demand system
#'@description The \dQuote{CESALM} class has the information for a CES demand system and
#'assumes that firms are playing a differentiated products Bertrand pricing game with unknown market elasticity.
#'@description  The \dQuote{CESNests} class has the information for a nested CES demand system.
#'@description Let k denote the number of products produced by all firms below.
#'
#'@section Objects from the Class:
#'For Logit, objects can be created by using the constructor function \code{\link{logit}}.
#'
#'For PriceLeadership, objects can be created by using the constructor function \code{\link{plm}}.
#'
#'For LogitALM, objects can be created by using the constructor function \code{\link{logit.alm}}.
#'
#'For LogitCap and LogitCapALM, objects can be created by using the constructor function \code{\link{logit.cap}} and \code{\link{logit.cap.alm}}.
#'
#'For LogitNests, objects can be created by using the constructor function \code{\link{logit.nests}}.
#'
#'For LogitNestsALM, objects can be created by using the constructor function \code{\link{logit.nests.alm}}.
#'
#'For CES, objects can be created by using the constructor function \code{\link{ces}}.
#'
#'For CESALM, objects can be created by using the constructor function \code{\link{ces.alm}}.
#'
#'For CESNests, objects can be created by using the constructor function \code{\link{ces.nests}}.
#'
#'@slot prices A length k vector of product prices.
#'@slot margins A length k vector of product margins, some of which may equal NA.
#'@slot normIndex An integer specifying the product index against which the mean values of all other products are normalized.
#'@slot shareInside The share of customers that purchase any of the products included in the `prices' vector.
#'@slot priceOutside The price of the outside good. Default is 0.
#'@slot slopes A list containing the coefficient on price (\sQuote{alpha}) and the vector of mean valuations (\sQuote{meanval}).
#'@slot mktElast A length 1 vector of market elasticities.
#'@slot priceStart A length-k vector of starting prices for the non-linear solver.
#'@slot insideSize A positive number equal to total pre-merger quantities (revenues for CES) for all products included in the simulation.
#'@slot mktSize A positive number equal to total quantities (revenues for CES) pre-merger for all products in the simulations
#'as well as the outside good.
#'@slot capacitiesPre A length k vector whose elements equal pre-merger product capacities. (LogitCap and LogitCapALM only)
#'@slot capacitiesPost A length k vector whose elements equal post-merger product capacities. (LogitCap and LogitCapALM only)
#'@slot nests A length k vector identifying the nest that each product belongs to. (LogitNests and CESNests Only)
#'@slot parmsStart A length k vector who elements equal an initial guess of the nesting parameter values. (LogitNests and CESNests Only)
#'@slot constraint A length 1 logical vector that equals TRUE if all nesting parameters are constrained to equal the same value
#'and FALSE otherwise. Default is TRUE. (LogitNests and CESNests Only)
#'@slot parmsStart A length 2 vector whose first element equals an initial guess of the price coefficient and whose second
#'element equals an initial guess of the outside share. The price
#'coefficient's initial value must be negative and the outside share's initial value must be between 0 and 1. (LogitALM and CESALM only)
#'@slot slopes A list containing the coefficient on the numeraire (`alpha'),  the coefficient on price (\sQuote{gamma}), and the vector of mean
#'valuations (\sQuote{meanval}) (CES only)
#'@slot priceOutside The price of the outside good. Default is 1. (CES only)

#'@section Extends:
#'Logit: Class \code{\linkS4class{Bertrand}}, directly.
#'Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 2.
#'
#'LogitCap: Class \code{\linkS4class{Logit}}, directly.
#'Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Logit}}, distance 2.
#'Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 3.
#'
#'#'LogitCapALM: Class \code{\linkS4class{LogitCap}}, directly.
#'Class \code{\linkS4class{Logit}}, by class \code{\linkS4class{LogitCap}}, distance 2.
#'Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Logit}}, distance 3.
#'Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 4.
#'
#'LogitNests: Class \code{\linkS4class{Logit}}, directly.
#'Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Logit}}, distance 2.
#'
#'LogitNestsALM: Class \code{\linkS4class{LogitNests}}, directly.
#'Class \code{\linkS4class{Logit}}, by class \code{\linkS4class{LogitNests}}, distance 2.
#'Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Logit}}, distance 3.
#'Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 4.
#'
#'LogitALM: Class \code{\linkS4class{Logit}}, directly.
#'Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Logit}}, distance 2.
#'Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 3.
#'
#'CES: Class \code{\linkS4class{Logit}}, directly.
#'Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Logit}}, distance 2.
#'Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 3.
#'
#'CESALM: Class \code{\linkS4class{CES}}, directly.
#'Class \code{\linkS4class{Logit}}, by class \code{\linkS4class{CES}}, distance 2.
#'Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Logit}}, distance 3.
#'Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 4.
#'
#'CESNests: Class \code{\linkS4class{CES}}, directly.
#'Class \code{\linkS4class{Logit}}, by class \code{\linkS4class{CES}}, distance 2.
#'Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Logit}}, distance 3.
#'Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 4.

#'@author Charles Taragin \email{ctaragin+antitrustr@gmail.com}
#'@examples
#'showClass("Logit")           # get a detailed description of the class
#'showClass("LogitCap")           # get a detailed description of the class
#'showClass("LogitNests")           # get a detailed description of the class
#'showClass("LogitNestsALM")           # get a detailed description of the class
#'showClass("LogitALM")           # get a detailed description of the class
#'showClass("CES")           # get a detailed description of the class
#'showClass("CESALM")           # get a detailed description of the class
#'showClass("CESNests")           # get a detailed description of the class
#'@include BertrandClasses.R

#'@rdname BertrandRUM-Classes
#'@export
setClass(
  Class   = "Logit",
  contains="Bertrand",
  representation=representation(

    prices           = "numeric",
    margins          = "numeric",
    priceStart       = "numeric",
    normIndex        = "vector",
    shareInside      = "numeric",
    priceOutside     = "numeric",
    mktElast         = "numeric",
    insideSize          = "numeric",
    mktSize             = "numeric"

  ),
  prototype=prototype(
    mktElast = NA_real_,
    insideSize  = NA_real_,
    mktSize = 1,
    priceStart  = numeric(),
    normIndex   = 1,
    shareInside = numeric(),
    priceOutside = 0,
    control.slopes = list(
      factr = 1e7
    )
  ),
  validity=function(object){



    margins <- object@margins
    nMargins <- length(margins[!is.na(margins)])

    nprods <- length(object@shares)



    if(
      nprods != length(margins) ||
      nprods != length(object@prices)){
      stop("'prices', 'margins' and 'shares' must all be vectors with the same length")}

    if(isTRUE(any(object@prices<0,na.rm=TRUE)))             stop("'prices' values must be positive")


    if(isTRUE(any(margins<0,na.rm=TRUE))) stop("'margins' values must be positive")

    if(nMargins == 0) stop("At least one margin must be supplied.")

    if(!(is.matrix(object@ownerPre))){
      ownerPre <- ownerToMatrix(object,TRUE)
    }
    else{ownerPre <- object@ownerPre}


    #isMargin    <- matrix(margins,nrow=nprods,ncol=nprods,byrow=TRUE)
    #isMargin[ownerPre==0]=0
    #isMargin    <- !is.na(rowSums(isMargin))

    #if(object@cls != "Auction2ndLogit" &&
    #    !any(isMargin)){ stop("Insufficient margin information to calibrate demand parameters.")}

    if(nprods != length(object@priceStart)){
      stop("'priceStart' must have the same length as 'shares'")}

    if(
      !(object@shareInside >=0 &&
        object@shareInside <=1) #||
      #!isTRUE(all.equal(object@shareInside,1,check.names=FALSE, tolerance=1e-3))
    ){
      stop("'shareInside' must be between 0 and 1")
    }

    if(
      !(all(object@shares >0) &&
        all(object@shares <=1))
    ){
      stop("elements of vector 'shares' must be between 0 and 1")
    }

    if(!(length(object@normIndex) == 1 &&
         object@normIndex %in% c(NA,1:nprods))){
      stop("'normIndex' must take on a value between 1 and ",nprods,
           " or NA")
    }

    if(length(object@priceOutside) != 1 || object@priceOutside<0
    ){stop("'priceOutside' must be a non-negative number")}

    if(!is.na(object@mktElast) && object@mktElast >0 ) stop("'mktElast' must be negative")
    #if(!is.na(object@mktElast) && !isTRUE(all.equal(sum(object@shares, na.rm=TRUE),1)) ) stop("`shares' must sum to 1 when 'mktElast' is supplied")

    if(#length(object@mktSize)!=1 ||
       any(is.na(object@mktSize) | object@mktSize<0)){
      stop("mktSize must be a positive number")}
    return(TRUE)

  })

#'@rdname BertrandRUM-Classes
#'@export
setClass(
  Class   = "LogitBLP",
  contains = "Logit",
  slots = list(
    nDraws = "numeric"
  ),
  prototype = prototype(
    nDraws = 1000
  ),
  validity = function(object){
    # nDraws must be a single positive numeric value
    if(length(object@nDraws) != 1 || is.na(object@nDraws) || !is.finite(object@nDraws) || object@nDraws <= 0){
      stop("'nDraws' must be a single positive number")
    }

    # If present, check slopes list types without requiring them at construction time
    if(!is.null(object@slopes)){
      if(!is.list(object@slopes)) stop("'slopes' must be a list when provided")
      if(!is.null(object@slopes$alpha) && !is.numeric(object@slopes$alpha)){
        stop("'slopes$alpha' must be numeric when provided")
      }
      if(!is.null(object@slopes$alphaMean) && !is.numeric(object@slopes$alphaMean)){
        stop("'slopes$alphaMean' must be numeric when provided")
      }
      if(!is.null(object@slopes$sigma) && (!is.numeric(object@slopes$sigma) || length(object@slopes$sigma)!=1)){
        stop("'slopes$sigma' must be a single numeric value when provided")
      }
      if(!is.null(object@slopes$piDemog) && !is.numeric(object@slopes$piDemog)){
        stop("'slopes$piDemog' must be numeric when provided")
      }
      if(!is.null(object@slopes$nDemog) && (length(object@slopes$nDemog)!=1 || is.na(object@slopes$nDemog) || object@slopes$nDemog < 0)){
        stop("'slopes$nDemog' must be a single non-negative number when provided")
      }
      if(!is.null(object@slopes$sigmaNest) && (!is.numeric(object@slopes$sigmaNest) || object@slopes$sigmaNest <= 0 || object@slopes$sigmaNest > 1)){
        stop("'slopes$sigmaNest' (nesting parameter) must be in (0,1] when provided")
      }
    }

    return(TRUE)
  }
)

#'@rdname BertrandRUM-Classes
#'@export
setClass(
  Class   = "PriceLeadership",
  contains = "Logit",
  slots = list(
    supermarkupPre = "numeric",   # Pre-merger supermarkup m above Bertrand prices
    supermarkupPost = "numeric",  # Post-merger supermarkup m above Bertrand prices
    timingParam = "numeric",      # Firm-specific timing/discount parameters delta in (0,1) - named vector
    coalitionPre = "ANY",         # Pre-merger coalition: vector of indices OR control matrix
    coalitionPost = "ANY",        # Post-merger coalition: vector of indices OR control matrix
    bindingFirm = "numeric",      # Index of firm with binding IC constraint (or NA if none)
    slackValues = "numeric"       # Slack function values for each coalition firm
  ),
  prototype = prototype(
    supermarkupPre = NA_real_,
    supermarkupPost = NA_real_,
    timingParam = numeric(),      # Empty named vector initially
    coalitionPre = numeric(),
    coalitionPost = numeric(),
    bindingFirm = NA_real_,
    slackValues = numeric()
  ),
  validity = function(object){
    
    # Call parent validity first
    # Parent (Logit) already checks: prices, margins, shares, normIndex, etc.
    
    nprods <- length(object@prices)
    
    # PriceLeadership-specific validations only
    
    # Pre-merger coalition (can be vector or matrix)
    if(length(object@coalitionPre) > 0){
      if(is.matrix(object@coalitionPre)){
        # Matrix form - check dimensions
        if(nrow(object@coalitionPre) != nprods || ncol(object@coalitionPre) != nprods){
          stop("'coalitionPre' matrix must be ", nprods, " x ", nprods)
        }
        # Extract coalition indices to validate
        coalitionIndices <- which(rowSums(object@coalitionPre) > 1)
      } else {
        # Vector form - validate indices
        coalitionIndices <- object@coalitionPre
        if(any(!(coalitionIndices %in% 1:nprods))){
          stop("'coalitionPre' must contain product indices between 1 and ", nprods)
        }
      }
      
      # Coalition must have at least 2 products
      if(length(coalitionIndices) < 2){
        stop("'coalitionPre' must contain at least 2 products for price leadership")
      }
      
      # Must have at least one fringe product
      if(length(coalitionIndices) >= nprods){
        stop("At least one fringe product is required for price leadership model")
      }
    }
    
    # Post-merger coalition (can be vector or matrix)
    if(length(object@coalitionPost) > 0){
      if(is.matrix(object@coalitionPost)){
        # Matrix form - check dimensions
        if(nrow(object@coalitionPost) != nprods || ncol(object@coalitionPost) != nprods){
          stop("'coalitionPost' matrix must be ", nprods, " x ", nprods)
        }
        # Extract coalition indices to validate
        coalitionIndices <- which(rowSums(object@coalitionPost) > 1)
      } else {
        # Vector form - validate indices
        coalitionIndices <- object@coalitionPost
        if(any(!(coalitionIndices %in% 1:nprods))){
          stop("'coalitionPost' must contain product indices between 1 and ", nprods)
        }
      }
      
      # Coalition must have at least 2 products
      if(length(coalitionIndices) < 2){
        stop("'coalitionPost' must contain at least 2 products for price leadership")
      }
      
      # Must have at least one fringe product
      if(length(coalitionIndices) >= nprods){
        stop("At least one fringe product is required for price leadership model")
      }
    }
    
    # Timing parameter must be in (0,1) if specified
    if(length(object@timingParam) > 0){
      timingValues <- object@timingParam[!is.na(object@timingParam)]
      if(length(timingValues) > 0 && isTRUE(any(timingValues <= 0 | timingValues >= 1))){
        stop("'timingParam' values must be between 0 and 1 (exclusive)")
      }
    }
    
    # Supermarkup should be non-negative if specified
    if(!is.na(object@supermarkupPre)){
      if(object@supermarkupPre < 0){
        stop("'supermarkupPre' must be non-negative")
      }
    }
    if(!is.na(object@supermarkupPost)){
      if(object@supermarkupPost < 0){
        stop("'supermarkupPost' must be non-negative")
      }
    }
    
    # bindingFirm must be valid product index if specified
    if(!is.na(object@bindingFirm)){
      if(length(object@bindingFirm) != 1 || !(object@bindingFirm %in% 1:nprods)){
        stop("'bindingFirm' must be a single product index between 1 and ", nprods, " or NA")
      }
    }
    
    return(TRUE)
  }
)

#'@rdname BertrandRUM-Classes
#'@export
setClass(
  Class   = "PriceLeadershipBLP",
  contains = "LogitBLP",
  slots = list(
    supermarkupPre = "numeric",   # Pre-merger supermarkup m above Bertrand prices
    supermarkupPost = "numeric",  # Post-merger supermarkup m above Bertrand prices
    timingParam = "numeric",      # Firm-specific timing/discount parameters delta in (0,1) - named vector
    coalitionPre = "ANY",         # Pre-merger coalition: vector of indices OR control matrix
    coalitionPost = "ANY",        # Post-merger coalition: vector of indices OR control matrix
    bindingFirm = "numeric",      # Index of firm with binding IC constraint (or NA if none)
    slackValues = "numeric"       # Slack function values for each coalition firm
  ),
  prototype = prototype(
    supermarkupPre = NA_real_,
    supermarkupPost = NA_real_,
    timingParam = numeric(),      # Empty named vector initially
    coalitionPre = numeric(),
    coalitionPost = numeric(),
    bindingFirm = NA_real_,
    slackValues = numeric()
  ),
  validity = function(object){
    
    # Call parent validity first (LogitBLP checks)
    # Parent (LogitBLP -> Logit) already checks: prices, shares, nDraws, etc.
    
    nprods <- length(object@prices)
    
    # PriceLeadershipBLP-specific validations
    
    # Validate coalitionPre is provided
    if(length(object@coalitionPre) == 0){
      stop("'coalitionPre' must be specified (vector of product indices in coordinating coalition)")
    }
    
    ## Validate that slopes is a list (BLP parameters will be checked more thoroughly)
    if(!is.null(object@slopes) && !is.list(object@slopes)){
      stop("'slopes' must be a list containing pre-calibrated BLP demand parameters")
    }
    
    # Note: Detailed BLP parameter validation (sigma, meanval, alpha) is deferred to calcSlopes
    # because slopes may be partially populated at object creation time
    
    # Set default values for optional BLP parameters if not provided
    if(!("sigmaNest" %in% names(object@slopes))){
      object@slopes$sigmaNest <- 1  # Default: no nesting (flat logit)
    }
    if(!("nDemog" %in% names(object@slopes))){
      object@slopes$nDemog <- 0
    }
    if(!("piDemog" %in% names(object@slopes))){
      object@slopes$piDemog <- numeric(0)
    }
    if(!("nDraws" %in% names(object@slopes))){
      if(length(object@nDraws) > 0 && !is.na(object@nDraws)){
        object@slopes$nDraws <- object@nDraws
      } else {
        object@slopes$nDraws <- 500  # Default
      }
    }
    
    # Pre-merger coalition (can be vector or matrix)
    if(length(object@coalitionPre) > 0){
      if(is.matrix(object@coalitionPre)){
        # Matrix form - check dimensions
        if(nrow(object@coalitionPre) != nprods || ncol(object@coalitionPre) != nprods){
          stop("'coalitionPre' matrix must be ", nprods, " x ", nprods)
        }
        # Extract coalition indices to validate
        coalitionIndices <- which(rowSums(object@coalitionPre) > 1)
      } else {
        # Vector form - validate indices
        coalitionIndices <- object@coalitionPre
        if(any(!(coalitionIndices %in% 1:nprods))){
          stop("'coalitionPre' must contain product indices between 1 and ", nprods)
        }
      }
      
      # Coalition must have at least 2 products
      if(length(coalitionIndices) < 2){
        stop("'coalitionPre' must contain at least 2 products for price leadership")
      }
      
      # Must have at least one fringe product
      if(length(coalitionIndices) >= nprods){
        stop("At least one fringe product is required for price leadership model")
      }
    }
    
    # Post-merger coalition (can be vector or matrix)
    if(length(object@coalitionPost) > 0){
      if(is.matrix(object@coalitionPost)){
        # Matrix form - check dimensions
        if(nrow(object@coalitionPost) != nprods || ncol(object@coalitionPost) != nprods){
          stop("'coalitionPost' matrix must be ", nprods, " x ", nprods)
        }
        # Extract coalition indices to validate
        coalitionIndices <- which(rowSums(object@coalitionPost) > 1)
      } else {
        # Vector form - validate indices
        coalitionIndices <- object@coalitionPost
        if(any(!(coalitionIndices %in% 1:nprods))){
          stop("'coalitionPost' must contain product indices between 1 and ", nprods)
        }
      }
      
      # Coalition must have at least 2 products
      if(length(coalitionIndices) < 2){
        stop("'coalitionPost' must contain at least 2 products for price leadership")
      }
      
      # Must have at least one fringe product
      if(length(coalitionIndices) >= nprods){
        stop("At least one fringe product is required for price leadership model")
      }
    }
    
    # Timing parameter must be in (0,1) if specified
    if(length(object@timingParam) > 0){
      timingValues <- object@timingParam[!is.na(object@timingParam)]
      if(length(timingValues) > 0 && isTRUE(any(timingValues <= 0 | timingValues >= 1))){
        stop("'timingParam' values must be between 0 and 1 (exclusive)")
      }
    }
    
    # Supermarkup should be non-negative if specified
    if(!is.na(object@supermarkupPre)){
      if(object@supermarkupPre < 0){
        stop("'supermarkupPre' must be non-negative")
      }
    }
    if(!is.na(object@supermarkupPost)){
      if(object@supermarkupPost < 0){
        stop("'supermarkupPost' must be non-negative")
      }
    }
    
    # bindingFirm must be valid product index if specified
    if(!is.na(object@bindingFirm)){
      if(length(object@bindingFirm) != 1 || !(object@bindingFirm %in% 1:nprods)){
        stop("'bindingFirm' must be a single product index between 1 and ", nprods, " or NA")
      }
    }
    
    return(TRUE)
  }
)

#'@rdname BertrandRUM-Classes
#'@export
setClass(
  Class   = "LogitCournot",
  contains="Logit"
)

#'@rdname BertrandRUM-Classes
#'@export
setClass(
Class   = "LogitCap",
contains="Logit",
representation=representation(
  capacitiesPre           = "numeric",
  capacitiesPost          = "numeric"

),
prototype=prototype(
control.slopes = list(
  reltol= .Machine$double.eps^0.5
)
),

validity=function(object){





  nprods <- length(object@shares)


  if(nprods != length(object@capacitiesPre)){
    stop("'shares', 'capacitiesPre' must all be vectors with the same length")}
  if(length(object@capacitiesPost) != length(object@capacitiesPre)){
    stop("'capacitiesPre', 'capacitiesPost', must be vectors with the same length")}


  if(any(is.na(object@capacitiesPre) |
         #!is.finite(object@capacitiesPre) |
         object@capacitiesPre<0 ,na.rm=TRUE)){stop("'capacitiesPre' values must be positive numbers")}


  if(any(is.na(object@capacitiesPost) |
         #!is.finite(object@capacitiesPost) |
         object@capacitiesPost<0 ,na.rm=TRUE)){stop("'capacitiesPost' values must be positive numbers")}

  if(is.na(object@insideSize)){stop("'insideSize' must equal the total pre-merger units sold in the market")}

  if(any(object@insideSize*object@shares > object@capacitiesPre)){warning("utilization is greater than capacity")}

  if(identical(object@insideSize*object@shares,object@capacitiesPre)){warning("utilization equal capacity for all products")}

  if(any(is.na(object@margins[object@insideSize*object@shares == object@capacitiesPre]))){
    stop("'margins' cannot equal NA for capacity constrained products")
  }

  return(TRUE)

})


#'@rdname BertrandRUM-Classes
#'@export
setClass(
  Class   = "LogitCapALM",
  contains="LogitCap",
  representation=representation(
    parmsStart="numeric"
  ),
  prototype=prototype(
    normIndex         =  NA,
    control.slopes = list(
      factr = 1e7
    )
  ),

  validity=function(object){



    nMargins  <- length(object@margins[!is.na(object@margins)])
    if(is.na(object@insideSize) || object@insideSize <= 0){stop("'insideSize' must be greater than or equal to 0")}

    if(nMargins<2 && is.na(object@mktElast)){stop("At least 2 elements of 'margins' must not be NA in order to calibrate demand parameters")}

    if(!isTRUE(all.equal(unname(as.vector(object@shareInside)),1))){
      stop("sum of 'shares' must equal 1")
    }


    if(length(object@parmsStart)!=2){
      stop("'parmsStart' must a vector of length 2")
    }
    return(TRUE)
  })




#'@rdname BertrandRUM-Classes
#'@export
setClass(
  Class   = "LogitNests",
  contains="Logit",

  representation=representation(
    nests="factor",
    parmsStart="numeric",
    constraint="logical"
  ),
  prototype=prototype(
    parmsStart      =  numeric(),
    control.slopes = list(
      factr = 1e7
    )
  ),

  validity=function(object){




    nprods    <- length(object@prices)
    nNestParm <- nlevels(object@nests) #calculate the number of nesting parameters

    ## Identify Singleton Nests
    nestCnt   <- tapply(object@prices,object@nests,length)
    nestCnt   <- nestCnt[object@nests]
    isSingleton <- nestCnt==1

    nNestParm <- nNestParm - sum(isSingleton) #singleton nests are not identified

    if(identical(nNestParm,1)) stop("'logit.nests', 'logit.nests.alm' may not be used for non-nested problems or problems with only singleton nests. Use 'logit', 'logit.alm' instead")

    if(nprods != length(object@nests)){
      stop("'nests' length must equal the number of products")}



    if(!object@constraint &&
       any(tapply(object@margins[!isSingleton],object@nests[!isSingleton],
                  function(x){if(all(is.na(x))){return(TRUE)} else{return(FALSE)}}
       )
       ,na.rm=TRUE)
    ){
      stop("when 'constraint' is FALSE, at least one product margin must be supplied for each non-singleton nest")
    }



    return(TRUE)
  }

)



#'@rdname BertrandRUM-Classes
#'@export
setClass(
  Class   = "LogitNestsALM",
  contains="LogitNests",
  prototype=prototype(
    normIndex         =  1
  ),
  validity=function(object){


    if(!isTRUE(all.equal(unname(as.vector(object@shareInside)),1))){
      stop("sum of 'shares' must equal 1")
    }


  }
)



#'@rdname BertrandRUM-Classes
#'@export
setClass(
  Class   = "LogitALM",
  contains="Logit",
  representation=representation(
    parmsStart="numeric"
  ),
  prototype=prototype(
    normIndex         =  NA,
    control.slopes = list(
      factr = 1e7
    )
  ),

  validity=function(object){



    nMargins  <- length(object@margins[!is.na(object@margins)])

    if(nMargins<2 && is.na(object@mktElast)){stop("At least 2 elements of 'margins' must not be NA in order to calibrate demand parameters")}

    if(!isTRUE(all.equal(unname(as.vector(object@shareInside)),1))){
      stop("sum of 'shares' must equal 1")
    }


    if(length(object@parmsStart)!=2){
      stop("'parmsStart' must a vector of length 2")
    }
  }
)


#'@rdname BertrandRUM-Classes
#'@export
setClass(
  Class   = "CES",
  contains="Logit",
  prototype=prototype(
    priceOutside=1
    
  ),
  validity=function(object){
    
   if(!is.na(object@mktElast) && object@mktElast > -1){
     stop("'mktElast' must be less than or equal to  -1")} 
  }
)


#'@rdname BertrandRUM-Classes
#'@export
setClass(
  Class   = "CESALM",
  contains="CES",
  representation=representation(
    parmsStart="numeric"
  ),
  prototype=prototype(
    normIndex         =  NA,
    control.slopes = list(
      factr = 1e7
    )
  ),

  validity=function(object){



    nMargins  <- length(object@margins[!is.na(object@margins)])

    if(nMargins<2 && is.na(object@mktElast)){stop("At least 2 elements of 'margins' must not be NA in order to calibrate demand parameters")}

    if(!isTRUE(all.equal(unname(as.vector(object@shareInside)),1))){
      stop("sum of 'shares' must equal 1")
    }

    if(length(object@parmsStart)!=2){
      stop("'parmsStart' must a vector of length 2")
    }
  }
)


#'@rdname BertrandRUM-Classes
#'@export
setClass(
Class   = "CESNests",
contains="CES",

representation=representation(
  nests="factor",
  parmsStart="numeric",
  constraint="logical"
),

prototype=prototype(
  parmsStart      =  numeric(),
  constraint=TRUE
),

validity=function(object){




  nprods    <- length(object@prices)
  nNestParm <- nlevels(object@nests) #calculate the number of nesting parameters
  nMargins  <- length(object@margins[!is.na(object@margins)])
  maxNests  <- nMargins - 1

  ## Identify Singleton Nests
  nestCnt   <- tapply(object@prices,object@nests,length)
  nestCnt   <- nestCnt[object@nests]
  isSingleton <- nestCnt==1

  nNestParm <- nNestParm - sum(isSingleton) #singleton nests are not identified

  if(identical(nNestParm,1)) stop("'ces.nests' cannot be used for non-nested problems or problems with only singleton nests. Use 'ces' instead")

  if(nprods != length(object@nests)){
    stop("'nests' length must equal the number of products")
  }

  if(object@constraint && length(object@parmsStart)!=2){
    stop("when 'constraint' is TRUE, 'parmsStart' must be a vector of length 2")
  }
  else if(!object@constraint && nNestParm + 1 != length(object@parmsStart)){
    stop("when 'constraint' is FALSE, 'parmsStart' must be a vector of length ",nNestParm + 1)
  }


  if(!object@constraint &&
     any(tapply(object@margins[!isSingleton],object@nests[!isSingleton],
                function(x){if(all(is.na(x))){return(TRUE)} else{return(FALSE)}}
     )
     ,na.rm=TRUE)
  ){
    stop("when 'constraint' is FALSE, at least one product margin must be supplied for each non-singleton nest")
  }


  if(nNestParm > nMargins){
    stop(paste(
      "Impossible to calibrate nest parameters with the number of margins supplied.\n",
      "The maximum number of nests supported by the supplied margin information is"
      ,maxNests,"."))
  }
}

)
