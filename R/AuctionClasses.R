#'@title Class \dQuote{Auction}
#'@name Auction-Classes

#'@aliases Auction2ndCap-class
#'Auction2ndLogit-class
#'Auction2ndLogitNests-class
#'Auction2ndLogitALM-class
#Auction2ndLogitNestsALM-class
#'@description The \dQuote{Auction2ndCap} class contains all the information needed to
#'calibrate a 2nd price auction with capacity constraints
#'@description The \dQuote{Auction2ndLogit} class contains all the information needed to
#'calibrate a Logit
#'demand system and perform a merger simulation analysis under the assumption that
#'firms are setting offers in a 2nd-score auction.
#'@description The \dQuote{Auction2ndLogitNests} class contains all the information needed to
#'calibrate a Nested Logit
#'demand system and perform a merger simulation analysis under the assumption that
#'firms are setting offers in a 2nd-score auction.
#'@description The \dQuote{Auction2ndLogitALM} class contains all the information needed to
#'calibrate a Logit
#'demand system with unobserved outside share and perform a merger simulation analysis under the assumption that
#'firms are setting offers in a 2nd-score auction.
#'@description Below, let k denote the number of firms.
#'
#'@section Objects from the Class:
#'Auction2ndCap: Objects can be created by using the constructor function \code{\link{auction2nd.cap}}.
#'
#'Auction2ndLogit: Objects can be created by using the constructor function \code{\link{auction2nd.logit}}.
#'
#'Auction2ndLogitNests: Objects can be created by using the constructor function \code{\link{auction2nd.logit.nests}}.
#'
#'Auction2ndLogitALM: Objects can be created by using the constructor function \code{\link{auction2nd.logit.alm}}.

#'@slot capacities A length k vector of firm capacities.
#'@slot margins A length k vector of product margins, some of which may
#'equal NA.
#'@slot prices A length k vector of product prices.
#'@slot reserve A length 1 vector equal to observed buyer's
#'reserve price. May equal NA.
#'@slot shareInside A length 1 vector equal to the
#'probability that a buyer does not select the outside option. May equal NA.
#'@slot sellerCostCDF A length 1 character vector equal to
#'the name of the function that calculates the Cumulative
#'Distribution (CDF) of SellerCosts.
#'@slot sellerCostCDFLowerTail A length 1 logical vector equal to
#'TRUE if the  probabilities are P[X <= x], otherwise, P[X > x].
#'@slot sellerCostPDF A function returning the Probability Density
#'of Seller Costs.
#'@slot sellerCostBounds The bounds on the seller's CDF.
#'@slot sellerCostParms The parameters of the seller's CDF.
#'@slot buyerValuation Buyer's self-supply cost.
#'@slot reservePre Buyer's optimal pre-merger reservation price.
#'@slot reservePost Buyer's optimal post-merger reservation
#'price.
#'@slot mcDelta A length k vector equal to the proportional
#'change in a firm's capacity following the merger.
#'@slot parmsStart A vector of starting values.

#'@section Extends:
#'Auction2ndCap: Class \code{\linkS4class{Antitrust}}, directly.
#'
#'Auction2ndLogit: Class \code{\linkS4class{Logit}}, directly.
#'Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Logit}}, distance 2.
#'Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 3.
#'
#'Auction2ndLogitALM: Class \code{\linkS4class{Auction2ndLogit}}, directly.
#'Class \code{\linkS4class{Logit}}, distance 2.
#'Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Logit}}, distance 3.
#'Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 4.

#'@author Charles Taragin \email{ctaragin+antitrustr@gmail.com}
#'@examples
#'showClass("Auction2ndCap")           # get a detailed description of the class
#'showClass("Auction2ndLogit")           # get a detailed description of the class
#'showClass("Auction2ndLogitALM")           # get a detailed description of the class
#'@keywords classes
#'@include BertrandRUMClasses.R
NULL

#'@rdname Auction-Classes
#'@export
setClass(
  Class   = "Auction2ndCap",
  contains="Antitrust",
  representation=representation(
    capacities       = "numeric",
    margins          = "numeric",
    prices           = "numeric",
    reserve          = "numeric",
    shareInside      = "numeric",
    sellerCostCDF    = "character",
    sellerCostCDFLowerTail    = "logical",
    sellerCostPDF    = "function",
    sellerCostBounds = "numeric",
    sellerCostParms  = "numeric",
    buyerValuation        = "numeric", #was buyerValuation
    reservePre       = "numeric",
    reservePost      = "numeric",
    mcDelta          = "numeric",
    parmsStart       = "numeric"

  ),
  prototype=prototype(
    reservePre      =  numeric(),
    reservePost     =  numeric(),
    buyerValuation  =  numeric(),
    sellerCostParms =  numeric(),
    sellerCostCDFLowerTail    = TRUE


  ),
  validity=function(object){

    nprods <- length(object@labels)

    cdf    <- object@sellerCostCDF
    if(is.na(object@reserve)){parmsStart <- object@parmsStart[-1]}
    else{parmsStart <- object@parmsStart}

    if(nprods != length(object@capacities) ||
       nprods != length(object@margins)    ||
       nprods != length(object@prices)     ||
       nprods != length(object@ownerPre)   ||
       nprods != length(object@ownerPost)  ||
       nprods != length(object@mcDelta)
    ){
      stop("'capacities', 'margins', 'prices', 'ownerPre', 'ownerPost', 'mcDelta' must all be the same length")}
    if(length(object@reserve) != 1 ||
       (!is.na(object@reserve) && object@reserve<0)){
      stop("'reserve'must be a length 1 vector whose value is greater than 0")}
    if(
      (!is.na(object@shareInside) &&
       (object@shareInside > 1 ||
        object@shareInside < 0)) ||
      length(object@shareInside) != 1){
      stop("'shareInside' must be a length 1 vector whose value is between 0 and 1")
    }
    if(any(object@capacities<0 | is.na(object@capacities),na.rm=TRUE)){
      stop("'capacities' cannot be negative or equal to NA")
    }
    if(any(object@prices<=0,na.rm=TRUE)){
      stop("'prices' must be positive")}
    if(any(object@margins > 1 | object@margins < 0,na.rm=TRUE)){
      stop("'margins' must be between 0 and 1")
    }


    if( identical(cdf,"punif")){
      if(length(parmsStart)!=2){
        if(is.na(object@reserve)){stop("parmsStart must be a length 3 vector whose first element is the starting value for 'reserve'")}
        else{stop("For the Uniform distribution, 'parmsStart' must be a numeric vector of length 2")}
      }
      if(parmsStart[1] >= parmsStart[2]){
        stop("The upper bound must be greater than the lower bound")
      }
    }

    else if( identical(cdf,"pexp")){
      if(length(parmsStart)!=1){
        if(is.na(object@reserve)){stop("parmsStart must be a length 2 vector whose first element is the starting value for 'reserve'")}
        else{stop("For the Exponential distribution, 'parmsStart' must be a numeric vector of length 1")}
      }
      if(parmsStart[1] <=0){
        stop("For the Exponential distribution, 'parmsStart' must be a numeric vector of length ",
             1 + is.na(object@reserve),"whose final element is greater than 0")
      }
    }

    else if( identical(cdf,"pweibull")){
      if(length(parmsStart)!=2){
        if(is.na(object@reserve)){stop("parmsStart must be a length 3 vector whose first element is the starting value for 'reserve'")}
        else{stop("For the Weibull distribution, 'parmsStart' must be a numeric vector of length 2")}
      }
      if(parmsStart[1] <=0  ||
         parmsStart[2] <=0 ){
        stop("For the Weibull distribution, 'parmsStart' must be a numeric vector of length ",
             2 + is.na(object@reserve)," whose final 2 elements must be greater than 0")
      }
    }

    else if( identical(cdf,"pgumbel")){
      if(length(parmsStart)!=2){
        if(is.na(object@reserve)){stop("parmsStart must be a length 3 vector whose first element is the starting value for 'reserve'")}
        else{stop("For the Gumbel distribution, 'parmsStart' must be a numeric vector of length 2")}
      }
      if(parmsStart[2] <=0){
        stop("For the Gumbel distribution, 'parmsStart' must be a numeric vector of length ",
             2 + is.na(object@reserve)," whose final element must be greater than 0")
      }
    }
    else if( identical(cdf,"pfrechet")){
      if(length(parmsStart)!=3){
        if(is.na(object@reserve)){stop("parmsStart must be a length 4 vector whose first element is the starting value for 'reserve'")}
        else{stop("For the Frechet distribution, 'parmsStart' must be a numeric vector of length 3")}
      }
      if(parmsStart[2] <=0  ||
         parmsStart[3] < 2  ){
        stop("For the Frechet distribution, 'parmsStart' must be a numeric vector of length ",
             2 + is.na(object@reserve)," whose next-to-last element must be positive and whose final element must be at least 2")
      }
    }


    data  <- sum(!is.na(object@reserve) , !is.na(object@shareInside) , sum(!is.na(object@prices)) , sum(!is.na(object@margins)))
    unknowns <-  is.na(object@reserve) + length(parmsStart)

    if(data < unknowns ){
      stop("Insufficient information to calibrate model parameters: ",unknowns, " unknowns but only ", data, " pieces of information")
    }

    return(TRUE)
  }
)

#'@rdname Auction-Classes
#'@export
setClass(
  Class   = "Auction2ndLogit",
  contains="Logit",
  
  prototype=prototype(
    priceStart  = numeric(),
    control.slopes = list(
      reltol= .Machine$double.eps^0.5
    )
  )
  
  )

#'@rdname Auction-Classes
#'@export
setClass(
  Class   = "Auction2ndLogitALM",
  contains="Auction2ndLogit",
  representation=representation(
    parmsStart="numeric"
  ),
  prototype=prototype(
    normIndex         =  NA,
    control.slopes = list(
      trace=FALSE,
      ftol = 1e-10
    )
  ),

  validity=function(object){



    nMargins  <- length(object@margins[!is.na(object@margins)])

    if(!is.na(object@mktElast) && all(is.na(object@prices))){stop("At least 1 price must be supplied")}

    if(nMargins<2 && is.na(object@mktElast)){stop("At least 2 elements of 'margins' must not be NA in order to calibrate demand parameters")}

    if(!isTRUE(all.equal(unname(as.vector(object@shareInside)),1))){
      stop("sum of 'shares' must equal 1")
    }

    if(length(object@parmsStart)!=2){
      stop("'parmsStart' must a vector of length 2")
    }
  }
)

#'@rdname Auction-Classes
#'@export
setClass(
  Class   = "Auction2ndLogitNests",
  contains="Auction2ndLogit",
  
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
    
    
    
    
    nprods    <- length(object@shares)
    nNestParm <- nlevels(object@nests) #calculate the number of nesting parameters
    
    ## Identify Singleton Nests
    nestCnt   <- tapply(object@shares,object@nests,length)
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

#' 
#' #'@rdname Auction-Classes
#' #'@export
#' setClass(
#'   Class   = "Auction2ndLogitNestsALM",
#'   contains="Auction2ndLogitNests"
#'   
#'   
#' )
