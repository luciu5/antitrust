#'@title \dQuote{Vertical} Classes
#'@name Vertical-Classes
#'@aliases VertBargBertLogit VertBarg2ndLogit VertBargBertLogitNests VertBarg2ndLogitNests

#'@description The \dQuote{Vertical} classes are building blocks used to create other classes
#'in this package. As such, it is most likely to be useful for developers
#'who wish to code their own calibration/simulation routines.

#'@description The \dQuote{VertBargBertLogit} class has the information for a Vertical Supply Chain with Logit demand and a downstream Nash-Bertrand Pricing game.

#'@description The \dQuote{VertBarg2ndLogit} class has the information fora Vertical Supply Chain with Logit demand and a downstream 2nd Score Auction.
#'@slot up an instance of \dQuote{Bargaining} class.
#'@slot down For  \dQuote{VertBargBertLogit}, an instance of  \dQuote{Logit} class.For  \dQuote{VertBarg2ndLogit}, an instance of  \dQuote{Auction2ndLogit} class.
#'@slot constrain  A length 1 character vector equal to "global", "pair", "wholesaler", or "retailer.

#'@author Charles Taragin \email{ctaragin+antitrustr@gmail.com}
#'@include BertrandRUMClasses.R AuctionClasses.R BargainingClasses.R
#'@keywords classes
#'
NULL

#'@rdname Vertical-Classes
#@export
setClass(
  
  Class = "VertBargBertLogit",
  representation=representation(
     up = "Bargaining",
     down = "Logit",
     supplyDown="character",
     isHorizontal="logical",
     isUpstream="logical",
     ownerDownPre="matrix",
     ownerDownPost="matrix",
     ownerDownLambdaPre="matrix",
     ownerDownLambdaPost="matrix",
     ownerUpLambdaPre="matrix",
     ownerUpLambdaPost="matrix",
     constrain        = "character",
     control.slopes =   "list",
     control.eq =   "list"
   ),
  prototype=prototype(
    control.slopes = list(
      trace=FALSE,
      ftol = 1e-10
    )
  ),
  validity = function(object){
    
    if( length(object@constrain) != 1 ||
        !object@constrain %in% c("global","pair","wholesaler","retailer")){
      stop("'constrain' must equal 'global','pair','wholesaler','retailer'")
    }
    if(length(object@supplyDown) != 1 ||
       !object@supplyDown %in% c("bertrand","2nd") ){
      stop("'supplyDown' must equal 'bertrand' or '2nd'")
    }
    
  }
 
)

#'@rdname Vertical-Classes
#@export
setClass(
  
  Class = "VertBarg2ndLogit",
  contains = "VertBargBertLogit"
  ,
  representation=representation(
    up = "Bargaining",
    down = "Auction2ndLogit"
  )
  
)


#'@rdname Vertical-Classes
#@export
setClass(
  
  Class = "VertBarg2ndLogitNests",
  contains = "VertBarg2ndLogit"
  ,
  representation=representation(
    up = "Bargaining",
    down = "Auction2ndLogitNests"
  )
  
)
