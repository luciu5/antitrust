#'@title \dQuote{Debt} Classes
#'@name Debt-Classes
#'@aliases LogitDebt

#'@description The \dQuote{Debt} classes are building blocks used to create other classes
#'in this package. As such, it is most likely to be useful for developers
#'who wish to code their own calibration/simulation routines.

#'@description The \dQuote{LogitDebt} class has the information for a multi-market Nash Bertrand Model with Logit Demand and Debt.

#'@slot up an instance of \dQuote{Bargaining} class.
#'@slot down For  \dQuote{VertBargBertLogit}, an instance of  \dQuote{Logit} class.For  \dQuote{VertBarg2ndLogit}, an instance of  \dQuote{Auction2ndLogit} class.
#'@slot constrain  A length 1 character vector equal to "global", "pair", "wholesaler", or "retailer.
#'@slot chain_level A length 1 character vector equal to "full", "wholesaler", or "retailer".

#'@author Charles Taragin \email{ctaragin+antitrustr@gmail.com}
#'@include AntitrustClasses.R 
#'@keywords classes
#'


#'@rdname Debt-Classes
#@export
# setClass(
#   Class   = "LogitDebt",
#   contains="Logit",
#   representation=representation(
#     marketID = "character",
#     productID =  "character",
#     debtPre  =   "numeric",
#     debtPost =   "numeric",
#     focal =      "numeric",
#     parmsStart = "numeric",
#     hFocalIntegral = "function",
#     gFocalIntegral = "function", 
#     tOtherIntegral = "function",
#     rOtherIntegral = "function"
#   ),
#   prototype=prototype(
#     normIndex         =  1,
#     focal             =  1,
#     labels= list(mkt=NULL,product=NULL),
#     control.slopes = list(
#       factr = 1e7
#     )
#   ),
#   
#   validity=function(object){
#     
#     nProducts =length(unique(object@productID))
#     nFirms=length(unique(object@ownerPre))
#     nMarkets=length(unique(object@marketID))
#     
#     if(length(object@debtPre) != nFirms || 
#        length(object@debtPre) != length(object@debtPost)){
#       stop("'debtPre' and 'debtPost' must have the same number of elements as columns in 'prices'")}
#     
#     nMargins  <- length(object@margins[!is.na(object@margins)])
#     
#     if(nMargins< nMarkets+2){stop("At least", nMarkets+2, "elements of 'margins' must not be NA in order to calibrate demand parameters")}
#     
#     if(!isTRUE(as.numeric(all.equal(tapply(object@shares, object@marketID,sum,na.rm=TRUE)),rep(1,nMarkets)))){
#       stop("sum of 'shares' must equal 1")
#     }
#     
#     
#     if(length(object@parmsStart)!=nMarkets+2){
#       stop("'parmsStart' must a vector of length ",nMarkets+2)
#     }
#   }
# )
