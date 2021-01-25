#' @title Methods for Calculating Diagnostics
#' @name Diagnostics-Methods
#' @docType methods
#'
#' @aliases calcDiagnostics
#' calcDiagnostics,ANY-method
#' calcDiagnostics,Bertrand-method
#' calcDiagnostics,Cournot-method
#' calcDiagnostics,VertBargBertLogit-method
#' 
#' @description Computes the percentage difference between predicted and observed pre-merger prices, shares,
#' margins and market elasticities (if supplied) . \sQuote{labels} is used to specify row labels.
#' @param object An instance of one of the classes listed above.
#' @param labels A length-k vector of product labels. Default is object@labels.
#'
#' @include CVMethods.R
#' @keywords methods
NULL

setGeneric (
  name= "calcDiagnostics",
  def=function(object,...){standardGeneric("calcDiagnostics")}
)

#'@rdname Diagnostics-Methods
#'@export
setMethod(
  f= "calcDiagnostics",
  signature= "Bertrand",
  definition=function(object,labels=object@labels){

    obsPrices <- object@prices
    obsShares <- object@shares
    obsShares <- obsShares/sum(obsShares,na.rm=TRUE)
    obsMargins <- object@margins
    obsElast <- object@mktElast

    prePrices <- unname(drop(object@pricePre))
    preMargins <- drop(calcMargins(object, preMerger=TRUE))
    preShares <- drop(calcShares(object, preMerger=TRUE))
    preShares <- drop(preShares/sum(preShares))
    preElast <- elast(object, preMerger=TRUE, market=TRUE)

    res <- data.frame(
      Prices= 1 - obsPrices/prePrices,
      Shares=1 - obsShares/preShares,
      Margins= 1 - obsMargins/preMargins,
      'Market Elasticity'= 1 - obsElast/preElast,
      check.names = FALSE
    )*100

    #rmThese <- colSums(abs(res),na.rm=TRUE)

    #res[-1,'Market Elasticity'] <- NA



    rownames(res) <- labels


    return(res)
  }
)


#'@rdname Diagnostics-Methods
#'@export
setMethod(
  f= "calcDiagnostics",
  signature= "VertBargBertLogit",
  definition=function(object,labels=object@down@labels){
    
    
    down <- object@down
    up <- object@up
    
    
    obsUpPrices <- up@prices
    obsUpMargins <- up@margins
    
    obsDownPrices <- down@prices
    obsDownMargins <- down@margins
    
    obsShares <- down@shares
    obsShares <- obsShares/sum(obsShares,na.rm=TRUE)
    obsElast <- down@mktElast
    
    
    preUpPrices <- unname(drop(up@pricePre))
    preDownPrices <- unname(drop(down@pricePre))
    
    preMargins <- calcMargins(object, preMerger=TRUE)
    preShares <- drop(calcShares(down, preMerger=TRUE))
    preShares <- drop(preShares/sum(preShares))
    preElast <- elast(object, preMerger=TRUE, market=TRUE)
    
    res <- data.frame(
      upPrices=1 - obsUpPrices/preUpPrices,
      downPrices= 1 - obsDownPrices/preDownPrices,
      Shares=1 - obsShares/preShares,
      upMargins= 1 - obsUpMargins/preMargins$up,
      downMargins= 1 - obsDownMargins/preMargins$down,
      'Market Elasticity'= 1 - obsElast/preElast,
      check.names = FALSE
    )*100
    
    #rmThese <- colSums(abs(res),na.rm=TRUE)
    
    #res[-1,'Market Elasticity'] <- NA
    
    
    
    rownames(res) <- labels
    
    return(res)
  }
)



#'@rdname Diagnostics-Methods
#'@export
setMethod(
  f= "calcDiagnostics",
  signature= "Cournot",
  definition=function(object){

    callNextMethod(object,label=object@labels[[1]])

  })
