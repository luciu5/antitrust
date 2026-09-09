#' @title Producer Surplus Methods
#' @name PS-methods
#' @docType methods
#' @aliases calcProducerSurplus calcProducerSurplus,ANY-method calcProducerSurplus,Bertrand-method calcProducerSurplus,Cournot-method calcProducerSurplus,VertBargBertLogit-method
#'
#' @description In the following methods, \code{calcProducerSurplus} computes the expected profits of each supplier
#' with the game depending on the class. The available classes are: Bertrand, Cournot, and Auction2ndCap.
#'
#' @param object An instance of one of the classes listed above.
#' @param preMerger If TRUE, returns pre-merger outcome. If FALSE, returns post-merger outcome. Default is TRUE.
#' @param exAnte If \sQuote{exAnte} equals TRUE then the
#' \emph{ex ante} expected result for each firm is produced, while FALSE produces the
#' expected result conditional on each firm winning the auction. Default is TRUE.
#' @include PriceDeltaMethods.R
#' @keywords methods
NULL


setGeneric (
  name= "calcProducerSurplus",
  def=function(object,...){standardGeneric("calcProducerSurplus")}
)


## compute producer surplus
#'@rdname PS-methods
#'@export
setMethod(
  f= "calcProducerSurplus",
  signature= "Bertrand",
  definition=function(object,preMerger=TRUE){


    margins <- calcMargins(object,preMerger,level=TRUE)


    output <- calcQuantities(object,preMerger)

    if (all(is.na(output))){
      warning("'insideSize' is missing; using normalized shares instead of quantities. Producer-surplus results are not in market units.")
      output <- calcShares(object,preMerger,revenue=FALSE)
    }

    ps <- margins * output
    names(ps) <- object@labels

    return(ps)
  }

)
#'@rdname PS-methods
#'@export
setMethod(
  f= "calcProducerSurplus",
  signature= "VertBargBertLogit",
  definition=function(object,preMerger=TRUE){

    mktSize <- object@down@mktSize

    margins <- calcMargins(object,preMerger=preMerger,level=TRUE)

    output <- calcShares(object,preMerger)

    if (is.na(mktSize)){
      warning("'insideSize' is missing; using normalized shares instead of quantities. Producer-surplus results are not in market units.")
      mktSize <- 1
    }

    psup <- margins$up * output * mktSize
    psdown <- margins$down * output * mktSize
    names(psup) <- names(psdown) <-  object@down@labels

    return(list(up=psup,down=psdown))
  }

)



#'@rdname PS-methods
#'@export
setMethod(
  f= "calcProducerSurplus",
  signature= "Auction2ndCap",
  definition=function(object,preMerger=TRUE,exAnte=TRUE){

    sellerCostBounds <-object@sellerCostBounds
    if(preMerger){r    <- object@reservePre}
    else{r    <- object@reservePost}


    if(preMerger) { capacities = object@capacities }
    else {          capacities = tapply(object@capacities*(1+object@mcDelta),object@ownerPost,sum) }

    totCap = sum(capacities)

    espIntegrand = function(c,t){
      sellerCostParms <- c(list(c),as.list(object@sellerCostParms),
                           lower.tail=as.list(object@sellerCostCDFLowerTail))
      Fc <- do.call(match.fun(object@sellerCostCDF),sellerCostParms)
      val <- (1-Fc)^(totCap-t)-(1-Fc)^totCap
    }


    retval <- sapply(
      capacities,
      function(t.i) {

        if( r < sellerCostBounds[2]) {
          retval <- integrate(espIntegrand,lower=sellerCostBounds[1],upper=r,
                              stop.on.error = FALSE,t=t.i)$value
        }
        else {
          retval <- integrate(espIntegrand,lower=sellerCostBounds[1],upper=sellerCostBounds[2],
                              stop.on.error = FALSE,t=t.i)$value
        }

        return(retval)
      })

    if(!preMerger){

      temp <- rep(NA, length(object@ownerPre))
      temp[object@ownerPre == object@ownerPost] <- retval
      retval <- temp

    }

    if(!exAnte){retval <- retval/calcShares(object,preMerger=preMerger,exAnte=TRUE)}

    return(retval)
  })


## compute producer surplus
#'@rdname PS-methods
#'@export
setMethod(
  f= "calcProducerSurplus",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){


    rev <-  calcRevenues(object, preMerger= preMerger)
    vc <- calcVC(object, preMerger= preMerger)

    ps <- rowSums(rev, na.rm=TRUE) - vc
    names(ps) <- object@labels[[1]]

    return(ps)
  }

)
