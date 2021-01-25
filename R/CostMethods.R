#' @title Methods for Calculating marginal and Variable Costs
#' @name Cost-Methods
#' @docType methods
#'
#' @aliases calcMC
#' calcMC,ANY-method
#' calcMC,Bertrand-method
#' calcMC,VertBargBertLogit-method
#' calcMC,Auction2ndLogit-method
#' calcMC,Cournot-method
#' calcMC,Auction2ndCap-method
#' calcdMC
#' calcdMC,ANY-method
#' calcdMC,Stackelberg-method
#' calcVC
#' calcVC,ANY-method
#' calcVC,Cournot-method
#'
#' @description
#' For Auction2ndCap, calcMC calculates (constant) marginal cost for each
#' product. For those classes that do not require prices, returns a
#' length-k vector of NAs when prices are not supplied.
#'
#' For Bertrand, calcMC computes either pre- or post-merger marginal costs. Marginal costs
#' are assumed to be constant. Post-merger marginal costs are equal to
#' pre-merger marginal costs multiplied by 1+\sQuote{mcDelta}, a length-k
#' vector of marginal cost changes. \sQuote{mcDelta} will typically be between 0 and 1.
#'
#' For Auction2ndLogit, calcMC computes constant marginal costs impied by the model.
#'
#' For Cournot, calcMC calculates marginal cost for each product.
#'
#' calcdMC computes the derivative of either pre- or post-merger marginal costs. The derivative of Marginal costs
#' is assumed to be constant. Post-merger marginal costs are equal to
#' pre-merger marginal costs multiplied by 1+\sQuote{mcDelta}, a length-k
#' vector of marginal cost changes. \sQuote{mcDelta} will typically be between 0 and 1.
#'
#' calcVC computes either pre- or post-merger variable costs. Variable costs
#' are assumed to be quadratic by default. Post-merger variable costs are equal to
#' pre-merger variable costs multiplied by 1+\sQuote{mcDelta}, a length-k
#' vector of marginal cost changes. \sQuote{mcDelta} will typically be between 0
#' and 1.
#'
#' @param object An instance of the respective class (see description for the classes)
#' @param  preMerger If TRUE, the pre-merger ownership structure is used. If FALSE, the post-merger ownership structure is used.
#' Default is TRUE.
#' @param exAnte If \sQuote{exAnte} equals TRUE then the
#' \emph{ex ante} expected result for each firm is produced, while FALSE produces the
#' expected result conditional on each firm winning the auction. Default is FALSE.
#' @param t The capacity profile of each supplier. Default is \sQuote{preMerger} capacities.
#'
#' @include PricesMethods.R
#' @keywords methods
NULL

setGeneric (
  name= "calcMC",
  def=function(object,...){standardGeneric("calcMC")}
)

setGeneric (
  name= "calcdMC",
  def=function(object,...){standardGeneric("calcdMC")}
)

setGeneric (
  name= "calcVC",
  def=function(object,...){standardGeneric("calcVC")}
)

## Create a method to recover marginal cost using
## demand parameters and supplied prices
#'@rdname Cost-Methods
#'@export
setMethod(
  f= "calcMC",
  signature= "Bertrand",
  definition= function(object,preMerger=TRUE){

    object@pricePre <- object@prices


    marginPre <- calcMargins(object,preMerger = TRUE)

    mc <- (1 - marginPre) * object@prices

    if(!preMerger){
      mc <- mc*(1+object@mcDelta)
    }

    mc <- as.vector(mc)

    names(mc) <- object@labels



    isNegMC <- mc < 0

    if( preMerger && any(isNegMC, na.rm = TRUE)){

      warning(paste("Negative marginal costs were calibrated for the following firms:", paste(object@labels[isNegMC], collapse=",")))
    }

    return(mc)
  }
)

#'@rdname Cost-Methods
#'@export
setMethod(
  f= "calcMC",
  signature= "VertBargBertLogit",
  definition= function(object,preMerger=TRUE){
    
    up <- object@up
    down <- object@down
    
  
      if(length(up@pricePre) == 0 ){
        priceUpPre <- up@prices
        object@up@pricePre <- up@prices
      }
      else{priceUpPre <- up@pricePre}
    
      if(length(down@pricePre) == 0 ){
        priceDownPre <- down@prices
        object@down@pricePre <- priceDownPre
      }
      else{priceDownPre <- object@down@pricePre}
    
    
    
    marginsPre <- calcMargins(object,preMerger = TRUE, level=TRUE)
    
   
    mcDown <- -(marginsPre$down - priceDownPre + priceUpPre)  
    mcUp <- -(marginsPre$up  - priceUpPre)  
    
    
    
    
    if(!preMerger){
      mcUp <- mcUp*(1+up@mcDelta)
      mcDown <- mcDown*(1+down@mcDelta)
    }
    
    mcUp <- as.vector(mcUp)
    mcDown <- as.vector(mcDown)
    
    names(mcUp) <- up@labels
    names(mcDown) <- down@labels
    
    
    
    isNegUpMC <- mcUp < 0
    isNegDownMC <- mcDown < 0
    
    if( preMerger && any(isNegUpMC, na.rm = TRUE)){
      
      warning(paste("Negative upstream marginal costs were calibrated for the following firms:", paste(up@labels[isNegUpMC & !is.na(isNegUpMC)], collapse=",")))
    }
    if( preMerger && any(isNegDownMC, na.rm = TRUE)){
      
      warning(paste("Negative downstream marginal costs were calibrated for the following firms:", paste(down@labels[isNegDownMC& !is.na(isNegDownMC)], collapse=",")))
    }
    
    return(list(up=mcUp,down=mcDown))
  }
)


#'@rdname Cost-Methods
#'@export
setMethod(
  f= "calcMC",
  signature= "Auction2ndCap",
  definition=function(object,t,preMerger=TRUE,exAnte=TRUE){


    cdfF <- match.fun(object@sellerCostCDF)
    pdfF <- object@sellerCostPDF
    sellerCostBounds <-object@sellerCostBounds



    if(preMerger) {
      capacities <- object@capacities
      r    <- object@reservePre
    }
    else {
      capacities <- tapply(object@capacities*(1+object@mcDelta),object@ownerPost,sum)
      r    <- object@reservePost
    }

    totCap <- sum(capacities)

    if(missing(t)){t <- capacities}



    ## The expected production cost
    ecIntegrand = function(c,t){
      sellerCostParms <- c(list(c),as.list(object@sellerCostParms))

      fc = do.call(pdfF,sellerCostParms)

      sellerCostParms <- c(sellerCostParms,
                           lower.tail=as.list(object@sellerCostCDFLowerTail))
      Fc = do.call(cdfF,sellerCostParms)

      retval = t*c*fc*(1-Fc)^(totCap-1)
      retval = ifelse(is.finite(retval),retval,0)

      return(retval)
    }

    result <- sapply(
      t,
      function(t.i) {
        if( r < sellerCostBounds[2]) {
          retval = integrate(ecIntegrand,lower=sellerCostBounds[1],upper=r, stop.on.error = FALSE,t=t.i)$value
        }
        else {
          retval = integrate(ecIntegrand,lower=sellerCostBounds[1],upper=sellerCostBounds[2], stop.on.error = FALSE,t=t.i)$value
        }

        return(retval)
      })


    if(!preMerger && length(t)>1){

      temp <- rep(NA, length(object@ownerPre))
      temp[object@ownerPre == object@ownerPost] <- result
      result <- temp

    }

    if(!exAnte){result <- result/calcShares(object,preMerger=preMerger,exAnte=TRUE)}

    return(result)

  }
)


#'@rdname Cost-Methods
#'@export
setMethod(
  f= "calcMC",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){

    if(preMerger){

      quantity  <- object@quantityPre
      mcfun <- object@mcfunPre
      cap <- object@capacitiesPre
    }
    else{

      quantity  <- object@quantityPost
      mcfun <- object@mcfunPost
      cap <- object@capacitiesPost
    }

    plantQuant <- rowSums(quantity,na.rm=TRUE)



    nplants <- nrow(quantity)

    mc <- rep(NA, nplants)

    for(f in 1:nplants){
      mc[f] <- mcfun[[f]](quantity[f,])
    }

    if(!preMerger){mc <- mc*(1 + object@mcDelta)}

    mc <- mc + 1/(100*(pmax(cap - plantQuant, 1e-16))) + 1/(100*(pmax(1e-16,plantQuant)))
    #mc <- ifelse(plantQuant <= cap & plantQuant >= 0 , mc, max(mc,na.rm=TRUE) * 1e3)


    names(mc) <- object@labels[[1]]

    return(mc)
  })


## Create a method to recover marginal cost using
## demand parameters and supplied prices
#'@rdname Cost-Methods
#'@export
setMethod(
  f= "calcMC",
  signature= "Auction2ndLogit",
  definition= function(object,preMerger=TRUE,exAnte=FALSE){

    prices <- object@prices

    marginPre <- calcMargins(object, preMerger = TRUE, level = TRUE)

    mc <-  prices - marginPre

    if(!preMerger){
      mc <- mc + object@mcDelta
    }

    if(exAnte){mc <- mc * calcShares(object,preMerger=preMerger)}

    names(mc) <- object@labels

    mc <- as.vector(mc)

    isNegMC <- mc < 0

    if( preMerger && any(isNegMC, na.rm = TRUE)){

      warning(paste("Negative marginal costs were calibrated for the following firms:", paste(object@labels[isNegMC], collapse=",")))

    }

    return(mc)
  }
)

#'@rdname Cost-Methods
#'@export
setMethod(
  f= "calcdMC",
  signature= "Stackelberg",
  definition=function(object,preMerger=TRUE){

    if(preMerger){

      quantity  <- object@quantityPre
      dmcfun <- object@dmcfunPre
    }
    else{

      quantity  <- object@quantityPost
      dmcfun <- object@dmcfunPost
    }




    nplants <- nrow(quantity)

    dmc <- rep(NA, nplants)

    for(f in 1:nplants){
      dmc[f] <- dmcfun[[f]](quantity[f,])
    }

    if(!preMerger){dmc <- dmc*(1 + object@mcDelta)}

    names(dmc) <- object@labels[[1]]

    return(dmc)
  })


#'@rdname Cost-Methods
#'@export
setMethod(
  f= "calcVC",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){

    if(preMerger){

      quantity  <- object@quantityPre
      vcfun <- object@vcfunPre
    }
    else{

      quantity  <- object@quantityPost
      vcfun <- object@vcfunPost
    }




    nplants <- nrow(quantity)

    vc <- rep(NA, nplants)

    for(f in 1:nplants){
      vc[f] <- vcfun[[f]](quantity[f,])
    }

    if(!preMerger){vc <- vc*(1 + object@mcDelta)}

    names(vc) <- object@labels[[1]]

    return(vc)
  })

