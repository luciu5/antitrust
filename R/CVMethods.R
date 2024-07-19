#' @title Methods For Calculating Compensating Variation (CV)
#' @name CV-Methods
#' @docType methods
#'
#' @aliases CV-methods
#' CV
#' CV,ANY-method
#' CV,AIDS-method
#' CV,CES-method
#' CV,CESNests-method
#' CV,Linear-method
#' CV,LogLin-method
#' CV,Logit-method
#' CV,LogitNests-method
#' CV,Auction2ndLogit-method
#' CV,VertBargBertLogit-method
#' CV,VertBarg2ndLogit-method
#' CV,Cournot-method
#'
#' @description Calculate the amount of money a consumer would need to
#' be paid to be just as well off as they were before the merger.
#'
#' @description
#' All the information needed to
#' compute CV is already available within the Logit, Nested Logit CES and Nested CES classes.
#' In CES and Nested CES, CV cannot be calculated if the outside share cannot be inferred.
#'
#' For AIDS, if the  \sQuote{insideSize}  slot to the \dQuote{AIDS} class equals NA, CV is calculated as a percentage of
#' total expenditure (revenues) on products included in the simulation. Otherwise CV is calculated in terms of dollars.
#' Pre-merger prices for all products in the market must be supplied in order for CV to be calculated.
#'
#' For Linear and LogLin, although no additional information is needed to calculate CV for
#' either the \dQuote{Linear} or \dQuote{LogLin} classes, The CV method will fail if
#' the appropriate restrictions on the demand parameters have not been imposed.
#'
#' @param object An instance of one of the classes listed above.
#'
#' @include CMCRMethods.R
#' @keywords methods
NULL

setGeneric (
  name= "CV",
  def=function(object,...){standardGeneric("CV")}
)

#'@rdname CV-Methods
#'@export
setMethod(
  f= "CV",
  signature= "Cournot",
  definition=function(object){

    demand <- object@demand
    slopes    <- object@slopes
    intercepts <- object@intercepts
    quantityPre <- colSums(object@quantityPre, na.rm=TRUE)
    quantityPost <- colSums(object@quantityPost, na.rm=TRUE)
    pricePre <- object@pricePre
    pricePost <- object@pricePost

    result <- ifelse(demand =="linear",
                     .5*(pricePost - pricePre)*(quantityPre - quantityPost)  ,
                     exp(intercepts)/(slopes + 1) * (quantityPre^(slopes + 1) - quantityPost^(slopes+1)) -  (quantityPre - quantityPost)* pricePre
    )

    result <- result  +  (pricePost - pricePre)*quantityPost
    names(result) <-  object@labels[[2]]
    return(result)
  })


#'@rdname CV-Methods
#'@export
setMethod(
  f= "CV",
  signature= "Linear",
  definition=function(object){

    slopes    <- object@slopes

    if(!isTRUE(all.equal(slopes,t(slopes)))){
      stop("price coefficient matrix must be symmetric in order to calculate compensating variation. Suggest setting 'symmetry=TRUE'")
    }

    intercept <- object@intercepts
    pricePre  <- object@pricePre
    pricePost <- object@pricePost

    result <- sum(intercept*(pricePost-pricePre)) + .5 * as.vector(pricePost%*%slopes%*%pricePost - pricePre%*%slopes%*%pricePre)

    return(result)
  })

#'@rdname CV-Methods
#'@export
setMethod(
  f= "CV",
  signature= "Logit",
  definition=function(object){

    alpha       <- object@slopes$alpha
    meanval     <- object@slopes$meanval
    subset <- object@subset
    mktSize <- object@mktSize

    # outVal <- ifelse(object@shareInside<1, 1, 0)
    outVal <- ifelse(is.na(object@normIndex), 1, 0)
    output <- ifelse(object@output,1,-1)

    VPre  <- sum(exp(meanval + (object@pricePre - object@priceOutside)*alpha))  + outVal
    VPost <- sum(exp(meanval + (object@pricePost - object@priceOutside)*alpha)[subset] ) + outVal

    result <- output*log(VPost/VPre)/alpha

    if(!is.na(mktSize)){ result <- result * mktSize}

    return(unname(result))

  })

#'@rdname CV-Methods
#'@export
setMethod(
  f= "CV",
  signature= "LogLin",
  definition=function(object){
    stop("CV method not currently available for 'LogLin' Class")

  })


#'@rdname CV-Methods
#'@export
setMethod(
  f= "CV",
  signature= "AIDS",
  definition=function(object){

    ## computes compensating variation using the closed-form solution found in LaFrance 2004, equation 10

    if(any(is.na(object@prices))){stop("Compensating Variation cannot be calculated without supplying values to 'prices'")}



    slopes <- object@slopes
    intercepts <- object@intercepts
    pricePre <- log(object@pricePre)
    pricePost <- log(object@pricePost)
    insideSizePre <- object@insideSize

    sharePost <- calcShares(object, preMerger = FALSE, revenue = TRUE)

    # source: epstein and rubinfeld Assumption C (appendix)
    insideSizeDelta <- sum( sharePost * object@priceDelta) * (object@mktElast + 1)

    insideSizePost <- insideSizePre*( 1 + insideSizeDelta)

    result <- sum(intercepts*(pricePost-pricePre)) + .5 * as.vector(t(pricePost)%*%slopes%*%pricePost - t(pricePre)%*%slopes%*%pricePre)


    if(is.na(insideSizePost)){
      warning("Slot 'insideSize' is missing. Calculating CV as a percentage change in (aggregate) income")
      return(result*100)}

    else{
      return(insideSizePost*(exp(result)-1))
    }

  }
)

#'@rdname CV-Methods
#'@export
setMethod(
  f= "CV",
  signature= "LogitNests",
  definition=function(object){

    nests       <- object@nests
    alpha       <- object@slopes$alpha
    sigma       <- object@slopes$sigma
    meanval     <- object@slopes$meanval
    mktSize <- object@mktSize

    outVal <- ifelse(is.na(object@normIndex), exp(alpha*object@priceOutside), 0)

    VPre  <- sum( tapply(exp((meanval + object@pricePre*alpha)  / sigma[nests]),nests,sum,na.rm=TRUE) ^ sigma ) + outVal
    VPost <- sum( tapply(exp((meanval + object@pricePost*alpha) / sigma[nests]),nests,sum,na.rm=TRUE) ^ sigma ) + outVal

  
    result <- log(VPost/VPre)/alpha
    if(!is.na(mktSize)){ result <- result * mktSize}
    names(result) <- NULL
    return(result)

  })

#'@rdname CV-Methods
#'@export
setMethod(
  f= "CV",
  signature= "Auction2ndLogit",
  definition=function(object){

    mktSize = object@mktSize
    alpha   = object@slopes$alpha
    meanvalPre = object@slopes$meanval
    idx <- object@normIndex
    subset <- object@subset
    output <- object@output
    outSign <- ifelse(output,-1,1)
    
    mcDelta <- object@mcDelta
    
    if (is.na(idx)) {
      outVal <- 1
      mcDeltaOut <- object@priceOutside
    }
    else {
      outVal <- 0  
      mcDeltaOut <- object@mcDelta[idx]
    }
    
    meanvalPost =  meanvalPre + alpha * (object@mcDelta - mcDeltaOut)
    
    
    marginPre <- calcMargins(object, preMerger = TRUE, exAnte= TRUE)
    marginPost <- calcMargins(object, preMerger = FALSE, exAnte= TRUE)

    ## Add the elimination of first best option
    VPre  <- sum(exp(meanvalPre),na.rm=TRUE)  + outVal
    VPost <- sum(exp(meanvalPost[subset]),na.rm=TRUE ) + outVal
    
    
    
    result <-   outSign * log(VPost/VPre)/alpha  
    result <- result + outSign*(sum(marginPost,na.rm=TRUE) - sum(marginPre, na.rm =TRUE))
    
    if(!is.na(mktSize)){result <- mktSize * result}

    return(outSign*result)
  })


#'@rdname CV-Methods
#'@export
setMethod(
  f= "CV",
  signature= "VertBargBertLogit",
  definition=function(object){
    
    down <- object@down
    logitCV <- selectMethod("CV",class(down))
    
    return(logitCV(down))
  }
  )


#'@rdname CV-Methods
#'@export
setMethod(
  f= "CV",
  signature= "VertBarg2ndLogit",
  definition=function(object){
    
    down <- object@down
    mktSize = down@mktSize

alpha = down@slopes$alpha
meanvalPre = down@slopes$meanval

priceDelta <- calcPriceDelta(object,levels=TRUE)

marginsPre <- calcMargins(object,preMerger=TRUE,level=TRUE)
marginsPost <- calcMargins(object,preMerger=FALSE,level=TRUE)

sharesPre <- calcShares(object, preMerger=TRUE, revenue=FALSE)
sharesPost <- calcShares(object, preMerger=FALSE, revenue=FALSE)

result <- sum(marginsPost$down*sharesPost,na.rm=TRUE) - sum(marginsPre$down*sharesPre,na.rm=TRUE)

idx <- down@normIndex
subset <- down@subset
if (is.na(idx)) {
  outVal <- 1
  mcDeltaOut <- down@priceOutside
}

else {
  outVal <- 0
  mcDeltaOut <- down@mcDelta[idx]
}


meanvalPost = meanvalPre +  alpha * (priceDelta$up + down@mcDelta - 
                                      mcDeltaOut)


VPre <- sum(exp(meanvalPre), na.rm = TRUE) + outVal
VPost <- sum(exp(meanvalPost[subset]), na.rm = TRUE) + outVal
result <- result + log(VPost/VPre)/alpha

if (!is.na(mktSize)) {
  result <- mktSize * result
}
return(result)
}
)


#'@rdname CV-Methods
#'@export
setMethod(
  f= "CV",
  signature= "CES",
  definition=function(object){

    alpha       <- object@slopes$alpha
    mktSize  <- object@mktSize

    if(is.null(alpha)) stop(" Sum of 'shares' must be less than 1  calculate Compensating Variation")

    gamma       <- object@slopes$gamma
    meanval     <- object@slopes$meanval
    output <- ifelse(object@output,1,-1)
    

    outVal <- ifelse(is.na(object@normIndex), 1, 0)

    VPre  <- sum(meanval * (object@pricePre / object@priceOutside)^(1-gamma),na.rm=TRUE) + outVal
    VPost <- sum(meanval * (object@pricePost/ object@priceOutside)^(1-gamma),na.rm=TRUE) + outVal

    result <- output*log(VPost/VPre) / ((1+alpha)*(1-gamma))

    result <- exp(result) - 1

    if(is.na(mktSize)){
      warning("'insideSize' is NA. Calculating CV as a percentage of (aggregate) expenditure")
      return(result*100)}

    else{

      return(mktSize*result)
    }


  })

#'@rdname CV-Methods
#'@export
setMethod(
  f= "CV",
  signature= "CESNests",
  definition=function(object){

    alpha       <- object@slopes$alpha
    mktSize  <- object@mktSize


    if(is.null(alpha)) stop("'shareInside' must be between 0 and 1 to  calculate Compensating Variation")

    nests       <- object@nests
    gamma       <- object@slopes$gamma
    sigma       <- object@slopes$sigma
    meanval     <- object@slopes$meanval
    shareInside <- object@shareInside



    VPre  <- sum(tapply(meanval *  object@pricePre^(1-sigma[nests]),nests,sum,na.rm=TRUE) ^((1-gamma)/(1-sigma)))
    VPost <- sum(tapply(meanval * object@pricePost^(1-sigma[nests]),nests,sum,na.rm=TRUE) ^((1-gamma)/(1-sigma)))

    ##tempPre  <- log( sum( tapply(meanval * object@pricePre^(1-sigma[nests]),nests,sum) ^((1-gamma)/(1-sigma)) ) )
    ##tempPre  <- (gamma/(gamma-1)) * log( sum( tapply( meanval^((1-gamma)/(1-sigma[nests])) * object@pricePre^(-gamma),nests,sum)^((gamma-1)/gamma)) ) - tempPre

    ##tempPost  <- log( sum( tapply(meanval * object@pricePost^(1-sigma[nests]),nests,sum) ^((1-gamma)/(1-sigma)) ) )
    ##tempPost  <- (gamma/(gamma-1)) * log( sum( tapply( meanval^((1-gamma)/(1-sigma[nests])) * object@pricePost^(-gamma),nests,sum)^((gamma-1)/gamma)) ) - tempPost


    result <- log(VPost/VPre) / ((1+alpha)*(1-gamma))
    names(result) <- NULL

    if(is.na(mktSize)){
      warning("'insideSize' is NA. Calculating CV as a percentage of (aggregate) expenditure")
      return(result*100)}

    else{
      totExp <- mktSize*(1+alpha)
      return(totExp*(exp(result)-1))
    }

  })


