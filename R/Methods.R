## Generate a bunch of generic functions

setGeneric (
 name= "ownerToMatrix",
 def=function(object,...){standardGeneric("ownerToMatrix")}
 )

setGeneric (
 name= "calcSlopes",
 def=function(object){standardGeneric("calcSlopes")}
 )

setGeneric (
 name= "calcQuantities",
 def=function(object,...){standardGeneric("calcQuantities")}
 )

setGeneric (
 name= "calcShares",
 def=function(object,...){standardGeneric("calcShares")}
 )


setGeneric (
name= "calcMargins",
 def=function(object,...){standardGeneric("calcMargins")}
 )

setGeneric (
 name= "calcPrices",
 def=function(object,...){standardGeneric("calcPrices")}
 )

setGeneric (
 name= "calcPriceDelta",
 def=function(object,...){standardGeneric("calcPriceDelta")}
 )

setGeneric (
 name= "elast",
 def=function(object,...){standardGeneric("elast")}
 )

setGeneric (
 name= "diversion",
 def=function(object,...){standardGeneric("diversion")}
 )

setGeneric (
 name= "CV",
 def=function(object,...){standardGeneric("CV")}
 )

setGeneric (
 name= "cmcr",
 def=function(object,...){standardGeneric("cmcr")}
 )


setGeneric (
 name= "getNestsParms",
 def=function(object,...){standardGeneric("getNestsParms")}
 )

setGeneric (name= "summary")


## Create some methods for the Antitrust Class
setMethod(
 f= "show",
 signature= "Antitrust",
 definition=function(object){

     pricePre  <- object@pricePre
     pricePost <- object@pricePost

     priceDelta <- pricePost/pricePre - 1
     names(priceDelta) <- object@labels
     print(priceDelta*100)

}
 )


setMethod(
 f= "summary",
 signature= "Antitrust",
 definition=function(object){

     if(hasMethod("calcQuantities",class(object))){
         outPre  <-  calcQuantities(object,TRUE)
         outPost <-  calcQuantities(object,FALSE)
     }

     else{
         outPre  <-  calcShares(object,TRUE) * 100
         outPost <-  calcShares(object,FALSE) * 100
     }


     outDelta <- (outPost/outPre - 1) * 100

     pricePre   <-  object@pricePre
     pricePost  <-  object@pricePost
     priceDelta <- (pricePost/pricePre - 1) * 100

     results <- data.frame(pricePre=pricePre,pricePost=pricePost,
                           priceDelta=priceDelta,outputPre=outPre,
                           outputPost=outPost,outputDelta=outDelta)

     rownames(results) <- object@labels

     sharesPost <- calcShares(object,FALSE)

     cat("\nMerger Simulation Results (Deltas are Percent Changes):\n\n")
     print(round(results,2))
     cat("\n\nShare-Weighted Price Change:",round(sum(sharesPost*priceDelta),2),sep="\t")
     cat("\nShare-Weighted CMCR:",round(sum(cmcr(object)*sharesPost),2),sep="\t")
     cat("\n\n")
     return(invisible(results))

 })


setMethod(
 f= "ownerToMatrix",
signature= "Antitrust",
definition=function(object,preMerger=TRUE){


    ## transform ownerPre and ownerPost vectors into matrices, when applicable

    if(preMerger) thisOwner <- object@ownerPre
    else{         thisOwner <- object@ownerPost}



    if(is.vector(thisOwner) || is.factor(thisOwner)){

        nprod <- length(object@shares)
        owners <- as.numeric(factor(thisOwner))
        thisOwner <- matrix(0,ncol=nprod,nrow=nprod)


        for( o in unique(owners)){
            thisOwner [owners == o, owners == o] = 1
        }


    }


    return(thisOwner)

}
)



setMethod(
 f= "diversion",
 signature= "Antitrust",
 definition=function(object,preMerger=TRUE){

     if(hasMethod("calcQuantities",class(object))){
         output  <-  calcQuantities(object,preMerger)
     }

     else{
         output  <-  calcShares(object,preMerger)
     }

    elasticity <- elast(object,preMerger)

    diversion <- -1 * t(elasticity) / diag(elasticity)
    diag(diversion) <- 1
    diversion <- diversion * tcrossprod(1/output,output)

    return(diversion)
}
 )


setMethod(
 f= "cmcr",
 signature= "Antitrust",
 definition=function(object){

    pricePre <- object@pricePre
    pricePre <- tcrossprod(1/pricePre,pricePre)

    elastPre <- elast(object,TRUE)
    divPre <- diversion(object,TRUE)

    Bpre =  -1 * divPre * pricePre * object@ownerPre;  diag(Bpre) = 1
    Bpost = -1 * divPre * pricePre * object@ownerPost; diag(Bpost) = 1

    marginPre <- -1 * as.vector(solve(Bpre)  %*% (1/diag(elastPre)))
    marginPost <-     as.vector(solve(Bpost) %*%  Bpre %*% marginPre)

    cmcr <- (marginPost - marginPre)/(1 - marginPre)
    names(cmcr) <- object@labels

    return(cmcr * 100)
}
 )
