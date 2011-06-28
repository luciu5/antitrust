setClass(
         Class = "LogLog",
         contains="Linear",
         representation=representation(
         priceStart = "vector"
         ),

          validity=function(object){

             ## Sanity Checks


             nprods <- length(object@prices)

             if(nprods != length(object@priceStart)){
                 stop("'priceStart' must have the same length as 'prices'")}
 })


setMethod(
 f= "calcSlopes",
 signature= "LogLog",
 definition=function(object){



     margins <- object@margins
     quantities <- object@quantities
     prices <- object@prices

     revenues <- prices * quantities

     nprods <- length(margins)

     diversion <- object@diversion * tcrossprod(quantities,1/quantities)
     diag(diversion) <- -1

    slopes <- matrix(margins * revenues,ncol=nprods, nrow=nprods,byrow=TRUE)
    slopes <- revenues / rowSums(slopes * diversion * object@ownerPre)
    slopes <- -t(slopes * diversion)

    dimnames(slopes) <- list(object@labels,object@labels)

    intercept <- as.vector(log(quantities) - slopes %*% log(prices))

     result <- cbind(Intercept=intercept,slopes)

     return(result)


 }
          )


setMethod(
 f= "calcPrices",
 signature= "LogLog",
 definition=function(object,preMerger=TRUE,...){
 if(preMerger){
     mcDelta <- rep(0,length(object@mcDelta))
     owner <- object@ownerPre
 }
 else{
     mcDelta <- object@mcDelta
     owner <- object@ownerPost
 }

 slopes    <- object@slopes[,-1]
 intercept <- object@slopes[,1]




 FOC <- function(price,object){

     require(nleqslv)

     quantity <- exp(intercept) * apply(price^slopes,1,prod) # log(price) can produce errors.
                                                             # this transformation avoidsusing logs
    #quantity <- exp(as.vector(intercept+slopes %*% log(price)))

     margin   <- 1 - (object@mc * (1 + mcDelta)) / price
     revenue  <- price*quantity

     thisFOC <- revenue + (t(slopes) * owner ) %*% (margin * revenue)

     return(as.vector(thisFOC))
 }

 minResult <- nleqslv(object@priceStart,FOC,object=object,...)

 if(minResult$termcd != 1){warning("'calcPrices' nonlinear solver may not have successfully converge. 'nleqslv' reports: '",minResult$message,"'")}

 priceEst        <- minResult$x
 names(priceEst) <- object@labels
 return(priceEst)

}
 )


setMethod(
 f= "calcQuantities",
 signature= "LogLog",
 definition=function(object,preMerger=TRUE,...){


     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

      slopes    <- object@slopes[,-1]
      intercept <- object@slopes[,1]

     quantities <- exp(as.vector(intercept+slopes %*% log(prices)))
     names(quantities) <- object@labels

     return(quantities)

}
 )




setMethod(
 f= "elast",
 signature= "LogLog",
 definition=function(object,margin=TRUE){

      elast    <- object@slopes[,-1]

      dimnames(elast) <- list(object@labels,object@labels)

      return(elast)

}
 )

setMethod(
 f= "CV",
 signature= "LogLog",
 definition=function(object){
     print("CV method not currently for 'LogLog' Class")
 })


loglog <- function(prices,quantities,margins,diversions=NULL,
                     ownerPre,ownerPost,
                     mcDelta=rep(0,length(prices)),
                     priceStart=prices,
                     labels=paste("Prod",1:length(prices),sep=""),...
                     ){




    shares=quantities/sum(quantities)


    if(is.null(diversions)){
        diversions <-  tcrossprod(1/(1-shares),shares)
        diag(diversions) <- 1
    }


    result <- new("LogLog",prices=prices, quantities=quantities,margins=margins,
                  shares=shares,mc=prices*(1-margins),mcDelta=mcDelta, priceStart=priceStart,
                  ownerPre=ownerPre,diversion=diversions, symmetry=FALSE,
                  ownerPost=ownerPost, labels=labels)


     ## Convert ownership vectors to ownership matrices
     result@ownerPre  <- ownerToMatrix(result,TRUE)
     result@ownerPost <- ownerToMatrix(result,FALSE)

     ## Calculate Demand Slope Coefficients
     result@slopes <- calcSlopes(result)

     ## Calculate pre and post merger equilibrium prices
     result@pricePre  <- calcPrices(result,TRUE,...)
     result@pricePost <- calcPrices(result,FALSE,...)


   return(result)

}

