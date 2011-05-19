## DEPRECATED: replaced by loglog.R

##source("pclinear.R")


setClass(
         Class = "PCLogLog",
         contains="PCLinear",
         representation=representation(
         priceStart = "vector"
         ),
          validity=function(object){

             ## Sanity Checks


             nprods <- length(object@shares)

             if(nprods != length(object@priceStart)){
                 stop("'priceStart' must have the same length as 'quantities'")}
         }
 )


setMethod(
 f= "calcSlopes",
 signature= "PCLogLog",
 definition=function(object){

     shares  <- object@shares
     margins <- object@margins
     quantities <- object@quantities
     prices <- object@prices

     revenues <- prices * quantities

     nprods <- length(shares)

     diversion <-  matrix(shares/(1-shares),ncol=nprods,nrow=nprods)

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
 signature= "PCLogLog",
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
     share  <- price*quantity/sum(price*quantity)

     thisFOC <- share + (t(slopes) * owner ) %*% (margin * share)

     return(as.vector(thisFOC))
 }

 price <- nleqslv(object@priceStart,FOC,object=object,...)$x

 return(price)

}
 )


setMethod(
 f= "calcQuantities",
 signature= "PCLogLog",
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
 signature= "PCLogLog",
 definition=function(object,margin=TRUE){

      elast    <- object@slopes[,-1]

      dimnames(elast) <- list(object@labels,object@labels)

      return(elast)

}
 )




pcloglog <- function(prices,quantities,margins,
                     shares=(prices*quantities)/sum(prices*quantities), #Revenue based shares
                     ownerPre,ownerPost,
                     mcDelta=rep(0,length(prices)),
                     priceStart=prices,
                     labels=paste("Prod",1:length(prices),sep=""),...
                     ){



     result <- new("PCLogLog",prices=prices, quantities=quantities,margins=margins,
                   shares=shares,mc=prices*(1-margins),mcDelta=mcDelta, priceStart=priceStart,
                   ownerPre=ownerPre,
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

