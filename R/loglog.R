#source("pcloglog.R")


setClass(
         Class = "LogLog",
         contains="PCLogLog",
         representation=representation(
         priceStart = "vector",
         diversion = "matrix"
         ),
          validity=function(object){

             ## Sanity Checks


             nprods <- length(object@prices)

             if(nprods != length(object@priceStart)){
                 stop("'priceStart' must have the same length as 'quantities'")}

             if(nprods != nrow(object@diversion) ||
                nprods != ncol(object@diversion)){
                 stop("'diversions' must be a square matrix whose dimension is the same as length 'price'")
             }

             if(!all(diag(object@diversion) == 1)){ stop("'diversion' diagonal elements must all equal 1")}
             if(any(abs(object@diversion) > 1)){ stop("'diversion' elements must be between -1 and 1")}
         }

 )


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





loglog <- function(prices,quantities,margins,diversions,
                     ownerPre,ownerPost,
                     mcDelta=rep(0,length(prices)),
                     priceStart=prices,
                     labels=paste("Prod",1:length(prices),sep=""),...
                     ){


    require(BB)

    shares=rep(NA,length(prices)) #not used with "LogLog",but required by "PCLogLog" class

    result <- new("LogLog",prices=prices, quantities=quantities,margins=margins,
                  shares=shares,mc=prices*(1-margins),mcDelta=mcDelta, priceStart=priceStart,
                  ownerPre=ownerPre,diversion=diversions,
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

