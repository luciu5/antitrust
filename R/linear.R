#setwd("h:/TaragIC/AntiTrustRPackage/antitrust/R")
#source("pclinear.R")






setClass(
         Class = "Linear",
         contains="PCLinear",
         representation=representation(
         diversion = "matrix"
         ),
         validity=function(object){

             ## Sanity Checks

             nprods <- length(object@shares) # count the number of products

              if(nprods != nrow(object@diversion) ||
                 nprods != ncol(object@diversion)){
                  stop("'diversions' must be a square matrix")
              }

             if(!all(diag(object@diversion) == 1)){ stop("'diversion' diagonal elements must all equal 1")}
             if(any(abs(object@diversion) > 1)){ stop("'diversion' elements must be between -1 and 1")}
         }
 )


setMethod(
 f= "calcSlopes",
 signature= "Linear",
 definition=function(object){

     margins    <- object@margins
     quantities <- object@quantities
     prices     <- object@prices
     diversion  <- object@diversion

     nprods <- length(margins)

     diag(diversion) <- -1


    slopes <- matrix(margins * prices,ncol=nprods, nrow=nprods,byrow=TRUE)
    slopes <- 1/rowSums(slopes * diversion * object@ownerPre) * quantities
    slopes <- -t(slopes * diversion)


     dimnames(slopes) <- list(object@labels,object@labels)

     intercept <- as.vector(quantities - slopes %*% prices)


     return(cbind(intercept,slopes))


 }
          )




linear <- function(prices,quantities,margins, diversions,
                     ownerPre,ownerPost,
                     mcDelta=rep(0,length(prices)),
                     labels=paste("Prod",1:length(prices),sep="")
                     ){

    shares=rep(NA,length(prices)) #not used with "Linear",but required by "PCLinear" class

     result <- new("Linear",prices=prices, quantities=quantities,margins=margins,
                   shares=shares,mc=prices*(1-margins),mcDelta=mcDelta,
                   ownerPre=ownerPre,diversion=diversions,
                   ownerPost=ownerPost, labels=labels)


     ## Convert ownership vectors to ownership matrices
     result@ownerPre  <- ownerToMatrix(result,TRUE)
     result@ownerPost <- ownerToMatrix(result,FALSE)

     ## Calculate Demand Slope Coefficients
     result@slopes <- calcSlopes(result)

      ## Calculate Pre and Post merger equilibrium prices
     result@pricePre  <- calcPrices(result,TRUE)
     result@pricePost <- calcPrices(result,FALSE)


   return(result)

}


