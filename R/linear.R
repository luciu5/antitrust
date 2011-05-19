#setwd("h:/TaragIC/AntiTrustRPackage/antitrust/R")
#source("pclinear.R")






setClass(
         Class = "Linear",
         contains="Antitrust",
         representation=representation(
         prices           = "vector",
         quantities       = "vector",
         margins          = "vector",
         mc               = "vector",
         pricePre         = "vector",
         pricePost        = "vector",
         diversion        = "matrix"
         ),
          prototype=prototype(

          pricePre      =  vector(),
          pricePost     =  vector()


        ),
         validity=function(object){

             ## Sanity Checks

             nprods <- length(object@shares) # count the number of products


             if(nprods != length(object@quantities) ||
                nprods != length(object@margins) ||
                nprods != length(object@prices)){
                 stop("'prices', 'quantities', 'margins', and 'shares' must all be vectors with the same length")}

             if(any(object@prices<0,na.rm=TRUE))             stop("'prices' values must be positive")
             if(any(object@quantities<0,na.rm=TRUE))          stop("'quantities' values must be positive")
             if(any(object@margins<0 | object@margins>1,na.rm=TRUE)) stop("'margins' values must be between 0 and 1")

             if(!all(diag(object@diversion) == 1)){ stop("'diversion' diagonal elements must all equal 1")}
             if(any(abs(object@diversion) > 1)){ stop("'diversion' elements must be between -1 and 1")}

              if(nprods != nrow(object@diversion) ||
                 nprods != ncol(object@diversion)){
                  stop("'diversions' must be a square matrix")
              }

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




setMethod(
 f= "calcPrices",
 signature= "Linear",
 definition=function(object,preMerger=TRUE){

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

     prices <- solve( slopes +  t(slopes*owner))
     prices <- prices %*% (crossprod(slopes * owner,object@mc * (1+mcDelta)) -  intercept)
     prices <- as.vector(prices)
     names(prices) <- object@labels

     return(prices)


 }
          )



setMethod(
 f= "calcQuantities",
 signature= "Linear",
 definition=function(object,preMerger=TRUE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

      slopes    <- object@slopes[,-1]
      intercept <- object@slopes[,1]


     quantities <- as.vector(intercept+slopes %*% prices)
     names(quantities) <- object@labels

     return(quantities)

}
 )




setMethod(
 f= "elast",
 signature= "Linear",
 definition=function(object,preMerger=TRUE){

       if(preMerger){ prices <- object@pricePre}
       else{          prices <- object@pricePost}

      slopes    <- object@slopes[,-1]


      quantities <-  calcQuantities(object,preMerger)

      elast <- slopes * tcrossprod(1/quantities,prices)
      dimnames(elast) <- list(object@labels,object@labels)

      return(elast)

}
 )




linear <- function(prices,quantities,margins, diversions=NULL,
                     ownerPre,ownerPost,
                     mcDelta=rep(0,length(prices)),
                     labels=paste("Prod",1:length(prices),sep="")
                     ){

    shares <- quantities/sum(quantities)

    if(is.null(diversions)){
        diversions <-  diversion <-  tcrossprod(1/(1-shares),shares)
        diag(diversion) <- 1
    }


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


