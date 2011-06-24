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
         diversion        = "matrix",
         symmetry        = "logical"
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

             if(length(object@symmetry)!=1){stop("'symmetry' must equal 'TRUE' or 'FALSE'")}

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
     ownerPre   <- object@ownerPre
     symmetry  <- object@symmetry

     nprods <- length(margins)

     diag(diversion) <- -1


      if(!symmetry){


          slopes <- matrix(margins * prices,ncol=nprods, nrow=nprods,byrow=TRUE)
          slopes <- 1/rowSums(slopes * diversion * ownerPre) * quantities
          slopes <- -t(slopes * diversion)


      }

     else{
         nMargins <-  length(margins[!is.na(margins)])
         revenues <- prices*quantities

          minD <- function(s){
              beta <- -s*diversion[1,]/diversion[,1] ## wlog, assume that b11 is unknown
              slopesCand <- matrix(beta,ncol=nprods,nrow=nprods,byrow=TRUE)
              slopesCand <- slopesCand*t(diversion)
              elast <- slopesCand * tcrossprod(1/quantities,prices)

              marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% revenues) / revenues

              measure <- sqrt(sum((margins - marginsCand)^2,na.rm=TRUE))/sqrt(nMargins)
              return(measure)
          }


         minSlope <- optimize(minD,c(-1e12,0))$minimum

         slopes <-  -minSlope*diversion[1,]/diversion[,1]
         slopes <- matrix(slopes,ncol=nprods,nrow=nprods,byrow=TRUE)
         slopes <- slopes*t(diversion)
     }




     dimnames(slopes) <- list(object@labels,object@labels)


     intercept <- as.vector(quantities - slopes %*% prices)

     if(!symmetry &&
        !isTRUE(all.equal(slopes[upper.tri(slopes)],slopes[lower.tri(slopes)]))){
         warn("Matrix of demand slopes coefficients is not symmetric. Demand parameters may not be consistent with utility maximization theory.")}

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
 f= "calcShares",
 signature= "Linear",
 definition=function(object,preMerger=TRUE){

     quantities <- calcQuantities(object,preMerger)

     return(quantities/sum(quantities))
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

setMethod(
 f= "CV",
 signature= "Linear",
 definition=function(object){

     slopes    <- object@slopes[,-1]

      if(!isTRUE(all.equal(slopes[upper.tri(slopes)],slopes[lower.tri(slopes)]))){
                  stop("price coefficient matrix must be symmetric in order to calculate compensating variation. Suggest setting 'symmetry=TRUE'")
              }

     intercept <- object@slopes[,1]
     pricePre  <- object@pricePre
     pricePost <- object@pricePost

     result <- sum(intercept*(pricePost-pricePre)) + .5 * as.vector(pricePost%*%slopes%*%pricePost - pricePre%*%slopes%*%pricePre)

     return(result)
 })


linear <- function(prices,quantities,margins, diversions=NULL, symmetry=TRUE,
                     ownerPre,ownerPost,
                     mcDelta=rep(0,length(prices)),
                     labels=paste("Prod",1:length(prices),sep="")
                     ){

    shares <- quantities/sum(quantities)

    if(is.null(diversions)){
        diversions <- tcrossprod(1/(1-shares),shares)
        diag(diversions) <- 1
    }


     result <- new("Linear",prices=prices, quantities=quantities,margins=margins,
                   shares=shares,mc=prices*(1-margins),mcDelta=mcDelta,
                   ownerPre=ownerPre,diversion=diversions, symmetry=symmetry,
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


