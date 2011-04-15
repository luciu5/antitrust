##source("pclinear.R")


setClass(
         Class   = "CES",
         contains="PCLinear",

         representation=representation(
         priceStart    = "vector",
         shareInside   = "numeric",
         totExp        = "numeric",
         normIndex     = "numeric"
         ),

         prototype=prototype(
         totExp      =  numeric()
         ),

         validity=function(object){

             ## Sanity Checks


             nprods <- length(object@shares)


             if(nprods != length(object@priceStart)){
                 stop("'priceStart' must have the same length as 'shares'")}

         }

         )


setMethod(
          f= "calcSlopes",
          signature= "CES",
          definition=function(object){

              ## Uncover Demand Coefficents


              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              margins      <-  object@margins
              prices       <-  object@prices
              quantities   <-  object@quantities
              idx          <-  object@normIndex
              shareInside  <-  object@shareInside

              ## uncover Numeraire Coefficients
              if(length(shareInside)>0) {alpha <- 1/shareInside -1}
              else{alpha <- NULL}



              ## Uncover price coefficient and mean valuation from margins and revenue shares


              nprods <- length(shares)


              nMargins <-  length(margins[!is.na(margins)])

              ## Minimize the distance between observed and predicted margins
              minD <- function(gamma){


                  elast <- (gamma - 1 ) * matrix(shares,ncol=nprods,nrow=nprods)
                  diag(elast) <- -1*gamma - diag(elast)

                  marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% shares) / shares

                  measure <- sqrt(sum((margins - marginsCand)^2,na.rm=TRUE))/sqrt(nMargins)

                  return(measure)
              }

              minGamma <- optimize(minD,c(1,1e12))$minimum


              meanval <- log(shares) - log(shares[idx]) + (minGamma - 1) * (log(prices) - log(prices[idx]))
              meanval <- exp(meanval)

              names(meanval)   <- object@labels

              object@slopes    <- list(alpha=alpha,gamma=minGamma,meanval=meanval)
              object@totExp    <- (1 + alpha) * sum(prices*quantities)

              return(object)
          }
          )




setMethod(
 f= "calcPrices",
 signature= "CES",
 definition=function(object,preMerger=TRUE,...){

     require(nleqslv) #needed to solve nonlinear system of firm FOC

     if(preMerger){
         mcDelta <- rep(0,length(object@mcDelta))
         owner <- object@ownerPre
     }

     else{
         mcDelta <- object@mcDelta
         owner <- object@ownerPost
     }


     nprods <- length(object@shares)

     gamma    <- object@slopes$gamma
     meanval  <- object@slopes$meanval


     ##Define system of FOC as a function of prices
     FOC <- function(priceCand){

         margins <- 1 - (object@mc * (1 + mcDelta)) / priceCand
         shares <- meanval*priceCand^(1-gamma)
         shares <- shares/sum(shares)

         elast <- (gamma - 1 ) * matrix(shares,ncol=nprods,nrow=nprods)
         diag(elast) <- -1*gamma - diag(elast)

         thisFOC <- shares + as.vector((elast * owner) %*% (margins * shares))

         return(thisFOC)
     }

     ## Find price changes that set FOCs equal to 0
     minResult <- nleqslv(object@priceStart,FOC,...)

     priceEst        <- minResult$x
     names(priceEst) <- object@labels

     return(priceEst)

 }
 )


setMethod(
 f= "calcShares",
 signature= "CES",
 definition=function(object,preMerger=TRUE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

     gamma    <- object@slopes$gamma
     meanval  <- object@slopes$meanval

     shares <- meanval*prices^(1-gamma)
     shares <- shares/sum(shares)

     names(shares) <- object@labels

     return(shares)

}
 )

setMethod(
 f= "calcQuantities",
 signature= "CES",
 definition=function(object,preMerger=TRUE){

     alpha       <- object@slopes$alpha
     totExp      <-  object@totExp

     if(is.null(alpha)) stop("'shareInside' must be provided in order to calculate Compensating Variation")

     quantities        <- calcShares(object,preMerger)*totExp / (1 + alpha)
     names(quantities) <- object@labels
     return(quantities)
 }

)

setMethod(
 f= "elast",
 signature= "CES",
 definition=function(object,preMerger=TRUE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

     gamma    <- object@slopes$gamma

     shares <-  calcShares(object,preMerger)
     nprods <-  length(shares)

     elast <- (gamma - 1 ) * matrix(shares,ncol=nprods,nrow=nprods,byrow=TRUE)
     diag(elast) <- -1*gamma - diag(elast)

     dimnames(elast) <- list(object@labels,object@labels)

      return(elast)

}
 )




setMethod(
 f= "diversion",
 signature= "CES",
 definition=function(object,preMerger=TRUE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

    elasticity <- elast(object,preMerger)
    shares     <- calcShares(object,preMerger)

    diversion <- -1 * t(elasticity) / diag(elasticity)
    diag(diversion) <- 1
    diversion <- diversion * tcrossprod(1/shares,shares) * tcrossprod(prices,1/prices)

    return(diversion)
}
 )


setMethod(
          f= "CV",
          signature= "CES",
          definition=function(object){

              alpha       <- object@slopes$alpha

              if(is.null(alpha)) stop("'shareInside' must be provided in order to calculate Compensating Variation")


              gamma       <- object@slopes$gamma
              meanval     <- object@slopes$meanval
              shareInside <- object@shareInside
              totExp      <- object@totExp
              price       <- object@prices
              quantities  <- object@quantities



              tempPre  <- sum(meanval * object@pricePre^(1-gamma))
              tempPost <- sum(meanval * object@pricePost^(1-gamma))

              CV <- exp(log(totExp) +  log(tempPre/tempPost)/((1+alpha) * (gamma-1))) - totExp

              return(CV)

 })



setMethod(
 f= "summary",
 signature= "CES",
 definition=function(object){

     alpha       <- object@slopes$alpha

     if(is.null(alpha)){
         pricePre  <- object@pricePre
         pricePost <- object@pricePost
         priceDelta <- (pricePost - pricePre)/pricePre *100
         sharesPre <-   calcShares(object,TRUE) *100
         sharesPost <-  calcShares(object,FALSE) *100
         sharesDelta <- (sharesPost - sharesPre)/sharesPre *100

         results <- data.frame(pricePre=pricePre,pricePost=pricePost,
                               priceDelta=priceDelta,sharesPre=sharesPre,
                               sharesPost=sharesPost,sharesDelta=sharesDelta)

         rownames(results) <- object@labels

         cat("\nMerger Simulation Results (Deltas are Percent Changes):\n\n")
         print(round(results,2))
     }
     else{   callNextMethod(object)}

}
 )



ces <- function(prices,quantities,margins,
                ownerPre,ownerPost,
                normIndex=1,
                mcDelta=rep(0,length(prices)),
                shareInside = NULL,
                priceStart = prices,
                labels=paste("Prod",1:length(prices),sep=""),
                ...
                ){



    revenues <- prices*quantities
    shares=(revenues)/sum(revenues) #Revenue based shares
    if(is.null(shareInside)){shareInside <- numeric()}

    ## Create CES  container to store relevant data
    result <- new("CES",prices=prices, quantities=quantities,margins=margins,
                  normIndex=normIndex,
                  shares=shares,mc=prices*(1-margins),mcDelta=mcDelta,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  priceStart=priceStart,
                  shareInside=shareInside,labels=labels)

    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate Demand Slope Coefficients
    result <- calcSlopes(result)

    ## Solve Non-Linear System for Price Changes
    result@pricePre  <- calcPrices(result,TRUE,...)
    result@pricePost <- calcPrices(result,FALSE,...)

    return(result)

}

