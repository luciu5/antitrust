setClass(
         Class   = "CES",
         contains="Logit",

         validity=function(object){

             ## Sanity Checks


             nprods <- length(object@prices)

             if(!(object@normIndex %in% 1:nprods)){
                 stop("'normIndex' must take on a value between 1 and ",nprods)
             }

         })


setMethod(
          f= "calcSlopes",
          signature= "CES",
          definition=function(object){

              ## Uncover Demand Coefficents


              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              margins      <-  object@margins
              prices       <-  object@prices
              idx          <-  object@normIndex
              shareInside  <-  object@shareInside

              ## uncover Numeraire Coefficients
              if(shareInside < 1 ) {alpha <- 1/shareInside - 1}
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


              return(object)
          }
          )




setMethod(
 f= "calcPrices",
 signature= "CES",
 definition=function(object,preMerger=TRUE,...){


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

      if(minResult$termcd != 1){warning("'calcPrices' nonlinear solver may not have successfully converged. 'nleqslv' reports: '",minResult$message,"'")}


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
          definition=function(object,revenueInside){

              alpha       <- object@slopes$alpha

              if(is.null(alpha)) stop("'shareInside' must be less than 1 to  calculate Compensating Variation")
              if(missing(revenueInside) || isTRUE(revenueInside<0)) stop("'revenueInside' must be greater than 0 to  calculate Compensating Variation")

              totExp    <- (1 + alpha) * revenueInside

              gamma       <- object@slopes$gamma
              meanval     <- object@slopes$meanval
              shareInside <- object@shareInside



              tempPre  <- sum(meanval * object@pricePre^(1-gamma))
              tempPost <- sum(meanval * object@pricePost^(1-gamma))

              CV <- exp(log(totExp) +  log(tempPre/tempPost)/((1+alpha) * (gamma-1))) - totExp

              return(CV)

 })

setMethod(
 f= "calcPriceDeltaHypoMon",
 signature= "CES",
 definition=function(object,prodIndex,...){

     nprods <- length(prodIndex)
     pricePreOld <- object@pricePre
     gamma    <- object@slopes$gamma
     meanval  <- object@slopes$meanval


     ##Define system of FOC as a function of priceDelta
     FOC <- function(priceCand){

         thisPrice <- pricePreOld
         thisPrice[prodIndex] <- priceCand

         margins <- 1 - object@mc / thisPrice
         shares <- meanval*thisPrice^(1-gamma)
         shares <- shares/sum(shares)

         shares <- shares[prodIndex]

         elast <- (gamma - 1 ) * matrix(shares,ncol=nprods,nrow=nprods)
         diag(elast) <- -1*gamma - diag(elast)

         thisFOC <- shares + as.vector(elast  %*% (margins[prodIndex] * shares))


         return(thisFOC)

     }




     ## Find price changes that set FOCs equal to 0
     minResult <- nleqslv(object@priceStart[prodIndex],FOC,...)

     if(minResult$termcd != 1){warning("'calcPriceDeltaHypoMon' nonlinear solver may not have successfully converged. 'nleqslv' reports: '",minResult$message,"'")}

     priceDelta <- minResult$x/pricePreOld[prodIndex] - 1

     names(priceDelta) <- object@labels[prodIndex]

     return(priceDelta)

 })


ces <- function(prices,shares,margins,
                ownerPre,ownerPost,
                shareInside = 1,
                normIndex=1,
                mcDelta=rep(0,length(prices)),
                priceStart = prices,
                labels=paste("Prod",1:length(prices),sep=""),
                ...
                ){


    ## Create CES  container to store relevant data
    result <- new("CES",prices=prices, shares=shares, margins=margins,
                  normIndex=normIndex,
                  mc=prices*(1-margins),mcDelta=mcDelta,
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

