##source("pclinear.R")


setClass(
         Class   = "Logit",
         contains="Antitrust",
         representation=representation(

         prices           = "vector",
         margins          = "vector",
         mc               = "vector",
         pricePre         = "vector",
         pricePost        = "vector",
         priceStart       = "vector",
         normIndex        = "numeric",
         shareInside     = "numeric"

         ),
          prototype=prototype(

          pricePre      =  vector(),
          pricePost     =  vector(),
          mc            =  vector()

        ),


         validity=function(object){

             ## Sanity Checks

             nprods <- length(object@shares)
             sumShares <- sum(object@shares)

             if(
                 nprods != length(object@margins) ||
                 nprods != length(object@prices)){
                 stop("'prices', 'margins' and 'shares' must all be vectors with the same length")}

             if(any(object@prices<0,na.rm=TRUE))             stop("'prices' values must be positive")


             if(any(object@margins<0 | object@margins>1,na.rm=TRUE)) stop("'margins' values must be between 0 and 1")

             if(!(is.matrix(object@ownerPre) || is(object@ownerPre,"Matrix"))){
                 ownerPre <- ownerToMatrix(object,TRUE)
             }
             else{ownerPre <- object@ownerPre}

             if(all(is.na(ownerPre %*% object@margins))) stop("Insufficient margin information to calibrate demand parameters.")

             if(nprods != length(object@priceStart)){
                 stop("'priceStart' must have the same length as 'shares'")}

              if(
                !(object@shareInside >0 &&
                  object@shareInside <=1)
                ){
                 stop("'shareInside' must be between 0 and 1")
         }

              if(
                 !(all(object@shares >0) &&
                   all(object@shares <=1))
                 ){
                 stop("elements of vector 'shares' must be between 0 and 1")
         }

             if(!(
                  (sumShares == 1  && object@normIndex %in% 1:nprods) ||
                   (sumShares < 1 && is.na(object@normIndex))
                  )){
                 stop("'normIndex' must take on a value between 1 and ",nprods,
                      ". If 'shares' sum to 1 , or NA if the sum of shares is less than 1")
             }


         })


setMethod(
          f= "calcSlopes",
          signature= "Logit",
          definition=function(object){

              ## Uncover Demand Coefficents


              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              margins      <-  object@margins
              prices       <-  object@prices
              idx          <-  object@normIndex

              if(is.na(idx)){
                  idxShare <- 1 - object@shareInside
                  idxPrice <- 0
              }
              else{
                  idxShare <- shares[idx]
                  idxPrice <- prices[idx]
               }

              ## Uncover price coefficient and mean valuation from margins and revenue shares


              nprods <- length(shares)


              nMargins <-  length(margins[!is.na(margins)])

              ## Minimize the distance between observed and predicted margins
              minD <- function(alpha){


                  elast <- -alpha *  matrix(prices * shares,ncol=nprods,nrow=nprods)
                  diag(elast) <- alpha*prices - diag(elast)

                  revenues <- shares * prices
                  marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% revenues) / revenues

                  measure <- sqrt(sum((margins - marginsCand)^2,na.rm=TRUE))/sqrt(nMargins)

                  return(measure)
              }

              minAlpha <- optimize(minD,c(-1e12,0))$minimum


              meanval <- log(shares) - log(idxShare) - minAlpha * (prices - idxPrice)

              names(meanval)   <- object@labels

              object@slopes    <- list(alpha=minAlpha,meanval=meanval)


              return(object)
          }
          )



## Create a function to recover marginal cost using
## demand parameters and supplied prices
setMethod(
          f= "calcMC",
          signature= "Logit",
          definition= function(object){

              object@pricePre <- object@prices


              marginPre <- calcMargins(object,TRUE)

              mc <- (1 - marginPre) * object@prices
              names(mc) <- object@labels

              return(mc)
          }
          )


setMethod(
 f= "calcPrices",
 signature= "Logit",
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
     shareInside <- object@shareInside

     alpha    <- object@slopes$alpha
     meanval  <- object@slopes$meanval

     isOutside   <- as.numeric(shareInside < 1)

     ##Define system of FOC as a function of prices
     FOC <- function(priceCand){

         margins <- 1 - (object@mc * (1 + mcDelta)) / priceCand
         shares <- exp(meanval + alpha*priceCand)
         shares <- shares/(isOutside + sum(shares))

         elast <- -alpha * matrix(priceCand*shares,ncol=nprods,nrow=nprods)
         diag(elast) <- alpha*priceCand - diag(elast)

         revenues <-  shares * priceCand
         thisFOC <- revenues + as.vector((elast * owner) %*% (margins * revenues))

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
 signature= "Logit",
 definition=function(object,preMerger=TRUE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

     alpha    <- object@slopes$alpha
     meanval  <- object@slopes$meanval

     shares <- exp(meanval + alpha*prices)
     shares <- shares/(as.numeric(object@shareInside<1) + sum(shares))

     names(shares) <- object@labels

     return(shares)

}
 )



setMethod(
 f= "elast",
 signature= "Logit",
 definition=function(object,preMerger=TRUE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

     alpha    <- object@slopes$alpha

     shares <-  calcShares(object,preMerger)
     nprods <-  length(shares)

     elast <- -alpha  * matrix(prices*shares,ncol=nprods,nrow=nprods,byrow=TRUE)
     diag(elast) <- alpha*prices - diag(elast)

     dimnames(elast) <- list(object@labels,object@labels)

      return(elast)

}
 )



setMethod(
          f= "CV",
          signature= "Logit",
          definition=function(object){

              alpha       <- object@slopes$alpha
              meanval     <- object@slopes$meanval


              tempPre  <- sum(exp(meanval + object@pricePre*alpha))
              tempPost <- sum(exp(meanval + object@pricePost*alpha))

              CV <- log(tempPre/tempPost)

              return(CV)

 })



setMethod(
 f= "calcPriceDeltaHypoMon",
 signature= "Logit",
 definition=function(object,prodIndex){

     nprods <- length(prodIndex)
     mc       <- object@mc[prodIndex]
     alpha    <- object@slopes$alpha
     meanval  <- object@slopes$meanval
     pricePre <- object@pricePre

     isOutside   <- as.numeric(object@shareInside < 1)

     calcMonopolySurplus <- function(priceCand){

         pricePre[prodIndex] <- priceCand #keep prices of products not included in HM fixed at premerger levels
         sharesCand <- exp(meanval + alpha*pricePre)
         sharesCand <- sharesCand/(isOutside + sum(sharesCand))

         surplus <- (priceCand-mc)*sharesCand[prodIndex]

         return(sum(surplus))
     }


     maxResult <- optim(object@prices[prodIndex],calcMonopolySurplus,
                              method = "L-BFGS-B",lower = 0,
                              control=list(fnscale=-1))

     pricesHM <- maxResult$par
     priceDelta <- pricesHM/object@pricePre[prodIndex] - 1

     names(priceDelta) <- object@labels[prodIndex]

     return(priceDelta)

 })




logit <- function(prices,shares,margins,
                ownerPre,ownerPost,
                normIndex=ifelse(sum(shares)<1,NA,1),
                mcDelta=rep(0,length(prices)),
                priceStart = prices,
                labels=paste("Prod",1:length(prices),sep=""),
                ...
                ){


    ## Create Logit  container to store relevant data
    result <- new("Logit",prices=prices, shares=shares,
                  margins=margins,
                  normIndex=normIndex,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  priceStart=priceStart,shareInside=sum(shares),
                  labels=labels)

    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate Demand Slope Coefficients
    result <- calcSlopes(result)

    ##Calculate constant marginal costs
    result@mc <- calcMC(result)

    ## Solve Non-Linear System for Price Changes
    result@pricePre  <- calcPrices(result,TRUE,...)
    result@pricePost <- calcPrices(result,FALSE,...)

    return(result)

}

