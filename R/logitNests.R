##source("logit.R")


setClass(
         Class   = "LogitNests",
         contains="Logit",

         representation=representation(
         nests="factor",
         parmsStart="vector"
         ),

         prototype=prototype(
         parmsStart      =  numeric()
         ),

         validity=function(object){

             ## Sanity Checks


             nprods    <- length(object@prices)
             nNestParm <- nlevels(object@nests) #calculate the number of nesting parameters
             nMargins  <- length(object@margins[!is.na(object@margins)])
             maxNests  <- nMargins - 1

             if(nNestParm==1) stop("'logit.nests' cannot be used for non-nested problems. Use 'logit' instead")

             if(nprods != length(object@nests)){
                 stop("'nests' length must equal the number of products")}


             if(!is.vector(object@parmsStart) || nNestParm + 1 != length(object@parmsStart)){
                 stop(paste("'parmsStart' must be a vector of length",nNestParm + 1))}

             if(nNestParm > nMargins){
                 stop(paste(
                            "Impossible to calibrate nest parameters with the number of margins supplied.\n",
                            "The maximum number of nests supported by the supplied margin information is"
                            ,maxNests,"."))
             }
         }

         )


setMethod(
          f= "calcSlopes",
          signature= "LogitNests",
          definition=function(object){

              ## Uncover Demand Coefficents


              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              margins      <-  object@margins
              prices       <-  object@prices
              idx          <-  object@normIndex

              nests        <- object@nests
              parmsStart   <- object@parmsStart

              nestIdx      <- which(levels(nests)==nests[idx]) # find index of nest whose parameter will be normalized to 1

              if(idx > 0){parmsStart <- parmsStart[-(nestIdx+1)]} # if inside shares all sum to 1, drop nest belonging to index

              ## Uncover price coefficient and mean valuation from margins and revenue shares


              nprods <- length(shares)

              sharesIn <- tapply(shares,nests,sum)[nests]

              sharesIn <- shares / sharesIn

              revenues <- prices * shares

              nMargins <-  length(margins[!is.na(margins)])




              ## Minimize the distance between observed and predicted margins
              minD <- function(theta){

                  alpha <- theta[1]

                  if(length(theta[-1]) < nlevels(nests)){
                      sigma           <- rep(1,length=nlevels(nests))
                      sigma[-nestIdx] <- theta[-1]
                  }

                  else{sigma <- theta[-1]}

                  elast <- diag((1/sigma-1)*alpha)
                  elast <- elast[nests,nests]
                  elast <- elast * matrix(sharesIn*prices,ncol=nprods,nrow=nprods)
                  elast <- -1*(elast + alpha * matrix(revenues,ncol=nprods,nrow=nprods))
                  diag(elast) <- diag(elast) + (1/sigma[nests])*alpha*prices
                  try(marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% revenues) / revenues,TRUE)
                  if(!exists("marginsCand")){marginsCand <- 1e3}

                  measure <- sqrt(sum((margins - marginsCand)^2,na.rm=TRUE))/sqrt(nMargins)

                  return(measure)
              }

              ## Constrained optimizer to look for solutions where alpha<0,  sigma > 0
              lowerB <- upperB <- rep(0,length(parmsStart))
              lowerB[1] <- -Inf

              upperB[-1] <- Inf

              minTheta <- optim(parmsStart,minD,method="L-BFGS-B",lower= lowerB,upper=upperB)

              if(minTheta$convergence != 0){
                  warning("'calcSlopes' nonlinear solver did not successfully converge. Reason: '",minTheta$message,"'")
              }



              minAlpha           <- minTheta$par[1]
              names(minAlpha)    <- "Alpha"

              if(length( minTheta$par[-1]) < nlevels(nests)){

                  minSigma           <- rep(1,length=nlevels(nests))
                  minSigma[-nestIdx] <- minTheta$par[-1]

              }

              else{
                  minSigma <- minTheta$par[-1]
                  }


              minSigmaOut        <- minSigma
              minSigma           <- minSigma[nests]
              names(minSigma)    <- as.character(nests)
              names(minSigmaOut) <- levels(nests)

              if(any(minSigma>1)){
                  warning("Some nesting parameters are greater than 1. These parameter values may not be consistent with economic theory")}


              if(idx==0){
                  idxShare <- 1 - object@shareInside
                  idxShareIn <- 1
                  idxPrice   <- 0
                  idxSigma   <- 1

              }
              else{
                  idxShare   <- shares[idx]
                  idxShareIn <- sharesIn[idx]
                  idxPrice   <- prices[idx]
                  idxSigma   <- minSigma[idx]
               }

              meanval <- log(shares) - log(idxShare) - minAlpha*(prices - idxPrice)

              meanval <- meanval + (minSigma-1)*log(sharesIn) -
                         (idxSigma-1)*log(idxShareIn)

              names(meanval)   <- object@labels

              object@slopes    <- list(alpha=minAlpha,sigma=minSigmaOut,meanval=meanval)

              return(object)
          }
          )




setMethod(
 f= "calcPrices",
 signature= "LogitNests",
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
     nests  <- object@nests

     alpha    <- object@slopes$alpha
     sigma    <- object@slopes$sigma
     meanval  <- object@slopes$meanval
     isOutside <- as.numeric(object@shareInside < 1)

     ##Define system of FOC as a function of prices
     FOC <- function(priceCand){

         margins      <- 1 - (object@mc * (1 + mcDelta)) / priceCand

         sharesIn     <- exp((meanval+alpha*priceCand)/sigma[nests])
         inclusiveValue <- log(tapply(sharesIn,nests,sum))
         sharesIn     <- sharesIn/sum(sharesIn)
         sharesAcross <-   exp(sigma*inclusiveValue)
         sharesAcross <- sharesAcross/(isOutside + sum(sharesAcross))

         shares       <- sharesIn * sharesAcross[nests]
         revenues <-  shares * priceCand

         elast <- diag((1/sigma-1)*alpha)
         elast <- elast[nests,nests]
         elast <- elast * matrix(sharesIn*priceCand,ncol=nprods,nrow=nprods)
         elast <- -1*(elast + alpha * matrix(revenues,ncol=nprods,nrow=nprods))
         diag(elast) <- diag(elast) + (1/sigma[nests])*alpha*priceCand



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
 signature= "LogitNests",
 definition=function(object,preMerger=TRUE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

     nests     <- object@nests
     alpha     <- object@slopes$alpha
     sigma     <- object@slopes$sigma
     meanval   <- object@slopes$meanval
     isOutside <- as.numeric(object@shareInside < 1)

     sharesIn     <- exp((meanval+alpha*prices)/sigma[nests])
     inclusiveValue <- log(tapply(sharesIn,nests,sum))
     sharesIn     <- sharesIn/(isOutside + sum(sharesIn))
     sharesAcross <-   exp(sigma*inclusiveValue)
     sharesAcross <- sharesAcross/(isOutside + sum(sharesAcross))

     shares       <- sharesIn * sharesAcross[nests]

     names(shares) <- object@labels

     return(shares)

}
 )



setMethod(
 f= "elast",
 signature= "LogitNests",
 definition=function(object,preMerger=TRUE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

     nests    <- object@nests
     alpha    <- object@slopes$alpha
     sigma    <- object@slopes$sigma
     meanval  <- object@slopes$meanval
     isOutside <- as.numeric(object@shareInside < 1)

     sharesIn     <- exp((meanval+alpha*prices)/sigma[nests])
     inclusiveValue <- log(tapply(sharesIn,nests,sum))
     sharesIn     <- sharesIn/(isOutside + sum(sharesIn))
     sharesAcross <-   exp(sigma*inclusiveValue)
     sharesAcross <- sharesAcross/(isOutside + sum(sharesAcross))

     shares       <- sharesIn * sharesAcross[nests]

     nprods <-  length(shares)

     elast <- diag((1/sigma-1)*alpha)
     elast <- elast[nests,nests]
     elast <- elast * matrix(sharesIn*prices,ncol=nprods,nrow=nprods,byrow=TRUE)
     elast <- -1*(elast + alpha * matrix(shares*prices,ncol=nprods,nrow=nprods,byrow=TRUE))
     diag(elast) <- diag(elast) + (1/sigma[nests])*alpha*prices

     dimnames(elast) <- list(object@labels,object@labels)

      return(elast)

}
 )





setMethod(
          f= "CV",
          signature= "LogitNests",
          definition=function(object){

              nests       <- object@nests
              alpha       <- object@slopes$alpha
              sigma       <- object@slopes$sigma
              meanval     <- object@slopes$meanval
              isOutside   <- as.numeric(object@shareInside < 1)


              tempPre  <- sum( tapply(exp((meanval + object@pricePre*alpha)  / sigma[nests]),nests,sum) ^ sigma ) + isOutside
              tempPost <- sum( tapply(exp((meanval + object@pricePost*alpha) / sigma[nests]),nests,sum) ^ sigma ) + isOutside



              CV <- log(tempPre/tempPost)

              names(CV) <- NULL

              return(CV)

 })

setMethod(
          f= "getNestsParms",
          signature= "LogitNests",
          definition=function(object){

              nestParms <- object@slopes$sigma


              return(nestParms)

          }
          )

setMethod(
 f= "summary",
 signature= "LogitNests",
 definition=function(object){

     callNextMethod(object)

     nestParms <- getNestsParms(object)

     cat("\nNesting Parameter Estimates:\n\n")
     print(round(nestParms,2))

     cat("\n\n")


 }


)


setMethod(
 f= "calcPriceDeltaHypoMon",
 signature= "LogitNests",
 definition=function(object,prodIndex){

     nprods <- length(prodIndex)
     mc       <- object@mc[prodIndex]
     nests     <- object@nests
     alpha     <- object@slopes$alpha
     sigma     <- object@slopes$sigma
     meanval  <- object@slopes$meanval
     pricePre <- object@pricePre

     isOutside   <- as.numeric(shareInside < 1)

     calcMonopolySurplus <- function(priceCand){

         pricePre[prodIndex] <- priceCand #keep prices of products not included in HM fixed at premerger levels

         sharesIn     <- exp((meanval+alpha*pricePre)/sigma[nests])
         inclusiveValue <- log(tapply(sharesIn,nests,sum))
         sharesIn     <- sharesIn/(isOutside + sum(sharesIn))
         sharesAcross <-   exp(sigma*inclusiveValue)
         sharesAcross <- sharesAcross/(isOutside + sum(sharesAcross))

         sharesCand       <- sharesIn * sharesAcross[nests]

         surplus <- (priceCand-mc)*sharesCand[prodIndex]

         return(sum(surplus))
     }


     minResult <- optim(object@prices[prodIndex],calcMonopolySurplus,
                              method = "L-BFGS-B",lower = 0,
                              control=list(fnscale=-1))

     pricesHM <- minResult$par
     priceDelta <- pricesHM/object@pricePre[prodIndex] - 1

     return(priceDelta)

 })



logit.nests <- function(prices,shares,margins,
                        ownerPre,ownerPost,
                        nests=rep(1,length(shares)),
                        shareInside = 1,
                        normIndex=ifelse(shareInside < 1,0,1),
                        mcDelta=rep(0,length(prices)),
                        priceStart = prices,
                        parmsStart=NULL,
                        labels=paste("Prod",1:length(prices),sep=""),
                        ...
                        ){




    shares <- shares * shareInside

    if(is.factor(nests)){nests <- nests[,drop=TRUE] }
    else{nests <- factor(nests)}


    if(is.null(parmsStart)){
        nNests <- nlevels(nests)
        parmsStart <- runif(nNests+1) # nesting parameter values are assumed to be between 0 and 1
        parmsStart[1] <- -1* parmsStart[1] # price coefficient is assumed to be negative
                            }

    ## Create LogitNests  container to store relevant data
    result <- new("LogitNests",prices=prices, margins=margins,
                  shares=shares,mc=prices*(1-margins),mcDelta=mcDelta,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  nests=nests,
                  normIndex=normIndex,
                  parmsStart=parmsStart,
                  priceStart=priceStart,shareInside=shareInside,
                  labels=labels)

    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)
    ##return(result)
    ## Calculate Demand Slope Coefficients
    result <- calcSlopes(result)

    ## Solve Non-Linear System for Price Changes
    result@pricePre  <- calcPrices(result,TRUE,...)
    result@pricePost <- calcPrices(result,FALSE,...)

    return(result)

}

