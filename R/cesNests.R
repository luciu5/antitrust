##source("ces.R")


setClass(
         Class   = "CESNests",
         contains="CES",

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

             if(nNestParm==1) stop("'ces.nests' cannot be used for non-nested problems. Use 'ces' instead")

             if(nprods != length(object@nests)){
                 stop("'nests' length must equal the number of products")}

             if(!(object@normIndex %in% seq(1,nprods)) ){
                 stop("'normIndex' value must be between 1 and the length of 'prices'")}


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
          signature= "CESNests",
          definition=function(object){

              ## Uncover Demand Coefficents


              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              margins      <-  object@margins
              prices       <-  object@prices
              quantities   <-  object@quantities
              idx          <-  object@normIndex
              shareInside  <-  object@shareInside
              nests        <- object@nests
              parmsStart   <- object@parmsStart

              ## uncover Numeraire Coefficients
              if(length(shareInside)>0) {alpha <- 1/shareInside -1}
              else{alpha <- NULL}



              ## Uncover price coefficient and mean valuation from margins and revenue shares


              nprods <- length(shares)

              sharesNests <- tapply(price * quantity,nests,sum)[nests]

              sharesNests <- (price * quantity) / sharesNests



              nMargins <-  length(margins[!is.na(margins)])

              ## Minimize the distance between observed and predicted margins
              minD <- function(theta){

                  gamma <- theta[1]
                  sigma <- theta[-1]



                  elast <- diag(sigma - gamma)

                  elast <- elast[nests,nests]
                  elast <- elast * matrix(sharesNests,ncol=nprods,nrow=nprods)
                  elast <- elast + (gamma-1) * matrix(shares,ncol=nprods,nrow=nprods)

                  diag(elast) <- diag(elast) - sigma[nests]

                  marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% shares) / shares

                  measure <- sqrt(sum((margins - marginsCand)^2,na.rm=TRUE))/sqrt(nMargins)

                  return(measure)
              }

              ## Constrain optimizer to look for solutions where sigma_i > gamma > 1 for all i
              constrA <- diag(nlevels(nests) + 1)
              constrA[-1,1] <- -1

              constrB <- rep(0,nlevels(nests) + 1)
              constrB[1] <- 1

              minTheta <- constrOptim(parmsStart,minD,grad=NULL,ui=constrA,ci=constrB)$par
              names(minTheta) <- c("Gamma",levels(nests))

              minGamma <- minTheta[1]
              minSigma <- minTheta[-1]
              minSigma <- minSigma[nests]

              meanval <- log(shares) - log(shares[idx]) + (minGamma - 1) * (log(prices) - log(prices[idx]))

              meanval <- meanval - (minSigma-minGamma)/(minSigma-1)*log(sharesNests) +
                         (minSigma[idx]-minGamma)/(minSigma[idx]-1)*log(sharesNests[idx])

              meanval <- exp( (minSigma-1)/(minGamma-1) * meanval )


              names(meanval)   <- object@labels

              object@slopes    <- list(alpha=alpha,gamma=minGamma,sigma=minTheta[-1],meanval=meanval)
              object@totExp    <- (1 + alpha) * sum(prices*quantities)

              return(object)
          }
          )




setMethod(
 f= "calcPrices",
 signature= "CESNests",
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

     gamma    <- object@slopes$gamma
     sigma    <- object@slopes$sigma
     meanval  <- object@slopes$meanval


     ##Define system of FOC as a function of prices
     FOC <- function(priceCand){

         margins      <- 1 - (object@mc * (1 + mcDelta)) / priceCand


         sharesIn     <- meanval*priceCand^(1-sigma[nests])
         sharesAcross <- tapply(sharesIn,nests,sum)
         sharesIn     <- sharesIn / sharesAcross[nests]
         sharesAcross <- sharesAcross^((1-gamma)/(1-sigma))
         sharesAcross <- sharesAcross / sum(sharesAcross)

         shares       <- sharesIn * sharesAcross[nests]



         elast <- diag(sigma - gamma)

         elast <- elast[nests,nests]
         elast <- elast * matrix(sharesIn,ncol=nprods,nrow=nprods)
         elast <- elast + (gamma-1) * matrix(shares,ncol=nprods,nrow=nprods)
         diag(elast) <- diag(elast) - sigma[nests]

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
 signature= "CESNests",
 definition=function(object,preMerger=TRUE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

     nests     <- object@nests
     gamma    <- object@slopes$gamma
     sigma    <- object@slopes$sigma
     meanval  <- object@slopes$meanval


     sharesIn     <- meanval*prices^(1-sigma[nests])
     sharesAcross <- tapply(sharesIn,nests,sum)
     sharesIn     <- sharesIn / sharesAcross[nests]
     sharesAcross <- sharesAcross^((1-gamma)/(1-sigma))
     sharesAcross <- sharesAcross / sum(sharesAcross)

     shares       <- sharesIn * sharesAcross[nests]

     names(shares) <- object@labels

     return(shares)

}
 )



setMethod(
 f= "elast",
 signature= "CESNests",
 definition=function(object,preMerger=TRUE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

     nests    <- object@nests
     gamma    <- object@slopes$gamma
     sigma    <- object@slopes$sigma
     meanval  <- object@slopes$meanval

     sharesIn     <- meanval*prices^(1-sigma[nests])
     sharesAcross <- tapply(sharesIn,nests,sum)
     sharesIn     <- sharesIn / sharesAcross[nests]
     sharesAcross <- sharesAcross^((1-gamma)/(1-sigma))
     sharesAcross <- sharesAcross / sum(sharesAcross)

     shares       <- sharesIn * sharesAcross[nests]

     nprods <-  length(shares)

     elast <- diag(sigma - gamma)
     elast <- elast[nests,nests]
     elast <- elast * matrix(sharesIn,ncol=nprods,nrow=nprods,byrow=TRUE)
     elast <- elast + (gamma-1) * matrix(shares,ncol=nprods,nrow=nprods,byrow=TRUE)
     diag(elast) <- diag(elast) - sigma[nests]

     dimnames(elast) <- list(object@labels,object@labels)

      return(elast)

}
 )





setMethod(
          f= "CV",
          signature= "CESNests",
          definition=function(object){

              alpha       <- object@slopes$alpha

              if(is.null(alpha)) stop("'shareInside' must be provided in order to calculate Compensating Variation")


              nests       <- object@nests
              gamma       <- object@slopes$gamma
              sigma       <- object@slopes$sigma
              meanval     <- object@slopes$meanval
              shareInside <- object@shareInside
              totExp      <- object@totExp
              price       <- object@prices
              quantities  <- object@quantities



              tempPre  <- log( sum( tapply(meanval * object@pricePre^(1-sigma[nests]),nests,sum) ^((1-gamma)/(1-sigma)) ) )
              tempPre  <- (gamma/(gamma-1)) * log( sum( tapply( meanval^((1-gamma)/(1-sigma[nests])) * object@pricePre^(-gamma),nests,sum)^((gamma-1)/gamma)) ) - tempPre

              tempPost  <- log( sum( tapply(meanval * object@pricePost^(1-sigma[nests]),nests,sum) ^((1-gamma)/(1-sigma)) ) )
              tempPost  <- (gamma/(gamma-1)) * log( sum( tapply( meanval^((1-gamma)/(1-sigma[nests])) * object@pricePost^(-gamma),nests,sum)^((gamma-1)/gamma)) ) - tempPost



              CV <- exp(log(totExp) + (tempPre - tempPost)/(1+alpha)) - totExp

              names(CV) <- NULL

              return(CV)

 })





setMethod(
 f= "summary",
 signature= "CESNests",
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

     cat("\nNesting Parameter Estimates:\n\n")
     print(round(object@slopes$sigma),2)

     cat("\n\n")
}
 )



ces.nests <- function(prices,quantities,margins,
                      ownerPre,ownerPost,
                      nests=rep(1,length(shares)),
                      normIndex=1,
                      mcDelta=rep(0,length(prices)),
                      shareInside = NULL,
                      priceStart = prices,
                      parmsStart=NULL,
                      labels=paste("Prod",1:length(prices),sep=""),
                      ...
                      ){



    revenues <- prices*quantities
    shares=(revenues)/sum(revenues) #Revenue based shares
    if(is.null(shareInside)){shareInside <- numeric()}

    if(is.factor(nests)){nests <- nests[,drop=TRUE] }
    else{nests <- factor(nests)}


    if(is.null(parmsStart)){
        nNests <- nlevels(nests)
        parmsStart <- cumsum(runif(nNests+1,1,1.5)) # parameter values are assumed to be greater than 1
                            }

    ## Create CESNests  container to store relevant data
    result <- new("CESNests",prices=prices, quantities=quantities,margins=margins,
                  shares=shares,mc=prices*(1-margins),mcDelta=mcDelta,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  nests=nests,
                  normIndex=normIndex,
                  parmsStart=parmsStart,
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

