
setClass(
         Class   = "Logit",
         contains="Antitrust",
         representation=representation(

         prices           = "numeric",
         margins          = "numeric",
         pricePre         = "numeric",
         pricePost        = "numeric",
         priceStart       = "numeric",
         normIndex        = "vector",
         shareInside     = "numeric"

         ),
          prototype=prototype(

          pricePre      =  numeric(),
          pricePost     =  numeric()

        ),


         validity=function(object){

             ## Sanity Checks

             margins <- object@margins

             nprods <- length(object@shares)
             sumShares <- sum(object@shares)

             if(
                 nprods != length(margins) ||
                 nprods != length(object@prices)){
                 stop("'prices', 'margins' and 'shares' must all be vectors with the same length")}

             if(any(object@prices<0,na.rm=TRUE))             stop("'prices' values must be positive")


             if(any(margins<0 | margins>1,na.rm=TRUE)) stop("'margins' values must be between 0 and 1")

             if(!(is.matrix(object@ownerPre) || is(object@ownerPre,"Matrix"))){
                 ownerPre <- ownerToMatrix(object,TRUE)
             }
             else{ownerPre <- object@ownerPre}


             margins[is.na(margins)]=0
             if(max(as.vector(ownerPre %*% margins))==0) stop("Insufficient margin information to calibrate demand parameters.")

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

             return(TRUE)

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
              revenues <- shares * prices


              ## Minimize the distance between observed and predicted margins
              minD <- function(alpha){

                  ## the following returns the elasticity TRANSPOSED
                  elast <- -alpha *  matrix(prices * shares,ncol=nprods,nrow=nprods)
                  diag(elast) <- alpha*prices + diag(elast)


                  marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% revenues) / revenues

                  measure <- sum((margins - marginsCand)^2,na.rm=TRUE)

                  return(measure)
              }

              minAlpha <- optimize(minD,c(-1e6,0))$minimum


              meanval <- log(shares) - log(idxShare) - minAlpha * (prices - idxPrice)

              names(meanval)   <- object@labels

              object@slopes    <- list(alpha=minAlpha,meanval=meanval)


              return(object)
          }
          )



setMethod(
 f= "calcPrices",
 signature= "Logit",
 definition=function(object,preMerger=TRUE,...){


     if(preMerger){owner <- object@ownerPre}
     else{owner <- object@ownerPost}

     mc <- calcMC(object,preMerger)


     ##Define system of FOC as a function of prices
     FOC <- function(priceCand){

         if(preMerger){ object@pricePre <- priceCand}
         else{          object@pricePost <- priceCand}


         margins   <- 1 - mc/priceCand
         revenues  <- calcShares(object,preMerger,revenue=TRUE)
         elasticities     <- elast(object,preMerger)

         thisFOC <- revenues + as.vector(t(elasticities * owner) %*% (margins * revenues))

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
 definition=function(object,preMerger=TRUE,revenue=FALSE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

     alpha    <- object@slopes$alpha
     meanval  <- object@slopes$meanval

     shares <- exp(meanval + alpha*prices)
     shares <- shares/(as.numeric(object@shareInside<1) + sum(shares))

     if(revenue){shares <- prices*shares/sum(prices*shares)}

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
     diag(elast) <- alpha*prices + diag(elast)

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


              VPre  <- sum(exp(meanval + object@pricePre*alpha))
              VPost <- sum(exp(meanval + object@pricePost*alpha))

              result <- log(VPost/VPre)/alpha

              return(result)

 })



setMethod(
 f= "calcPriceDeltaHypoMon",
 signature= "Logit",
 definition=function(object,prodIndex){


     mc       <- calcMC(object,TRUE)[prodIndex]
     pricePre <- object@pricePre

     calcMonopolySurplus <- function(priceCand){

         pricePre[prodIndex] <- priceCand #keep prices of products not included in HM fixed at premerger levels
         object@pricePre <- pricePre
         sharesCand <- calcShares(object,TRUE,revenue=FALSE)

         surplus <- (priceCand-mc)*sharesCand[prodIndex]

         return(sum(surplus))
     }


     maxResult <- optim(object@prices[prodIndex],calcMonopolySurplus,
                              method = "L-BFGS-B",lower = 0,
                              control=list(fnscale=-1))

     pricesHM <- maxResult$par
     priceDelta <- pricesHM/pricePre[prodIndex] - 1

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
                  mcDelta=mcDelta,
                  priceStart=priceStart,shareInside=sum(shares),
                  labels=labels)

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

