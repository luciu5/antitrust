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
          pricePost     =  vector()

        ),


         validity=function(object){

             ## Sanity Checks

             nprods <- length(object@shares)

             if(
                 nprods != length(object@margins) ||
                 nprods != length(object@prices)){
                 stop("'prices', 'margins' and 'shares' must all be vectors with the same length")}

             if(any(object@prices<0,na.rm=TRUE))             stop("'prices' values must be positive")


             if(any(object@margins<0 | object@margins>1,na.rm=TRUE)) stop("'margins' values must be between 0 and 1")

             if(nprods != length(object@priceStart)){
                 stop("'priceStart' must have the same length as 'shares'")}

              if(
                !(object@shareInside >0 &&
                  object@shareInside <=1)
                ){
                 stop("'shareInside' must be between 0 and 1")
         }

             if(!(object@normIndex %in% 0:nprods) ||
                (object@shareInside==1 && object@normIndex==0)){
                 stop("'normIndex' must take on a value between 1 and ",nprods, " if 'shares' sum to 1 , or '0' if the sum of shares is less than 1")
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

              if(idx==0){
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




setMethod(
 f= "calcPrices",
 signature= "Logit",
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





logit <- function(prices,shares,margins,
                ownerPre,ownerPost,
                normIndex=ifelse(sum(shares)<1,0,1),
                mcDelta=rep(0,length(prices)),
                priceStart = prices,
                labels=paste("Prod",1:length(prices),sep=""),
                ...
                ){


    ## Create Logit  container to store relevant data
    result <- new("Logit",prices=prices, shares=shares,
                  margins=margins,
                  normIndex=normIndex,
                  mc=prices*(1-margins),mcDelta=mcDelta,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
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

