
setClass(
         Class   = "LogitDistance",
         contains="Logit",

         representation=representation(
         distances="numeric"
         ),

         validity=function(object){

             ## Sanity Checks


             nprods    <- length(object@prices)

             if(nprods != length(object@distances)){
                 stop("'distances' length must equal the number of products")}

             if(any(is.na(object@distances))){
                 stop("'distances' elements cannot be missing")}


         }

         )


setMethod(
          f= "calcSlopes",
          signature= "LogitDistance",
          definition=function(object){
                ## Uncover Demand Coefficents


              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              margins      <-  object@margins
              prices       <-  object@prices
              idx          <-  object@normIndex
              distances     <-  object@distances

              if(is.na(idx)){
                  idxShare    <- 1 - object@shareInside
                  idxPrice    <- 0
                  idxDistance <- 0
              }
              else{
                  idxShare <- shares[idx]
                  idxPrice <- prices[idx]
                  idxDistance <- distances[idx]
               }

              ## Uncover price coefficient and mean valuation from margins and revenue shares


              nprods <- length(shares)
              revenues <- shares * prices


              ## Minimize the distance between observed and predicted margins
              minD <- function(theta){

                  alpha   <- theta[1]
                  tau    <- theta[2]

                  if(!is.na(idx)){
                      meanval <- rep(0,nprods)
                      meanval[-idx] <- theta[-(1:2)]
                  }
                  else{
                      meanval <- theta[-(1:2)]
                  }


                  sharesCand <- exp(meanval + tau*distances + alpha*prices)
                  sharesCand <- sharesCand/(as.numeric(object@shareInside<1) + sum(sharesCand))

                  ## the following returns the elasticity TRANSPOSED
                  elast <- -alpha * matrix(prices * sharesCand,ncol=nprods,nrow=nprods)
                  diag(elast) <- alpha*prices + diag(elast)


                  marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% revenues) / revenues

                  measure <- c(margins - marginsCand,
                               log(share) - log(idxShare) - meanval - tau*(distances - idxDistance) - alpha*(prices - idxPrice) )
                  measure <- sum(measure^2,na.rm=TRUE)

                  return(measure)
              }

              lowerB      <- upperB <- rep(Inf,2 + nprods - object@shareInside==1) # alpha,tau, n meanval, less 1 if no outside good
              lowerB      <- -lowerB
              upperB[1:2] <- 0

              minTheta <- optim(parmsStart,minD,method="L-BFGS-B",lower= lowerB,upper=upperB)

               if(minTheta$convergence != 0){
                   warning("'calcSlopes' nonlinear solver did not successfully converge. Reason: '",minTheta$message,"'")
               }



              meanval <- rep(0,nprods)
              if(!is.na(idx)){
                      meanval[-idx] <- minTheta$par[-(1:2)]
                  }
              else{
                  meanval <- minTheta$par[-(1:2)]
              }

              names(meanval)   <- object@labels

              object@slopes    <- list(alpha=minTheta$par[1],tau=minTheta$par[2], meanval=meanval)


              return(object)

              }

          )

setMethod(
 f= "calcShares",
 signature= "LogitDistance",
 definition=function(object,preMerger=TRUE,revenue=FALSE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

     alpha    <- object@slopes$alpha
     tau      <- object@slopes$tau
     meanval  <- object@slopes$meanval

     shares <- exp(meanval + tau*object@distances + alpha*prices)
     shares <- shares/(as.numeric(object@shareInside<1) + sum(shares))

     if(revenue){shares <- prices*shares/sum(prices*shares)}

     names(shares) <- object@labels

     return(shares)
 }
          )




setMethod(
          f= "CV",
          signature= "LogitDistance",
          definition=function(object){

              alpha       <- object@slopes$alpha
              meanval     <- object@slopes$meanval
              tau         <- object@slopes$tau

              VPre  <- sum(exp(meanval +  tau*object@distances + object@pricePre*alpha))
              VPost <- sum(exp(meanval +  tau*object@distances + object@pricePost*alpha))

              result <- log(VPost/VPre)/alpha

              return(result)

 })


logit.nests <- function(prices,shares,margins,
                        ownerPre,ownerPost,
                        nests=rep(1,length(shares)),
                        normIndex=ifelse(sum(shares) < 1,NA,1),
                        mcDelta=rep(0,length(prices)),
                        priceStart = prices,
                        constraint = TRUE,
                        parmsStart,
                        labels=paste("Prod",1:length(prices),sep=""),
                        ...
                        ){


    if(is.factor(nests)){nests <- nests[,drop=TRUE] }
    else{nests <- factor(nests)}


    if(missing(parmsStart)){

        nNests <- nlevels(nests)
        parmsStart <- runif(nNests+1) # nesting parameter values are assumed to be between 0 and 1
        parmsStart[1] <- -1* parmsStart[1] # price coefficient is assumed to be negative

        if(constraint){parmsStart <- parmsStart[1:2]}
                   }

    ## Create LogitDistance  container to store relevant data
    result <- new("LogitDistance",prices=prices, margins=margins,
                  shares=shares,mcDelta=mcDelta,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  nests=nests,
                  normIndex=normIndex,
                  parmsStart=parmsStart,
                  constraint=constraint,
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

