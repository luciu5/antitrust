
setClass(
         Class   = "LogitNests",
         contains="Logit",

         representation=representation(
         nests="factor",
         parmsStart="numeric",
         constraint="logical"
         ),

         prototype=prototype(
         parmsStart      =  numeric()
         ),

         validity=function(object){




             nprods    <- length(object@prices)
             nNestParm <- nlevels(object@nests) #calculate the number of nesting parameters
             nMargins  <- length(object@margins[!is.na(object@margins)])
             maxNests  <- nMargins - 1

             ## Identify Singleton Nests
             nestCnt   <- tapply(object@prices,object@nests,length)
             nestCnt   <- nestCnt[object@nests]
             isSingleton <- nestCnt==1

             if(nNestParm==1) stop("'logit.nests' cannot be used for non-nested problems. Use 'logit' instead")

             if(nprods != length(object@nests)){
                 stop("'nests' length must equal the number of products")}


             if(object@constraint && length(object@parmsStart)!=2){
                 stop("when 'constraint' is TRUE, 'parmsStart' must be a vector of length 2")
                 }
             else if(!object@constraint && nNestParm + 1 != length(object@parmsStart)){
                 stop("when 'constraint' is FALSE, 'parmsStart' must be a vector of length ",nNestParm + 1)
             }


             if(!object@constraint &&
                any(tapply(object@margins[!isSingleton],object@nests[!isSingleton],
                           function(x){if(all(is.na(x))){return(TRUE)} else{return(FALSE)}}
                           )
                    )
                ){
                 stop("when 'constraint' is FALSE, at least one product margin must be supplied for each non-singleton nest")
             }

             if(nNestParm > nMargins){
                 stop(paste(
                            "Impossible to calibrate nest parameters with the number of margins supplied.\n",
                            "The maximum number of nests supported by the supplied margin information is"
                            ,maxNests))
             }


             return(TRUE)
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
              naMargins    <- !is.na(object@margins)

              prices       <-  object@prices
              idx          <-  object@normIndex

              parmsStart   <- object@parmsStart
              nests        <- object@nests
              nestCnt      <- tapply(prices,nests,length)
              constraint   <- object@constraint

              isSingletonNest <- nestCnt==1

              if(any(isSingletonNest)){
                  warning("Some nests contain only one product; their nesting parameters are not identified.
Normalizing these parameters to 1.")

              }




              if(!constraint){
                  parmsStart   <- parmsStart[c(TRUE,!isSingletonNest)] #always retain first element; this is
                                                                          # the initial value for price coefficient
              }
               ## Uncover price coefficient and mean valuation from margins and revenue shares


              nprods <- length(shares)

              sharesNests <- tapply(shares,nests,sum)[nests]

              sharesNests <- shares / sharesNests

              revenues <- prices * shares


              ## create index variables, contingent on whether an outside good is defined
              if(is.na(idx)){
                  idxShare <- 1 - object@shareInside
                  idxShareIn <- 1
                  idxPrice   <- 0


              }

              else{


                  idxShare   <- shares[idx]
                  idxShareIn <- sharesNests[idx]
                  idxPrice   <- prices[idx]

               }




              ## Minimize the distance between observed and predicted margins
              minD <- function(theta){

                  alpha <- theta[1]
                  sigma <- as.numeric(isSingletonNest)
                  sigma[!isSingletonNest] <- theta[-1]

                  ## The following returns the elasticity TRANSPOSED
                  elasticity <- diag((1/sigma-1)*alpha)
                  elasticity <- elasticity[nests,nests]
                  elasticity <- elasticity * matrix(sharesNests*prices,ncol=nprods,nrow=nprods)
                  elasticity <- -1*(elasticity + alpha * matrix(shares*prices,ncol=nprods,nrow=nprods))
                  diag(elasticity) <- diag(elasticity) + (1/sigma[nests])*alpha*prices



                  #marginsCand <- -1 * as.vector(solve(elasticity * ownerPre) %*% revenues) / revenues
                  #marginsCand[is.nan(marginsCand)] <- 1e30

                  #measure <- sum((margins[naMargins] - marginsCand[naMargins])^2)

                  measure <- revenues * diag(ownerPre) + as.vector((elasticity * ownerPre) %*% (margins * revenues))
                  #measure[is.nan(measure)] <- 1e30

                  measure <- sum(measure[naMargins]^2,na.rm=TRUE)

                  return(measure)
              }

              ##  Constrained optimizer to look for solutions where alpha<0,  1 > sigma > 0.
              ##  sigma > 1 or sigma < 0 imply complements
              lowerB <- upperB <- rep(0,length(parmsStart))
              lowerB[1] <- -Inf

              upperB[-1] <- 1

              minTheta <- optim(parmsStart,minD,method="L-BFGS-B",lower= lowerB,upper=upperB)


              if(minTheta$convergence != 0){
                  warning("'calcSlopes' nonlinear solver did not successfully converge. Reason: '",minTheta$message,"'")
              }



              minAlpha           <- minTheta$par[1]
              names(minAlpha)    <- "Alpha"


              minSigma <-  as.numeric(isSingletonNest)
              minSigma[!isSingletonNest] <- minTheta$par[-1]


              minSigmaOut        <- minSigma
              minSigma           <- minSigma[nests]
              names(minSigmaOut) <- levels(nests)


              ##normalize outside good nesting parameter to 1
              if(is.na(idx)){
                  idxSigma   <- minSigma[idx]
              }
              else{idxSigma <- 1}


              meanval <-
                  log(shares) - log(idxShare) -
                      minAlpha*(prices - idxPrice) +
                          (minSigma-1)*log(sharesNests) -
                              (idxSigma-1)*log(idxShareIn)

              names(meanval)   <- object@labels

              object@slopes    <- list(alpha=minAlpha,sigma=minSigmaOut,meanval=meanval)

              return(object)
          }
          )

setMethod(
 f= "calcShares",
 signature= "LogitNests",
 definition=function(object,preMerger=TRUE,revenue=FALSE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

     nests     <- object@nests
     alpha     <- object@slopes$alpha
     sigma     <- object@slopes$sigma
     meanval   <- object@slopes$meanval
     isOutside <- as.numeric(object@shareInside < 1)

     sharesIn     <- exp((meanval+alpha*prices)/sigma[nests])

     inclusiveValue <- log(tapply(sharesIn,nests,sum))
     sharesAcross <-   exp(sigma*inclusiveValue)
     sharesAcross <- sharesAcross/(isOutside + sum(sharesAcross))


     sharesIn     <- sharesIn/(isOutside + sum(sharesIn))


     shares       <- sharesIn * sharesAcross[nests]

     if(revenue){shares <- prices*shares/sum(prices*shares,1-sum(shares))}

     names(shares) <- object@labels

     return(as.vector(shares))

}
 )



setMethod(
 f= "elast",
 signature= "LogitNests",
 definition=function(object,preMerger=TRUE,market=FALSE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}


     nests    <- object@nests
     alpha    <- object@slopes$alpha
     sigma    <- object@slopes$sigma
     meanval  <- object@slopes$meanval

     shares <- calcShares(object,preMerger,revenue=FALSE)

     if(market){

         elast <- alpha * sum(shares*prices) * (1 - sum(shares))
         names(elast) <- NULL
         }

     else{

         sharesNests <- shares/tapply(shares,nests,sum)[nests]


         nprods <-  length(shares)

         elast <- diag((1/sigma-1)*alpha)
         elast <- elast[nests,nests]
         elast <- elast * matrix(sharesNests*prices,ncol=nprods,nrow=nprods,byrow=TRUE)
         elast <- -1*(elast + alpha * matrix(shares*prices,ncol=nprods,nrow=nprods,byrow=TRUE))
         diag(elast) <- diag(elast) + (1/sigma[nests])*alpha*prices

         dimnames(elast) <- list(object@labels,object@labels)

         }
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



              VPre  <- sum( tapply(exp((meanval + object@pricePre*alpha)  / sigma[nests]),nests,sum) ^ sigma )
              VPost <- sum( tapply(exp((meanval + object@pricePost*alpha) / sigma[nests]),nests,sum) ^ sigma )



              result <- log(VPost/VPre)/alpha
              names(result) <- NULL
              return(result)

 })



logit.nests <- function(prices,shares,margins,
                        ownerPre,ownerPost,
                        nests=rep(1,length(shares)),
                        normIndex=ifelse(sum(shares) < 1,NA,1),
                        mcDelta=rep(0,length(prices)),
                        priceStart = prices,
                        isMax=FALSE,
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

    ## Create LogitNests  container to store relevant data
    result <- new("LogitNests",prices=prices, margins=margins,
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
    result@pricePre  <- calcPrices(result,preMerger=TRUE,isMax=isMax,...)
    result@pricePost <- calcPrices(result,preMerger=FALSE,isMax=isMax,...)


    return(result)

}

