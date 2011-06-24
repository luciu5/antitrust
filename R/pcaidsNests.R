#source("pcaids.R")


setClass(
         Class = "PCAIDSNests",
         contains="PCAIDS",
         representation=
         representation(
                        margins="vector",
                        nests="factor",
                        nestsParms="vector"),

         validity=function(object){

             nprods <- length(object@shares)

             if(!is.vector(object@margins) || any(object@margins<=0 | object@margins>=1,na.rm=TRUE)){
                 stop("'margins' must be a numeric vector whose values are between 0 and 1")}

             if( nprods != length(object@margins)){
                 stop("'margins' length must equal the number of products")}


             if(nprods != length(object@nests)){
                 stop("'nests' length must equal the number of products")}


             ## Test to see if enough margin info has been supplied to identify all nesting parameters
             nNestParm <- nlevels(object@nests)
             nNestParm <- nNestParm*(nNestParm -1)/2 #calculate the number of nesting parameters
             nMargins  <- length(object@margins[!is.na(object@margins)])

             maxNests <- floor((sqrt(8 * nMargins + 1) + 1)/2) # compute the maximum number of nests that may be used
                                                               # given available margins

             if(!is.vector(object@nestsParms) || nNestParm != length(object@nestsParms)){
                 stop(paste("'nestsParmStart' must be a vector of length",nNestParm))}

             if(nNestParm > nMargins){
                 stop(paste(
                            "Impossible to calibrate nest parameters with the number of margins supplied.\n",
                            "The maximum number of nests supported by the supplied margin information is"
                            , maxNests,"."))
             }
         }

         )





setMethod(
          f= "getNestsParms",
          signature= "PCAIDSNests",
          definition=function(object){

              nests <- object@nests

              nNests <- nlevels(nests)

              labels <- levels(nests)

              nestWeights <- diag(nNests)
              nestWeights[upper.tri(nestWeights)] <- nestWeights[lower.tri(nestWeights)] <- object@nestsParms

              dimnames(nestWeights) <- list(labels,labels)

              return(nestWeights)
          }
          )

setMethod(
          f= "calcSlopes",
          signature= "PCAIDSNests",
          definition=function(object){


              ## Uncover linear demand slopes from shares, knownElast, mktElast, margins, and nesting structure
              ## Since demand is linear IN LOG PRICE, model assumes that slopes remain unchanged following merger


              knownIndx  <- object@knownElastIndex
              knownElast <- object@knownElast
              mktElast   <- object@mktElast
              ownerPre   <- object@ownerPre
              shares <-  object@shares
              margins <- object@margins
              nests <- object@nests


              shareKnown <- shares[knownIndx]
              nprods <- length(shares)


              nMargins <-  length(margins[!is.na(margins)])


              nNests <- nlevels(nests)
              nests <- as.numeric(nests)

              bKnown <- shareKnown * (knownElast + 1 - shareKnown * (mktElast + 1))


              calcB <- function(n){
                  nestWeights <- diag(nNests)
                  nestWeights[upper.tri(nestWeights)] <- nestWeights[lower.tri(nestWeights)] <- n

                  nestWeights <- nestWeights[nests,nests]

                  sumWeights <- sum(shares * nestWeights[,knownIndx], na.rm=TRUE) - shareKnown

                  beta <- shares/shareKnown
                  beta <- beta * (rowSums(t(shares * t(nestWeights)),na.rm=TRUE) - shares)
                  beta <- beta / sumWeights

                  b    <- beta * bKnown


                  B <- -bKnown * (tcrossprod(shares) * nestWeights)/(shareKnown*sumWeights)


                  diag(B) <- b

                  return(B)

              }

              minD <- function(theseNests){

                  Bcand <- calcB(theseNests)

                  elast <- t(Bcand/shares) + shares * (mktElast + 1)
                  diag(elast) <- diag(elast) - 1

                  marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% shares) / shares

                  measure <- sqrt(sum((margins - marginsCand)^2,na.rm=TRUE))/sqrt(nMargins)

                  return(measure)
              }

              minNests <- optim(object@nestsParms,minD,method ="L-BFGS-B",lower=0,upper=1)$par

              B <- calcB(minNests)

              dimnames(B) <- list(object@labels,object@labels)



              object@nestsParms <- minNests
              object@slopes    <- B


              return(object)
          }
          )



setMethod(
 f= "summary",
 signature= "PCAIDSNests",
 definition=function(object){

     callNextMethod(object)

     cat("\nNesting Parameter Estimates:\n\n")
     print(round(getNestsParms(object),2))

     cat("\n\n")


 }


)


pcaids.nests <- function(shares,margins,knownElast,mktElast=-1,ownerPre,ownerPost,
                         nests=rep(1,length(shares)),
                         knownElastIndex=1,
                         mcDelta=rep(0, length(shares)),
                         priceDeltaStart=runif(length(shares)),
                         nestsParmStart=NULL,
                         labels=paste("Prod",1:length(shares),sep=""),
                         ...){

    if(is.factor(nests)){nests <- nests[,drop=TRUE] }
    else{nests <- factor(nests)}

    if(is.null(nestsParmStart)){
        nNests <- nlevels(nests)
        nestsParmStart <- runif(nNests*(nNests -1)/2)
                            }

    ## Create PCAIDS Nests  container to store relevant data
    result <- new("PCAIDSNests",shares=shares,margins=margins,mcDelta=mcDelta
                  ,knownElast=knownElast,mktElast=mktElast,nests=nests,
                  nestsParms=nestsParmStart,
                  ownerPre=ownerPre,ownerPost=ownerPost,knownElastIndex=knownElastIndex,
                  priceDeltaStart=priceDeltaStart,labels=labels)

    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate Demand Slope Coefficients
    result <- calcSlopes(result)

    ## Solve Non-Linear System for Price Changes
    result@priceDelta <- calcPriceDelta(result,...)

    return(result)

}

