setClass(
         Class   = "LogitALM",
         contains="Logit",

         prototype=prototype(
         mktElast      =  numeric()
         ),

         validity=function(object){

             ## Sanity Checks

             if(length(object@mktElast)!=1 ||
                object@mktElast>=0){
                 stop("'mktElast' must be negative")
             }
         }
         )

setMethod(
          f= "calcSlopes",
          signature= "LogitALM",
          definition=function(object){

              ## Uncover Demand Coefficents

              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              margins      <-  object@margins
              prices       <-  object@prices
              idx          <-  object@normIndex
              mktElast     <-  object@mktElast

              avgPrice     <- sum(price*shares)


              nMargins <-  length(margins[!is.na(margins)])

              minD <- function(alpha){


                  probs <- shares*(1 - mktElast/(alpha*avgPrice))

                  elast <- -alpha *  matrix(prices * probs,ncol=nprods,nrow=nprods)
                  diag(elast) <- alpha*prices - diag(elast)

                  revenues <- probs * prices
                  marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% revenues) / revenues

                  measure <- sqrt(sum((margins - marginsCand)^2,na.rm=TRUE))/sqrt(nMargins)

                  return(measure)
              }

              minAlpha <- optimize(minD,c(-1e12,0))$minimum

              probs   <- shares * (1-mktElast/(minalpha*avgPrice))
              idxProb <- probs[idx]
              idxPrice <- prices[idx]

              meanval <- log(probs) - log(idxProb) - minAlpha * (prices - idxPrice)

              names(meanval)   <- object@labels

              object@slopes    <- list(alpha=minAlpha,meanval=meanval)

              return(object)

          }

          )



logit.alm <- function(prices,shares,margins,
                mktElast=1,
                ownerPre,ownerPost,
                normIndex=1,
                mcDelta=rep(0,length(prices)),
                priceStart = prices,
                labels=paste("Prod",1:length(prices),sep=""),
                ...
                ){


    ## Create Logit  container to store relevant data
    result <- new("LogitALM",prices=prices, shares=shares,
                  margins=margins,
                  mktElast=mktElast,
                  normIndex=normIndex,
                  mc=prices*(1-margins),mcDelta=mcDelta,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  priceStart=priceStart,shareInside=1,
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
