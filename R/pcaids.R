setClass(
         Class = "PCAIDS",
         contains="AIDS",
         representation=representation(

         knownElast      = "numeric",
         knownElastIndex = "numeric"
         ),

         validity=function(object){




             nprods <- length(object@shares)

             if(!(object@knownElastIndex %in% seq(1,nprods)) ){
                 stop("'knownElastIndex' value must be between 1 and the length of 'shares'")}
             if(nprods != length(object@mcDelta)){
                 stop("'mcDelta' must have the same length as 'shares'")}
             if(object@knownElast>0 || object@mktElast > -1 ){
                 stop("All elasticities must be negative. 'mktElast' must be less than -1")}
             if(object@knownElast > object@mktElast){
                 stop("'mktElast' must be greater than 'knownElast'")}

         }

         )








setMethod(
 f= "calcSlopes",
 signature= "PCAIDS",
 definition=function(object){


     ## Uncover linear demand slopes from shares, knownElast and mktElast
     ## Since demand is linear IN LOG PRICE, model assumes that slopes remain unchanged following merger

     shares    <- object@shares
     diversion <- object@diversion
     labels    <- object@labels

     idx        <- object@knownElastIndex
     shareKnown <- shares[idx]

     nprod <- length(shares)

     bKnown <- shareKnown * (object@knownElast + 1 - shareKnown * (object@mktElast + 1))

     b <- bKnown*diversion[idx,]/diversion[,idx]

     B <- -diversion * b

     dimnames(B) <- list(labels,labels)
     object@slopes <- B
     object@intercepts <- as.vector(shares - B%*%object@prices)
     names(object@intercepts) <- object@labels

     return(object)
 }

      )









pcaids <- function(shares,knownElast,mktElast=-1,
                   prices,diversions,
                   ownerPre,ownerPost,
                   knownElastIndex=1,
                   mcDelta=rep(0, length(shares)),
                   priceStart=runif(length(shares)),
                   isMax=FALSE,
                   labels=paste("Prod",1:length(shares),sep=""),
                   ...){



    if(missing(prices)){ prices <- rep(NA,length(shares))}

    if(missing(diversions)){
        diversions <- tcrossprod(1/(1-shares),shares)
        diag(diversions) <- -1
    }

  ## Create PCAIDS container to store relevant data
    result <- new("PCAIDS",shares=shares,prices=prices,
                   quantities=shares, margins=shares,mcDelta=mcDelta,
                  knownElast=knownElast,mktElast=mktElast,
                  ownerPre=ownerPre,ownerPost=ownerPost,
                  knownElastIndex=knownElastIndex,
                  diversion=diversions,
                  priceStart=priceStart,labels=labels)

    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate Demand Slope Coefficients
    result <- calcSlopes(result)

    ## Solve Non-Linear System for Price Changes
    result@priceDelta <- calcPriceDelta(result,isMax=isMax,...)


    ## Calculate Pre and Post merger equilibrium prices
    ## These are equal to NA in pcaids
    result@pricePre  <- calcPrices(result,TRUE)
    result@pricePost <- calcPrices(result,FALSE)


    return(result)
}



