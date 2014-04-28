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
             if(object@knownElast>0 || object@mktElast > 0 ){
                 stop("'mktElast', 'knownElast' must be non-positive")}
             if(abs(object@knownElast) < abs(object@mktElast) ){
                  stop("'mktElast' must be less  than 'knownElast' in absolute value")}
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

     nprod    <- length(shares)
     
     idx      <- object@knownElastIndex
     
     shareKnown <- shares[idx]
     
     minD <- function(s){
    
        m1 <- s - shareKnown * (object@knownElast + 1 - shareKnown * (object@mktElast + 1))
        b <- s*diversion[idx,]/diversion[,idx]
        B  <- -diversion * b
        m2=rowSums(t(B)/diag(B)) # row sums of diversion matrix must sum to 0
        m3=B[upper.tri(B)] - t(B)[upper.tri(B)] # slope coefficients must be symmetry
        
        measure=c(m1,m2,m3)
        
       return(sum(measure^2))
     }
     
     
     bKnown <- optimize(minD,c(-1e12,0))$minimum
     
     b <- bKnown*diversion[idx,]/diversion[,idx]
     
     B <- -diversion * b
     
     if(!isTRUE(all.equal(B[upper.tri(B)],t(B)[upper.tri(B)]))){
       stop("Matrix of calibrated slope parameters is not symmetric. User-supplied 'diversions' may be inconsistent with model.")
     }
     
     if(!isTRUE(all.equal(sum(t(B)/diag(B)),0))){
       stop("Matrix of calibrated slope parameters yields a diversion matrix whose rows do not sum to 0. User-supplied 'diversions' may be inconsistent with model.")
     }
     
    
     dimnames(B) <- list(labels,labels)
     object@slopes <- B
     object@intercepts <- as.vector(shares - B%*%log(object@prices))
     names(object@intercepts) <- object@labels

     return(object)
 }

      )









pcaids <- function(shares,knownElast,mktElast=-1,
                   prices,diversions,
                   ownerPre,ownerPost,
                   knownElastIndex=1,
                   mcDelta=rep(0, length(shares)),
                   subset=rep(TRUE, length(shares)),
                   priceStart=runif(length(shares)),
                   isMax=FALSE,
                   labels=paste("Prod",1:length(shares),sep=""),
                   ...){



    if(missing(prices)){ prices <- rep(NA_real_,length(shares))}

    if(missing(diversions)){
        diversions <- tcrossprod(1/(1-shares),shares)
        diag(diversions) <- -1
    }

  ## Create PCAIDS container to store relevant data
    result <- new("PCAIDS",shares=shares,prices=prices,
                   quantities=shares, margins=shares,mcDelta=mcDelta,
                  subset=subset,
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
    result@priceDelta <- calcPriceDelta(result,isMax=isMax,subset=subset,...)


    ## Calculate Pre and Post merger equilibrium prices
    ## These are equal to NA in pcaids
    result@pricePre  <- calcPrices(result,TRUE)
    result@pricePost <- calcPrices(result,FALSE)


    return(result)
}



