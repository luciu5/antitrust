setClass(
         Class = "PCAIDS",
         contains="AIDS",
         representation=representation(

         knownElast      = "numeric",
         knownElastIndex = "numeric"
         ),

         validity=function(object){




             nprods <- length(object@shares)

             if(length(object@knownElastIndex) != 1 ){stop("'knownElastIndex' must be length 1")}
             if(length(object@knownElast) != 1 ){stop("'knownElast' must be length 1")}
             if(length(object@mktElast) != 1 ){stop("'mktElast' must be length 1")}
             
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


     minD <- function(bknown){

       #enforce symmetry
       betas <- -diversion[idx,]/diversion[,idx] * bknown
       

       B = t(diversion * betas)



       m1 = bknown - shareKnown * (object@knownElast + 1 - shareKnown * (object@mktElast + 1))
       
       m2 <- diversion/t(diversion) - tcrossprod(1/diag(B), diag(B)) 
       m2 <-  m2[upper.tri(m2)]
       m2 <- m2[is.finite(m2) & m2 != 0]
       
       measure=c(m1,m2)

       return(sum(measure^2))
     }


     bestBeta <- optimize(minD,
                     upper=0,lower=-1e6)


     betas <- -diversion[idx,]/diversion[,idx] * bestBeta$minimum
     
     
     B = t(diversion * betas)
     

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
                   control.slopes,
                   control.equ,
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

    if(!missing(control.slopes)){
      result@control.slopes <- control.slopes
    }
    if(!missing(control.equ)){
      result@control.equ <- control.equ
    }
    
    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate Demand Slope Coefficients
    result <- calcSlopes(result)

    ## Solve Non-Linear System for Price Changes
    result@priceDelta <- calcPriceDelta(result,isMax=isMax,subset=subset,...)

    ## Calculate marginal cost
    result@mcPre <-  calcMC(result,TRUE)
    result@mcPost <- calcMC(result,FALSE)
    

    ## Calculate Pre and Post merger equilibrium prices
    ## These are equal to NA in pcaids
    result@pricePre  <- calcPrices(result,TRUE)
    result@pricePost <- calcPrices(result,FALSE)


    return(result)
}



