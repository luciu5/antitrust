setClass(
         Class   = "Auction2ndLogit",
         contains="Logit",
         representation=representation(
           
           meanvalPost           = "numeric"
         ),
         prototype=prototype(
           priceStart  = numeric()
         ),
         validity=function(object){
           
           idx <- object@normIndex
           
           if(any(is.na(object@prices) & 
                  !is.na(object@mcDelta) &
                  object@mcDelta!=0)){
             stop("'prices'  cannot equal NA when 'mcDelta' supplied")
           }
           
           if(is.na(object@prices[idx]) &&
              any(!is.na(object@prices[-idx]))){
             stop("'prices' for index 'normIndex' cannot be NA")
           }
           
         }
        )


setMethod(
          f= "calcSlopes",
          signature= "Auction2ndLogit",
          definition=function(object){

              ## Uncover Demand Coefficents


              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              margins      <-  object@margins
              prices       <-  object@prices
              idx          <-  object@normIndex
              
              margins <- margins*prices

             
              ## Uncover price coefficient and mean valuation from margins and revenue shares


              nprods <- length(shares)

              
              ## identify which products have enough margin information
              ##  to impute  margins
              isMargin    <- matrix(margins,nrow=nprods,ncol=nprods,byrow=TRUE)
              isMargin[ownerPre==0]=0
              isMargin    <- !is.na(rowSums(isMargin))
              
              
              margshares <- margins * shares
              margshares[is.na(margshares)] = 0             
              
              firmShares <- drop(ownerPre %*% shares)
              firmMargins <- drop(ownerPre %*% (margshares))
              
              
              ## Minimize the distance between observed and predicted  ex Ante margins
              minD <- function(alpha){

                  measure <- margins + log((1 - firmShares + shares)/(1-firmShares))/( alpha * shares)
                  measure <- unique(measure[isMargin])
                  measure <- sum((measure)^2,na.rm=TRUE)

                  return(measure)
              }

              minAlpha <- optimize(minD,c(-1e6,0),
                                   tol=object@control.slopes$reltol)$minimum

              ## calculate costs conditional on a product winning
              marginsPre <- - log((1 - firmShares + shares)/(1-firmShares))/(minAlpha * shares)
              mcPre <- prices - marginsPre 
              
              if(is.na(idx)){
                idxShare <- 1 - object@shareInside
                idxPrice <- object@priceOutside
              }
              else{
                idxShare <- shares[idx]
                idxPrice <- mcPre[idx]
              }

              mcPre[is.na(mcPre)] <- 0             

              meanval <- log(shares) - log(idxShare) - minAlpha * (mcPre - idxPrice)
              
              
              names(meanval)   <- object@labels

              object@slopes    <- list(alpha=minAlpha,meanval=meanval)


              return(object)
          }
          )




## compute margins
setMethod(
  f= "calcMargins",
  signature= "Auction2ndLogit",
  definition=function(object,preMerger=TRUE,exAnte=FALSE,subset){
    
    
    nprods <- length(object@shares)
    
    if(missing(subset)){
      subset <- rep(TRUE,nprods)
    }
    
    if(!is.logical(subset) || length(subset) != nprods ){stop("'subset' must be a logical vector the same length as 'shares'")}
    
    margins <- rep(NA,nprods)
    
    if( preMerger) {
      owner  <- object@ownerPre}
    else{
      owner  <- object@ownerPost}
    
    
      owner <- owner[subset,subset]
    
     
    alpha <- object@slopes$alpha
    shares <- calcShares(object,preMerger=preMerger,revenue=FALSE)
    shares <- shares[subset]
    firmShares <- drop(owner %*% shares)
    margins[subset] <-  -log((1 - firmShares + shares)/(1-firmShares))/(alpha * shares) 
    
    if(exAnte){ margins <-  margins * shares}
    
    names(margins) <- object@labels
    
    return(as.vector(margins))
  }
  
)


## Create a method to recover marginal cost using
## demand parameters and supplied prices
setMethod(
  f= "calcMC",
  signature= "Auction2ndLogit",
  definition= function(object,preMerger=TRUE,exAnte=FALSE){
    
    prices <- object@prices
    ownerPre <- object@ownerPre
    shares <- object@shares
    alpha <- object@slopes$alpha
    
    firmShares <- drop(ownerPre %*% shares)
    
    marginPre <- - log((1 - firmShares + shares)/(1-firmShares))/(alpha * shares)
    
    mc <-  prices - marginPre
    
    if(!preMerger){
      mc <- mc*(1+object@mcDelta)
    }
    
    if(exAnte){mc <- mc*shares}
    
    names(mc) <- object@labels
    
    return(as.vector(mc))
  }
)

setMethod(
  f= "calcShares",
  signature= "Auction2ndLogit",
  definition=function(object,preMerger=TRUE,revenue=FALSE,subset){

    nprods <- length(object@shares)
    
    if(missing(subset)){
      subset <- rep(TRUE,nprods)
    }
    
    if(!is.logical(subset) || length(subset) != nprods ){stop("'subset' must be a logical vector the same length as 'shares'")}
    
    mc <- rep(NA,nprods)
    
    if( preMerger) {
      prices <- object@pricePre
      mc[subset] <- object@mcPre[subset]
      }
    else{
      prices <- object@pricePost
      mc[subset] <- object@mcPost[subset]
      }
    
    
    mc <- ifelse(is.na(mc) & subset, 0, mc)
    
    alpha    <- object@slopes$alpha
    meanval  <- object@slopes$meanval
    
    outVal <- ifelse(isTRUE(all.equal(object@shareInside,1)), 0, exp(alpha*object@priceOutside))
    
    shares <- exp(meanval + alpha*mc)
    shares <- shares/(outVal + sum(shares,na.rm=TRUE))
    
    if(revenue){
      res <- rep(NA,nprods)
      res[subset] <- prices[subset]*shares[subset]/sum(prices[subset]*shares[subset],object@priceOutside*(1-sum(shares[subset])))
      shares <- res
      }
    
    names(shares) <- object@labels
    
    return(shares)
    
    
    
    }
    )

setMethod(
 f= "calcPrices",
 signature= "Auction2ndLogit",
 definition=function(object,preMerger=TRUE,exAnte=FALSE,subset){

   nprods <- length(object@shares)
   if(missing(subset)){subset <- rep(TRUE,nprods)}   

   if(preMerger){mc <- object@mcPre}
   else{mc <- object@mcPost}
   
   margins <- calcMargins(object,preMerger,exAnte=FALSE,subset=subset)
   
   prices <- margins + mc 
   
   if(exAnte){ prices <- prices * calcShares(object,preMerger=TRUE,preMerger)}
   
   return(prices)
  }
 )




setMethod(
 f= "calcPricesHypoMon",
 signature= "Auction2ndLogit",
 definition=function(object,prodIndex){


     ownerMon <- object@ownerPre
     ownerMon[prodIndex,] <- 0
     ownerMon[,prodIndex] <- 0
     ownerMon <- ownerMon[prodIndex,prodIndex] <- 1
     
     object@ownerPre <- ownerMon
  
     pricesHM <- calcPrices(object,preMerger=TRUE)
     pricesHM <- pricesHM[prodIndex]
     names(pricesHM) <- object@labels[prodIndex]

     return(pricesHM)

 })


setMethod(
  f= "calcPriceDelta",
  signature= "Auction2ndLogit",
  definition=function(object,exAnte=FALSE,levels=TRUE){
    
    subset <- object@subset
    
    result <- calcPrices(object, preMerger=FALSE,exAnte=exAnte) -
              calcPrices(object, preMerger=TRUE,exAnte=exAnte)
    
    if(levels){
      marginDelta <- calcMargins(object, preMerger=FALSE,exAnte=exAnte) - 
                     calcMargins(object, preMerger=TRUE,exAnte=exAnte)
      result[is.na(result) & subset] <- marginDelta[is.na(result) & subset]
    }
    
    return(result)
  }
)
    

setMethod(
  f= "summary",
  signature= "Auction2ndLogit",
  definition=function(object,levels=TRUE,revenue=FALSE,...){
  
    callNextMethod(object,levels=levels,revenue=revenue,...)
  }
    
)  

## CMCR does not make sense in 2nd score auction
## delete
setMethod(
  f= "cmcr",
  signature= "Auction2ndLogit",
  definition=function(object){
   
    stop("'cmcr' is not defined for a 2nd-score auction")
  }
)


setMethod(
  f= "upp",
  signature= "Auction2ndLogit",
  definition=function(object){
    
    stop("'upp' is not defined for a 2nd-score auction")
  }
)




setMethod(
  f= "CV",
  signature= "Auction2ndLogit",
  definition=function(object){
    
    sharesPost <- calcShares(object, preMerger=FALSE)
    
    if(all(object@mcDelta==0,na.rm=TRUE) && any(is.na(object@pricePre))){
      
      priceDelta <- calcPriceDelta(object,levels=TRUE)
      
      result <- sum(priceDelta*sharesPost, na.rm = TRUE)
      return(result)
    }
    else{
    callNextMethod(object)  
    }
    
   
    
  })


auction2nd.logit <- function(prices,shares,margins,
                  ownerPre,ownerPost,
                  normIndex=ifelse(isTRUE(all.equal(sum(shares),1)),NA,1),
                  mcDelta=rep(0,length(prices)),
                  subset=rep(TRUE,length(prices)),
                  mcOutside = 0,
                  control.slopes,
                  labels=paste("Prod",1:length(prices),sep="")
                  ){


    ## Create Auction2ndLogit  container to store relevant data
    result <- new("Auction2ndLogit",prices=prices, shares=shares,
                  margins=margins,
                  normIndex=normIndex,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  mcDelta=mcDelta,
                  subset=subset,
                  priceOutside=mcOutside,
                  shareInside=ifelse(isTRUE(all.equal(sum(shares),1)),1,sum(shares)),
                  priceStart=rep(0,length(shares)),
                  labels=labels)

    if(!missing(control.slopes)){
      result@control.slopes <- control.slopes
    }
    
    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate Demand Slope Coefficients
    result <- calcSlopes(result)

    ## Calculate marginal cost
    result@mcPre <-  calcMC(result,TRUE)
    result@mcPost <- calcMC(result,FALSE)

    ## Solve Non-Linear System for Price Changes
    result@pricePre  <- calcPrices(result,preMerger=TRUE)
    result@pricePost <- calcPrices(result,preMerger=FALSE,subset=subset)

    return(result)

}

