
setClass(
         Class   = "Auction2ndLogit",
         contains="Logit"
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
              ##  to impute model margins
              isMargin    <- matrix(margins,nrow=nprods,ncol=nprods,byrow=TRUE)
              isMargin[ownerPre==0]=0
              isMargin    <- !is.na(rowSums(isMargin))


              firmShares <- drop(ownerPre %*% shares)
              
              ## Minimize the distance between observed and predicted margins
              minD <- function(alpha){

                 
                  firmShares <- firmShares[isMargin]
                  margins    <- margins[isMargin]
                  
                  measure <- margins + log(1/(1-firmShares))/alpha
                  measure <- sum((measure/prices)^2,na.rm=TRUE)

                  return(measure)
              }

              minAlpha <- optimize(minD,c(-1e6,0))$minimum

              marginsPre <- - log(1/(1-firmShares))/minAlpha
              mcPre <- prices - marginsPre 
              
              if(is.na(idx)){
                idxShare <- 1 - object@shareInside
                idxPrice <- object@priceOutside
              }
              else{
                idxShare <- shares[idx]
                idxPrice <- mcPre[idx]
              }
              

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
  definition=function(object,preMerger=TRUE,subset){
    
    
    nprods <- length(object@shares)
    
    if(missing(subset)){
      subset <- object@subset
    }
    
    if(!is.logical(subset) || length(subset) != nprods ){stop("'subset' must be a logical vector the same length as 'shares'")}
    
    margins <- rep(NA,nprods)
    
    if( preMerger) {
      owner  <- object@ownerPre}
    else{
      owner  <- object@ownerPost}
    
    
      owner <- owner[subset,subset]
    
     
    alpha <- object@slopes$alpha
    shares <- calcShares(object,preMerger,revenue=FALSE)
    shares <- shares[subset]
    firmShares <- drop(owner %*% shares)
    margins[subset] <-  -log(1/(1-firmShares))/alpha 
    
    names(margins) <- object@labels
    
    return(as.vector(margins))
  }
  
)


## Create a method to recover marginal cost using
## demand parameters and supplied prices
setMethod(
  f= "calcMC",
  signature= "Auction2ndLogit",
  definition= function(object,preMerger=TRUE){
    
    prices <- object@prices
    ownerPre <- object@ownerPre
    shares <- object@shares
    alpha <- object@slopes$alpha
    
    firmShares <- drop(ownerPre %*% shares)
    
    marginPre <- - log(1/(1-firmShares))/alpha
    
    mc <-  prices - marginPre
    
    if(!preMerger){
      mc <- mc*(1+object@mcDelta)
    }
    
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
      subset <- object@subset
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
    
    
    alpha    <- object@slopes$alpha
    meanval  <- object@slopes$meanval
    
    outVal <- ifelse(object@shareInside<1, exp(alpha*object@priceOutside), 0)
    
    shares <- exp(meanval + alpha*mc)
    shares <- shares/(outVal+ sum(shares,na.rm=TRUE))
    
    if(revenue){shares <- prices*shares/sum(prices*shares,object@priceOutside*(1-sum(shares,na.rm=TRUE)),na.rm=TRUE)}
    
    names(shares) <- object@labels
    
    return(shares)
    
    
    
    }
    )

setMethod(
 f= "calcPrices",
 signature= "Auction2ndLogit",
 definition=function(object,preMerger=TRUE,subset){


   if(preMerger){mc <- object@mcPre}
   else{mc <- object@mcPost}
   
   margins <- calcMargins(object,preMerger,subset)
   
   prices <- margins + mc 
   
   return(prices)
  }
 )




setMethod(
 f= "calcPricesHypoMon",
 signature= "Auction2ndLogit",
 definition=function(object,prodIndex){


     mc       <- object@mcPre[prodIndex]
     pricePre <- object@pricePre

     calcMonopolySurplus <- function(priceCand){

         pricePre[prodIndex] <- priceCand #keep prices of products not included in HM fixed at premerger levels
         object@pricePre     <- pricePre
         sharesCand          <- calcShares(object,TRUE,revenue=FALSE)

         surplus             <- (priceCand-mc)*sharesCand[prodIndex]

         return(sum(surplus,na.rm=TRUE))
     }


     maxResult <- optim(object@prices[prodIndex],calcMonopolySurplus,
                              method = "L-BFGS-B",lower = 0,
                              control=list(fnscale=-1))

     pricesHM <- maxResult$par

     #priceDelta <- pricesHM/pricePre[prodIndex] - 1
     #names(priceDelta) <- object@labels[prodIndex]
     names(pricesHM) <- object@labels[prodIndex]

     return(pricesHM)

 })


## account for the fact that calcMargins method for 2nd price 
## is in dollars rather than %
setMethod(
  f= "cmcr",
  signature= "Auction2ndLogit",
  definition=function(object){
   
    isParty <- rowSums( abs(object@ownerPost - object@ownerPre) ) > 0
    
    ##Compute pre-merger margins
    marginPre  <- calcMargins(object,TRUE)
    
    
    ##compute post-merger margins evaluated at pre-merger prices
    object@ownerPre <- object@ownerPost
    marginPost <- calcMargins(object,TRUE)
    
    cmcr <- (marginPost - marginPre)/(object@pricePre - marginPre)
    names(cmcr) <- object@labels
    
    cmcr <- cmcr[isParty]
    return(cmcr * 100)
  }
)

setMethod(
  f= "CV",
  signature= "Auction2ndLogit",
  definition=function(object){
    
    priceDelta <-  object@pricePost - object@pricePre 
    sharesPost <- calcShares(object)
    
    result <- sum(priceDelta*sharesPost, na.rm = TRUE)
    
    return(result)
    
  })


auction2nd.logit <- function(prices,shares,margins,
                  ownerPre,ownerPost,
                  normIndex=ifelse(sum(shares)<1,NA,1),
                  mcDelta=rep(0,length(prices)),
                  subset=rep(TRUE,length(prices)),
                  mcOutside = 0,
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
                  shareInside=sum(shares),
                  priceStart=rep(0,length(shares)),
                  labels=labels)

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

