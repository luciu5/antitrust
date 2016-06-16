
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

              if(is.na(idx)){
                  idxShare <- 1 - object@shareInside
                  idxPrice <- object@priceOutside
              }
              else{
                  idxShare <- shares[idx]
                  idxPrice <- prices[idx]
               }

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
                  
                  measure <- margins - log(1/(1-firmShares))/alpha
                  measure <- sum(measure^2,na.rm=TRUE)

                  return(measure)
              }

              minAlpha <- optimize(minD,c(-1e6,0))$minimum

              marginsPre <- - log(1/(1-firmShares))/minAlpha
              mcPre <- prices - marginsPre 

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
      subset <- rep(TRUE,nprods)
    }
    
    if(!is.logical(subset) || length(subset) != nprods ){stop("'subset' must be a logical vector the same length as 'shares'")}
    
    mc <- margins <- rep(NA,nprods)
    
    if( preMerger) {
      mc[subset] <- object@mcPre[subset] 
      object@pricePre <- mc
      owner  <- object@ownerPre}
    else{
      mc[subset] <- object@mcPost[subset]
      object@pricePost <- mc
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
    
    
    marginPre <- calcMargins(object,TRUE)
    
    mc <-  prices - marginPre
    
    if(!preMerger){
      mc <- mc*(1+object@mcDelta)
    }
    
    names(mc) <- object@labels
    
    return(as.vector(mc))
  }
)


setMethod(
 f= "calcPrices",
 signature= "Auction2ndLogit",
 definition=function(object,preMerger=TRUE,subset){


   if(preMerger){mc <- object@mcPre}
   else{mc <- object@mcPost}
   
   margins <- calcMargins(object,PreMerger,subset)
   
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




Auction2ndLogit <- function(prices,shares,margins,
                  ownerPre,ownerPost,
                  normIndex=ifelse(sum(shares)<1,NA,1),
                  mcDelta=rep(0,length(prices)),
                  subset=rep(TRUE,length(prices)),
                  priceOutside = 0,
                  priceStart = prices,
                  isMax=FALSE,
                  labels=paste("Prod",1:length(prices),sep=""),
                  ...
                  ){


    ## Create Auction2ndLogit  container to store relevant data
    result <- new("Auction2ndLogit",prices=prices, shares=shares,
                  margins=margins,
                  normIndex=normIndex,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  mcDelta=mcDelta,
                  subset=subset,
                  priceOutside=priceOutside,
                  priceStart=priceStart,shareInside=sum(shares),
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
    result@pricePre  <- calcPrices(result,preMerger=TRUE,isMax=isMax,...)
    result@pricePost <- calcPrices(result,preMerger=FALSE,isMax=isMax,subset=subset,...)

    return(result)

}

