setClass(
         Class   = "Auction2ndLogit",
         contains="Logit",
         prototype=prototype(
           priceStart  = numeric()
         )
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
              
  

             
              ## Uncover price coefficient and mean valuation from margins and revenue shares


              nprods <- length(shares)

              
             
              firmShares <- drop(ownerPre %*% shares)
              
              
              ## Minimize the distance between observed and predicted  ex Ante margins
              minD <- function(alpha){

                  measure <- margins + log(1 / (1-firmShares))/( alpha * firmShares)
                 
                  measure <- sum((measure)^2,na.rm=TRUE)

                  return(measure)
              }

              minAlpha <- optimize(minD,c(-1e6,0),
                                   tol=object@control.slopes$reltol)$minimum

              ## calculate costs conditional on a product winning
              marginsPre <- - log(1 /(1-firmShares))/(minAlpha * firmShares)
              
              if(is.na(idx)){
                idxShare <- 1 - object@shareInside
                
              }
              else{
                idxShare <- shares[idx]
               
              }

             
              meanval <- log(shares) - log(idxShare)
              
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
    margins[subset] <-  -log(1/(1-firmShares))/(alpha * firmShares) 
    
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
    
    marginPre <- calcMargins(object, preMerger = TRUE)
    
    mc <-  prices - marginPre
    
    if(!preMerger){
      mc <- mc + object@mcDelta
    }
    
    if(exAnte){mc <- mc * calcShares(object,preMerger=preMerger)}
    
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
    
    
    idx <- object@normIndex
    alpha    <- object@slopes$alpha
    meanval  <- object@slopes$meanval
    
    
    
    if(is.na(idx)){
      outVal <- 1
      mcDeltaOut <- object@priceOutside
    }
    
    else{
      outVal <- 0
      mcDeltaOut <- object@mcDelta[idx]
    }
  
    
    if( preMerger) {
      prices <- object@pricePre
      
      }
    else{
      prices <- object@pricePost
      meanval <- meanval + alpha * (object@mcDelta - mcDeltaOut)
      }
    
    
    
    shares <- exp(meanval)
    shares <- shares/(outVal + sum(shares,na.rm=TRUE))
    
    if(revenue){
      res <- rep(NA,nprods)
      res[subset] <- prices[subset]*shares[subset]/sum(prices[subset]*shares[subset],mcDeltaOut*(1-sum(shares[subset])))
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

   if(preMerger){
     owner <- object@ownerPre
     mc <- object@mcPre
     }
   else{
     owner <- object@ownerPost
     mc <- object@mcPost}
   
   margins <- calcMargins(object,preMerger,exAnte=FALSE,subset=subset)
   
   prices <- margins + mc 
   
   if(exAnte){ 
     
     
     prices <- prices * calcShares(object, preMerger=preMerger,revenue=FALSE)
     }
   
   
   names(prices) <- object@labels
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
    
    mcDelta <- object@mcDelta
    
    if(exAnte){sharesPost <- calcShares(object, preMerger=FALSE)}
    else{sharesPost <- rep(1,length(subset))}
    
    if(levels){
      result <- calcMargins(object, preMerger=FALSE,exAnte=exAnte,subset=subset) + mcDelta*sharesPost -
        calcMargins(object, preMerger=TRUE,exAnte=exAnte)
    }
    else{result <- result/calcPrices(object,preMerger = TRUE, exAnte = exAnte )}
    
    names(result) <- object@labels
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
   
    stop("'cmcr' is currently not available")
  }
)


setMethod(
  f= "upp",
  signature= "Auction2ndLogit",
  definition=function(object){
    
    shares <- calcShares(object,preMerger=TRUE,revenue=FALSE)
    diversion <- t(tcrossprod(shares,1/(1-shares)))
    diag(diversion) <- -1
    
    mcDelta <- object@mcDelta
    
    margins <- calcMargins(object,preMerger=TRUE)
    
    isParty <- abs(object@ownerPost - object@ownerPre)
    
    gross <- margins * shares * diversion * isParty
    gross <- sum(gross) 
    
    isParty <- rowSums(isParty) > 0 
    dpp <-  sum(mcDelta * diversion[,isParty], na.rm=TRUE)


    result <- gross + dpp
    
    
    return(result)
    
  }
)




setMethod(
  f= "CV",
  signature= "Auction2ndLogit",
  definition=function(object){
    
  meanvalPre <-   meanvalPost <- object@slopes$meanval
  alpha <- object@slopes$alpha
  sigma <- -1/alpha
  
  mcDelta <- object@mcDelta
  mcDelta[is.na(mcDelta)] <- 0
  
  idx <- object@normIndex
  if(is.na(idx)){
    mcDeltaOut <- object@priceOutside
  }
  
  else{
    mcDeltaOut <- object@mcDelta[idx]
  }

  meanvalPost <-   meanvalPre + alpha*(mcDelta - mcDeltaOut)
  
  sharePre <- calcShares(object,preMerger=TRUE)
  sharePost <- calcShares(object,preMerger=FALSE)
  
  firmSharePre <- drop(object@ownerPre %*% sharePre)  
  firmSharePost <- drop(object@ownerPost %*% sharePost)
  
  incValPre <- log(sum(exp(meanvalPre)))
  incValPost <- log(sum(exp(meanvalPost)))
  
  firmIncVal <- function(x, preMerger = TRUE){
    x <- x == 0
    
    if(preMerger){return(log(sum(exp(meanvalPre[x]))))}
    else{return(log(sum(exp(meanvalPost[x]))))}
    
    }
  
  csPre <- apply(object@ownerPre, 1,firmIncVal,preMerger=TRUE )
  csPost <- apply(object@ownerPost, 1,firmIncVal,preMerger=FALSE )  
  
  csPre <- (sigma/firmSharePre) * (csPre - (1-firmSharePre)*incValPre)
  
  csPost <- (sigma/firmSharePost) * (csPost - (1-firmSharePost)*incValPost)
  
  result <-  mcDeltaOut + sum(csPre*sharePre) - sum(csPost*sharePost)
  
  return(result)
  })


auction2nd.logit <- function(prices,shares,margins,
                  ownerPre,ownerPost,
                  normIndex=ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE)),1, NA),
                  mcDelta=rep(0,length(prices)),
                  subset=rep(TRUE,length(prices)),
                  mcDeltaOutside=0,
                  control.slopes,
                  labels=paste("Prod",1:length(prices),sep="")
                  ){

if(missing(prices)){prices <- rep(NA_integer_, length(shares))}
    ## Create Auction2ndLogit  container to store relevant data
    result <- new("Auction2ndLogit",prices=prices, shares=shares,
                  margins=margins,
                  normIndex=normIndex,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  mcDelta=mcDelta,
                  subset=subset,
                  priceOutside=mcDeltaOutside,
                  shareInside=ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE)),1,sum(shares)),
                  priceStart=rep(0,length(shares)),
                  labels=labels,
                  cls = "Auction2ndLogit")

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

