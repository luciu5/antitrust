setClass(
    Class   = "CES",
    contains="Logit",
    prototype=prototype(
    priceOutside=0
    )
    )


setMethod(
          f= "calcSlopes",
          signature= "CES",
          definition=function(object){

              ## Uncover Demand Coefficents


              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              margins      <-  object@margins
              prices       <-  object@prices
              idx          <-  object@normIndex
              shareInside  <-  object@shareInside
              insideSize   <-  object@insideSize

              ## uncover Numeraire Coefficients
              if(shareInside <= 1 && shareInside>0) {alpha <- 1/shareInside - 1}
              else{alpha <- NULL}

              ## if sum of shares is less than 1, add numeraire
               if(is.na(idx)){
                  idxShare <- 1 - sum(shares)
                  idxPrice <- object@priceOutside
              }
              else{
                  idxShare <- shares[idx]
                  idxPrice <- prices[idx]
               }


              ## Uncover price coefficient and mean valuation from margins and revenue shares


              nprods <- length(shares)

              ## identify which products have enough margin information
              ##  to impute Bertrand margins
              isMargin    <- matrix(margins,nrow=nprods,ncol=nprods,byrow=TRUE)
              isMargin[ownerPre==0]=0
              isMargin    <- !is.na(rowSums(isMargin))





              ## Minimize the distance between observed and predicted margins
              minD <- function(gamma){


                  elasticity <- (gamma - 1 ) * matrix(shares,ncol=nprods,nrow=nprods)
                  diag(elasticity) <- -gamma + diag(elasticity)

                  elasticity <- elasticity[isMargin,isMargin]
                  shares     <- shares[isMargin]
                  ownerPre   <- ownerPre[isMargin,isMargin]
                  margins    <- margins[isMargin]

                  #marginsCand <- -1 * as.vector(MASS::ginv(elasticity * ownerPre) %*% (shares * diag(ownerPre))) / shares
                  #measure <- sum((margins - marginsCand)^2,na.rm=TRUE)
                   FOC <- (shares * diag(ownerPre)) + (elasticity * ownerPre) %*% (shares * margins)
                   measure<-sum(FOC^2,na.rm=TRUE)

                  return(measure)
              }

              minGamma <- optimize(minD,c(1,1e6),
                                   tol=object@control.slopes$reltol)$minimum


              meanval <- log(shares) - log(idxShare) + (minGamma - 1) * (log(prices) - log(idxPrice))
              meanval <- exp(meanval)

              names(meanval)   <- object@labels

              object@slopes    <- list(alpha=alpha,gamma=minGamma,meanval=meanval)
              object@priceOutside <- idxPrice
              object@mktSize <- insideSize*(1+alpha)


              return(object)
          }
          )

## compute product revenues
setMethod(
  f= "calcRevenues",
  signature= "CES",
  definition=function(object,preMerger=TRUE, market =FALSE){
    
    shares <- calcShares(object, preMerger = preMerger, revenue=TRUE)
    
    mktSize <- object@mktSize
    
    
    res <- shares * mktSize
    
    if(market){return(sum(res))}
    
    else{return(res)}
    
  })


setMethod(
  f= "calcQuantities",
  signature= "CES",
  definition=function(object,preMerger=TRUE, market=FALSE){
    
    if(preMerger){ prices <- object@pricePre}
    else{          prices <- object@pricePost}
    
    mktSize <- object@mktSize
    
    shares <- calcShares(object, preMerger= preMerger, revenue = TRUE) / prices
    
    if(market) shares <- sum(shares,na.rm=TRUE)
    
    return(mktSize*shares )
    
    
  })

setMethod(
 f= "calcShares",
 signature= "CES",
 definition=function(object,preMerger=TRUE,revenue=FALSE){




     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}


     gamma    <- object@slopes$gamma
     meanval  <- object@slopes$meanval
     priceOutside <- object@priceOutside

     #outVal <- ifelse(object@shareInside<1, object@priceOutside^(1-gamma), 0)
     outVal <- ifelse(is.na(object@normIndex), priceOutside^(1-gamma), 0)
     
     shares <- meanval*prices^(1-gamma)
     shares <- shares/(sum(shares,na.rm=TRUE) + outVal)

     ##transform revenue shares to quantity shares
     if(!revenue){shares <- (shares/prices)/sum((1-sum(shares,na.rm=TRUE))/priceOutside,shares/prices,na.rm=TRUE)}

     names(shares) <- object@labels

     return(as.vector(shares))

}
 )





setMethod(
 f= "elast",
 signature= "CES",
 definition=function(object,preMerger=TRUE,market=FALSE){

   
     gamma    <- object@slopes$gamma

     shares <-  calcShares(object,preMerger,revenue=TRUE)


      if(market){

        if(preMerger){ prices <- object@pricePre}
        else{          prices <- object@pricePost}
        
          avgPrice <- sum(prices*shares)/sum(shares)
          
          elast <- ( 1 - gamma  )  * (1 - sum(shares)) * avgPrice 
          
         }

     else{

         nprods <-  length(shares)
         elast <- (gamma - 1 ) * matrix(shares,ncol=nprods,nrow=nprods,byrow=TRUE)
         diag(elast) <- -gamma + diag(elast)

         dimnames(elast) <- list(object@labels,object@labels)
     }
      return(elast)

}
 )



setMethod(
          f= "CV",
          signature= "CES",
          definition=function(object){

              alpha       <- object@slopes$alpha
              mktSize  <- object@mktSize 
  
             if(is.null(alpha)) stop(" Sum of 'shares' must be less than 1  calculate Compensating Variation")

              gamma       <- object@slopes$gamma
              meanval     <- object@slopes$meanval
              shareInside <- object@shareInside


              outVal <- ifelse(shareInside<1, 1, 0)

              VPre  <- sum(meanval * (object@pricePre / object@priceOutside)^(1-gamma),na.rm=TRUE) + outVal
              VPost <- sum(meanval * (object@pricePost/ object@priceOutside)^(1-gamma),na.rm=TRUE) + outVal

              result <- log(VPost/VPre) / ((1+alpha)*(1-gamma))
              
              result <- exp(result) - 1

              if(is.na(mktSize)){
                  warning("'insideSize' is NA. Calculating CV as a percentage of (aggregate) expenditure")
                  return(result*100)}

              else{
                  
                  return(mktSize*result)
              }


 })


ces <- function(prices,shares,margins,
                ownerPre,ownerPost,
                normIndex=ifelse(sum(shares)<1,NA,1),
                insideSize = NA_real_,
                mcDelta=rep(0,length(prices)),
                subset=rep(TRUE,length(prices)),
                priceOutside=1,
                priceStart = prices,
                isMax=FALSE,
                control.slopes,
                control.equ,
                labels=paste("Prod",1:length(prices),sep=""),
                ...
                ){




    ## Create CES  container to store relevant data
    result <- new("CES",prices=prices, shares=shares, margins=margins,
                  normIndex=normIndex,
                  mcDelta=mcDelta,
                  insideSize = insideSize,
                  subset=subset,
                  priceOutside=priceOutside,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  priceStart=priceStart,
                  shareInside= sum(shares),labels=labels)

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

    ## Calculate marginal cost
    result@mcPre <-  calcMC(result,TRUE)
    result@mcPost <- calcMC(result,FALSE)


    ## Solve Non-Linear System for Price Changes
    result@pricePre  <- calcPrices(result,preMerger=TRUE,isMax=isMax,...)
    result@pricePost <- calcPrices(result,preMerger=FALSE,isMax=isMax,subset=subset,...)

    return(result)

}

