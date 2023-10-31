#' @title Methods for Calculating Diagnostics
#' @name Margins-Methods
#' @docType methods
#'
#' @aliases calcMargins
#' calcMargins,ANY-method
#' calcMargins,AIDS-method
#' calcMargins,Bertrand-method
#' calcMargins,VertBargBertLogit-method
#' calcMargins,LogitCap-method
#' calcMargins,LogitDebt-method
#' calcMargins,Auction2ndLogit-method
#' calcMargins,Auction2ndLogitNests-method
#' calcMargins,Cournot-method
#' calcMargins,BargainingLogit-method
#' calcMargins,Bargaining2ndLogit-method
#'
#' @description Computes equilibrium product margins assuming that firms are playing a
#' Nash-Bertrand, Cournot, 2nd Score Auction, or Bargaining game. For "LogitCap", assumes firms are
#' playing a Nash-Bertrand or Cournot game with capacity constraints.
#'
#' @param object An instance of one of the classes listed above.
#' @param preMerger If TRUE, returns pre-merger outcome. If
#' FALSE, returns post-merger outcome.  Default is TRUE.
#' @param level IF TRUE, return margins in dollars. If FALSE, returns 
#' margins in proportions. Default for most classes is FALSE.
#' @param exAnte If \sQuote{exAnte} equals TRUE then the
#' \emph{ex ante} expected result for each firm is produced, while FALSE produces the
#' expected result conditional on each firm winning the auction. Default is FALSE.
#'
#' @include CostMethods.R
#' @keywords methods
NULL

setGeneric (
  name= "calcMargins",
  def=function(object,...){standardGeneric("calcMargins")}
)

## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Bertrand",
  definition=function(object,preMerger=TRUE, level=FALSE){



    if( preMerger) {

      prices <- object@pricePre
      owner  <- object@ownerPre
      
    }

    else{
      prices <- object@pricePost
      owner  <- object@ownerPost

    }

    revenue<- calcShares(object,preMerger,revenue=TRUE)
    
    elast <-  elast(object,preMerger)
    
    margins <-  try(-1 * as.vector(solve(t(elast)*owner) %*% (revenue * diag(owner))) / revenue,silent=TRUE)
    if(any(class(margins) == "try-error")){margins <- -1 * as.vector(MASS::ginv(t(elast)*owner) %*% (revenue * diag(owner))) / revenue}
    
    
    
    if(level) {margins <- margins * prices }
    
    names(margins) <- object@labels

    return(as.vector(margins))
  }

)


## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Bargaining2ndLogit",
  definition=function(object,preMerger=TRUE,exAnte=FALSE,level=TRUE){
  
    if( preMerger) {
      
      barg <- object@bargpowerPre
      
      
    }
    
    else{
      
      barg <- object@bargpowerPost
      
    }
      
    (1-barg)*callNextMethod()
    
  })

## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "BargainingLogit",
  definition=function(object,preMerger=TRUE, level=FALSE){
    
    
    alpha <- object@slopes$alpha
    
    if( preMerger) {
      
      prices <- object@pricePre
      owner  <- object@ownerPre
      barg <- object@bargpowerPre
      
      
    }
    
    else{
      prices <- object@pricePost
      owner  <- object@ownerPost
      barg <- object@bargpowerPost
      
    }
    
    
    barg <- barg/(1-barg) #relative bargaining
    
    nprods <- length(prices)
    
    shares <- calcShares(object,preMerger,revenue=FALSE)
    
    div <- shares/(1-shares)
    
    
    #diag(owner) <- -1*diag(owner)
    
    margins <- -owner * shares
    diag(margins) <- diag(owner) +  diag(margins)
    margins <- solve(t(margins))
    
    margins <-  as.vector(margins %*% ((log(1-shares)*diag(owner))/(alpha*(barg*div - 
                                                               log(1-shares)))))
    
    
    if(!level) {margins <- margins / prices }
    
    names(margins) <- object@labels
    
    return(as.vector(margins))
  }
  
)


#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "VertBargBertLogit",
  definition=function(object,preMerger=TRUE, level=FALSE){
    
    #check if class is 2nd score Logit
    is2nd <- any(grepl("2nd",class(object)))
    
    up <- object@up
    down <- object@down
    alpha <- down@slopes$alpha
    
    chain_level <- object@chain_level
    
    ## Set downstream margin equal to pre-merger value when chain_level=="wholesaler
    preMerger_chain <- preMerger
    if(chain_level=="wholesaler") {
      preMerger_chain <- TRUE
       if(!preMerger) down@ownerPre <- ownerToMatrix(down,preMerger= TRUE)
    }
  
    #is2nd <- grepl("2nd",class(object))
  
    ownerUp <- ownerToMatrix(up,preMerger = preMerger)
    ownerDown <- ownerToMatrix(down,preMerger= preMerger)
    
    nprods <- nrow(ownerDown)
    
    if( preMerger) {
      bargparm <- up@bargpowerPre
       
      #if(length(up@pricePre) == 0 ){
      #priceUp <- up@prices
      #}
      
      #else{
        priceUp <- up@pricePre
       # }
      #if(length(down@pricePre) == 0 ){
       # priceDown <- down@prices
      #  down@pricePre <- priceDown
      #}
      #else{
        priceDown <- down@pricePre
        
       # }
     
      ownerDownLambda <- object@ownerDownLambdaPre
      ownerUpLambda <- object@ownerUpLambdaPre
      ownerVDown <- object@ownerDownPre
      
      
      #down@mcPre <- down@mcPre + priceUp
      down@ownerPre <- ownerDown
      
    }
    
    else{
      bargparm <- up@bargpowerPost
      
      
      priceUp <- up@pricePost
      priceDown <- down@pricePost
      
      ownerDownLambda <- object@ownerDownLambdaPost
      ownerUpLambda <- object@ownerUpLambdaPost
      ownerVDown <- object@ownerDownPost
      
      down@ownerPost <- ownerDown
      
    }
    
    shareDown <- calcShares(object,preMerger=preMerger, revenue=FALSE)
    elast <-  -alpha*tcrossprod(shareDown)
    diag(elast) <- alpha*shareDown + diag(elast)
    elast.inv <- try(solve(ownerDown * elast),silent=TRUE)
    if(any(class(elast.inv) == "try-error")){elast.inv <- MASS::ginv(ownerDown * elast)}
    
    if(is2nd){
      alpha <- down@slopes$alpha
      down@slopes$meanval <- down@slopes$meanval + alpha *(priceUp - down@priceOutside) 
    }
    
    marginsDown <- calcMargins(down, preMerger=preMerger_chain,level=TRUE)
    
    div <- tcrossprod(1/(1-shareDown),shareDown)*shareDown
    diag(div) <- -shareDown
    div <- as.vector(div)
    
    
    #marginsDown <-  marginsDown - elast.inv %*% ( (ownerVDown * elast) %*% (priceCandUp-mcUp) )
    
    #marginsUp <-  as.vector(solve(ownerUpLambda * div) %*% (((ownerDownLambda * div) %*% (marginsDown)))) 
    
    if(chain_level == "retailer" && !preMerger){ 
      #marginsUp <- calcMargins(object,preMerger=TRUE,level=TRUE)$up
      marginsUp <- up@pricePre - up@mcPost
      }
    else{
    marginsUpPart <-  try(solve(ownerUpLambda * div) %*% (ownerDownLambda * div) ,silent=TRUE)
    if(any(class(marginsUpPart) == "try-error")){
      marginsUpPart <-  MASS::ginv(ownerUpLambda * div) %*% (ownerDownLambda * div) 
    }
    
    marginsUp <- try(solve(diag(nprods) + (marginsUpPart %*% elast.inv %*%  (ownerVDown * elast))),silent=TRUE)
    
    if(any(class(marginsUp) == "try-error")){
    marginsUp <- MASS::ginv(diag(nprods) + (marginsUpPart %*% elast.inv %*%  (ownerVDown * elast)))
    }
    marginsUp <- drop(marginsUp %*% marginsUpPart %*% marginsDown)
    
    }

    if(!is2nd){
      
      if(chain_level == "wholesaler" && !preMerger){ 
        #marginsDown <- calcMargins(object,preMerger=TRUE,level=TRUE)$down
        marginsDown <- down@pricePre - down@mcPost - priceUp
        } 
      else{ marginsDown <-  marginsDown - elast.inv %*% ( (ownerVDown * elast) %*% marginsUp )}
      
      }
  
    if(!level) {
      marginsUp <- marginsUp/priceUp
      marginsDown <- marginsDown/priceDown
    }
     
      
      
    
    names(marginsUp) <- up@labels
    names(marginsDown) <- down@labels
    
    return(list(up=as.vector(marginsUp), down = as.vector(marginsDown))
           )
  }
  
)



#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Auction2ndCap",
  definition=function(object,preMerger=TRUE,exAnte=TRUE,level=FALSE){

    result <- calcProducerSurplus(object,preMerger=preMerger,exAnte=exAnte)
    
    if(!level) result <- result / calcPrices(object,preMerger=preMerger,exAnte=exAnte)
    
    return(result)
  })

## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE, level=FALSE){


    if(preMerger){
      prices <- object@pricePre
    }
    else{prices <- object@pricePost}

    mc <- calcMC(object,preMerger = preMerger)
    prices <- matrix(prices, ncol=length(prices), nrow=length(mc),byrow=TRUE)



    margin <- prices - mc
    if(!level){margin <- margin/prices}

    dimnames(margin) <- object@labels
    return(margin)
  }

)


#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "AIDS",
  definition=function(object,preMerger=TRUE,level=FALSE){

    priceDelta <- object@priceDelta
    ownerPre   <- object@ownerPre
    shares     <- calcShares(object,TRUE)

    elastPre <-  t(elast(object,TRUE))
    
    elastInv <- try(solve(elastPre * ownerPre),silent=TRUE)
    if(any(class(elastInv)=="try-error")){elastInv <- MASS::ginv(elastPre * ownerPre)}
    
    marginPre <-  -1 * as.vector(elastInv %*% (shares * diag(ownerPre))) / shares

    if(preMerger){
      names(marginPre) <- object@labels
      
      if(level) {marginPre <- marginPre * object@pricePre}
      return(marginPre)}

    else{

      marginPost <- 1 - ((1 + object@mcDelta) * (1 - marginPre) / (priceDelta + 1) )
      names(marginPost) <- object@labels
      if(level) {marginPost <- marginPost * object@pricePost}
      return(marginPost)
    }

  }
)

## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "LogitCap",
  definition=function(object,preMerger=TRUE,level=FALSE){

    margins <- object@margins #capacity-constrained margins not identified -- use supplied margins

    if( preMerger) {
      capacities <- object@capacitiesPre
      prices <- object@pricePre

    }
    else{

      capacities <- object@capacitiesPost
      prices <- object@pricePost
    }



    quantities <- calcQuantities(object, preMerger=TRUE)
    constrained <-  abs(capacities - quantities) < 1e-5

    owner  <- object@ownerPre
    revenue<- calcShares(object,preMerger,revenue=TRUE)[!constrained]
    elast <-  elast(object,preMerger)
    margins[!constrained] <-  -1 * as.vector(MASS::ginv(t(elast*owner)[!constrained,!constrained]) %*% revenue) / revenue



    names(margins) <- object@labels

    if(level) {margins <- margins * prices} 
    return(as.vector(margins))
  }

)


## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "LogitDebt",
  definition=function(object,preMerger=TRUE,level=FALSE){
    
    
    if( preMerger) {
      firmVec <- object@ownerPre # 1:nFirms
      debt    <- object@debtPre
      debtDelta <- rep(0,length(debt))
      prices <- object@pricePre
      owner <- object@ownerPre
      mc <- object@mcPost
    }
    else{
      firmVec <- object@ownerPost #2:nFirms
      debt    <- object@debtPost
      prices <- object@pricePre
      owner <- object@ownerPost
      mc <- object@mcPost
    }
    
    
    nFirms <- object@nFirms 
    nMarkets <- object@nMarkets
    slopes <- object@slopes
    mval <- object@mval
    M <- object@mktSize
    shareOutParm <- object@shareOutParm
    outPrice=object@priceOutside
    focal=object@focal 
    z_crit=mkt$z_crit

    hFocalIntegral <- object@hFocalIntegral
    gFocalIntegral <- object@gFocalIntegral
    tOtherIntegral <- object@tOtherIntegral
    rOtherIntegral <- object@rOtherIntegral
    
    
    lowerB <- rep(0,nMarkets-1)
    upperB <- rep(1,nMarkets-1)
    
    
    
    bDens <- function(u,...) {dbeta(u, shape1 = shareOutParm[1],shape2 = shareOutParm[2],...)}
    bProb <- function(u,...) {pbeta(u, shape1 = shareOutParm[1],lower.tail=TRUE,shape2 = shareOutParm[2],...)}
    
    integrateDerivOtherFun <- function(u,f,pre=preMerger){s0=s0Focal(u,f,pre); u*(1-u)*bProb(s0)}
    integrateInsideOtherFun <- function(u,f,pre=preMerger){s0=s0Focal(u,f,pre); (1-u)*bProb(s0)}
    
    s0Focal <- function(sOut=rep(NA,nMarkets),i,pre=preMerger,n=nMarkets,l=focal){
      ## allow firm index f to be a vector of length>1.
      ## choose firm 2's product in market 1 to be the reference good
      if(!pre && any(i %in% 1:2)){i=2:1}
      if(n==1){ res <- 1- sum(debt[i])/sum(profitsCond[,i]); return(res)}
      res <- (sum(profitsCond[,i,drop=FALSE]) - sum(debt[i]))/profitsCond[l,i[1],drop=TRUE] - 
        sum(as.vector(sOut) * profitsCond[-l,i,drop=FALSE])/profitsCond[l,i[1],drop=TRUE]
      
      
      
      return(res)
    }
    
   
    
    if(nMarkets>1) theseStates <-  expand.grid(lapply(1:(nMarkets-1),function(x){0:1}))
    
    
    
   
    margin <- prices - mc
    
    sCond <- exp(mval + slopes * prices)
    sOutBE <-exp(z_crit + slopes * outPrice)
    sOutBE <- sOutBE/(sOutBE+rowSums(sOutBE))
    sCond <- sCond/rowSums(sCond)
    profitsCond <- M*margin*sCond
    
    
    
    if(BE){
      sharesCand <- sCond*(1-sOutBE)
      if(sOnly){return(as.vector(t(sharesCand)))}
      
      
      if(!preMerger){
        if(nMarkets>1){sharesCand[,1:2] <- rowSums(sharesCand[,1:2])}
        else{sharesCand[,1:2] <- sum(sharesCand[,1:2])}
      }
      
      marginsCand <- -1/(slopes*(1-sharesCand))
      res <- as.vector(t(margin - marginsCand))
      
      return(res)
      
      
    }
    
    
    
    
    if(nMarkets ==1){hFun <- sapply(firmVec,function(q){s0_focal=s0Focal(NA,q);hFocalIntegral(s0_focal)})  }
    else if(nMarkets>=2 && shareOutParm[1] ==1 && shareOutParm[2] ==1){
      hFun <- sapply(firmVec,function(q,pre=preMerger) {
        s0_focal <- apply(theseStates,1,s0Focal,i=q)
        betaStates <- hFocalIntegral(s0_focal)
        st=ncol(theseStates)
        while(st>=1){
          lastState <- theseStates[,st]
          betaStates <- betaStates[lastState==1]-betaStates[lastState==0]
          theseStates <- unique(theseStates[,-st,drop=FALSE])
          st <- st-1
        }
        #betaStates <- rev((-1)^(1:length(betaStates)))*betaStates
        #betaStates <- sum(betaStates)
        thisProfits <- rowSums(profitsCond[,q,drop=FALSE])
        
        res <- (-1)^(nMarkets+1)*(prod(thisProfits[-1])/(thisProfits[1]^(nMarkets-1)))*betaStates
        return(res)
      })
    }
    
    else{
      hFun <- 
        sapply(firmVec,function(g) {cubintegrate(function(x,f=g){
          x <- as.matrix(x)
          matrix(apply(x, 2, function(s){
            s0_focal=s0Focal(s,f)
            res <- hFocalIntegral(s0_focal)*prod(bDens(s))
            return(res)
          }),ncol=col(x)
          )}
          ,lower=lowerB,upper=upperB,nVec=nVec,method= method_foc,relTol = relTol_foc, absTol = absTol_foc,maxEval = maxEval_foc)$integral
        })
      
    }
    
    if(nMarkets ==1){gFun <- sapply(firmVec,function(q){s0_focal=s0Focal(NA,q);gFocalIntegral(s0_focal)})  }
    else if(nMarkets>=2 && shareOutParm[1] ==1 && shareOutParm[2] ==1){
      gFun <- sapply(firmVec,function(q,pre=preMerger) {
        
        s0_focal <- apply(theseStates,1,s0Focal,i=q)
        betaStates <- gFocalIntegral(s0_focal)
        st=ncol(theseStates)
        while(st>=1){
          lastState <- theseStates[,st]
          betaStates <- betaStates[lastState==1]-betaStates[lastState==0]
          theseStates <- unique(theseStates[,-st,drop=FALSE])
          st <- st-1
        }
        
        if(!pre && q==2){q <- 2:1}
        thisProfits <- rowSums(profitsCond[,q,drop=FALSE])
        res <- (-1)^(nMarkets+1)*(prod(thisProfits[-1])/(thisProfits[1]^(nMarkets-1)))*betaStates
        return(res)
      })
    }
    
    else{
      gFun <- 
        sapply(firmVec,function(g) {cubintegrate(function(x,f=g){
          x <- as.matrix(x)
          matrix(apply(x, 2, function(s){
            s0_focal=s0Focal(s,f)
            
            gFocalIntegral(s0_focal)*prod(bDens(s))}),ncol=col(x)
          )}
          ,lower=lowerB,upper=upperB,nVec=nVec,method= method_foc,relTol = relTol_foc,absTol = absTol_foc, maxEval = maxEval_foc)$integral
        })
    }
    
    if(nMarkets >1){
      if(shareOutParm[1] ==1 && shareOutParm[2] ==1){
        tFun <- sapply(firmVec,function(q,pre=preMerger) {
          
          
          theseStates <- filter(theseStates,Var1==1)
          s0_focal <- apply(theseStates,1,s0Focal,i=q)
          betaStates <- tOtherIntegral(s0_focal)
          st=ncol(theseStates)
          while(st>=2){
            lastState <- theseStates[,st]
            betaStates <- betaStates[lastState==1]-betaStates[lastState==0]
            theseStates <- unique(theseStates[,-st,drop=FALSE])
            st <- st-1
          }
          
          if(!pre && q==2){q <- 2:1}
          thisProfits <- rowSums(profitsCond[,q,drop=FALSE])
          if(nMarkets==2){res <- betaStates}
          else{res <- (-1)^(nMarkets+1)*(prod(thisProfits[-(1:2)])/(thisProfits[1]^(nMarkets-2)))*betaStates}
          
          res <- res +.5 * thisProfits[2]/thisProfits[1]
          return(res)
        })
        
        tFun <- matrix(rep(tFun,nMarkets-1),byrow = TRUE,nrow = nMarkets-1)
        
        rFun <- sapply(firmVec,function(q,pre=preMerger) {
          
          
          theseStates <- filter(theseStates,Var1==1)
          s0_focal <- apply(theseStates,1,s0Focal,i=q)
          betaStates <- rOtherIntegral(s0_focal)
          st=ncol(theseStates)
          while(st>=2){
            lastState <- theseStates[,st]
            betaStates <- betaStates[lastState==1]-betaStates[lastState==0]
            theseStates <- unique(theseStates[,-st,drop=FALSE])
            st <- st-1
          }
          
          if(!pre && q==2){q <- 2:1}
          thisProfits <- rowSums(profitsCond[,q,drop=FALSE])
          if(nMarkets==2){res <- betaStates}
          else{res <- (-1)^(nMarkets+1)*(prod(thisProfits[-(1:2)])/(thisProfits[1]^(nMarkets-2)))*betaStates}
          
          res <- res + (2/3)*thisProfits[2]/thisProfits[1]
          return(res)
        })
        rFun <- matrix(rep(rFun,nMarkets-1),byrow = TRUE,nrow = nMarkets-1)
      }
      
      else{
        tFun <- 
          sapply(firmVec,function(g) {cubintegrate(function(x,f=g){
            x <- as.matrix(x)
            matrix(apply(x, 2, function(s){
              
              integrateDerivOtherFun(s,f)*prod(bDens(s))}),ncol=col(x)
            )}
            ,fDim=1,lower=lowerB,upper=upperB,nVec=nVec,method= method_foc,relTol = relTol_foc, absTol = absTol_foc,maxEval = maxEval_foc)$integral
          })
        
        rFun <- 
          sapply(firmVec,function(g) {cubintegrate(function(x,f=g){
            x <- as.matrix(x)
            matrix(apply(x, 2, function(s){
              
              integrateInsideOtherFun(s,f)*prod(bDens(s))}),ncol=col(x)
            )}
            ,fDim=1,lower=lowerB,upper=upperB,nVec=nVec,method= method_foc,relTol = relTol_foc,absTol = absTol_foc, maxEval = maxEval_foc)$integral
          })
        
      }
      
      #shareOutOther <- exp(log(abs(tFun))-log(abs(rFun)))
      shareOutOther <- tFun/rFun
      shareOutOther[shareOutOther>1] <- 1
      shareOutOther[shareOutOther<0] <- 0
    }
    
    
    
    #replace hFun/gFun with log difference because of stability concerns
    #shareOutFocal <- matrix(exp(log(abs(hFun))-log(abs(gFun))),nrow=1)
    shareOutFocal <- matrix(hFun/gFun,nrow=1)
    shareOutFocal[shareOutFocal > 1] <- 1
    shareOutFocal[shareOutFocal < 0] <- 0
    if(nMarkets==2){shareOutOther <- matrix(shareOutOther,nrow=nMarkets - 1)}
    
    
    if(!preMerger){
      shareOutFocal <- c(shareOutFocal[1],shareOutFocal)
      if(nMarkets>1) shareOutOther <- cbind(shareOutOther[,1],shareOutOther)
    }
    
    SharesOutInt <- matrix(NA,nrow=nMarkets,ncol=nFirms)
    SharesOutInt[focal,]=shareOutFocal
    if(nMarkets>1) SharesOutInt[-focal,]=shareOutOther
    
    
    if(sOnly){
      return(as.vector(t(sCond*(1-SharesOutInt))))
    }
    if(!preMerger){
      if(nMarkets>1)  sCond[,1:2] <- rowSums(sCond[,1:2])
      else{ sCond[,1:2] <- sum(sCond[,1:2])}
    }
    
    margins <- -1/(slopes*(1 - sCond*(1-SharesOutInt) ))
    
    
    
    names(margins) <- object@labels
    
    if(level) {margins <- margins * prices} 
    return(as.vector(margins))
  }
  
)



## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Auction2ndLogit",
  definition=function(object,preMerger=TRUE,exAnte=FALSE,level=TRUE){


    nprods <- length(object@shares)
   
   
    margins <- rep(NA,nprods)

    if( preMerger) {
      owner  <- object@ownerPre
      subset <- rep(TRUE,nprods)
      prices <- object@pricePre
    }
    else{
    subset <- object@subset  
    owner  <- object@ownerPost
    prices <- object@pricePost
    }

    owner <- owner[subset,subset]
    alpha <- object@slopes$alpha
    shares <- calcShares(object,preMerger=preMerger,revenue=FALSE)
    shares <- shares[subset]
    firmShares <- drop(owner %*% shares)
    margins[subset] <-  log(1-firmShares)/(alpha * firmShares)

    if(exAnte){ margins[subset] <-  margins[subset] * shares}
    if(!level){margins <- margins/prices}
    names(margins) <- object@labels

    return(as.vector(margins))
  }

)


#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Auction2ndLogitNests",
  definition=function(object,preMerger=TRUE,exAnte=FALSE,level=FALSE){
    
    
    nprods <- length(object@shares)
    
    
    if( preMerger) {
      owner  <- object@ownerPre
      subset <- rep(TRUE, nprods)
      prices <- object@pricePre}
    else{
      owner  <- object@ownerPost
      subset <- object@subset
      prices <- object@pricePost}
    
    
    nests <- object@nests
    nests <- droplevels(nests[subset])
    nestMat <- tcrossprod(model.matrix(~-1+nests))
    
    margins <- rep(NA,nprods)
    
    
    owner <- owner[subset,subset]
    
    
    alpha <- object@slopes$alpha
    sigma <- object@slopes$sigma
    
    shares <- calcShares(object,preMerger=preMerger,revenue=FALSE)
    
    shares <- shares[subset]
    
    sharesAcross <- as.vector(tapply(shares,nests,sum))[nests]
    sharesIn <- shares/sharesAcross
    
    firmShares <- drop(owner %*% shares)
    
    dupCnt <- rowSums(owner*nestMat) #only include the values in a given nest once
    
    ownerValue <-   1 - (1 - ((owner*nestMat)%*%sharesIn))^sigma[nests]
    ownerValue <-  drop(1 - owner %*%(ownerValue * sharesAcross/dupCnt))
    
    margins[subset] <-  log(ownerValue)/(alpha*firmShares) 
    
    
    if(exAnte){ margins[subset] <-  margins[subset] * shares}
    if(!level){margins <- margins/prices}
    names(margins) <- object@labels
    
    return(as.vector(margins))
  }
  
)
