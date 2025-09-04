#' @title Methods for Calculating Demand Parameters
#' @name Params-Methods
#' @docType methods

#' @aliases calcSlopes
#' calcSlopes,ANY-method
#' calcSlopes,AIDS-method
#' calcSlopes,CES-method
#' calcSlopes,CESNests-method
#' calcSlopes,Linear-method
#' calcSlopes,LogLin-method
#' calcSlopes,Logit-method
#' calcSlopes,LogitALM-method
#' calcSlopes,CESALM-method
#' calcSlopes,LogitCap-method
#' calcSlopes,LogitCapALM-method
#' calcSlopes,LogitNests-method
#' calcSlopes,LogitNestsALM-method
#' calcSlopes,PCAIDS-method
#' calcSlopes,PCAIDSNests-method
#' calcSlopes,Auction2ndLogit-method
#' calcSlopes,Auction2ndLogitNests-method
#' calcSlopes,Auction2ndLogitALM-method
#' calcSlopes,Cournot-method
#' calcSlopes,Stackelberg-method
#' calcSlopes,VertBargBertLogit-method
#' calcSlopes,BargainingLogit-method
#' calcSlopes,Bargaining2ndLogit-method
#' getParms
#' getParms,ANY-method
#' getParms,Bertrand-method
#' getParms,VertBargBertLogit-method
#' getNestsParms
#' getNestsParms,PCAIDSNests-method
#'
#' @description The calcSlopes methods calculate demand parameters assuming that firms are playing
#' a differentitated product Nash-Bertrand pricing game or
#' (as in the case of the Cournot and Stackelberg classes), a Cournot game.
#' @description getNestsParms returns a matrix containing the calibrated nesting parameters.
#' @description getParms returns a list of model-specific demand parameters.
#' \sQuote{digits} specifies the number of significant digit to return (default 10).
#'
#' @param object An instance of the respective class (see description for the classes)
#' @param digits Number of significant digits to report. Default is 2.
#'
#' @include MarginsMethods.R VerticalClasses.R
#' @keywords methods
NULL

setGeneric (
  name= "calcSlopes",
  def=function(object,...){standardGeneric("calcSlopes")}
)
setGeneric (
  name= "getParms",
  def=function(object,...){standardGeneric("getParms")}
)
setGeneric (
  name= "getNestsParms",
  def=function(object,...){standardGeneric("getNestsParms")}
)


## Create a method to recover marginal cost using
## demand parameters and supplied prices
#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "Cournot",
  definition=function(object){

    prices <- object@prices
    quantities <- object@quantities
    quantities[is.na(quantities)] <- 0
    margins <- object@margins
    mktElast <- object@mktElast
    output <- object@output
    outSign <- ifelse(output,-1,1)

    cap <- object@capacitiesPre

    if(output) {mc <- t(t(1 - margins) * prices)}
    else{mc <- t(t(1 + margins) * prices) }

    products <- object@productsPre
    demand <- object@demand
    owner <- object@ownerPre
    mcfunPre <- object@mcfunPre
    nprods <- ncol(quantities)
    nplants <- nrow(quantities)

    noCosts <- length(mcfunPre) == 0
    isLinearD <- demand=="linear"
    isLinearC <- object@cost=="linear"

    quantTot <- colSums(quantities, na.rm = TRUE)
    quantPlants <- rowSums(quantities, na.rm = TRUE)
    quantOwner <- owner %*% quantities

    isConstrained <- quantPlants >= cap

    if(!noCosts){
      mcPre <- sapply(1:nplants, function(i){object@mcfunPre[[i]](quantities[i,])})
    }

    sharesOwner <- t(t(quantOwner)/quantTot)

    minDemand <- function(theta){

      if(noCosts){

        thiscap <- theta[1:nplants]
        theta <- theta[-(1:nplants)]
        mcPre <- ifelse(isLinearC, quantPlants/thiscap, thiscap)

      }

      thisints <- theta[1:nprods]
      thisslopes <- theta[-(1:nprods)]

      thisprices <- ifelse(isLinearD, thisints + thisslopes*quantTot,
                           exp(thisints)*quantTot^thisslopes)

      thisPartial <- ifelse(isLinearD,
                            thisslopes,
                            exp(thisints)*thisslopes*quantTot^(thisslopes - 1))


      thisFOC <- (t(quantities) * thisPartial ) %*% owner  + thisprices
      thisFOC <- t(thisFOC)/mcPre - 1
      thisFOC <- thisFOC[!isConstrained,]

      dist <- c(thisFOC,thisprices/prices -1 , (1/mktElast)/(thisPartial*quantTot/prices) - 1 )

      if(noCosts){ dist <- c(dist, mcPre/mc - 1)}

      return(sum((dist*10)^2,na.rm=TRUE))
    }


    margGuess <- margins
    margGuess[is.na(margGuess)] <- -t(t(sharesOwner)/mktElast)[is.na(margGuess)]

    bStart      =   ifelse(isLinearD,
                           colMeans(outSign*(prices*margGuess)/(sharesOwner*quantTot),na.rm=TRUE),
                           colMeans(outSign*margGuess/sharesOwner,na.rm=TRUE))
    intStart    =   ifelse(isLinearD,
                           prices - bStart*quantTot,
                           log(prices/(quantTot^bStart)))
    intStart    =   -1*outSign*abs(intStart)

    parmStart   =   c( intStart,bStart)

    if(output){
     lowerB <- c(rep(0, nprods), rep(-Inf, nprods))
      upperB <- c(rep(Inf,  nprods), rep(0, nprods))
    }
    else{ lowerB <- c(rep(-Inf,  nprods), rep(0, nprods))
          upperB <- c(rep(0, nprods), rep(Inf, nprods))
    }
    
    if(noCosts){

      margStart     <- matrix(rowMeans(outSign*(sharesOwner*quantTot)/(prices/bStart),na.rm=TRUE), ncol=nprods, nrow=nplants)
      margStart[,!isLinearD]     <-  rowMeans(outSign*sharesOwner*bStart,na.rm=TRUE)
      
      mcStart  <- abs(prices*(margStart + 1))
      capStart <- ifelse(isLinearC, quantPlants/mcStart, mcStart)
      parmStart <- c(capStart,parmStart)

      lowerB <- c(rep(0, nplants), lowerB)
      upperB <- c(rep(Inf,nplants), upperB)
    }



    bestParms=optim(parmStart,minDemand, method="L-BFGS-B", lower=lowerB, upper= upperB)$par

    if(isTRUE(all.equal(bestParms[1:nplants],rep(0, nplants),check.names=FALSE))){warning("Some plant-level cost parameters are close to 0.")}
    if(isTRUE(all.equal(bestParms[-(1:nplants)],rep(0, nprods),check.names=FALSE))){warning("Some demand parameters are close to 0.")}



    ## if no marginal cost functions are supplied
    ## assume that plant i's marginal cost is
    ## q_i/k_i, where k_i is calculated from FOC

    if(noCosts){

      mcparm <- bestParms[1:nplants]
      bestParms <- bestParms[-(1:nplants)]


      mcdef <- ifelse(isLinearC,"function(q,mcparm = %f){ val <- sum(q, na.rm=TRUE) / mcparm; return(val)}",
                      "function(q,mcparm = %f){ val <- mcparm; return(val)}")
      mcdef <- sprintf(mcdef,mcparm)
      mcdef <- lapply(mcdef, function(x){eval(parse(text=x ))})

      object@mcfunPre <- mcdef
      names(object@mcfunPre) <- object@labels[[1]]

      vcdef <- ifelse(isLinearC,"function(q,mcparm = %f){  val <-  sum(q, na.rm=TRUE)^2 / (mcparm * 2); return(val)}",
                      "function(q,mcparm = %f){  val <-  sum(q, na.rm=TRUE) * mcparm; return(val)}")
      vcdef <- sprintf(vcdef,mcparm)
      vcdef <- lapply(vcdef, function(x){eval(parse(text=x ))})

      object@vcfunPre <- vcdef
      names(object@vcfunPre) <- object@labels[[1]]

    }
    if(length(object@mcfunPost)==0){
      object@mcfunPost <- object@mcfunPre
      object@vcfunPost <- object@vcfunPre}

    intercepts = bestParms[1:nprods]
    slopes = bestParms[-(1:nprods)]


    object@intercepts <- intercepts
    object@slopes <-     slopes


    return(object)

  })


#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "Stackelberg",
  definition=function(object){


    prices <- object@prices
    quantities <- object@quantities
    quantities[is.na(quantities)] <- 0
    margins <- object@margins

    mc <- t(t(1 - margins) * prices)

    products <- object@productsPre
    demand <- object@demand
    owner <- object@ownerPre
    vcfunPre <- object@vcfunPre
    nprods <- ncol(quantities)
    nplants <- nrow(quantities)

    noCosts <- length(vcfunPre) == 0
    isLeader <- object@isLeaderPre
    isLinearD <- demand=="linear"
    isLinearC <- object@cost=="linear"

    quantTot <- colSums(quantities, na.rm = TRUE)
    quantPlants <- rowSums(quantities, na.rm = TRUE)
    quantOwner <- owner %*% quantities



    sharesOwner <- t(t(quantOwner)/quantTot)

    ## if variable cost function is missing
    ## assume vc_i = q_i^2/k_i for plant i
    ## adding create functions for mc and dmc
    if(!noCosts){
      mcPre <- sapply(1:nplants, function(i){object@mcfunPre[[i]](quantities[i,])})
      dmcPre <- sapply(1:nplants, function(i){object@dmcfunPre[[i]](quantities[i,])})
    }

    minParms <- function(theta){

      if(noCosts){

        thiscap <- theta[1:nplants]
        theta <- theta[-(1:nplants)]
        mcPre <- ifelse(isLinearC, quantPlants/thiscap, thiscap)
        dmcPre <- ifelse(isLinearC, 1/thiscap, 0)
      }

      thisints <- theta[1:nprods]
      thisslopes <- theta[-(1:nprods)]

      thisprices <- ifelse(isLinearD, thisints + thisslopes*quantTot,
                           exp(thisints)*quantTot^thisslopes)

      thisPartial <- ifelse(isLinearD,
                            thisslopes,
                            exp(thisints)*thisslopes*quantTot^(thisslopes - 1))

      dthisPartial <- ifelse(isLinearD,
                             0,
                             exp(thisints)*thisslopes*(thisslopes - 1)*quantTot^(thisslopes - 2))


      demPass <- dthisPartial * t(!isLeader * quantOwner)
      thisPass <- -t((thisPartial + demPass)/
                       (2*thisPartial  + t(t(demPass) - dmcPre)))

      thisPass[isLeader] <- 0

      thisFOC <- (t(quantities) * thisPartial + t(isLeader * quantities) * thisPartial * colSums(thisPass)) %*% owner  + thisprices
      thisFOC <- t(thisFOC) -   mcPre


      dist <- c(thisFOC,thisprices - prices)

      if(noCosts){ dist <- c(dist, mcPre - mc)}

      return(sum(dist^2,na.rm=TRUE))
    }


    bStart      =   ifelse(isLinearD,
                           colMeans(-(prices*margins)/(sharesOwner*quantTot),na.rm=TRUE),
                           colMeans(-margins/sharesOwner,na.rm=TRUE))
    intStart    =   ifelse(isLinearD,
                           prices - bStart*quantTot,
                           log(prices/(quantTot^bStart)))
    intStart    =   abs(intStart)

    parmStart   =   c( intStart,bStart)

    if(noCosts){

      if(isLinearD){margStart <- rowMeans(-(sharesOwner*quantTot)/(prices/bStart),na.rm=TRUE)}
      else{margStart <- rowMeans(-sharesOwner*bStart,na.rm=TRUE)}

      mcStart  <- abs(prices*(margStart - 1))
      capStart <- ifelse(isLinearC, quantPlants/mcStart, mcStart)
      parmStart <- c(capStart,parmStart)
    }

    bestParms=optim(parmStart,minParms)$par

    ## if no variable cost functions are supplied
    ## assume that plant i's varialbe cost is
    ## q_i^2/(2*k_i), where k_i is calculated from FOC

    if(noCosts){

      mcparm <- bestParms[1:nplants]
      bestParms <- bestParms[-(1:nplants)]


      dmcdef <- ifelse(isLinearC,"function(q,mcparm = %f){ val <-   1/mcparm; return(val)}",
                       "function(q,mcparm = %f){ val <-   0; return(val)}")
      dmcdef <- sprintf(dmcdef,mcparm)
      dmcdef <- lapply(dmcdef, function(x){eval(parse(text=x ))})

      object@dmcfunPre <- dmcdef
      names(object@dmcfunPre) <- object@labels[[1]]

      mcdef <- ifelse(isLinearC,"function(q,mcparm = %f){ val <- sum(q, na.rm=TRUE) / mcparm; return(val)}",
                      "function(q,mcparm = %f){ val <- mcparm; return(val)}")
      mcdef <- sprintf(mcdef,mcparm)
      mcdef <- lapply(mcdef, function(x){eval(parse(text=x ))})

      object@mcfunPre <- mcdef
      names(object@mcfunPre) <- object@labels[[1]]

      vcdef <- ifelse(isLinearC,"function(q,mcparm = %f){  val <-  sum(q, na.rm=TRUE)^2 / (mcparm * 2); return(val)}",
                      "function(q,mcparm = %f){  val <-  sum(q, na.rm=TRUE) * mcparm; return(val)}")
      vcdef <- sprintf(vcdef,mcparm)
      vcdef <- lapply(vcdef, function(x){eval(parse(text=x ))})

      object@vcfunPre <- vcdef
      names(object@vcfunPre) <- object@labels[[1]]


    }


    intercepts = bestParms[1:nprods]
    slopes = bestParms[-(1:nprods)]





    if(length(object@vcfunPost)==0){
      object@dmcfunPost <- object@dmcfunPre
      object@mcfunPost <- object@mcfunPre
      object@vcfunPost <- object@vcfunPre}

    object@intercepts <- intercepts
    object@slopes <-     slopes


    return(object)

  })

#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "Linear",
  definition=function(object){

    margins    <- object@margins
    quantities <- object@quantities
    prices     <- object@prices
    diversion  <- object@diversion
    ownerPre   <- object@ownerPre
    symmetry  <- object@symmetry

    nprod <- length(margins)



    if(!symmetry){


      slopes <- matrix(margins * prices,ncol=nprod, nrow=nprod,byrow=TRUE)
      slopes <- diag(ownerPre)/rowSums(slopes * diversion * ownerPre) * quantities
      slopes <- -t(slopes * diversion)


    }

    else{

      existMargins <- which(!is.na(margins))

      revenues <- prices*quantities
      k <- existMargins[1] ## choose a diagonal demand parameter corresponding to a provided margin


      minD <- function(betas){

        #enforce symmetry

        B=diag(betas[1:nprod])
        B[upper.tri(B,diag=FALSE)] <- betas[-(1:nprod)]
        B=t(B)
        B[upper.tri(B,diag=FALSE)] <- betas[-(1:nprod)]

        elast <- t(B * tcrossprod(1/quantities,prices))

        marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% (revenues * diag(ownerPre))) / revenues


        m1 <- margins - marginsCand
        m2 <- as.vector(diversion +  t(B)/diag(B)) #measure distance between observed and predicted diversion


        measure=c(m1,m2)

        return(sum(measure^2,na.rm=TRUE))
      }


      ## Create starting values for optimizer

      bKnown      =  -quantities[k]/(prices[k]*margins[k])
      bStart      =   bKnown*diversion[k,]/diversion[,k]

      ## change starting guess to ensure that it satisfies constraints
      mltplyr <- 1.01 # increase starting guess by 1%
      isneg <- as.vector(t(diversion) %*% bStart < 0)

      while(any(isneg)){

        bStart[!isneg] <- bStart[!isneg] * mltplyr
        mltplyr <- mltplyr + .01 # decrement by 1%

        isneg <- as.vector(t(diversion) %*% bStart < 0)

        if(any(is.na(isneg))){
          stop("'calcSlopes' cannot find initial values that satisfy symmetry constraints using supplied data. Consider setting 'symmetry' equal to FALSE."
          )}
      }

      bStart      =  -diversion*bStart
      parmStart   =   c(diag(bStart),bStart[upper.tri(bStart,diag=FALSE)])



      ## constrain diagonal elements so that D'b >=0
      ## constrain off-diagonal elements to be non-negative.

      ui          =  diag(length(parmStart))
      ui[1:nprod,1:nprod] = t(diversion)

      ci = rep(0,length(parmStart))



      bestParms=constrOptim(parmStart,minD,grad=NULL,ui=ui,ci=ci,
                            control=object@control.slopes)

      slopes = diag(bestParms$par[1:nprod])

      slopes[upper.tri(slopes,diag=FALSE)] <- bestParms$par[-(1:nprod)]
      slopes=t(slopes)
      slopes[upper.tri(slopes,diag=FALSE)] <- bestParms$par[-(1:nprod)]


    }


    dimnames(slopes) <- list(object@labels,object@labels)


    intercept <- as.vector(quantities - slopes %*% prices)
    names(intercept) <- object@labels

    if(!symmetry &&
       !isTRUE(all.equal(slopes,t(slopes)))){
      warning("Matrix of demand slopes coefficients is not symmetric. Demand parameters may not be consistent with utility maximization theory.")}

    if(any(intercept<0))   {warning(  "Some demand intercepts are negative")}
    if(any(diag(slopes)>0)){warning(  "Some own-slope coefficients are positive")}

    object@slopes <- slopes
    object@intercepts <- intercept

    return(object)


  }
)


#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "Logit",
  definition=function(object){

    ## Uncover Demand Coefficents


    ownerPre     <-  object@ownerPre
    shares       <-  object@shares
    margins      <-  object@margins
    prices       <-  object@prices
    idx          <-  object@normIndex
    mktElast     <-  object@mktElast
    shareInside  <-  object@shareInside
    diversion    <-  object@diversion
    output       <-  object@output 
    outSign <- ifelse(output,-1,1)
   
    margins <- margins*prices
    
    if(is.na(idx)){
      idxShare <- 1 - shareInside
      idxPrice <- object@priceOutside
    }
    else{
      idxShare <- shares[idx]
      idxPrice <- prices[idx]
    }

    ## Choose starting parameter values
    notMissing <- which(!is.na(margins))[1]

    parmStart <- outSign/(margins[notMissing]*(1 - shares[notMissing]))
    mvalStart <-  log(shares) - log(idxShare) - parmStart * (prices - idxPrice)
    if(!is.na(idx)) mvalStart <-  mvalStart[-idx]
    parmStart <- c(parmStart, mvalStart)

    ## Uncover price coefficient and mean valuation from margins and revenue shares


    nprods <- length(shares)

    avgPrice <- sum(shares * prices, na.rm=TRUE) / sum(shares)

    ## identify which products have enough margin information
    ##  to impute Bertrand margins
    #isMargin    <- matrix(margins,nrow=nprods,ncol=nprods,byrow=TRUE)
    #isMargin[ownerPre==0]=0
    #isMargin    <- !is.na(rowSums(isMargin))


    ## Minimize the distance between observed and predicted margins
    minD <- function(theta){
      
      alpha <- theta[1]
      
      if(!is.na(idx)){
        meanval <- rep(0,nprods)
        meanval[-idx] <- theta[-1]
      }
      else{meanval <- theta[-1]}

      probs <- shares
      
      predshares <- exp(meanval + alpha*(prices-idxPrice))
      predshares <- predshares/(is.na(idx) + sum(predshares) )

      preddiversion <-tcrossprod( 1/(1-predshares),predshares)
      diag(preddiversion) <- -1
      
      if(!is.na(mktElast)){
        shareInside <-   1 - mktElast/( alpha * avgPrice )
        probs <- probs/sum(probs,na.rm=TRUE) * shareInside

      }

      
    
     marginsCand <- outSign/(alpha*(1-as.numeric(crossprod(ownerPre,predshares))))

      m1 <- (margins - marginsCand)
      m2 <- log(predshares/probs)
      m3 <- drop(diversion - preddiversion)
      measure <- sum((c(m1, m2, m3))^2,na.rm=TRUE)

      return(measure)
    }
# 
#     alphaBounds <- c(-1e6,0)
#     if(!is.na(mktElast)){ alphaBounds[2] <- mktElast/ avgPrice}
# 
#     minAlpha <- optimize(minD, alphaBounds,
#                          tol=object@control.slopes$reltol)$minimum
# 
#     if(!is.na(mktElast)){
# 
# 
#       object@shareInside <-    1 - mktElast/(minAlpha * avgPrice )
#       idxShare <-  1 - object@shareInside
# 
#     }

    
    
    ##  Constrained optimizer to look for solutions where alpha<0,  
    lowerB <- upperB <- c(1e6,rep(12,length(parmStart)-1))
    lowerB <- lowerB * -1
    
    if(output){
      upperB[1] <- -1e-5
    }
    else{
      lowerB[1] <- 1e-5
    }
    
    minTheta <- optim(parmStart,minD,method="L-BFGS-B",
                      lower= lowerB,upper=upperB,
                      control=object@control.slopes)
    
    
    if(minTheta$convergence != 0){
      warning("'calcSlopes' nonlinear solver may not have successfully converge. Reason: '",minTheta$message,"'")
    }
    
    
    # ui=diag(length(parmStart))
    # #ui[1,1] <- -1
    # if(!is.na(idx)){ui[-1,1] <- prices[-idx] - idxPrice}
    # else{ui[-1,1] <- prices - idxPrice}
    # ci_hi=rep(log(.9999/(1-.9999)), length(mvalStart)) 
    # ci_low=rep(log((1-.9999)/.9999), length(mvalStart)) 
    # 
    # ui=rbind(-ui,ui[-1,])
    # ci=c(0,-ci_hi,ci_low)
    # 
    # minTheta <- constrOptim(parmStart,minD,grad=NULL,ui=ui,ci=ci)
    
    minAlpha           <- minTheta$par[1]
    names(minAlpha)    <- "alpha"
    
    if(is.na(idx)) meanval <-  minTheta$par[-1]
    else{
      meanval <- rep(0,nprods)
      meanval[-idx] <- minTheta$par[-1]
    }
    
    names(meanval)   <- object@labels



    object@slopes    <- list(alpha=minAlpha,meanval=meanval)
    object@priceOutside <- idxPrice
    object@mktSize <- object@insideSize / sum(shares)


    return(object)
  }
)


#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "LogitCournot",
  definition=function(object){
    
    ## Uncover Demand Coefficents
    
    
    ownerPre     <-  object@ownerPre
    shares       <-  object@shares
    margins      <-  object@margins
    
    prices       <-  object@prices
    idx          <-  object@normIndex
    mktElast     <-  object@mktElast
    shareInside  <-  object@shareInside
    diversion    <-  object@diversion
    output       <-  object@output
    outSign <- ifelse(output,-1,1)
    margins      <- margins*prices 
    
    if(is.na(idx)){
      idxShare <- 1 - shareInside
      idxPrice <- object@priceOutside
    }
    else{
      idxShare <- shares[idx]
      idxPrice <- prices[idx]
    }
    
    ## Choose starting parameter values
    notMissing <- which(!is.na(margins))[1]
    
    parmStart <- outSign/(margins[notMissing])*(1 + shares[notMissing]/idxShare)
    mvalStart <-  log(shares) - log(idxShare) - parmStart * (prices - idxPrice)
    if(!is.na(idx)) mvalStart <-  mvalStart[-idx]
    parmStart <- c(parmStart, mvalStart)
    
    ## Uncover price coefficient and mean valuation from margins and revenue shares
    
    
    nprods <- length(shares)
    
    avgPrice <- sum(shares * prices, na.rm=TRUE) / sum(shares)
    
    
    ## Minimize the distance between observed and predicted margins
    minD <- function(theta){
      
      alpha <- theta[1]
      
      if(!is.na(idx)){
        meanval <- rep(0,nprods)
        meanval[-idx] <- theta[-1]
      }
      else{meanval <- theta[-1]}
      
      probs <- shares
      
      predshares <- exp(meanval + alpha*(prices-idxPrice))
      predshares <- predshares/(is.na(idx) + sum(predshares) )
      
      preddiversion <-tcrossprod( 1/(1-predshares),predshares)
      diag(preddiversion) <- -1
      
      if(!is.na(mktElast)){
        shareInside <-   1 - mktElast/( alpha * avgPrice )
        probs <- probs/sum(probs,na.rm=TRUE) * shareInside
        
      }
      
      predsharesFirm <- as.numeric(ownerPre %*% predshares)
      
      if(is.na(idx)){idxPredShare <- 1-sum(predshares)}
      else{idxPredShare <- predshares[idx]}
      marginsCand <- (1+predsharesFirm/idxPredShare)/alpha
      
      marginsCand <- outSign*(marginsCand)
      
      m1 <- (margins - marginsCand)
      m2 <- log(predshares/probs)
      m3 <- drop(log(diversion/ preddiversion))
      measure <- sum((c(m1, m2, m3))^2,na.rm=TRUE)
      
      return(measure)
    }
    
    
    
    ##  Constrained optimizer to look for solutions where alpha<0,  
    lowerB <- upperB <- c(1e6,rep(12,length(parmStart)-1))
    lowerB <- lowerB * -1
    
    if(output){
      upperB[1] <- -1e-5
    }
    else{
      lowerB[1] <- 1e-5
    }
    
    minTheta <- optim(parmStart,minD,method="L-BFGS-B",
                      lower= lowerB,upper=upperB,
                      control=object@control.slopes)
    
    
    if(minTheta$convergence != 0){
      warning("'calcSlopes' nonlinear solver may not have successfully converge. Reason: '",minTheta$message,"'")
    }
    
    
    minAlpha           <- minTheta$par[1]
    names(minAlpha)    <- "alpha"
    
    if(is.na(idx)) meanval <-  minTheta$par[-1]
    else{
      meanval <- rep(0,nprods)
      meanval[-idx] <- minTheta$par[-1]
    }
    
    names(meanval)   <- object@labels
    
    
    
    object@slopes    <- list(alpha=minAlpha,meanval=meanval)
    object@priceOutside <- idxPrice
    object@mktSize <- object@insideSize / sum(shares)
    
    
    return(object)
  }
)


#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "LogLin",
  definition=function(object){



    margins <- object@margins
    quantities <- object@quantities
    prices <- object@prices
    ownerPre <- object@ownerPre

    revenues <- prices * quantities

    nprods <- length(margins)

    diversion <- object@diversion * tcrossprod(quantities,1/quantities)

    slopes <- matrix(margins * revenues,ncol=nprods, nrow=nprods,byrow=TRUE)

    slopes <- (revenues * diag(ownerPre)) / rowSums(slopes * diversion * ownerPre)
    slopes <- -t(slopes * diversion)

    dimnames(slopes) <- list(object@labels,object@labels)

    intercept <- as.vector(log(quantities) - slopes %*% log(prices))
    names(intercept) <-  object@labels

    object@slopes <- slopes
    object@intercepts <- intercept

    return(object)


  }
)


#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "AIDS",
  definition=function(object){


    ## Uncover linear demand slopes
    shares     <- object@shares
    prices     <- object@prices
    margins    <- object@margins
    diversion  <- object@diversion
    labels     <- object@labels
    ownerPre   <- object@ownerPre
    parmStart  <- object@parmStart
    mktElast   <- ifelse(length(object@mktElast) == 0, NA, object@mktElast)
    nprod=length(shares)

    cancalibrate <- apply(diversion,1,function(x){!any(x==0)})
    idx <- which.max(ifelse(cancalibrate,shares, -1))


    if(any(is.na(parmStart))){
      parmStart[2] <-  -1.2
      parmStart[1] <-  (1 - shares[idx] + parmStart[2] * (1 - shares[idx])) * shares[idx] - .1
      if(parmStart[1] >= 0){parmStart[1] <- -.5}
    }
    minD <- function(s){

      #enforce symmetry
      thismktElast = s[2]
      betas  =   s[1]

      betas <- -diversion[idx,]/diversion[,idx] * betas

      B = t(diversion * betas)
      #diag(B)=  betas - rowSums(B) #enforce homogeneity of degree zero

      elast <- t(B/shares) + shares * (thismktElast + 1) #Caution: returns TRANSPOSED elasticity matrix
      diag(elast) <- diag(elast) - 1

      marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% (shares * diag(ownerPre))) / shares


      m1 <- margins - marginsCand
      m2 <- diversion/t(diversion) - tcrossprod(1/betas,betas)
      m2 <-  m2[upper.tri(m2)]
      m2 <- m2[is.finite(m2) & m2 != 0]
      m3 <- thismktElast - mktElast
      #m2 <- as.vector(diversion +  t(B)/diag(B)) #measure distance between observed and predicted diversion


      measure=c(m1,m2,m3)


      return(sum(measure^2,na.rm=TRUE))
    }



    ui = -diag(2)
    ui[2,1] = -1/shares[idx]
    ui[2,2] =  1 - shares[idx]
    ui <- rbind(ui,c(0,-1))
    ci = rep(0,3)
    ci[2] = -(1 - shares[idx])

    bestParms=constrOptim(parmStart,f=minD, ui=ui,ci=ci, grad=NULL,
                          control=object@control.slopes)


    betas <- -diversion[idx,]/diversion[,idx]*bestParms$par[1]

    B = t(diversion * betas)

    dimnames(B) <- list(object@labels,object@labels)

    object@slopes <- B

    if(abs(bestParms$par[2])>5){warning("'mktElast' estimate is large: ",bestParms$par[2])}
    if(isTRUE(all.equal(bestParms$par[2] , object@parmStart[2]))){warning("'mktElast' estimate is not identified" )}
    object@mktElast <- bestParms$par[2]
    object@intercepts <- as.vector(shares - B%*%log(prices))
    names(object@intercepts) <- object@labels


    return(object)

  }

)

#'@rdname Params-Methods
#'@export
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


#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "PCAIDSNests",
  definition=function(object){


    ## Uncover linear demand slopes from shares, knownElast, mktElast, margins, and nesting structure



    knownIndx  <- object@knownElastIndex
    knownElast <- object@knownElast
    mktElast   <- object@mktElast
    ownerPre   <- object@ownerPre
    shares <-  object@shares
    margins <- object@margins
    nests <- object@nests


    shareKnown <- shares[knownIndx]
    nprods <- length(shares)


    nNests <- nlevels(nests)
    nests <- as.numeric(nests)

    bKnown <- shareKnown * (knownElast + 1 - shareKnown * (mktElast + 1))


    calcB <- function(n){
      nestWeights <- diag(nNests)
      nestWeights[upper.tri(nestWeights)] <- nestWeights[lower.tri(nestWeights)] <- n

      nestWeights <- nestWeights[nests,nests]

      sumWeights <- sum(shares * nestWeights[,knownIndx], na.rm=TRUE) - shareKnown

      beta <- shares/shareKnown
      beta <- beta * (rowSums(t(shares * t(nestWeights)),na.rm=TRUE) - shares)
      beta <- beta / sumWeights

      b    <- beta * bKnown


      B <- -bKnown * (tcrossprod(shares) * nestWeights)/(shareKnown*sumWeights)


      diag(B) <- b

      return(B)

    }

    minD <- function(theseNests){

      Bcand <- calcB(theseNests)

      elast <- t(Bcand/shares) + shares * (mktElast + 1)
      diag(elast) <- diag(elast) - 1

      marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% (shares * diag(ownerPre))) / shares

      measure <- sum((margins - marginsCand)^2,na.rm=TRUE)

      return(measure)
    }

    minNests <- optim(object@nestsParms,minD,method ="L-BFGS-B",
                      lower=0,upper=1,
                      control=object@control.slopes)$par

    B <- calcB(minNests)

    object@nestsParms <- minNests


    dimnames(B) <- list(object@labels,object@labels)
    object@slopes <- B
    object@intercepts <- as.vector(shares - B%*%log(object@prices))
    names(object@intercepts) <- object@labels




    return(object)
  }
)


#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "LogitCap",
  definition=function(object){

    ## Uncover Demand Coefficents


    ownerPre     <-  object@ownerPre
    shares       <-  object@shares
    margins      <-  object@margins
    prices       <-  object@prices
    capacities  <-  object@capacitiesPre/object@insideSize
    idx          <-  object@normIndex
    output      <- object@output
    
  
    
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
    revenues <- shares * prices

    ##create a matrix of 1s and 0s where the i,jth element equals 1 if product i is NOT producing at capacity
    notBinds <- matrix(as.numeric(capacities > shares),ncol=nprods,nrow=nprods,byrow=TRUE)
    ## create a TRUE/FALSE vector equal to TRUE if a single product firm is capacity constrained
    singleConstrained <- rowSums( object@ownerPre>0) == 1 & capacities == shares

    ## Minimize the distance between observed and predicted margins
    minD <- function(alpha){

      ## the following returns the elasticity TRANSPOSED
      elast <- -alpha *  matrix(prices * shares,ncol=nprods,nrow=nprods)
      diag(elast) <- alpha*prices + diag(elast)


      FOC <- revenues * diag(ownerPre) + (elast * ownerPre * notBinds) %*% (margins * revenues)

      ## omit the FOCs of single product, capacity constrained firms
      measure <- sum(as.vector(FOC[!singleConstrained])^2,na.rm=TRUE)

      return(measure)
    }

    
    if(output){bounds <- c(-1e6,0)}
    else{bounds <- c(0,1e6)}
    
    minAlpha <- optimize(minD,bounds,
                         tol=object@control.slopes$reltol)$minimum


    meanval <- log(shares) - log(idxShare) - minAlpha * (prices - idxPrice)

    names(meanval)   <- object@labels

    object@slopes    <- list(alpha=minAlpha,meanval=meanval)
    object@mktSize <- object@insideSize / sum(shares)



    return(object)
  }
)


#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "LogitNests",
  definition=function(object){

    ## Uncover Demand Coefficents


    ownerPre     <-  object@ownerPre
    shares       <-  object@shares
    nprods      <-   length(shares)

    diversion <- object@diversion
    margins      <-  object@margins

    
    prices       <-  object@prices
    idx          <-  object@normIndex

    parmsStart   <- object@parmsStart
    nests        <- object@nests
    nestCnt      <- tapply(prices,nests,length)
    constraint   <- object@constraint
    output       <- object@output
    
    isSingletonNest <- nestCnt==1

    
    sharesNests <- tapply(shares,nests,sum)[nests]
    sharesNests <- shares / sharesNests
    
    if(any(isSingletonNest)){
      warning("Some nests contain only one product; their nesting parameters are not identified.
              Normalizing these parameters to 1.")

    }


    ## create index variables, contingent on whether an outside good is defined
    isnotIdx <- rep(TRUE,nprods)
    
    if(is.na(idx)){
      idxShare <- 1 - object@shareInside
      idxShareIn <- 1
      idxPrice   <- object@priceOutside
      idxSigma   <- 1
      
     
      
    }
    
    else{
      
      
      idxShare   <- shares[idx]
      idxShareIn <- sharesNests[idx]
      idxPrice   <- prices[idx]
      
      isnotIdx[idx] <- FALSE
     
      
    }
    
   
    if(!constraint){
      parmsStart   <- parmsStart[c(TRUE,!isSingletonNest)] #always retain first element; this is
      # the initial value for price coefficient
    }
    ## Uncover price coefficient and mean valuation from margins and revenue shares


    ## Choose starting parameter values according to flat logit
    notMissing <- which(!is.na(margins))[1]
    
    mvalStart <- log(shares) - log(idxShare) - parmsStart[1] * (prices - idxPrice)
   mvalStart <- mvalStart[isnotIdx]
   
    parmsStart <- c(parmsStart[1], mvalStart ,parmsStart[-1] )
    



    





    ## Minimize the distance between observed and predicted margins
    minD <- function(theta){

      alpha <- theta[1]
      theta <- theta[-1]
      meanval <- rep(0,nprods)
      meanvalIdx <- 1:(nprods - !is.na(idx))
      meanval[isnotIdx] <- theta[meanvalIdx]
      theta <-   theta[-meanvalIdx]
      sigma <- as.numeric(isSingletonNest)
      sigma[!isSingletonNest] <- theta

      
      outVal    <- ifelse(is.na(idx), exp(alpha*object@priceOutside), 0)
      
      predSharesIn     <- exp((meanval+alpha*prices)/sigma[nests])
      
      inclusiveValue <- log(tapply(predSharesIn,nests,sum,na.rm=TRUE))
      
      predSharesAcross <-   exp(sigma*inclusiveValue)
      predSharesAcross <- predSharesAcross/(outVal + sum(predSharesAcross,na.rm=TRUE))
      
      
      predSharesIn     <- predSharesIn/exp(inclusiveValue)[nests]
      
      predSharesAcross <- predSharesAcross[nests]
      predshares       <- predSharesIn * predSharesAcross
      
      revenues <- prices * predshares
      
      ## The following returns the transpose of the elasticity matrix
      elasticity <- diag((1/sigma-1)*alpha)
      elasticity <- elasticity[nests,nests]
      elasticity <- elasticity * matrix( predSharesAcross*prices,ncol=nprods,nrow=nprods)
      elasticity <- -1*(elasticity + alpha * matrix(predshares*prices,ncol=nprods,nrow=nprods))
      diag(elasticity) <- diag(elasticity) + (1/sigma[nests])*alpha*prices

      preddiversion <- -elasticity/diag(elasticity)*tcrossprod(1/predshares,predshares)
      diag(preddiversion) <- -1

      elastInv <- try(solve(elasticity * ownerPre),silent=TRUE)
      if(any(class(elastInv)=="try-error")){elastInv <- MASS::ginv(elasticity * ownerPre)}
      
      
      marginsCand <- -1 * as.vector(elastInv %*% (revenues * diag(ownerPre))) / revenues
      m1 <- margins - marginsCand
      m2 <- predshares - shares
      m3 <- drop(diversion - preddiversion)
      measure <- sum((c(m1, m2, m3)*100)^2,na.rm=TRUE)
      

      return(measure)
    }

    ##  Constrained optimizer to look for solutions where alpha<0,  1 > sigma > 0.
    ##  sigma > 1 or sigma < 0 imply complements
    lowerB <- upperB <- rep(-Inf,length(parmsStart))
    upperB <- upperB*-1
    
    if(!output){
      lowerB[1] <- 0
    }
    else{
      upperB[1] <- 0
    }
    
    upperB[-(1:((nprods - !is.na(idx))+1))] <- 1
    lowerB[-(1:((nprods - !is.na(idx))+1))] <- 1e-3
    
    minTheta <- optim(parmsStart,minD,method="L-BFGS-B",
                      lower= lowerB,upper=upperB,
                      control=object@control.slopes)


    if(minTheta$convergence != 0){
      warning("'calcSlopes' nonlinear solver did not successfully converge. Reason: '",minTheta$message,"'")
    }
    
    



    minAlpha           <- minTheta$par[1]
    names(minAlpha)    <- "alpha"


    minSigma <-  as.numeric(isSingletonNest)
    meanval <- rep(0,nprods)
  
    meanval[isnotIdx] <-  minTheta$par[2:((nprods - !is.na(idx)) + 1)]
    minSigma[!isSingletonNest] <- minTheta$par[-(1:((nprods - !is.na(idx)) + 1))]


    minSigmaOut        <- minSigma
    minSigma           <- minSigma[nests]
    names(minSigmaOut) <- levels(nests)

    if(any(minSigmaOut %in%(c(1e-3,1)))) warning("nesting parameter 'sigma' at boundary constraint.")



    names(meanval)   <- object@labels

    object@slopes    <- list(alpha=minAlpha,sigma=minSigmaOut,meanval=meanval)

    return(object)
  }
)



#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "Auction2ndLogitNests",
  definition=function(object){
    
    ## Uncover Demand Coefficents
    
    
    ownerPre     <-  object@ownerPre
    shares       <-  object@shares
    nprods      <-   length(shares)
    
    diversion <- object@diversion
    margins      <-  object@margins
    output    <-  object@output
    
    prices <- object@prices
    
    idx          <-  object@normIndex
    
    parmsStart   <- object@parmsStart
    nests        <- object@nests
    nestCnt      <- tapply(shares,nests,length)
    constraint   <- object@constraint
    
    nestMat <- tcrossprod(model.matrix(~-1+nests)) 
    
    
    dupCnt <- rowSums(ownerPre*nestMat) #only include the values in a given nest once
    
    isSingletonNest <- nestCnt==1
    
    
    sharesBetween <- sharesNests <- as.vector(tapply(shares,nests,sum)[nests])
    sharesNests <- shares / sharesNests
    
    ownerShares <- drop(ownerPre %*% shares)
    
    if(any(isSingletonNest)){
      warning("Some nests contain only one product; their nesting parameters are not identified.
              Normalizing these parameters to 1.")
      
    }
    
    
    ## create index variables, contingent on whether an outside good is defined
    isnotIdx <- rep(TRUE,nprods)
    
    if(is.na(idx)){
      idxShare <- 1 - object@shareInside
      idxShareIn <- 1
      idxPrice   <- object@priceOutside
      idxSigma   <- 1
      
      
      
    }
    
    else{
      
      
      idxShare   <- shares[idx]
      idxShareIn <- sharesNests[idx]
      idxPrice   <- prices[idx]
      
      isnotIdx[idx] <- FALSE
      
      
    }
    
    
    if(!constraint){
      parmsStart   <- parmsStart[c(TRUE,!isSingletonNest)] #always retain first element; this is
      # the initial value for price coefficient
    }
    
    
    
    
    ## Minimize the distance between observed and predicted margins,
    ## diversions, shares
    minD <- function(theta){
      
      alpha <- theta[1]
      theta <- theta[-1]
      sigma <- as.numeric(isSingletonNest)
      sigma[!isSingletonNest] <- theta
      
      
      ownerValue <-   1 - (1 - ((ownerPre*nestMat)%*%sharesNests))^sigma[nests]
      ownerValue <-  drop(1 - ownerPre %*%(ownerValue * sharesBetween/dupCnt))
      
      marginsCand <-  log(ownerValue)/(alpha*ownerShares) 
   
      ## The following returns the transpose of the elasticity matrix
      elasticity <- diag((1/sigma-1)*alpha)
      elasticity <- elasticity[nests,nests]
      elasticity <- elasticity * matrix( shares*prices,ncol=nprods,nrow=nprods)
      elasticity <- -1*(elasticity + alpha * matrix(shares*prices,ncol=nprods,nrow=nprods))
      diag(elasticity) <- diag(elasticity) + (1/sigma[nests])*alpha*prices
      
      preddiversion <- -elasticity/diag(elasticity)*tcrossprod(1/shares,shares)
      diag(preddiversion) <- -1
      
      
      m1 <- marginsCand/margins-1
      m3 <- drop(diversion - preddiversion)
      measure <- sum((c(m1,  m3)*100)^2,na.rm=TRUE)
      
      
      return(measure)
    }
    
    ##  Constrained optimizer to look for solutions where alpha<0,  1 > sigma > 0.
    ##  sigma > 1 or sigma < 0 imply complements
    
    
    if(output){
      lowerB <- c(-Inf,0)
      upperB <- c(0,1)
    }
    
    else{
      lowerB <- c(0,0)
      upperB <- c(Inf,1)
    }
    
    minTheta <- optim(parmsStart,minD,method="L-BFGS-B",
                      lower= lowerB,upper=upperB,
                      control=object@control.slopes)
    
    
    if(minTheta$convergence != 0){
      warning("'calcSlopes' nonlinear solver did not successfully converge. Reason: '",minTheta$message,"'")
    }
    
    
    
    minAlpha           <- minTheta$par[1]
    names(minAlpha)    <- "alpha"
    
    meanval <- log(shares) - log(idxShare)
    
    minSigma <-  as.numeric(isSingletonNest)
    
    minSigma[!isSingletonNest] <- minTheta$par[-1]
    
    
    names(minSigma) <- levels(nests)
    
    if(is.na(idx)) minSigmaOut <- 0
    else {minSigmaOut        <- minSigma[idx]}
    
  
    meanvalDown=minSigma*(log(shares)-log(idxShare))
    
    names(meanval)   <- object@labels
    
    object@slopes    <- list(alpha=minAlpha,sigma=minSigma,meanval=meanval)
    
    return(object)
  }
)



#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "LogitCapALM",
  definition=function(object){

    ## Uncover Demand Coefficents


    ownerPre     <-  object@ownerPre
    shares       <-  object@shares
    margins      <-  object@margins
    prices       <-  object@prices
    insideSize   <-  object@insideSize
    capacities  <-  object@capacitiesPre/insideSize
    mktElast     <-  object@mktElast
    priceOutside <- object@priceOutside

    avgPrice <- sum(shares*prices)

    nprods <- length(object@shares)


    ## Uncover price coefficient and mean valuation from margins and revenue shares



    revenues <- shares * prices

    ##create a matrix of 1s and 0s where the i,jth element equals 1 if product i is NOT producing at capacity
    notBinds <- matrix(as.numeric(capacities > shares),ncol=nprods,nrow=nprods,byrow=TRUE)
    ## create a TRUE/FALSE vector equal to TRUE if a single product firm is capacity constrained
    singleConstrained <- rowSums( ownerPre>0) == 1 & capacities == shares



    minD <- function(theta){

      alpha <- theta[1]
      sOut  <- theta[2]

      probs <- shares * (1 - sOut)
      elast <- -alpha *  matrix(prices * probs,ncol=nprods,nrow=nprods)
      diag(elast) <- alpha*prices + diag(elast)


      revenues <- probs * prices
      
      elastInv <- try(solve(elast * ownerPre * notBinds),silent=TRUE)
      if(any(class(elastInv)=="try-error")){elastInv <- MASS::ginv(elast * ownerPre *notBinds)}
      
      marginsCand <- -1 * as.vector(elastInv %*% (revenues * diag(ownerPre))) / revenues

      m1 <- margins - marginsCand
      m1 <- m1[!singleConstrained]
      m2 <- mktElast/(avgPrice * alpha ) - sOut
      measure <- sum(c(m1,m2)^2,na.rm=TRUE)

      #elast      <-   elast[isMargin,isMargin]
      #revenues   <-   revenues[isMargin]
      #ownerPre   <-   ownerPre[isMargin,isMargin]
      #margins    <-   margins[isMargin]

      #marginsCand <- -1 * as.vector(MASS::ginv(elasticity * ownerPre) %*% (revenues * diag(ownerPre))) / revenues
      #measure <- sum((margins - marginsCand)^2,na.rm=TRUE)

      #measure <- revenues * diag(ownerPre) + as.vector((elast * ownerPre) %*% (margins * revenues))
      #measure <- sum(measure^2,na.rm=TRUE)

      return(measure)
    }

    ## Constrain optimizer to look  alpha <0,  0 < sOut < 1
    lowerB <- c(-Inf,0)
    upperB <- c(-1e-10,.99999)


    if(!is.na(mktElast)){
      upperB[1] <- mktElast/avgPrice
    }

    minTheta <- optim(object@parmsStart,minD,
                      method="L-BFGS-B",
                      lower= lowerB,upper=upperB,
                      control=object@control.slopes)$par

    if(isTRUE(all.equal(minTheta[2],lowerB[2],check.names=FALSE))){

      warning("Estimated outside share is close to 0. Normalizing relative to largest good.")

      idx <- which.max(shares)
      shares[idx]
      priceOutside <- prices[idx]
      minTheta[2] <- 0
      object@normIndex <- idx

      meanval <- log(shares)  - log(shares[idx]) - minTheta[1] * (prices - priceOutside)

    }
    else{meanval <- log(shares * (1 - minTheta[2])) - log(minTheta[2]) - minTheta[1] * (prices - priceOutside)}

    if(isTRUE(all.equal(minTheta[2],upperB[2],check.names=FALSE))){stop("Estimated outside share is close to 1.")}



    names(meanval)   <- object@labels


    object@slopes      <- list(alpha=minTheta[1],meanval=meanval)
    object@shareInside <- 1-minTheta[2]
    object@priceOutside <- priceOutside
    object@mktSize <-  insideSize/object@shareInside

    return(object)


  }
)


#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "LogitNestsALM",
  definition=function(object){

    ## Uncover Demand Coefficents

    ownerPre     <-  object@ownerPre
    shares       <-  object@shares

    margins      <-  object@margins
    prices       <-  object@prices

    nprods       <-  length(shares)

    parmsStart   <- object@parmsStart
    nests        <- object@nests
    nestCnt      <- tapply(prices,nests,length)
    constraint   <- object@constraint
    output       <- object@output
    
    isSingletonNest <- nestCnt==1



    if(any(isSingletonNest)){
      warning("Some nests contain only one product; their nesting parameters are not identified.
              Normalizing these parameters to 1.")

    }


    if(!constraint){
      parmsStart   <- parmsStart[c(TRUE,TRUE,!isSingletonNest)] #always retain first two elements; these are
      # the initial value for price coefficient, outside share
    }

    ## Uncover price coefficient and mean valuation from margins and revenue shares


    nprods <- length(shares)

    sharesNests <- tapply(shares,nests,sum)[nests]

    sharesNests <- shares / sharesNests



    minD <- function(theta){

      alpha <- theta[1]
      sOut  <- theta[2]
      sigma <- as.numeric(isSingletonNest)
      sigma[!isSingletonNest] <- theta[-c(1,2)]

      probs <- shares * (1 - sOut)
      elast <- diag((1/sigma-1)*alpha)
      elast <- elast[nests,nests]
      elast <- elast * matrix(sharesNests*prices,ncol=nprods,nrow=nprods)
      elast <- -1*(elast + alpha * matrix(probs*prices,ncol=nprods,nrow=nprods))
      diag(elast) <- diag(elast) + (1/sigma[nests])*alpha*prices


      revenues <- probs * prices
      marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% (revenues * diag(ownerPre))) / revenues

      measure <- marginsCand-margins
      measure <- sum((measure)^2,na.rm=TRUE)

      return(measure)
    }

    ## Constrain optimizer to look  alpha <0, 0 < sOut < 1, sigma

    lowerB <- upperB <- rep(0,length(parmsStart))
    
    upperB[-1] <- 1

    if(output){lowerB[1] <- -Inf}
    else{upperB[1] <- Inf}
    
    minTheta <- optim(parmsStart,minD,method="L-BFGS-B",
                      lower= lowerB,upper=upperB,
                      control=object@control.slopes)


    minAlpha           <- minTheta$par[1]
    names(minAlpha)    <- "alpha"

    shareOut           <- minTheta$par[2]


    minSigma <-  as.numeric(isSingletonNest)
    minSigma[!isSingletonNest] <- minTheta$par[-c(1,2)]


    minSigmaOut        <- minSigma
    minSigma           <- minSigma[nests]
    names(minSigmaOut) <- levels(nests)


    meanval <-
      log(shares * (1 - shareOut)) - log(shareOut) -
      minAlpha*(prices - object@priceOutside) +
      (minSigma-1)*log(sharesNests)


    names(meanval)   <- object@labels


    object@slopes      <- list(alpha=minAlpha,sigma=minSigmaOut,meanval=meanval)
    object@shareInside <- 1-shareOut

    return(object)

  }

)


#'@rdname Params-Methods
#'@export
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
    mktElast     <-  object@mktElast
    output       <-  object@output
    outSign <- ifelse(output,1,-1)

    avgPrice <- weighted.mean(prices,shares)
    
    ## Uncover price coefficient and mean valuation from margins and revenue shares


    nprods <- length(shares)



    firmShares <- drop(ownerPre %*% shares)

    if(is.na(idx)){
      idxShare <- 1 - object@shareInside
      
    }
    else{
      idxShare <- shares[idx]
      
    }
    
    shareOut <-  1 - object@shareInside

    ## Minimize the distance between observed and predicted  ex Ante margins
    minD <- function(alpha){

      

      m1 <-  mktElast / (alpha * avgPrice) - shareOut
      
      m2 <- 1 - outSign*log(1-firmShares)/( alpha * firmShares)/margins

      measure <- sum(c(m1,m2)^2,na.rm=TRUE)

      return(measure)
    }

    if(output){bounds <- c(-1e6,0)}
    else{bounds <- c(0,1e6)}
    minAlpha <- optimize(minD,bounds,
                         tol=object@control.slopes$reltol)$minimum

    

    meanval <- log(shares) - log(idxShare)

    names(meanval)   <- object@labels

    object@slopes    <- list(alpha=minAlpha,meanval=meanval)
    object@mktSize   <- object@insideSize/object@shareInside
    

    return(object)
  }
)

#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "Auction2ndLogitALM",
  definition=function(object){

    ## Uncover Demand Coefficents

    ownerPre     <-  object@ownerPre
    shares       <-  object@shares
    margins      <-  object@margins
    mktElast     <-  object@mktElast
    prices       <-  object@prices
   output        <- object@output
   outSign <- ifelse(output,1,-1)
   
    avgPrice     <- sum(shares * prices,na.rm=TRUE)/sum(shares[!is.na(prices)])


    nprods <- length(object@shares)


    minD <- function(theta){

      alpha <- theta[1]
      sOut  <- theta[2]

      probs <- shares * (1 - sOut)

      firmShares <- drop(ownerPre %*% probs)


      m1 <- 1 - outSign*(log((1-firmShares))/( alpha * firmShares))/margins
      m2 <-  mktElast / (alpha * avgPrice) - sOut

      measure <- sum(c(m1 , m2)^2,na.rm=TRUE)

      return(measure)
    }

    ## Constrain optimizer to look  alpha <0,  0 < sOut < 1

    lowerB <- c(-Inf,1e-9)
    upperB <- c(-1e-10,.9999999999)

    if(!output){
      lowerB[1] <- 1e-10
      upperB[1] <- Inf
    }
    if(!is.na(mktElast)){upperB[1] <- mktElast/avgPrice}

    # minTheta <- optim(object@parmsStart,minD,
    #                   method="L-BFGS-B",
    #                   lower= lowerB,upper=upperB,
    #                   control=object@control.slopes)
    # 
    minTheta <- BBoptim(object@parmsStart,minD,
                      method="L-BFGS-B",
                      lower= lowerB,upper=upperB,quiet=TRUE,
                      control=object@control.slopes)

    if(minTheta$convergence != 0){
      warning("'calcSlopes' nonlinear solver may not have successfully converged. Reason: '",minTheta$message,"'")
    }
      minTheta <- minTheta$par
    
    if(isTRUE(all.equal(minTheta[2],lowerB[2],check.names=FALSE))){warning("Estimated outside share is close to 0. Normalizing relative to largest good.")
      idx <- which.max(shares)
      meanval <- log(shares) - log(shares[idx])
      minTheta[2] <- 0
      object@normIndex <- idx

    }
    else{ meanval <- log(shares * (1 - minTheta[2])) - log(minTheta[2]) }
    if(isTRUE(all.equal(minTheta[2],upperB[2],check.names=FALSE))){stop("Estimated outside share is close to 1.")}




    names(meanval)   <- object@labels


    object@slopes      <- list(alpha=minTheta[1],meanval=meanval)
    object@shareInside <- 1-minTheta[2]
    object@mktSize <- object@insideSize/object@shareInside


    return(object)

  }

)

#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "LogitALM",
  definition=function(object){

    ## Uncover Demand Coefficents

    ownerPre     <-  object@ownerPre
    shares       <-  object@shares
    margins      <-  object@margins
    prices       <-  object@prices
    mktElast     <-  object@mktElast
    priceOutside <- object@priceOutside
    output <- object@output
    outSign <- ifelse(output,-1,1)

    avgPrice <- sum(shares*prices)

    nprods <- length(object@shares)
    
    #margins <- margins*prices

    ##identify which products have enough margin information
    ##  to impute Bertrand margins
    #isMargin    <- matrix(margins,nrow=nprods,ncol=nprods,byrow=TRUE)
    #isMargin[ownerPre==0]=0
    #isMargin    <- !is.na(rowSums(isMargin))

    minD <- function(theta){

      alpha <- theta[1]
      sOut  <- theta[2]

      probs <- shares * (1 - sOut)
      elast <- -alpha *  matrix(prices * probs,ncol=nprods,nrow=nprods)
      diag(elast) <- alpha*prices + diag(elast)

      revenues <- probs * prices
      
      elastInv <- try(solve(elast * ownerPre),silent=TRUE)
      if(any(class(elastInv)=="try-error")){elastInv <- MASS::ginv(elast * ownerPre)}
      
      #marginsCand <- outSign/(alpha*(1-as.numeric(crossprod(ownerPre,predshares))))
      marginsCand <-  outSign*elastInv
      
      m1 <- margins - marginsCand
      m2 <- mktElast/(avgPrice * alpha ) - sOut
      measure <- sum((c(m1,m2)*100)^2,na.rm=TRUE)

      #elast      <-   elast[isMargin,isMargin]
      #revenues   <-   revenues[isMargin]
      #ownerPre   <-   ownerPre[isMargin,isMargin]
      #margins    <-   margins[isMargin]

      #marginsCand <- -1 * as.vector(MASS::ginv(elasticity * ownerPre) %*% (revenues * diag(ownerPre))) / revenues
      #measure <- sum((margins - marginsCand)^2,na.rm=TRUE)

      #measure <- revenues * diag(ownerPre) + as.vector((elast * ownerPre) %*% (margins * revenues))
      #measure <- sum(measure^2,na.rm=TRUE)

      return(measure)
    }

    ## Constrain optimizer to look  alpha <0,  0 < sOut < 1
    lowerB <- c(-Inf,0)
    upperB <- c(-1e-10,.99999)


    if(!is.na(mktElast)){
      upperB[1] <- mktElast/avgPrice
    }

    minTheta <- optim(object@parmsStart,minD,
                      method="L-BFGS-B",
                      lower= lowerB,upper=upperB,
                      control=object@control.slopes)$par

    if(isTRUE(all.equal(minTheta[2],lowerB[2],check.names=FALSE))){

      warning("Estimated outside share is close to 0. Normalizing relative to largest good.")

      idx <- which.max(shares)
      shares[idx]
      priceOutside <- prices[idx]
      minTheta[2] <- 0
      object@normIndex <- idx

      meanval <- log(shares)  - log(shares[idx]) - minTheta[1] * (prices - priceOutside)

    }
    else{meanval <- log(shares * (1 - minTheta[2])) - log(minTheta[2]) - minTheta[1] * (prices - priceOutside)}

    if(isTRUE(all.equal(minTheta[2],upperB[2],check.names=FALSE))){stop("Estimated outside share is close to 1.")}



    names(meanval)   <- object@labels


    object@slopes      <- list(alpha=minTheta[1],meanval=meanval)
    object@shareInside <- 1-minTheta[2]
    object@priceOutside <- priceOutside
    object@mktSize <-  object@insideSize/object@shareInside

    return(object)

  }

)


#'@rdname Params-Methods
#'@export
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
    diversion    <-  object@diversion
    mktElast     <-  object@mktElast 
    output       <-  object@output
    
    outSign <- ifelse(output,1,-1)
    shareOut <- 1 - shareInside

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

    ## Choose starting paramter values
    notMissing <- which(!is.na(margins))[1]
    
    parmStart <- (shares[notMissing] - outSign/margins[notMissing])/(shares[notMissing] - 1 ) 
    parmStart <- c(parmStart, exp(log(shares) - log(idxShare) - (parmStart - 1) * (log(prices) - log(idxPrice))))
    

    
    ## Uncover price coefficient and mean valuation from margins and revenue shares


    nprods <- length(shares)

   

    ## Minimize the distance between observed and predicted margins
    minD <- function(theta){

      gamma <- theta[1]
      meanval <- theta[-1]

      predshares <- meanval * (prices/idxPrice)^(1-gamma)
      predshares <- predshares/(is.na(idx) + sum(predshares) )
         
      preddiversion <-tcrossprod( 1/(1-predshares),predshares)
      diag(preddiversion) <- -1
       
      
      elasticity <- (gamma - 1 ) * matrix(predshares,ncol=nprods,nrow=nprods)
      diag(elasticity) <- -gamma + diag(elasticity)

      elastInv <- try(solve(elasticity * ownerPre),silent=TRUE)
      if(any(class(elastInv)=="try-error")){elastInv <- MASS::ginv(elasticity * ownerPre)}
      
      
      marginsCand <- -1 * outSign * as.vector(elastInv %*% (predshares * diag(ownerPre))) / predshares
      #measure <- sum((margins - marginsCand)^2,na.rm=TRUE)
      #FOC <- (shares * diag(ownerPre)) + (elasticity * ownerPre) %*% (shares * margins)
      
      m1 <- margins - marginsCand
      m2 <- predshares - shares
      m3 <- drop(diversion - preddiversion)
      m4 <- (mktElast + 1)/(1-gamma) - shareOut
      measure <- sum((c(m1, m2, m3, m4)*100)^2,na.rm=TRUE)
      
      #measure<-sum(FOC^2,na.rm=TRUE)

      return(measure)
    }

    
    ##  Constrained optimizer to look for solutions where gamma>1
    
    
    lowerB <- upperB <- rep(Inf,length(parmStart))
    lowerB <- lowerB * -1
    
    if(output){lowerB[1] <- 1}
    else{upperB[1] <- 1 }
    
    minTheta <- optim(parmStart,minD,method="L-BFGS-B",
                      lower= lowerB,upper=upperB,
                      control=object@control.slopes)


    if(minTheta$convergence != 0){
      warning("'calcSlopes' nonlinear solver did not successfully converge. Reason: '",minTheta$message,"'")
    }



    minGamma           <- minTheta$par[1]
    names(minGamma)    <- "Gamma"

    meanval <-  minTheta$par[-1]

    if(!is.na(idx)) meanval <- meanval/meanval[idx]
    
    names(meanval)   <- object@labels

    
    object@slopes    <- list(alpha=alpha,gamma=minGamma,meanval=meanval)
    object@priceOutside <- idxPrice
    object@mktSize <- insideSize*(1+alpha)


    return(object)
  }
)


#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "CESALM",
  definition=function(object){

    ## Uncover Demand Coefficents

    ownerPre     <-  object@ownerPre
    shares       <-  object@shares
    margins      <-  object@margins
    prices       <-  object@prices
    mktElast    <-   object@mktElast
    priceOutside <- object@priceOutside

    avgPrice <- sum(prices*shares)/sum(shares)

    nprods <- length(object@shares)

    ##identify which products have enough margin information
    ##  to impute Bertrand margins
    #isMargin    <- matrix(margins,nrow=nprods,ncol=nprods,byrow=TRUE)
    #isMargin[ownerPre==0]=0
    #isMargin    <- !is.na(rowSums(isMargin))

    minD <- function(theta){

      gamma <- theta[1]
      sOut  <- theta[2]


      probs <- shares * (1 - sOut)

      elasticity <- (gamma - 1 ) * matrix(probs,ncol=nprods,nrow=nprods)
      diag(elasticity) <- -gamma + diag(elasticity)


      elastInv <- try(solve(elasticity * ownerPre),silent=TRUE)
      if(any(class(elastInv)=="try-error")){elastInv <- MASS::ginv(elasticity * ownerPre)}
      
      
      marginsCand <- -1 * as.vector(elastInv %*% (probs * diag(ownerPre))) / probs

      m1 <- margins - marginsCand
      m2 <- (mktElast + 1)/(1-gamma) - sOut
     
      measure <- sum(c(m1 , m2)^2,na.rm=TRUE)


      return(measure)
    }

    ## Constrain optimizer to look  gamma > 1,  0 < sOut < 1
    lowerB <- c(1,0)
    upperB <- c(Inf,.99999)



    minGamma <- optim(object@parmsStart,minD,
                      method="L-BFGS-B",
                      lower= lowerB,upper=upperB,
                      control=object@control.slopes)$par

    if(isTRUE(all.equal(minGamma[2],lowerB[2],check.names=FALSE))){warning("Estimated outside share is close to 0. Normalizing relative to largest good.")

      idx <- which.max(shares)
      object@normIndex <- idx
      priceOutside <- priceOutside[idx]
      minGamma[2] <- 0

      meanval <- log(shares) - log(shares[idx]) + (minGamma[1] - 1) * (log(prices) - log(priceOutside))
    }
    else{ meanval <- log(shares * (1 - minGamma[2])) - log(minGamma[2]) + (minGamma[1] - 1) * (log(prices) - log(object@priceOutside))}
    if(isTRUE(all.equal(minGamma[2],upperB[2],check.names=FALSE))){stop("Estimated outside share is close to 1.")}


    meanval <- exp(meanval)



    names(meanval)   <- object@labels


    object@slopes      <- list(alpha=1/(1 - minGamma[2]) - 1  ,gamma=minGamma[1],meanval=meanval)
    object@shareInside <- 1-minGamma[2]
    object@priceOutside <- priceOutside
    object@mktSize <- object@insideSize*(1+object@slopes$alpha)
    return(object)

  }

)

#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "CESNests",
  definition=function(object){

    ## Uncover Demand Coefficents


    ownerPre     <-  object@ownerPre
    shares       <-  object@shares
    margins      <-  object@margins
    prices       <-  object@prices
    idx          <-  object@normIndex
    shareInside  <-  object@shareInside
    nests        <- object@nests
    parmsStart   <- object@parmsStart
    constraint   <- object@constraint
    insideSize   <- object@insideSize

    nestCnt      <- tapply(prices,nests,length)


    isSingletonNest <- nestCnt==1

    if(any(isSingletonNest)){
      warning("Some nests contain only one product; their nesting parameters are not identified.
              Normalizing these parameters to 1.")

    }



    if(!constraint){
      parmsStart   <- parmsStart[c(TRUE,!isSingletonNest)] #always retain first element; this is
      # the initial value for price coefficient
    }

    ## Uncover price coefficient and mean valuation from margins and revenue shares


    nprods <- length(shares)

    ## identify which products have enough margin
    ## information to impute Bertrand margins
    isMargin    <- matrix(margins,nrow=nprods,ncol=nprods,byrow=TRUE)
    isMargin[ownerPre==0]=0
    isMargin    <- !is.na(rowSums(isMargin))

    sharesNests <- tapply(shares,nests,sum)[nests]

    sharesNests <- shares / sharesNests


    ## back out the parameter on the numeraire, when appropriate
    if(shareInside<1) {alpha <- 1/shareInside -1}
    else{ alpha <- NULL}






    ## Estimate parameters by
    ## Minimizing the distance between observed and predicted margins
    minD <- function(theta){

      gamma  <- theta[1]

      sigma <- as.numeric(!isSingletonNest) # normalize singleton nest parms to 0
      sigma[!isSingletonNest] <- theta[-1]

      elast <- diag(sigma - gamma)

      elast <- elast[nests,nests]
      elast <- elast * matrix(sharesNests,ncol=nprods,nrow=nprods)
      elast <- elast + (gamma-1) * matrix(shares,ncol=nprods,nrow=nprods)

      diag(elast) <- diag(elast) - sigma[nests]

      #marginsCand <- -1 * as.vector(MASS::ginv(elast * ownerPre) %*% (shares * diag(ownerPre))) / shares
      #measure <- sum((margins - marginsCand)^2,na.rm=TRUE)


      elast      <- elast[isMargin,isMargin]
      shares     <- shares[isMargin]
      ownerPre   <- ownerPre[isMargin,isMargin]
      margins    <- margins[isMargin]

      FOC <- (shares * diag(ownerPre)) + (elast * ownerPre) %*% (shares * margins)
      measure<-sum(FOC^2,na.rm=TRUE)

      return(measure)
    }

    ## Constrain optimizer to look for solutions where sigma_i > gamma > 1 for all i
    constrA <- diag(length(parmsStart))
    constrA[-1,1] <- -1

    constrB <- rep(0,length(parmsStart))
    constrB[1] <- 1

    minTheta <- constrOptim(parmsStart,minD,grad=NULL,ui=constrA,ci=constrB,
                            control=object@control.slopes)


    if(minTheta$convergence != 0){
      warning("'calcSlopes' nonlinear solver did not successfully converge. Reason: '",minTheta$message,"'")
    }





    minGamma <- minTheta$par[1]
    names(minGamma) <- "Gamma"

    minSigma <-  as.numeric(!isSingletonNest)
    minSigma[!isSingletonNest] <- minTheta$par[-1]


    minSigmaOut        <- minSigma
    minSigma           <- minSigma[nests]
    names(minSigmaOut)    <- levels(nests)


    if(is.na(idx)){
      idxShare      <- 1 - sum(shares)
      idxShareNests <- 1
      idxPrice      <- object@priceOutside
      idxSigma      <- 0
    }

    else{

      idxShare      <- shares[idx]
      idxShareNests <- sharesNests[idx]
      idxPrice      <- prices[idx]
      idxSigma      <- minSigma[idx]
    }


    meanval <-
      log(shares) - log(idxShare) + (minGamma - 1) *
      (log(prices) - log(idxPrice)) -
      (minSigma-minGamma)/(minSigma-1)*log(sharesNests) +
      (idxSigma - minGamma)/(idxSigma-1)*log(idxShareNests)

    meanval <- exp( (minSigma-1)/(minGamma-1) * meanval )


    names(meanval)   <- object@labels

    object@slopes    <- list(alpha=alpha,gamma=minGamma,sigma=minSigmaOut,meanval=meanval)
    object@mktSize <- insideSize*(1+alpha)

    return(object)
  }
)



#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "BargainingLogit",
  definition=function(object){
    
    ## Uncover Demand Coefficents
    
    
    ownerPre     <-  object@ownerPre
    shares       <-  object@shares
    margins      <-  object@margins
    prices       <-  object@prices
    idx          <-  object@normIndex
    mktElast     <-  object@mktElast
    shareInside  <-  object@shareInside
    diversion    <-  object@diversion
    barg         <-  object@bargpowerPre
    output       <-  object@output
    outSign <- ifelse(output,-1,1)
    
    
    if(is.na(idx)){
      idxShare <- 1 - shareInside
      idxPrice <- object@priceOutside
    }
    else{
      idxShare <- shares[idx]
      idxPrice <- prices[idx]
    }
    
    ## Choose starting parameter values
    notMissing <- which(!is.na(margins))[1]
    
    ##Start at 50/50 Bargaining
    parmStart <- -1*outSign*log(1- shares[notMissing])/(margins[notMissing]*prices[notMissing]*(1 - shares[notMissing])*(shares[notMissing]/(1- shares[notMissing]) - log(1- shares[notMissing])))
   
    mvalStart <-  log(shares) - log(idxShare) - parmStart * (prices - idxPrice)
    if(!is.na(idx)) mvalStart <-  mvalStart[-idx]
    parmStart <- c(parmStart, mvalStart)
    
    ## if any bargaining parameters are missing, set starting bargaining parameter to 0.5
    if(any(is.na(barg))){parmStart <- c(parmStart,0.5)}
    
    
    nprods <- length(shares)
    
    avgPrice <- sum(shares * prices, na.rm=TRUE) / sum(shares)
    
 
    
    
    nParm <- length(parmStart)
    ## identify which products have enough margin information
    ##  to impute Bertrand margins
    #isMargin    <- matrix(margins,nrow=nprods,ncol=nprods,byrow=TRUE)
    #isMargin[ownerPre==0]=0
    #isMargin    <- !is.na(rowSums(isMargin))
    
    
    ## Minimize the distance between observed and predicted margins
    minD <- function(theta){
      
      alpha <- theta[1]
      
      if(any(is.na(barg))){
        thisBarg <- theta[nParm]
        theta <- theta[-nParm]
        barg[is.na(barg)] <- thisBarg
        }
      
      if(!is.na(idx)){
        meanval <- rep(0,nprods)
        meanval[-idx] <- theta[-1]
      }
      else{meanval <- theta[-1]}
      
      barg <- barg/(1-barg)
      
      probs <- shares
      
      predshares <- exp(meanval + alpha*(prices-idxPrice))
      predshares <- predshares/(is.na(idx) + sum(predshares) )
      
      preddiversion <-predshares/(1-predshares)
      
      
      if(!is.na(mktElast)){
        shareInside <-   1 - mktElast/( alpha * avgPrice )
        probs <- probs/sum(probs,na.rm=TRUE) * shareInside
        
      }
      
      ownerPreInv <- ownerPre
      #diag(ownerPreInv) <- -1*diag(ownerPreInv)
      ownerPreInv <- -1*t(ownerPreInv * predshares)
      diag(ownerPreInv) <- diag(ownerPre) + diag(ownerPreInv)
      
      tmp <- try(solve(ownerPreInv),silent=TRUE)
      if(any(class(tmp)=="try-error")){ownerPreInv=MASS::ginv(ownerPreInv)}
      else{ ownerPreInv <- tmp}
      
      marginsCand <-  ownerPreInv %*% ((log(1-predshares)*diag(ownerPre))/(-1*outSign*alpha*(barg*predshares/(1-predshares) - 
                                                                   log(1-predshares))))
      marginsCand <-  as.vector(marginsCand)
      m1 <- margins - marginsCand/prices
      m2 <- (predshares - probs)
      measure <- sum((c(m1,m2)*100)^2,na.rm=TRUE)
      
      return(measure)
    }
 
    
    
    
    
    ##  Constrained optimizer to look for solutions where alpha<0,  
    lowerB <- upperB <- c(1e6,rep(12,length(parmStart)-1))
    lowerB <- lowerB * -1
    
    if(output){
      upperB[1] <- -1e-5
    }
    else{
      lowerB[1] <- 1e-5
    }
    
    minTheta <- optim(parmStart,minD,method="L-BFGS-B",
                      lower= lowerB,upper=upperB,
                      control=object@control.slopes)
    
    
    if(minTheta$convergence != 0){
      warning("'calcSlopes' nonlinear solver may not have successfully converge. Reason: '",minTheta$message,"'")
    }
    
    
    # ui=diag(length(parmStart))
    # #ui[1,1] <- -1
    # if(!is.na(idx)){ui[-1,1] <- prices[-idx] - idxPrice}
    # else{ui[-1,1] <- prices - idxPrice}
    # ci_hi=rep(log(.9999/(1-.9999)), length(mvalStart))
    # ci_low=rep(log((1-.9999)/.9999), length(mvalStart))
    # 
    # ui=rbind(-ui,ui[-1,])
    # ci=c(0,-ci_hi,ci_low)
    # 
    # minTheta <- constrOptim(parmStart,minD,grad=NULL,ui=ui,ci=ci)
    # 
    minAlpha           <- minTheta$par[1]
    names(minAlpha)    <- "alpha"
    
    if(any(is.na(barg))){
    minBarg <- minTheta$par[nParm]
    minTheta$par <- minTheta$par[-nParm]
    object@bargpowerPre[is.na(barg)] <- object@bargpowerPost[is.na(object@bargpowerPost)] <- minBarg
    
    }
    
    if(is.na(idx)) meanval <-  minTheta$par[-1]
    else{
      meanval <- rep(0,nprods)
      meanval[-idx] <- minTheta$par[-1]
    }
    
    names(meanval)   <- object@labels
    
    
    
    
    object@slopes    <- list(alpha=minAlpha,meanval=meanval)
    if(any(is.na(barg))){ object@slopes$barg <- minBarg }
    
    object@priceOutside <- idxPrice
    object@mktSize <- object@insideSize / sum(shares)
    
    
    return(object)
  }
)



#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "Bargaining2ndLogit",
  definition=function(object){
    
    ## Uncover Demand Coefficents
    
    
    ownerPre     <-  object@ownerPre
    shares       <-  object@shares
    margins      <-  object@margins
    prices       <-  object@prices
    idx          <-  object@normIndex
    mktElast     <-  object@mktElast
    barg         <-  object@bargpowerPre
    output       <-  object@output
    outSign <- ifelse(output,-1,1)
    
    
    
    avgPrice <- weighted.mean(prices,shares)
    
    ## Uncover price coefficient and mean valuation from margins and revenue shares
    
    
    nprods <- length(shares)
    
    
    
    firmShares <- drop(ownerPre %*% shares)
    
    if(is.na(idx)){
      idxShare <- 1 - object@shareInside
      
    }
    else{
      idxShare <- shares[idx]
      
    }
    
    shareOut <-  1 - object@shareInside
    
    ## Minimize the distance between observed and predicted  ex Ante margins
    minD <- function(alpha){
      
      
      
      m1 <-  mktElast / (alpha * avgPrice) - shareOut
      
      m2 <- 1 - (1-barg)*(log(1-firmShares)/( -1*outSign*alpha * firmShares))/margins
      
      measure <- sum(c(m1,m2)^2,na.rm=TRUE)
      
      return(measure)
    }
    
    minAlpha <- optimize(minD,c(-1e6,0),
                         tol=object@control.slopes$reltol)$minimum
    
  
    
    
    meanval <- log(shares) - log(idxShare)
    
    names(meanval)   <- object@labels
    
    object@slopes    <- list(alpha=minAlpha,meanval=meanval)
    object@mktSize   <- object@insideSize/object@shareInside
    
    
    return(object)
  }
)


#'@rdname Params-Methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "VertBargBertLogit",
  definition=function(object){
   
    
    
    constrain <- object@constrain
    
    is2nd <- any(grepl("2nd",class(object)))
    
    up <- object@up
    down <- object@down
    
    sigma <- 0.1
    
    if(grepl("Nest",class(object))){
      nests <- down@nests
      sigma <- rep(0.1,nlevels(nests))
      names(sigma) <- levels(nests)
    }
    
    constrain <- object@constrain
    

    owner.up.pre <- up@ownerPre
    owner.down.pre <- down@ownerPre
    
    owner.up.post <- up@ownerPost
    owner.down.post <- down@ownerPost
  
    
    pricesUp    <- up@prices
    marginsUp   <- up@margins
    
    marginsDown <- down@margins
    idx         <- down@normIndex
    sharesDown  <- down@shares
    
    
    pricesDown <- down@prices
    down@pricePre <- pricesDown
    
    if(is.na(idx)){
      idxShare <- 1 - down@shareInside
      idxPrice <- down@priceOutside
    }
    else{
      idxShare <- sharesDown[idx]
      idxPrice <- pricesDown[idx]
    }
    
    
    id <- data.frame(up.firm=owner.up.pre,
                     down.firm=owner.down.pre
                     )
    
    
   
    
    nprods <- nrow(id)
    
    if(constrain == "pair"){
      id <- with(id,interaction(up.firm,down.firm))
    }
    else if(constrain =="wholesaler"){
      id <- with(id,up.firm)
    }
    else if(constrain =="retailer"){
      id <- with(id,down.firm)
    }
    else{ id <- rep(1,nprods)}
    
    id <- factor(id)
   
    marginsUp <- marginsUp*pricesUp
    if (!is2nd) marginsDown <- marginsDown*pricesDown
     
    #set starting value for bargaining parameter equal to 0.5
    bStart <- rep(0.5,nlevels(id))
    # set starting value for alpha equal to single product 
    # unintegrated firm
    firstAvail <- which(!is.na(marginsDown))[1]
    alphaStart <- -1/(marginsDown[firstAvail]*(1 - sharesDown[firstAvail]))
    
    parmStart <- c(alphaStart,bStart)
    
   
    
    
    
   
    div <- tcrossprod(1/(1-sharesDown),sharesDown)*sharesDown
    diag(div) <- -sharesDown
    div <- as.vector(div)
    
    
    vertFirms <- intersect(owner.up.pre,owner.down.pre)
    
    ownerDownMat <-  ownerToMatrix(down, preMerger=TRUE)
    ownerBargUpVert<- ownerToMatrix(up, preMerger=TRUE)
    
    ownerDownMatVertical <- matrix(0,nrow=nprods,ncol=nprods)
    
    minD <- function(theta){
      
      alpha <- theta[1]
      b <- theta[-1]
      b <- b[as.numeric(id)]
      
  
      b[owner.up.pre == owner.down.pre] <- 1 
      
      for( v in vertFirms){
        
        vertrows <- owner.up.pre != v  & owner.down.pre == v
        ownerBargUpVert[vertrows, owner.up.pre == v] <- -(1-b[vertrows])/b[vertrows]
      }
      
      ## set integrated margin disagreement payoff to 0,
      ## constrain upstream integrated margin to zero
      
      for(n in which(owner.up.pre==owner.down.pre)){
        ownerBargUpVert[n,-n] <- ownerBargUpVert[-n,n] <- 0
      }
      
      ownerBargDownVert  <-  ownerDownMat  * (1-b)/b
      
      for( v in vertFirms){
        
        vertrows <-  owner.up.pre == v  & owner.down.pre != v
        
        ## only change downstream matrix when firms are playing Bertrand
        if(!is2nd){ownerDownMatVertical[owner.down.pre == v, vertrows] <- 1}
        #ownerDownMatVertical[owner.down.pre == v, !vertrows] <- 0
        
        
        ownerBargDownVert [vertrows, owner.down.pre == v] <- -1
        
      }
      
      #ownerDownMatVertical[!owner.down.pre %in% vertFirms, ] <- 0
      
      down@ownerPre <- ownerDownMat
      
      if(is2nd) mval <- log(sharesDown) - log(idxShare) - alpha*(pricesUp - idxPrice)
      else{mval <- log(sharesDown) - log(idxShare) - alpha*(pricesDown - idxPrice)}
      
      down@slopes <- list(alpha = alpha,
                          meanval = mval,
                          sigma=sigma
                          )
      
      marginsCandDown <- calcMargins(down, preMerger= TRUE,level=TRUE)
       
      shareCandDown   <- calcShares(down,preMerger=TRUE,revenue=FALSE)

        if(!is2nd){
        elast <-  -alpha*tcrossprod(sharesDown)
        diag(elast) <- alpha*sharesDown + diag(elast)
        elast.inv <- try(solve(ownerDownMat * elast),silent=TRUE)
        if(any(class(elast.inv) == "try-error")){elast.inv <- MASS::ginv(ownerDownMat * elast)}
        
        marginsCandDown <- marginsCandDown - elast.inv %*% ( (ownerDownMatVertical * elast) %*% (marginsUp) )
      }
      
      depVar <- as.vector((ownerBargUpVert  * div) %*% marginsUp)
      regressor <- as.vector( ( ownerBargDownVert  * div) %*% marginsCandDown)
      
      err <- c(depVar - regressor, marginsDown - marginsCandDown 
               , (sharesDown - shareCandDown)
               )
      return(sum((err)^2,na.rm = TRUE))
    }
    
    #optmethod <- "L-BFGS-B"
    #if(length(bStart) ==1) optmethod <- "Brent"
    lowerB <- rep(.01,length(parmStart))
    lowerB[1] <- -1e9
    upperB <- rep(.99, length(parmStart))
    upperB[1] <- -1e-9
    
    #thetaOpt <- optim(parmStart,minD,method=optmethod,lower = lowerB,upper = upperB)
    thetaOpt <- BBoptim(parmStart,minD,lower = lowerB,upper = upperB,control = object@control.slopes,quiet=TRUE)
    
    if(thetaOpt$convergence !=0){
      warning("Calibration routine may not have converged. Optimizer Reports:\n\t",thetaOpt$message)}
    
    ## Pre-merger bargaining parameter
    alphaOpt <- thetaOpt$par[1]
    
    if(!is2nd) {mvalOpt <- log(sharesDown) - log(idxShare) - alphaOpt*(pricesDown - idxPrice)}
    else{mvalOpt <- log(sharesDown) - log(idxShare)- alphaOpt*(pricesUp - idxPrice)}
    
    bOpt     <- thetaOpt$par[-1]
    bargparmPre <- bargparmPost <-  bOpt[as.numeric(id)]
    bargparmPre[owner.up.pre  == owner.down.pre ] <- 1 
    names(bargparmPre) <- down@labels
    
    ## Post-merger bargaining parameter
    
    #owner.up.pre <- up@ownerPost
    #owner.down <- down@ownerPost
    bargparmPost[owner.up.post  == owner.down.post ] <- 1 
    names(bargparmPost) <- down@labels
      
      
    down@slopes <- list(alpha=alphaOpt,meanval=mvalOpt,sigma=sigma)  
    down@mktSize <- down@insideSize/down@shareInside
    object@down <- down
    
    up@bargpowerPre <- bargparmPre
    up@bargpowerPost <- bargparmPost
    object@up <- up
  
    
    
    
    object <- ownerToMatrix(object, preMerger=TRUE) #create ownership matrices
    object <- ownerToMatrix(object, preMerger=FALSE) #create ownership matrices
    
    return(object)
    
  }
)



#'@rdname Params-Methods
#'@export
setMethod(
  f= "getParms",
  signature= "Bertrand",
  definition=function(object,digits=10){
    if(is.list(object@slopes)){
      result <- lapply(object@slopes,round,digits=digits)
    }
    else{

      result <-  list(slopes = round(object@slopes,digits),
                      intercepts =  round(object@intercepts,digits)
      )

    }

    result$mc <- round(calcMC(object, preMerger=TRUE),digits)

    return(result)

  })

#'@rdname Params-Methods
#'@export
setMethod(
  f= "getParms",
  signature= "VertBargBertLogit",
  definition=function(object,digits=10){
    up <- object@up
    down <- object@down
    
    if(is.list(down@slopes)){
      result <- lapply(down@slopes,round,digits=digits)
    }
    else{
      
      result <-  list(slopes = round(object@slopes,digits),
                      intercepts =  round(object@intercepts,digits)
      )
      
    }
    
    mcPre <- calcMC(object, preMerger=TRUE)
    result$mcUpPre <- round(mcPre$up,digits)
    result$mcDownPre <- round(mcPre$down,digits)
    result$bargpower <- round(up@bargpowerPre,digits)
    
    return(result)
    
  })


#'@rdname Params-Methods
#'@export
setMethod(
  f= "getNestsParms",
  signature= "PCAIDSNests",
  definition=function(object){

    nests <- object@nests

    nNests <- nlevels(nests)

    labels <- levels(nests)

    nestWeights <- diag(nNests)
    nestWeights[upper.tri(nestWeights)] <- nestWeights[lower.tri(nestWeights)] <- object@nestsParms

    dimnames(nestWeights) <- list(labels,labels)

    return(nestWeights)
  }
)
