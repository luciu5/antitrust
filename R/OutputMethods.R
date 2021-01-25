#' @title Output Methods
#' @name Output-Methods
#' @docType methods
#'
#' @aliases calcQuantities
#' calcQuantities,ANY-method
#' calcQuantities,Logit-method
#' calcQuantities,CES-method
#' calcQuantities,Linear-method
#' calcQuantities,LogLin-method
#' calcQuantities,LogitCap-method
#' calcQuantities,Logit-method
#' calcQuantities,Cournot-method
#' calcQuantities,Stackelberg-method
#' calcQuantities,AIDS-method
#' calcShares
#' calcShares,ANY-method
#' calcShares,AIDS-method
#' calcShares,CES-method
#' calcShares,CESNests-method
#' calcShares,Linear-method
#' calcShares,Logit-method
#' calcShares,LogitNests-method
#' calcShares,Auction2ndLogit-method
#' calcShares,Auction2ndLogitNests-method
#' calcShares,Cournot-method
#' calcRevenues
#' calcRevenues,ANY-method
#' calcRevenues,Bertrand-method
#' calcRevenues,CES-method
#' calcRevenues,AIDS-method
#' calcRevenues,Cournot-method
#' calcRevenues,VertBargBertLogit-method
#'
#' @description This section contains three types of methods: calcShares, calcQuantities, and calcRevenues.
#'  calcShares computes equilibrium product shares assuming that firms are playing a
#'  Nash-Bertrand or Cournot  game. \sQuote{revenue} takes
#'  on a value of TRUE or FALSE, where TRUE calculates revenue shares,
#'  while FALSE calculates quantity shares.
#'
#'  calcQuantities computes equilibrium product quantities assuming that firms are playing a
#'  Nash-Bertrand, 2nd Score Auction, or Cournot game. Setting `market' to TRUE returns total market quantity.
#'
#'  calcRevenues computes equilibrium product revenues assuming that firms are playing a
#'  Nash-Bertrand, 2nd Score Auction, or Cournot game. Setting `market' to TRUE returns total market revenue.
#'
#' @param object An instance of one of the classes listed above.
#' @param preMerger If TRUE, returns pre-merger outcome. If
#' FALSE, returns post-merger outcome.  Default is TRUE.
#' @param revenue If TRUE, returns revenues. If FALSE,
#' returns quantities. Default is TRUE.
#' @param market If TRUE, reports market-level summary.
#' Otherwise reports product/plant level summary. Default is FALSE.
#' @param exAnte If \sQuote{exAnte} equals TRUE then the
#' \emph{ex ante} expected result for each firm is produced, while FALSE produces the
#' expected result conditional on each firm winning the auction. Default is FALSE.
#' @param ... Additional arguments to pass to \code{calcQuantities}.
#'
#' @include PSMethods.R
#' @keywords methods
NULL

setGeneric (
  name= "calcQuantities",
  def=function(object,...){standardGeneric("calcQuantities")}
)
setGeneric (
  name= "calcRevenues",
  def=function(object,...){standardGeneric("calcRevenues")}
)
setGeneric (
  name= "calcShares",
  def=function(object,...){standardGeneric("calcShares")}
)


#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcQuantities",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE,market=FALSE,...){

    slopes <- object@slopes
    intercepts <- object@intercepts
    quantityStart <- object@quantityStart
    quantityStart[is.na(quantityStart)] = 0

    if(preMerger){
      if(market) return(sum(object@quantityPre, na.rm=TRUE))
      owner  <- object@ownerPre
      products  <- object@productsPre
      cap <- object@capacitiesPre
    }
    else{
      if(market) return(sum(object@quantityPost, na.rm=TRUE))
      owner <-  object@ownerPost
      products <-  object@productsPost
      cap <- object@capacitiesPost
    }

    nprods <- ncol(products)
    #isProducts <- rowSums(products) > 0
    products <- as.vector(products)

    isConstrained <- is.finite(cap)

    FOC <- function(quantCand){

      #quantCand <- quantCand^2 # constrain positive

      thisQuant <- rep(0,length(products))
      thisQuant[products] = quantCand
      thisQuant = matrix(thisQuant,ncol=nprods)

      if(preMerger){ object@quantityPre  <- thisQuant}
      else{          object@quantityPost <- thisQuant}

      thisPrice <- calcPrices(object, preMerger= preMerger)

      thisMC <- calcMC(object, preMerger= preMerger)


      mktQuant <- colSums(thisQuant, na.rm=TRUE)
      plantQuant <- rowSums(thisQuant, na.rm=TRUE)

      thisPartial <- ifelse(object@demand=="linear",
                            slopes,
                            exp(intercepts)*slopes*mktQuant^(slopes - 1))


      thisFOC <- (t(thisQuant) * thisPartial) %*% owner + thisPrice
      thisFOC <- t(thisFOC)/thisMC - 1
      #thisFOC <- t(t(thisFOC)/thisPrice) # rescale
      #thisCons <- (plantQuant - cap)/cap # rescale
      #thisFOC[isConstrained,] <- thisFOC[isConstrained,] +
      #   thisCons[isConstrained] +
      #  sqrt(thisFOC[isConstrained,]^2 +
      #         thisCons[isConstrained]^2)
      thisFOC <- as.vector(thisFOC)[products]
      return(thisFOC)
    }


    #quantityStart <- sqrt(object@quantityStart[products]) #constrain positive
    quantityStart <- ifelse(quantityStart >= cap, cap-1, quantityStart)
    quantityStart <- quantityStart[products]


    ## Find price changes that set FOCs equal to 0
    minResult <- BB::BBsolve( quantityStart,FOC, quiet=TRUE,control=object@control.equ)

    if(minResult$convergence != 0){warning("'calcQuantities' nonlinear solver may not have successfully converged. 'BBsolve' reports: '",minResult$message,"'")}

    quantEst        <- rep(0, length(products))
    quantEst[products] <- minResult$par#^2
    quantEst <- matrix(quantEst,ncol = nprods)

    dimnames(quantEst) <- object@labels

    return(quantEst)
  })


#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcQuantities",
  signature= "Stackelberg",
  definition=function(object,preMerger=TRUE,market=FALSE,...){

    slopes <- object@slopes
    intercepts <- object@intercepts
    isLinearD <- object@demand=="linear"


    if(preMerger){
      if(market) return(sum(object@quantityPre,na.rm=TRUE))
      owner  <- object@ownerPre
      products  <- object@productsPre
      isLeader <- object@isLeaderPre
    }
    else{

      if(market) return(sum(object@quantityPost,na.rm=TRUE))
      owner <-  object@ownerPost
      products <-  object@productsPost
      isLeader <- object@isLeaderPost
    }

    nprods <- ncol(products)
    isProducts <- rowSums(products) > 0
    products <- as.vector(products)

    FOC <- function(quantCand){

      quantCand <- quantCand^2 # constrain positive

      allquant <- rep(0,length(products))
      allquant[products] <- quantCand
      quantCand <- matrix(allquant,ncol=nprods)

      if(preMerger){ object@quantityPre  <- quantCand}
      else{          object@quantityPost <- quantCand}

      thisPrice <- calcPrices(object, preMerger= preMerger)

      thisMC <- calcMC(object, preMerger= preMerger)
      thisdMC <- calcdMC(object, preMerger= preMerger)

      mktQuant <- colSums(quantCand, na.rm=TRUE)
      ownerQuant <- owner %*% quantCand

      thisPartial <- ifelse(isLinearD,
                            slopes,
                            exp(intercepts)*slopes*mktQuant^(slopes - 1))

      dthisPartial <- ifelse(isLinearD,
                             0,
                             exp(intercepts)*slopes*(slopes - 1)*mktQuant^(slopes - 2))


      demPass <- dthisPartial * t(!isLeader * ownerQuant)
      thisPass <- -t((thisPartial + demPass)/
                       (2*thisPartial  + t(t(demPass) - thisdMC)))


      thisPass[isLeader | !isProducts] <- 0


      thisFOC <- (t(quantCand) * thisPartial  + t(isLeader * quantCand) * thisPartial*colSums(thisPass)) %*% owner + thisPrice
      thisFOC <- t(thisFOC) - thisMC

      thisFOC <- thisFOC[isProducts,]
      return(as.vector(thisFOC))
    }


    quantityStart <- sqrt(object@quantityStart[products]) #constrain positive

    ## Find price changes that set FOCs equal to 0
    minResult <- BBsolve( quantityStart,FOC, quiet=TRUE,control=object@control.equ,...)

    if(minResult$convergence != 0){warning("'calcQuantities' nonlinear solver may not have successfully converged. 'BBsolve' reports: '",minResult$message,"'")}

    quantEst        <- rep(NA, length(products))
    quantEst[products] <- minResult$par^2
    quantEst <- matrix(quantEst,ncol = nprods)

    dimnames(quantEst) <- object@labels

    return(quantEst)
  })

#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcQuantities",
  signature= "Linear",
  definition=function(object,preMerger=TRUE, market = FALSE){

    if(preMerger){ prices <- object@pricePre}
    else{          prices <- object@pricePost}

    slopes    <- object@slopes
    intercept <- object@intercepts


    quantities <- as.vector(intercept+slopes %*% prices)
    names(quantities) <- object@labels

    if (market) quantities <- sum(quantities, na.rm=TRUE)

    return(quantities)

  }
)

#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcQuantities",
  signature= "Logit",
  definition=function(object,preMerger=TRUE, market=FALSE){

    mktSize <- object@mktSize

    shares <- calcShares(object, preMerger= preMerger, revenue = FALSE)
    if(market) shares <- sum(shares,na.rm=TRUE)
    return(mktSize*shares)


  })

#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcQuantities",
  signature= "LogLin",
  definition=function(object,preMerger=TRUE,market=FALSE,...){


    if(preMerger){ prices <- object@pricePre}
    else{          prices <- object@pricePost}

    slopes    <- object@slopes
    intercept <- object@intercepts
    quantities <- exp(intercept) * apply(prices^t(slopes),2,prod)
    names(quantities) <- object@labels

    if (market) quantities <- sum(quantities, na.rm=TRUE)

    return(quantities)

  }
)

#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcQuantities",
  signature= "AIDS",
  definition=function(object,preMerger=TRUE, market=FALSE){

    if(preMerger){ prices <- object@pricePre}
    else{          prices <- object@pricePost}

    insideSize <- object@insideSize

    shares <- calcShares(object, preMerger = preMerger, revenue=TRUE)

    if(!preMerger){

      # source: epstein and rubinfeld Assumption C (appendix)
      insideSizeDelta <- sum( shares * object@priceDelta) * (object@mktElast + 1)

      insideSize <- insideSize*( 1 + insideSizeDelta)

    }

    shares <- calcShares(object, preMerger= preMerger, revenue = TRUE)

    res <- insideSize*shares / prices
    if(market) res <- sum(res)
    return(res )


  })

#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcQuantities",
  signature= "CES",
  definition=function(object,preMerger=TRUE, market=FALSE){

    mktSizeRev <- object@mktSize
    quantOutPre <- mktSizeRev / object@priceOutside
    
    mktSizeQuant <- sum(calcShares(object, preMerger= FALSE, revenue = TRUE) * mktSizeRev / object@pricePre)
    mkSizeQuant <- sum(mktSizeQuant,quantOutPre)
    

    shares <- calcShares(object, preMerger= preMerger, revenue = FALSE)

    if(market) shares <- sum(shares,na.rm=TRUE)

    return(mktSizeQuant*shares )


  })

## compute product revenues
#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcRevenues",
  signature= "Bertrand",
  definition=function(object,preMerger=TRUE, market = FALSE){


    if( preMerger) { prices <- object@pricePre}
    else{prices <- object@pricePost}

    quantities <- calcQuantities(object, preMerger)

    res <- prices * quantities
    if(market){return( sum(res) )}
    else{return(res)}



  })

## compute product revenues
#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcRevenues",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE, market = FALSE){


    if( preMerger) { quantities <- object@quantityPre}
    else{quantities <- object@quantityPost}

    prices <- calcPrices(object, preMerger)

    res <- t(prices * t(quantities))

    if(market){return( sum(res) )}
    else{return(res)}



  })


## compute product revenues
#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcRevenues",
  signature= "VertBargBertLogit",
  definition=function(object,preMerger=TRUE, market = FALSE){
    
    result <- calcRevenues(object@down, preMerger=preMerger, market=market)
    return(result)
    
  })

#'@rdname Output-Methods
#'@export

setMethod(
  f= "calcQuantities",
  signature= "AIDS",
  definition=function(object,preMerger=TRUE, market=FALSE){

    if(preMerger){ prices <- object@pricePre}
    else{          prices <- object@pricePost}

    insideSize <- object@insideSize

    shares <- calcShares(object, preMerger = preMerger, revenue=TRUE)

    if(!preMerger){

      # source: epstein and rubinfeld Assumption C (appendix)
      insideSizeDelta <- sum( shares * object@priceDelta) * (object@mktElast + 1)

      insideSize <- insideSize*( 1 + insideSizeDelta)

    }

    shares <- calcShares(object, preMerger= preMerger, revenue = TRUE)

    res <- insideSize*shares / prices
    if(market) res <- sum(res)
    return(res )


  })

## compute product revenues
#'@rdname Output-Methods
#'@export
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

#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcShares",
  signature= "Auction2ndCap",
  definition=function(object,preMerger=TRUE,exAnte=TRUE){


    ownerPre  <- object@ownerPre
    ownerPost  <- object@ownerPost

    shareInside <- cdfG(object,preMerger=preMerger)

    if(preMerger){
      capacities <- object@capacities/sum(object@capacities)
      names(capacities) <- object@labels
      if(exAnte){return((shareInside*capacities))}
      else{return(capacities)}
    }

    else{

      capacities <- object@capacities*(1+object@mcDelta)
      result <- rep(NA, length(object@ownerPre))
      names(result) <- object@labels
      result[object@ownerPre == ownerPost] <- tapply(capacities,ownerPost,sum)

      if(exAnte){return((shareInside*result))}
      else{return(result)}

    }


  }
)

#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcShares",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE,revenue=FALSE){

    if(preMerger) quantities <- object@quantityPre
    else{ quantities <- object@quantityPost}

    if (revenue){
      if(preMerger){ prices <- object@pricePre}
      else{          prices <- object@pricePost}

      totrev <- rowSums(prices*t(quantities), na.rm = TRUE)
      return(t(prices*t(quantities)/totrev))
    }

    else{
      totquant <- colSums(quantities,na.rm=TRUE)
      return(t(t(quantities)/totquant))}
  }
)

#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcShares",
  signature= "Linear",
  definition=function(object,preMerger=TRUE,revenue=FALSE){

    quantities <- calcQuantities(object,preMerger)

    if (revenue){
      if(preMerger){ prices <- object@pricePre}
      else{          prices <- object@pricePost}

      return(prices*quantities/sum(prices*quantities))
    }

    else{return(quantities/sum(quantities))}
  }
)




#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcQuantities",
  signature= "VertBargBertLogit",
  definition=function(object,preMerger=TRUE, market=FALSE){
    
    down <- object@down
    
    result <- calcQuantities(down, preMerger=preMerger,market=market)
    
    return(result)
  })



#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcShares",
  signature= "VertBargBertLogit",
  definition=function(object,preMerger=TRUE,revenue=FALSE){
    
    down <- object@down 
   
    result <- calcShares(down, preMerger=preMerger,revenue=revenue)
    return(result)
    
    }
)



#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcShares",
  signature= "VertBarg2ndLogit",
  definition=function(object,preMerger=TRUE,revenue=FALSE){
  
    
    down <- object@down 
    up <- object@up
    
    priceOutside <- down@priceOutside
    
    if(preMerger){ upPrice <- up@pricePre}
    else{ upPrice <- up@pricePost}
    
   
    meanval <- down@slopes$meanval
    alpha <- down@slopes$alpha
    
    down@slopes$meanval <- meanval + alpha *(upPrice - priceOutside) 
    result <- calcShares(down, preMerger=preMerger,revenue=revenue)
    return(result)
    
    
  }
)

#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcShares",
  signature= "Logit",
  definition=function(object,preMerger=TRUE,revenue=FALSE){


    if(preMerger){ prices <- object@pricePre}
    else{          prices <- object@pricePost}


    alpha    <- object@slopes$alpha
    meanval  <- object@slopes$meanval

    #outVal <- ifelse(object@shareInside<1, 1, 0)
    outVal <- ifelse(is.na(object@normIndex), 1, 0)

    shares <- exp(meanval + alpha*(prices - object@priceOutside))
    shares <- shares/(outVal + sum(shares,na.rm=TRUE))

    if(revenue){shares <- prices*shares/sum(prices*shares,object@priceOutside*(1-sum(shares,na.rm=TRUE)),na.rm=TRUE)}

    names(shares) <- object@labels

    return(shares)

  }
)

#'@rdname Output-Methods
#'@export
setMethod(
f= "calcShares",
signature= "AIDS",
definition=function(object,preMerger=TRUE,revenue=TRUE){

  if(!revenue &&
     any(is.na(object@prices))
  ){
    warning("'prices' contains missing values. Some results are missing")
  }

  prices <- calcPrices(object,preMerger)
  shares <- object@shares


  if(!preMerger){
    shares <-  shares + as.vector(object@slopes %*% log(object@priceDelta + 1))
  }

  if(!revenue){shares <- (shares/prices)/sum(shares/prices)}

  names(shares) <- object@labels
  return(shares)
}
)

#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcShares",
  signature= "LogitNests",
  definition=function(object,preMerger=TRUE,revenue=FALSE){

    if(preMerger){ prices <- object@pricePre}
    else{          prices <- object@pricePost}

    nests     <- object@nests
    alpha     <- object@slopes$alpha
    sigma     <- object@slopes$sigma
    meanval   <- object@slopes$meanval
    isOutside <- as.numeric(object@shareInside < 1)
    outVal    <- ifelse(object@shareInside<1, exp(alpha*object@priceOutside), 0)

    sharesIn     <- exp((meanval+alpha*prices)/sigma[nests])

    inclusiveValue <- log(tapply(sharesIn,nests,sum,na.rm=TRUE))
    sharesAcross <-   exp(sigma*inclusiveValue)
    sharesAcross <- sharesAcross/(outVal + sum(sharesAcross,na.rm=TRUE))


    sharesIn     <- sharesIn/exp(inclusiveValue)[nests]


    shares       <- sharesIn * sharesAcross[nests]

    if(revenue){shares <- prices*shares/sum(prices*shares,object@priceOutside*(1-sum(shares,na.rm=TRUE)),na.rm=TRUE)}

    names(shares) <- object@labels

    return(as.vector(shares))

  }
)

#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcShares",
  signature= "Auction2ndLogit",
  definition=function(object,preMerger=TRUE,revenue=FALSE){

    nprods <- length(object@shares)

    if(preMerger){ subset <- rep(TRUE,nprods) }
    else{subset <- object@subset}


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
    shares[!subset] <- NA
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



#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcShares",
  signature= "Auction2ndLogitNests",
  definition=function(object,preMerger=TRUE,revenue=FALSE){
    
    nprods <- length(object@shares)
    
    if(preMerger){ subset <- rep(TRUE,nprods) }
    else{subset <- object@subset}
    
    
    idx <- object@normIndex
    alpha    <- object@slopes$alpha
    meanval  <- object@slopes$meanval
    sigma    <- object@slopes$sigma
    
    nests <- object@nests
    
    
    
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
    
    
    
    
    sharesBetween <- sharesWithin <- ifelse(subset,exp(meanval/sigma[nests]),NA)
    sharesBetween <-  as.vector(tapply(sharesBetween,nests,sum,na.rm=TRUE))
    sharesWithin <- sharesWithin/sharesBetween[nests]
    
    
    sharesBetween <- sharesBetween^sigma
    sharesBetween <- sharesBetween/(exp(alpha*mcDeltaOut) + sum(sharesBetween,na.rm=TRUE))
    sharesBetween <- sharesBetween[nests]
    
    shares <- sharesWithin*sharesBetween
    
    
    if(revenue){
      res <- rep(NA,nprods)
      res[subset] <- prices[subset]*shares[subset]/sum(prices[subset]*shares[subset],mcDeltaOut*(1-sum(shares[subset])))
      shares <- res
    }
    
    names(shares) <- object@labels
    
    return(shares)
    
    
    
  }
)

#'@rdname Output-Methods
#'@export
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

#'@rdname Output-Methods
#'@export
setMethod(
  f= "calcShares",
  signature= "CESNests",
  definition=function(object,preMerger=TRUE,revenue=FALSE){



    if(preMerger){ prices <- object@pricePre}
    else{          prices <- object@pricePost}


    nests     <- object@nests
    gamma    <- object@slopes$gamma
    sigma    <- object@slopes$sigma
    meanval  <- object@slopes$meanval

    outVal <- ifelse(sum(object@shares)<1, object@priceOutside^(1-gamma), 0)

    sharesIn     <- meanval*prices^(1-sigma[nests])
    sharesAcross <- tapply(sharesIn,nests,sum,na.rm=TRUE)
    sharesIn     <- sharesIn / sharesAcross[nests]
    sharesAcross <- sharesAcross^((1-gamma)/(1-sigma))
    sharesAcross <- sharesAcross / (sum(sharesAcross,na.rm=TRUE) + outVal)

    shares       <- sharesIn * sharesAcross[nests]

    ##transform revenue shares to quantity shares
    if(!revenue){shares <- (shares/prices)/sum((1-sum(shares,na.rm=TRUE))/object@priceOutside,shares/prices,na.rm=TRUE)}

    names(shares) <- object@labels

    return(as.vector(shares))

  }
)


