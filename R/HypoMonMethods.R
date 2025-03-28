#' @title Methods For Implementing The Hypothetical Monopolist Test
#' @name defineMarketTools-methods
#' @docType methods
#' @aliases HypoMonTest
#' HypoMonTest,ANY-method
#' HypoMonTest,Bertrand-method
#' HypoMonTest,VertBargBertLogit-method
#' calcPricesHypoMon
#' calcPricesHypoMon,ANY-method
#' calcPricesHypoMon,AIDS-method
#' calcPricesHypoMon,Linear-method
#' calcPricesHypoMon,LogLin-method
#' calcPricesHypoMon,Logit-method
#' calcPricesHypoMon,LogitCap-method
#' calcPricesHypoMon,Auction2ndLogit-method
#' calcPricesHypoMon,Cournot-method
#' calcPriceDeltaHypoMon
#' calcPriceDeltaHypoMon,ANY-method
#' calcPriceDeltaHypoMon,AIDS-method
#' calcPriceDeltaHypoMon,Bertrand-method
#' calcPriceDeltaHypoMon,Cournot-method
#' diversionHypoMon
#' diversionHypoMon,ANY-method
#' diversionHypoMon,AIDS-method
#' diversionHypoMon,Bertrand-method
#'
#' @description An Implementation of the Hypothetical Monopolist Test described in either the 2023 or 2010 Merger Guidelines.
#' @description \code{\link{HypoMonTest}} implements the Hypothetical Monopolist Test for a given \sQuote{ssnip}.
#' @description \code{calcPricesHypoMon} computes prices for a subset of firms under the control of a hypothetical monopolist
#' under the specified demand function or auction.
#' @description \code{\link{diversionHypoMon}} calculates the matrix of revenue diversions between all products included in the
#' merger simulation, \emph{irrespective of whether or not they are also included in \sQuote{prodIndex}}.
#' @description \code{\link{calcPriceDeltaHypoMon}} computes the proportional difference in product prices between the
#' prices of products in \sQuote{prodIndex} (i.e. prices set by the
#' Hypothetical Monopolist) and prices set in the pre-merger equilibrium.
#' \sQuote{...} may be used to pass arguments to the optimizer.
#'
#' @param object An instance of one of the classes listed above.
#' @param prodIndex A vector of product indices that are to be placed under the control of the Hypothetical Monopolist.
#' @param plantIndex A vector of plant indices that are to be placed under the control of the Hypothetical Monopolist (Cournot).
#' @param ssnip A number between 0 and 1 that equals the threshold for a \dQuote{Small but Significant and
#' Non-transitory Increase in Price} (SSNIP). Default is .05, or 5\%.
#' @param hmg Either the string "2023" (default) or "2010". Implements the Hypothetical Monopolist Test described in either the 2010 or 2023 Merger Guidelines.
#' @param ... Pass options to the optimizer used to solve for equilibrium prices.
#'
#' @details
#' Let k denote the number of products produced by all firms playing the Bertrand pricing game above.
#'
#' @details
#' \code{HypoMonTest} is an implementation of the 2023 Merger Guidelines Hypothetical Monopolist Test
#' on the products indexed by \sQuote{prodIndex} for a \sQuote{ssnip}. The
#' Hypothetical Monopolist Test described in the 2023 Merger Guidelines determines whether a profit-maximizing
#' Hypothetical Monopolist who controls the products indexed by
#' \sQuote{prodIndex} would increase the price of at least one of the products in \sQuote{prodIndex} by a
#' small, significant, and non-transitory amount (i.e. impose a SSNIP). Setting \sQuote{hmg} to "2010" implements the
#' Hypothetical Monopolist Test described in the 2010 Merger Guidelines, which requires the Hypothetical Monopolist to 
#' increase the price of one of the merging parties' products in \sQuote{prodIndex} by a SSNIP. 
#'
#' @details
#' \code{calcPriceDeltaHypoMon} calculates the price changes relative to (predicted) pre-merger prices that a
#' Hypothetical Monopolist would impose on the products indexed by \sQuote{prodIndex}, holding the prices of products not
#' controlled by the Hypothetical Monopolist fixed at pre-merger levels. With the exception of \sQuote{AIDS}, the
#' \code{calcPriceDeltaHypoMon} for all the classes listed above calls \code{calcPricesHypoMon} to compute price
#' levels. \code{calcPriceDeltaHypoMon} is in turn called by \code{HypoMonTest}.
#'
#' @details
#' \code{diversionHypoMon} calculates the matrix of revenue diversions between all products included in the merger simulation,
#' \emph{irrespective} of whether or not they are also included in
#' \sQuote{prodIndex}. This matrix is useful for diagnosing whether or not a
#' product not included in \sQuote{prodIndex} may have a higher revenue diversion
#' either to or from a product included in \sQuote{prodIndex}. Note that the \sQuote{AIDS}
#' \code{diversionHypoMon} method does not contain the \sQuote{prodIndex}
#' argument, as AIDS revenue diversions are only a function of demand parameters.
#'
#' @return
#' \code{HypoMonTest} returns TRUE if a profit-maximizing Hypothetical Monopolist who controls the products indexed by
#' \sQuote{prodIndex} would increase the price of at least one of the merging
#' parties' products in \sQuote{prodIndex} by a \sQuote{ssnip}, and
#' FALSE otherwise. \code{HypoMonTest} returns an error if \sQuote{prodIndex}
#' does not contain at least one of the merging parties products.
#'
#' @return
#' \code{calcPriceDeltaHypoMon} returns a vector of proportional price changes for
#' all products placed under the control of the Hypothetical
#' Monopolist (i.e. all products indexed by \sQuote{prodIndex}).
#' @return \code{calcPricesHypoMon} is identical, but for price levels.
#' @return \code{diversionHypoMon} returns a k x k matrix of diversions,
#' where element i,j is the diversion from product i to product j.
#'
#' @references U.S. Department of Justice and Federal Trade Commission,
#' \emph{Horizontal Merger Guidelines}. Washington DC: U.S. Department of Justice, 2010.
#' \url{https://www.justice.gov/atr/horizontal-merger-guidelines-08192010} (accessed May 5, 2021).
#'
#' @include PlotMethods.R
#' @keywords methods
NULL

setGeneric (
  name= "HypoMonTest",
  def=function(object,...){standardGeneric("HypoMonTest")}
)

setGeneric (
  name= "calcPricesHypoMon",
  def=function(object,...){standardGeneric("calcPricesHypoMon")}
)


setGeneric (
  name= "diversionHypoMon",
  def=function(object,...){standardGeneric("diversionHypoMon")}
)


setGeneric (
  name= "calcPriceDeltaHypoMon",
  def=function(object,...){standardGeneric("calcPriceDeltaHypoMon")}
)


#'@rdname defineMarketTools-methods
#'@export
setMethod(
  f= "HypoMonTest",
  signature= "Bertrand",
  definition=function(object,prodIndex,ssnip=ifelse(object@output,.05,-.05),hmg=c("2023","2010"),...){

    hmg <- match.arg(hmg)
    ownerPre <- object@ownerPre
    nprods   <- ncol(ownerPre)
    pricesDelta <- rep(0,nprods)
    output <- object@output

    if(missing(prodIndex) || any(prodIndex>nprods | prodIndex <1 ) ){
      stop("'prodIndex' must be a vector of product indices between 1 and ",nprods)
    }

    if(length(ssnip)>1 || abs(ssnip)>1 ){stop("absolute value of 'ssnip' must be less than 1")}

    if(hmg=="2010"){ isParty <- rowSums( abs(object@ownerPost - ownerPre) )>0} #identify which products belong to the merging parties}
    else if(hmg=="2023"){isParty <- rep(TRUE,nprods)}

    if(identical(length(intersect(which(isParty),prodIndex)),0)){
      stop("'prodIndex' does not contain any of the merging parties' products. Add at least one of the following indices: ",
           paste(which(isParty),collapse=","))
    }



    pricesDelta[prodIndex] <-  calcPriceDeltaHypoMon(object,prodIndex,...)


    
    result <- ifelse(output,max(pricesDelta[isParty]) > ssnip,
                     min(pricesDelta[isParty]) < ssnip)
                     

    return( result)
  }

)


setMethod(
  f= "HypoMonTest",
  signature= "VertBargBertLogit",
  
  definition=function(object,prodIndex,ssnip,hmg=c("2023","2010"),...){
    
    if(missing(ssnip)){
      ssnip <- ifelse(object@down@output,.05,-.05)}
    
    hmg=match.arg(hmg)
    down <- object@down
    down@ownerPre <- ownerToMatrix(down,preMerger=TRUE)
    down@ownerPost <- ownerToMatrix(down,preMerger=FALSE)
    down@pricePre <- calcPrices(down,preMerger=TRUE)
    
    HypoMonTest(object=down,prodIndex=prodIndex,ssnip=ssnip,hmg=hmg,...)
    
  })

#'@rdname defineMarketTools-methods
#'@export
setMethod(
  f= "HypoMonTest",
  signature= "Cournot",
  definition=function(object,plantIndex,prodIndex,ssnip=ifelse(object@output,.05,-.05),hmg=c("2023","2010"),...){
    
    hmg <- match.arg(hmg)
    ownerPre <- object@ownerPre
    nprods   <- ncol(ownerPre)
    nplants <- nrow(ownerPre)
    output <- object@output

    
    if(missing(plantIndex) || any(plantIndex>nplants | plantIndex <1 ) ){
      stop("'plantIndex' must be a vector of plant indices between 1 and ",nplants)
    }
    
    if(missing(prodIndex) || length(prodIndex) != 1 || any(prodIndex>nprods | prodIndex <1 ) ){
      stop("'prodIndex' must be between between 1 and ",nprods)
    }
    if(length(ssnip)>1 || abs(ssnip)>1 ){stop("absolute value of 'ssnip' must be a number between 0 and 1")}
    
    if(hmg=="2010"){     isParty <- rowSums( abs(object@ownerPost - object@ownerPre) )>0} #identify which plants belong to the merging parties
    else if(hmg=="2023"){isParty <- rep(TRUE,nplants)}
    
    if(identical(length(intersect(which(isParty),plantIndex)),0)){
      stop("'plantIndex' does not contain any of the merging parties' plants. Add at least one of the following indices: ",
           paste(which(isParty),collapse=","))
    }
    
    
    
    pricesDelta <-  calcPriceDeltaHypoMon(object,prodIndex=prodIndex,plantIndex=plantIndex,...)
    
    

    result <- ifelse(output,max(pricesDelta) > ssnip,
                     min(pricesDelta) < ssnip)
    
    return( result)
  }
  
)


#'@rdname defineMarketTools-methods
#'@export
setMethod(
  f= "calcPricesHypoMon",
  signature= "Cournot",
  definition=function(object,plantIndex,prodIndex){
    
    
    nhypoplants <- length(plantIndex)
    
    nprods <- length(prodIndex)
    
    intercept <- object@intercepts[prodIndex]
    slopes <- object@slopes[prodIndex]
    quantityPre <- as.vector(object@quantityPre[plantIndex,prodIndex])
    quantFixed <- colSums(object@quantityPre[-plantIndex,prodIndex,drop=FALSE])
    demand <- object@demand[prodIndex]
    

    ## how to deal with multiple products?
    #stop("A work in progress!! May not properly handle multiple products")

    calcMonopolySurplus <- function(quantCand){

      
      quantCand <- matrix(quantCand,ncol=nprods, nrow=nhypoplants)
      object@quantityPre[plantIndex,prodIndex] <- quantCand
      mktQuant <- quantFixed + colSums(quantCand, na.rm = TRUE)

      priceCand <- ifelse(demand == "linear",
                          intercept + slopes * mktQuant,
                          exp(intercept)*mktQuant^slopes)

      vcCand <- calcVC(object, preMerger=TRUE)
      vcCand <- vcCand[plantIndex]

      revCand <-  colSums(priceCand*t(quantCand), na.rm=TRUE)


      surplus <- sum(revCand - vcCand, na.rm =TRUE)

      return(sum(surplus))
    }

    if( nhypoplants > 1){

      maxResult <- optim(quantityPre,
                         calcMonopolySurplus,
                         method="L-BFGS-B",
                         lower = rep(0,nhypoplants),
                         control = list(fnscale=-1)
      )

      quantitiesHM <- maxResult$par
    }


    else{

      upperB <- sum(quantityPre,na.rm=TRUE)
      maxResult <- optimize(calcMonopolySurplus,c(0, upperB),maximum = TRUE)
      quantitiesHM <- maxResult$maximum
    }


    quantitiesHM <- matrix(quantitiesHM, nrow=nhypoplants,ncol=nprods)
    mktQuant <- quantFixed + colSums(quantitiesHM)
    priceHM <- ifelse(demand == "linear",
                        intercept + slopes * mktQuant,
                        exp(intercept)*mktQuant^slopes)
    
    
    names(priceHM) <- object@labels[[2]][prodIndex]
    return(priceHM)


  })


#'@rdname defineMarketTools-methods
#'@export
setMethod(
  f= "calcPricesHypoMon",
  signature= "Linear",
  definition=function(object,prodIndex){

    nprods <- length(prodIndex)
    intercept <- object@intercepts
    slopes <- object@slopes
    mc <- object@mcPre[prodIndex]
    pricePre <- object@pricePre

    calcMonopolySurplus <- function(priceCand){


      pricePre[prodIndex] <- priceCand
      quantityCand <- intercept + as.vector(slopes %*% pricePre)

      surplus <- (priceCand-mc)*quantityCand[prodIndex]

      return(sum(surplus))
    }

    ##Find starting value that always meets boundary conditions
    ##Note: if nprods=1, need to use a more accurate optimizer.

    if(nprods > 1){

      if(det(slopes)!=0){startParm <- as.vector(solve(slopes) %*% (1 - intercept ))}
      else{startParm <- rep(0,nprods)}


      priceConstr <- pricePre
      priceConstr[prodIndex] <- 0

      maxResult <- constrOptim(startParm[prodIndex],calcMonopolySurplus,
                               grad=NULL,
                               ui=slopes[prodIndex,prodIndex],
                               ci=-intercept[prodIndex] - as.vector(slopes %*% priceConstr)[prodIndex],
                               control=list(fnscale=-1))

      pricesHM <- maxResult$par
    }


    else{

      upperB <- -(intercept[prodIndex] + sum(pricePre[-prodIndex]*slopes[prodIndex,-prodIndex]))/slopes[prodIndex,prodIndex]

      maxResult <- optimize(calcMonopolySurplus,c(0,upperB),maximum = TRUE)
      pricesHM <- maxResult$maximum
    }

    #priceDelta <- pricesHM/pricePre[prodIndex] - 1
    #names(priceDelta) <- object@labels[prodIndex]
    names(pricesHM) <- object@labels[prodIndex]

    return(pricesHM)


  })


#'@rdname defineMarketTools-methods
#'@export
setMethod(
  f= "calcPricesHypoMon",
  signature= "Logit",
  definition=function(object,prodIndex){


    mc       <- object@mcPre[prodIndex]
    pricePre <- object@pricePre
    output <- object@output
    outSign <- ifelse(output,1,-1)

    calcMonopolySurplus <- function(priceCand){

      pricePre[prodIndex] <- priceCand #keep prices of products not included in HM fixed at premerger levels
      object@pricePre     <- pricePre
      sharesCand          <- calcShares(object,TRUE,revenue=FALSE)

      surplus             <- outSign*(priceCand-mc)*sharesCand[prodIndex]

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


#'@rdname defineMarketTools-methods
#'@export
setMethod(
  f= "calcPricesHypoMon",
  signature= "LogLin",
  definition=function(object,prodIndex){


    mc <- object@mcPre[prodIndex]
    pricePre <- object@pricePre
    output <- object@output
    outSign <- ifelse(output,1,-1)

    calcMonopolySurplus <- function(priceCand){

      pricePre[prodIndex] <- priceCand
      object@pricePre     <- pricePre
      quantityCand        <- calcQuantities(object,TRUE)


      surplus             <- outSign*(priceCand-mc)*quantityCand[prodIndex]


      return(sum(surplus))
    }


    minResult <- optim(object@prices[prodIndex],calcMonopolySurplus,
                       method = "L-BFGS-B",lower = 0,
                       control=list(fnscale=-1))

    pricesHM <- minResult$par

    #priceDelta <- pricesHM/pricePre[prodIndex] - 1
    #names(priceDelta) <- object@labels[prodIndex]
    names(pricesHM) <- object@labels[prodIndex]

    return(pricesHM)

  })


#'@rdname defineMarketTools-methods
#'@export
setMethod(
  f= "calcPricesHypoMon",
  signature= "AIDS",
  definition=function(object,prodIndex,...){


    priceDeltaHM <- calcPriceDeltaHypoMon(object,prodIndex,...)

    prices <- object@prices[prodIndex] * (1 + priceDeltaHM)


    return(prices)

  })


#'@rdname defineMarketTools-methods
#'@export
setMethod(
  f= "calcPricesHypoMon",
  signature= "LogitCap",
  definition=function(object,prodIndex,...){


    mc       <- object@mcPre[prodIndex]
    capacities <- object@capacitiesPre[prodIndex]
    pricePre <- object@pricePre

    FOC <- function(priceCand){

      thisPrice <- pricePre
      thisPrice[prodIndex] <- priceCand

      object@pricePre <- thisPrice

      output <- object@output
      outSign <- ifelse(output,1,-1)
      
      margins          <- outSign*(1 - mc/priceCand)
      
      quantities       <- calcQuantities(object,preMerger=TRUE)[prodIndex]
      revenues         <- quantities * priceCand
      elasticities     <- elast(object,preMerger=TRUE)[prodIndex,prodIndex]

      thisFOC <- revenues + as.vector(t(elasticities) %*% (margins * revenues))
      constraint <- ifelse(!is.finite(capacities),0, quantities - capacities)

      measure <- thisFOC + constraint + sqrt(thisFOC^2 + constraint^2)

      return(measure)
    }



    ## Find price changes that set FOCs equal to 0
    minResult <- BBsolve(object@priceStart[prodIndex],FOC,quiet=TRUE,control=object@control.equ,...)

    if(minResult$convergence != 0){warning("'calcPricesHypoMon' nonlinear solver may not have successfully converged. 'BBSolve' reports: '",minResult$message,"'")}


    pricesHM <- minResult$par
    #priceDelta <- pricesHM/pricePre[prodIndex] - 1
    #names(priceDelta) <- object@labels[prodIndex]
    names(priceHM) <- object@labels[prodIndex]

    return(priceHM)

  })


#'@rdname defineMarketTools-methods
#'@export
setMethod(
  f= "calcPricesHypoMon",
  signature= "Auction2ndLogit",
  definition=function(object,prodIndex){


    ownerMon <- object@ownerPre
    ownerMon[prodIndex,] <- 0
    ownerMon[,prodIndex] <- 0
    ownerMon[prodIndex,prodIndex] <- 1

    object@ownerPre <- ownerMon

    pricesHM <- calcPrices(object,preMerger=TRUE)
    pricesHM <- pricesHM[prodIndex]
    names(pricesHM) <- object@labels[prodIndex]

    return(pricesHM)

  })

#'@rdname defineMarketTools-methods
#'@export
setMethod(
  f= "diversionHypoMon",
  signature= "Bertrand",
  definition=function(object,prodIndex,...){

    object@pricePre[prodIndex] <- calcPricesHypoMon(object,prodIndex,...)

    return(diversion(object,preMerger=TRUE,revenue=TRUE))



  }
)


#'@rdname defineMarketTools-methods
#'@export
setMethod(
  f= "diversionHypoMon",
  signature= "AIDS",
  definition=function(object){

    return(diversion(object,revenue=FALSE))

  })


## Use the Hypothetical Monopolist Test to determine whether a candidate market satisfies a SSNIP.
#'@rdname defineMarketTools-methods
#'@export
setMethod(
  f= "calcPriceDeltaHypoMon",
  signature= "Bertrand",
  definition=function(object,prodIndex,...){


    pricesHM <-  calcPricesHypoMon(object,prodIndex,...)

    pricesDelta <- pricesHM/object@pricePre[prodIndex] - 1

    return(pricesDelta)

  })


#'@rdname defineMarketTools-methods
#'@export
setMethod(
  f= "calcPriceDeltaHypoMon",
  signature= "Cournot",
  definition=function(object,prodIndex,plantIndex,...){
    
    
    pricesHM <-  calcPricesHypoMon(object,prodIndex=prodIndex,plantIndex=plantIndex,...)
    
    pricesDelta <- pricesHM/object@pricePre[prodIndex] - 1
    
    return(pricesDelta)
    
  })



#'@rdname defineMarketTools-methods
#'@export
setMethod(
  f= "calcPriceDeltaHypoMon",
  signature= "AIDS",
  definition=function(object,prodIndex,...){

    priceDeltaOld <- object@priceDelta

    ##Define system of FOC as a function of priceDelta
    FOC <- function(priceDelta){

      priceCand <- priceDeltaOld
      priceCand[prodIndex] <- priceDelta
      object@priceDelta <- exp(priceCand)-1

      shareCand <-  calcShares(object,FALSE)
      elastCand <-  elast(object,FALSE)
      marginCand <- calcMargins(object,FALSE)

      elastCand <-   elastCand[prodIndex,prodIndex]
      shareCand <-   shareCand[prodIndex]
      marginCand <-  marginCand[prodIndex]

      thisFOC <- shareCand + as.vector(t(elastCand) %*% (shareCand*marginCand))
      return(thisFOC)

    }



    ## Find price changes that set FOCs equal to 0
    minResult <- BBsolve(object@priceStart[prodIndex],FOC,quiet=TRUE,...)

    if(minResult$convergence != 0){warning("'calcPricesHypoMon' nonlinear solver may not have successfully converged. 'BBsolve' reports: '",minResult$message,"'")}


    deltaPrice <- (exp(minResult$par)-1)

    names(deltaPrice) <- object@labels[prodIndex]

    return(deltaPrice[prodIndex])

  })

