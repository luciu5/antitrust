#' @title Methods For Calculating Price Delta
#' @name PriceDelta-Methods
#' @docType methods

#' @aliases calcPriceDelta
#' calcPriceDelta,ANY-method
#' calcPriceDelta,Antitrust-method
#' calcPriceDelta,AIDS-method
#' calcPriceDelta,Auction2ndLogit-method
#' calcPriceDelta,Cournot-method
#' calcPriceDelta,VertBargBertLogit-method
#'
#' @description For Antitrust, the method computes equilibrium price changes
#' due to a merger assuming that firms are playing a
#' Nash-Bertrand or Cournot game. This is a wrapper method for computing
#' the difference between pre- and post-merger equilbrium prices.
#' @description For AIDS, the method computes equilibrium price changes
#' due to a merger assuming that firms are playing a
#' Nash-Bertrand or Cournot game and LA-AIDS. This method calls a non-linear
#' equations solver to find a sequence of price changes that satisfy
#' the Bertrand FOCs.
#' @param object An instance of one of the classes listed above.
#' @param levels If TRUE, report results in levels. If FALSE, report results in percents. Default is FALSE.
#' @param market If TRUE, calculates (post-merger) share-weighted average of metric. Default is FALSE.
#' @param party If TRUE, calculates (post-merger) share-weighted average of metric for merging parties. Default is FALSE.
#' @param index If "paasche",calculates market-wide price changes using post-merger predicted shares. If  "laspeyres", 
#' calculates price index using pre-merger shares. Default is "paasche".
#' @param isMax If TRUE, uses numerical derivatives to determine if
#' equilibrium price vector is a local maximum. Default is FALSE.
#' @param subset A vector of length k where each element equals TRUE if
#' the product indexed by that element should be included in the
#' post-merger simulation and FALSE if it should be excluded.Default is a
#' length k vector of TRUE.
#' @param exAnte If \sQuote{exAnte} equals TRUE then the
#' \emph{ex ante} expected result for each firm is produced, while FALSE produces the
#' expected result conditional on each firm winning the auction. Default is FALSE, unless \sQuote{market} is TRUE.
#' @param ... Additional values that may be used to change the default values of the non-linear
#' equation solver.
#'
#' @include ParamsMethods.R
#' @keywords methods
NULL

setGeneric (
  name= "calcPriceDelta",
  def=function(object,levels = FALSE, market = FALSE, party = FALSE, isMax=FALSE,index=c("paasche","laspeyres"),...){standardGeneric("calcPriceDelta")},signature="object")

## Method to compute price changes
#'@rdname PriceDelta-Methods
#'@export
setMethod(
  f= "calcPriceDelta",
  signature= "Antitrust",
  definition=function(object, levels = FALSE, market = FALSE, party = FALSE, index=c("paasche","laspeyres"), ...  ){

    index <- match.arg(index)
    
    pricePre  <- object@pricePre
    pricePost <- object@pricePost

    
    if(levels){priceDelta <- pricePost - pricePre}
    else{priceDelta <- pricePost/pricePre - 1}
    #names(priceDelta) <- object@labels

    if(market || party){
      
      sharesPre <- calcShares(object, preMerger=TRUE,revenue=FALSE,...)
      sharesPre <- sharesPre/sum(sharesPre,na.rm=TRUE)
      
      sharesPost <- calcShares(object, preMerger=FALSE,revenue=FALSE,...)
      sharesPost <- sharesPost/sum(sharesPost,revenue=FALSE,na.rm=TRUE)
      
      
      
      if(party){
        ownerPre <- ownerToMatrix(object,preMerger=TRUE)
        ownerPost <-ownerToMatrix(object,preMerger=FALSE)
        
        isParty <- rowSums( abs(ownerPost - ownerPre))>0
        sharesPre <- sharesPre[isParty]
        sharesPost <- sharesPost[isParty]
        pricePre <- pricePre[isParty]
        pricePost <- pricePost[isParty]
        
      }
      
      if(index=="paasche")  priceDelta <- sum(sharesPost*pricePost)/sum(sharesPost*pricePre) - 1
      else if (index=="laspeyres")  priceDelta <- sum(sharesPre*pricePost)/sum(sharesPre*pricePre) - 1
      
      
      
      
    }

    return(priceDelta)

  }
)

#'@rdname PriceDelta-Methods
#'@export
setMethod(
  f= "calcPriceDelta",
  signature= "Cournot",
  definition=function(object, levels = FALSE, market=TRUE, ...  ){
    
    callNextMethod()
    
  }
)

#'@rdname PriceDelta-Methods
#'@export
setMethod(
  f= "calcPriceDelta",
  signature= "VertBargBertLogit",
  definition=function(object, levels = FALSE, market = FALSE, ...  ){
    
    up <- object@up
    down <- object@down
    
    chain_level <- object@chain_level

    marginsPre <- calcMargins(object,preMerger=TRUE,level=TRUE)
    marginsPost <- calcMargins(object,preMerger=FALSE,level=TRUE)
    
    sharesPre <- calcShares(object, preMerger=TRUE,revenue=FALSE)
    sharesPost <- calcShares(object, preMerger=FALSE,revenue=FALSE)
    
      upMCPre=up@mcPre
      upMCPre=ifelse(is.na(upMCPre),0,upMCPre)
      downMCPre=down@mcPre
      downMCPre=ifelse(is.na(downMCPre),0,downMCPre)
      upPricePre <- up@pricePre  
   
      upMCPost=up@mcPost
      upMCPost=ifelse(is.na(upMCPost),0,upMCPost)
      downMCPost=down@mcPost
      downMCPost=ifelse(is.na(downMCPost),0,downMCPost)
      upPricePost <-  up@pricePost
    
        if(!market){
        mcDeltaUp <- upMCPost  - upMCPre
        ## assume 0 marginal cost changes if unkown
        #mcDeltaUp <- ifelse(is.na(mcDeltaUp),0,mcDeltaUp)
        
        mcDeltaDown <- (downMCPost- downMCPre )
        #mcDeltaDown <- ifelse(is.na(mcDeltaDown),0,mcDeltaDown)
        mcDeltaDown <- mcDeltaDown + upPricePost - upPricePre
        
   
       
    
      upDelta <- marginsPost$up - marginsPre$up + mcDeltaUp
      downDelta <- marginsPost$down - marginsPre$down + mcDeltaDown
      
      if(chain_level =="retailer") upDelta <-rep(0,length(upDelta)) 
      else if(chain_level =="wholesaler") downDelta <-rep(0,length(downDelta))
      upPricePre <- up@pricePre
      downPricePre <- down@pricePre
    }
      
      
    else{
      
      
      mcDeltaUp <- upMCPost *sharesPost - upMCPre*sharesPre
      mcDeltaDown <- (downMCPost+upPricePost)*sharesPost - (downMCPre + upPricePre)*sharesPre
      
      ## assume 0 marginal cost changes if unkown
      #mcDeltaUp <- ifelse(is.na(mcDeltaUp),0,mcDeltaUp)
      #mcDeltaDown <- ifelse(is.na(mcDeltaDown),0,mcDeltaDown)
      
      upDelta <- marginsPost$up*sharesPost - marginsPre$up*sharesPre + mcDeltaUp
      downDelta <- marginsPost$down*sharesPost - marginsPre$down*sharesPre + mcDeltaDown
      
      upPricePre <- up@pricePre*sharesPre
      downPricePre <- down@pricePre*sharesPre
      
      
      upDelta <- sum(upDelta,na.rm=TRUE)
      downDelta <- sum(downDelta,na.rm=TRUE)
      
      if(chain_level =="retailer") upDelta <-0
      else if(chain_level =="wholesaler") downDelta <- 0
      
      upPricePre <- sum(upPricePre,na.rm=TRUE)
      downPricePre <- sum(downPricePre,na.rm=TRUE)
    }
   
    if(!levels){
      upDelta <- upDelta/upPricePre
      downDelta <- downDelta/downPricePre
    }
    priceDelta <- list(up = upDelta,
                       down= downDelta)
    
    return(priceDelta)
    
  }
)


#'@rdname PriceDelta-Methods
#'@export
setMethod(
  f= "calcPriceDelta",
  signature= "AIDS",
  definition=function(object,levels=FALSE,market=FALSE,party=FALSE,isMax=FALSE, index=c("paasche","laspeyres"),subset,...){

    index <- match.arg(index)
    
  
    priceDelta <- object@priceDelta 

      
  if (market || party){
    
    if(!is.null(object@pricePre) && all(!is.na(object@pricePre))){
      result <- callNextMethod()
      return(result)
    }
    
    if(party){
      if(index=="paasche") shares <-  calcShares(object, preMerger = FALSE)
      else{shares <-  calcShares(object, preMerger = TRUE)}
      
      
      isParty <- rowSums( abs(object@ownerPost - object@ownerPre))>0
      priceDelta <- priceDelta[isParty]
      shares <- shares[isParty]
    }
    

     
      return(sum( priceDelta * shares,na.rm=TRUE))
           }

    ownerPost <- object@ownerPost

    nprods <- length(object@shares)
    if(missing(subset)){subset <- rep(TRUE,nprods)}

    if(!is.logical(subset) || length(subset) != nprods ){stop("'subset' must be a logical vector the same length as 'shares'")}


    ##Define system of FOC as a function of priceDelta
    FOC <- function(priceDelta){

      object@priceDelta <- exp(priceDelta)-1

      sharePost <-  calcShares(object,FALSE)
      elastPost <-  t(elast(object,FALSE))
      marginPost <- calcMargins(object,FALSE)


      thisFOC <- sharePost*diag(ownerPost) + as.vector((elastPost*ownerPost) %*% (sharePost*marginPost))
      thisFOC[!subset] <- sharePost[!subset]
      return(thisFOC)

    }




    ## Find price changes that set FOCs equal to 0
    minResult <- BBsolve(object@priceStart,FOC,quiet=TRUE,control=object@control.equ)


    if(minResult$convergence != 0){warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBsolve' reports: '",minResult$message,"'")}

    if(isMax){

      hess <- genD(FOC,minResult$par) #compute the numerical approximation of the FOC hessian at optimium
      hess <- hess$D[,1:hess$p]


      if(any(eigen(hess)$values>0)){warning("Hessian of first-order conditions is not positive definite. Price vector may not maximize profits. Consider rerunning 'calcPrices' using different starting values")}
    }

    deltaPrice <- exp(minResult$par)-1
    names(deltaPrice) <- object@labels

    if(levels){deltaPrice <- calcPrices(object,FALSE) - calcPrices(object,TRUE)}

    if(market){
      sharePost <-  calcShares(object,FALSE,...)
      sharePost <- sharePost/sum(sharePost, na.rm=TRUE)

      deltaPrice <- sum(deltaPrice*sharePost,na.rm=TRUE)

    }



    return(deltaPrice)
  }
)


#'@rdname PriceDelta-Methods
#'@export
setMethod(
  f= "calcPriceDelta",
  signature= "Auction2ndLogit",
  definition=function(object,levels=TRUE, market=FALSE,party=FALSE,exAnte=ifelse(market,TRUE,FALSE),...){

    output <- object@output
    
    if(all(!is.na(object@pricePre)) || !levels){
      result <- callNextMethod()
      return(result)
    }
    
    subset <- object@subset

    mcDelta <- object@mcDelta


    
    if(exAnte || market || party){
      sharesPost <- calcShares(object, preMerger=FALSE)
      mcDelta <- mcDelta*sharesPost
    }
   
    result <- calcMargins(object, preMerger=FALSE,exAnte=exAnte) + mcDelta -
      calcMargins(object, preMerger=TRUE,exAnte=exAnte)

    if(!output) result <- -1*result
    
    if(market) result <- sum(result,na.rm=TRUE)
    if(party){
      isParty <- rowSums( abs(object@ownerPost - object@ownerPre))>0
      result <- sum(result[isParty],na.rm=TRUE)
    }
    
 

    if(!market && !party) names(result) <- object@labels
    
    return(result)
  }
)
