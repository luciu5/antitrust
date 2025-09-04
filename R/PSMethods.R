#' @title Producer Surplus Methods
#' @name PS-methods
#' @docType methods
#' @aliases calcProducerSurplus
#' calcProducerSurplus,ANY-method
#' calcProducerSurplus,Bertrand-method
#' calcProducerSurplus,Cournot-method
#' calcProducerSurplus,VertBargBertLogit-method
#' calcProducerSurplusGrimTrigger
#'
#' @description In the following methods, \code{calcProducerSurplus} computes the expected profits of each supplier
#' with the game depending on the class. The available classes are: Bertrand, Cournot, and Auction2ndCap.
#' @description \code{calcProducerSurplusGrimTrigger} is a method that
#' may be used to explore how a merger affects firms' incentives to collude.
#'
#' @details calcProducerSurplusGrimTrigger calculates \sQuote{preMerger} product
#' producer surplus (as well as other statistics -- see below), under the
#' assumption that firms are playing an N-player iterated Prisoner's Dilemma where in each period a
#' coalition of firms decides whether to \emph{cooperate} with
#' one another by setting the joint surplus maximizing price on some
#' \sQuote{coalition} of their products, or \emph{defect} from the coalition by setting all of their products'
#' prices to optimally undercut the prices of the coalition's products. Moreover, firms are assumed
#' to play Grim Trigger strategies where each firm cooperates in the
#' current period so long as \emph{every} firm in the coalition cooperated last period and
#' defects otherwise. product level \sQuote{discount} rates are then employed to determine
#' whether a firm's discounted surplus from remaining in the coalition are greater than
#' its surplus from optimally undercutting the coalition prices' for one
#' period plus its discounted surplus when all firms set Nash-Bertrand prices in all subsequent periods.
#'
#' @return \code{calcProducerSurplusGrimTrigger} returns a data frame with rows
#' equal to the number of products produced by any firm participating in
#' the coalition and the following 5 columns
#'
#'     \item{Discount:}{The user-supplied discount rate}
#'     \item{Coord:}{Single period producer surplus from coordinating}
#'     \item{Defect:}{Single period producer surplus from defecting}
#'     \item{Punish:}{Single period producer surplus from punishing using Bertrand price}
#'     \item{IC:}{TRUE if the discounted producer surplus from coordinating across all firm products are
#'     greater than the surplus from defecting across all firm products for one period and receiving
#'     discounted Bertrand surplus for all subsequent periods under Grim Trigger.}
#'
#'
#' @param object An instance of one of the classes listed above.
#' @param preMerger If TRUE, returns pre-merger outcome. If FALSE, returns post-merger outcome. Default is TRUE.
#' @param exAnte If \sQuote{exAnte} equals TRUE then the
#' \emph{ex ante} expected result for each firm is produced, while FALSE produces the
#' expected result conditional on each firm winning the auction. Default is TRUE.
#' @param coalition A length c vector of integers indicating the index of the products participating in the coalition.
#' @param discount A length k vector of values between 0 and 1 that
#' represent the product-specific discount rate for all products
#' produced by firms particiapting in the coalition. NAs are allowed.
#' @param isCollusion TRUE recalculates demand and cost parameters under
#' the assumption that the coalition specified in \sQuote{coalition} is
#' operating pre-merger. FALSE (the default) uses demand
#' and cost parameters calculated from the \sQuote{ownerPre} matrix.
#' @param ... Additional argument to pass to calcPrices (for \code{calcProducerSurplusGrimTrigger})
#' @include PriceDeltaMethods.R
#' @keywords methods
NULL


setGeneric (
  name= "calcProducerSurplus",
  def=function(object,...){standardGeneric("calcProducerSurplus")}
)

setGeneric (

  name= "calcProducerSurplusGrimTrigger",
  def=function(object,...){standardGeneric("calcProducerSurplusGrimTrigger")}
)


## compute producer surplus
#'@rdname PS-methods
#'@export
setMethod(
  f= "calcProducerSurplus",
  signature= "Bertrand",
  definition=function(object,preMerger=TRUE){


    margins <- calcMargins(object,preMerger,level=TRUE)


    output <- calcQuantities(object,preMerger)

    if (all(is.na(output))){
      warning("'calcQuantities' yielded all NAs. Using 'calcShares' instead")
      output <- calcShares(object,preMerger,revenue=FALSE)
    }

    ps <- margins * output
    names(ps) <- object@labels

    return(ps)
  }

)
#'@rdname PS-methods
#'@export
setMethod(
  f= "calcProducerSurplus",
  signature= "VertBargBertLogit",
  definition=function(object,preMerger=TRUE){
    
    mktSize <- object@down@mktSize
    
    margins <- calcMargins(object,preMerger=preMerger,level=TRUE)
    
    output <- calcShares(object,preMerger)
    
    if (is.na(mktSize)){
      warning("'calcQuantities' yielded all NAs. Using 'calcShares' instead")
      mktSize <- 1
    }
    
    psup <- margins$up * output * mktSize
    psdown <- margins$down * output * mktSize
    names(psup) <- names(psdown) <-  object@down@labels
    
    return(list(up=psup,down=psdown))
  }
  
)



#'@rdname PS-methods
#'@export
setMethod(
  f= "calcProducerSurplus",
  signature= "Auction2ndCap",
  definition=function(object,preMerger=TRUE,exAnte=TRUE){

    sellerCostBounds <-object@sellerCostBounds
    if(preMerger){r    <- object@reservePre}
    else{r    <- object@reservePost}


    if(preMerger) { capacities = object@capacities }
    else {          capacities = tapply(object@capacities*(1+object@mcDelta),object@ownerPost,sum) }

    totCap = sum(capacities)

    espIntegrand = function(c,t){
      sellerCostParms <- c(list(c),as.list(object@sellerCostParms),
                           lower.tail=as.list(object@sellerCostCDFLowerTail))
      Fc <- do.call(match.fun(object@sellerCostCDF),sellerCostParms)
      val <- (1-Fc)^(totCap-t)-(1-Fc)^totCap
    }


    retval <- sapply(
      capacities,
      function(t.i) {

        if( r < sellerCostBounds[2]) {
          retval <- integrate(espIntegrand,lower=sellerCostBounds[1],upper=r,
                              stop.on.error = FALSE,t=t.i)$value
        }
        else {
          retval <- integrate(espIntegrand,lower=sellerCostBounds[1],upper=sellerCostBounds[2],
                              stop.on.error = FALSE,t=t.i)$value
        }

        return(retval)
      })

    if(!preMerger){

      temp <- rep(NA, length(object@ownerPre))
      temp[object@ownerPre == object@ownerPost] <- retval
      retval <- temp

    }

    if(!exAnte){retval <- retval/calcShares(object,preMerger=preMerger,exAnte=TRUE)}

    return(retval)
  })


## compute producer surplus
#'@rdname PS-methods
#'@export
setMethod(
  f= "calcProducerSurplus",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){


    rev <-  calcRevenues(object, preMerger= preMerger)
    vc <- calcVC(object, preMerger= preMerger)

    ps <- rowSums(rev, na.rm=TRUE) - vc
    names(ps) <- object@labels[[1]]

    return(ps)
  }

)


## Coordinated effects using Grim Trigger strategies
## and optimal defection
#'@rdname PS-methods
#'@export
setMethod(
  f= "calcProducerSurplusGrimTrigger",
  signature= "Bertrand",
  definition= function(object,coalition,discount,preMerger=TRUE,isCollusion=FALSE,...){

    subset <- object@subset
    nprod  <- length(object@labels)
    if(!is.numeric(coalition) ||
       length(coalition) > nprod ||
       !all(coalition %in% 1:nprod)){
      stop ("'coalition' must be a vector of product indices no greater than the number of products")
    }
    if(any(discount<=0 | discount>=1,na.rm=TRUE)){
      stop ("'discount' must be a vector of values between 0 and 1 or NA")
    }


    if(preMerger){owner <- object@ownerPre}
    else{owner <- object@ownerPost}
    anyProds <- as.logical(apply(owner[coalition,]>0,2,max)) # TRUE if a product is made by firm participating in coalition

    if(any(is.na(discount[anyProds]))){
      stop("'discount' must include the discount factors for all products produced by firms with products in the coalition")}


    ## re-calibrate demand and cost parameters under the assumption
    ## that firms are currently colluding
    if(isCollusion){

      ownerPre <- object@ownerPre
      object@ownerPre[coalition,coalition] <- 1
      ## Calculate Demand Slope Coefficients
      object <- calcSlopes(object)
      ## Calculate marginal cost
      object@mcPre <-  calcMC(object,TRUE)
      object@mcPost <- calcMC(object,FALSE)
      object@ownerPre <- ownerPre
    }



    psPunish <- calcProducerSurplus(object,preMerger)
    owner    <- ownerToVec(object,preMerger)

    if(preMerger){
      ownerPre <- object@ownerPre
      object@ownerPre[coalition,coalition] <- 1
      ownerCoalition <- object@ownerPre
      object@pricePre <- calcPrices(object,preMerger,...)
      psCoord  <- calcProducerSurplus(object,preMerger)

      psDefect <- psPunish

      ##compute the producer surplus from defecting
      for(c in unique(owner[coalition])){
        thisOwner <- c==owner
        object@ownerPre <- ownerCoalition
        object@ownerPre[thisOwner,] <- ownerPre[thisOwner,]
        object@pricePre <- calcPrices(object,preMerger,...)
        psDefect[thisOwner] <- calcProducerSurplus(object,preMerger)[thisOwner]
      }


      ## Determine if the firm finds it profitable to cooperate or defect under
      ## GRIM TRIGGER
      IC <- as.vector(ownerPre %*% (psCoord/(1-discount))) >
        as.vector(ownerPre %*% (psDefect + (psPunish*discount/(1-discount))))

    }
    else{

      ownerPost <- object@ownerPost
      object@ownerPost[coalition,coalition] <- 1
      ownerCoalition <- object@ownerPost
      object@pricePost <- calcPrices(object,preMerger,subset=subset,...)
      psCoord  <- calcProducerSurplus(object,preMerger)

      psDefect <- psPunish

      ##compute the producer surplus from defecting
      for(c in unique(owner[coalition])){
        thisOwner <- c==owner
        object@ownerPost <- ownerCoalition
        object@ownerPost[thisOwner,] <- ownerPost[thisOwner,]
        object@pricePost <- calcPrices(object,preMerger,subset=subset,...)
        psDefect[thisOwner] <- calcProducerSurplus(object,preMerger)[thisOwner]
      }


      IC <- as.vector(ownerPost %*% (psCoord/(1-discount))) >=
        as.vector(ownerPost %*% (psDefect + (psPunish*discount/(1-discount))))


    }


    result <- data.frame(Coalition=coalition,Discount=discount,Coord=psCoord,Defect=psDefect,Punish=psPunish,IC=IC)
    rownames(result) <- object@labels

    result <- result[anyProds,]



    return(result)
  }
)
