#' @title \dQuote{Calculating Prices} Methods
#' @name Prices-Methods
#' @docType methods
#' @aliases calcPrices,ANY-method
#' calcPrices,Auction2ndCap-method
#' calcPrices,Cournot-method
#' calcPrices,Linear-method
#' calcPrices,Logit-method
#' calcPrices,LogLin-method
#' calcPrices,AIDS-method
#' calcPrices,LogitCap-method
#' calcPrices,Auction2ndLogit-method
#' calcPrices,Auction2ndLogit-method
#' calcPrices

#' @param object An instance of the respective class (see description for the classes)
#' @param  preMerger If TRUE, the pre-merger ownership structure is used. If FALSE, the post-merger ownership structure is used.
#' Default is TRUE.
#' @param exAnte If \sQuote{exAnte} equals TRUE then the
#'\emph{ex ante} expected result for each firm is produced, while FALSE produces the
#' expected result conditional on each firm winning the auction. Default is FALSE.
#' @param subset A vector of length k where each element equals TRUE if
#' the product indexed by that element should be included in the
#' post-merger simulation and FALSE if it should be excluded. Default is a
#' length k vector of TRUE.
#' @param isMax If TRUE, a check is run to determine if the calculated equilibrium price vector locally maximizes profits.
#' Default is FALSE.
#'
#' @param ... For Logit, additional values that may be used to change the
#' default values of \code{\link[BB]{BBsolve}}, the non-linear equation solver.
#'
#' For others, additional values that may be used to change the default values of \code{\link{constrOptim}}, the non-linear
#' equation solver used to enforce non-negative equilibrium quantities.

#' @description For Auction2ndCap, the calcPrices method computes the expected price that the buyer pays,
#' conditional on the buyer purchasing from a particular firm.
#' @description For Logit, the calcPrices method computes either pre-merger or post-merger equilibrium prices under the assumptions
#' that consumer demand is Logit and firms play a differentiated product Bertrand Nash pricing game.
#' @description For LogitCap, the calcPrices method computes either pre-merger or post-merger equilibrium shares under the assumptions that
#' consumer demand is Logit and firms play a differentiated product Bertrand Nash pricing game with capacity constraints.
#' @description For Logit, the calcPrices method computes either pre-merger or post-merger equilibrium prices under the assumptions
#' that consumer demand is Logit and firms play a differentiated product Bertrand Nash pricing game.
#' @description For LogLin, the calcPrices method computes either pre-merger or post-merger equilibrium prices under the assumptions
#' that consumer demand is Log-Linear and firms play a differentiated product Bertrand Nash pricing game.
#' @description For AIDS, the calcPrices method computes either pre-merger or post-merger equilibrium prices under the assumptions
#' that consumer demand is AIDS and firms play a differentiated product Bertrand Nash pricing game.
#' It returns a length-k vector of NAs if the user did not supply prices.
#'
#' @include AuctionCapMethods.R VerticalClasses.R
#' @keywords methods
NULL

setGeneric (
  name= "calcPrices",
  def=function(object,...){standardGeneric("calcPrices")}
)

#'@rdname Prices-Methods
#'@export
setMethod(
  f= "calcPrices",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){

    if(preMerger){

      quantities <- object@quantityPre
    }
    else{
      quantities <- object@quantityPost
    }

    intercepts <- object@intercepts
    slopes     <- object@slopes


    mktQuant <- colSums(quantities, na.rm=TRUE)

    prices <- ifelse(object@demand == "linear",
                     intercepts + slopes * mktQuant,
                     exp(intercepts) * mktQuant^slopes
    )



    names(prices) <- object@labels[[2]]
    return(prices)


  })

#'@rdname Prices-Methods
#'@export
setMethod(
  f= "calcPrices",
  signature= "Auction2ndCap",
  definition=function(object,preMerger=TRUE,exAnte=TRUE){


    val <- calcProducerSurplus(object,preMerger=preMerger,exAnte=exAnte) + calcMC(object,preMerger=preMerger,exAnte=exAnte)

    names(val) <- object@labels
    return(val)
  })

#'@rdname Prices-Methods
#'@export
setMethod(
  f= "calcPrices",
  signature= "Logit",
  definition=function(object,preMerger=TRUE,isMax=FALSE,subset,...){


    priceStart <- object@priceStart

    if(preMerger){
      owner <- object@ownerPre
      mc    <- object@mcPre
    }
    else{
      owner <- object@ownerPost
      mc    <- object@mcPost
    }

    nprods <- length(object@shares)
    if(missing(subset)){
      subset <- rep(TRUE,nprods)
    }

    if(!is.logical(subset) || length(subset) != nprods ){stop("'subset' must be a logical vector the same length as 'shares'")}

    if(any(!subset)){
      owner <- owner[subset,subset]
      mc    <- mc[subset]
      priceStart <- priceStart[subset]
    }

    priceEst <- rep(NA,nprods)



    ##Define system of FOC as a function of prices
    FOC <- function(priceCand){

      if(preMerger){ object@pricePre[subset] <- priceCand}
      else{          object@pricePost[subset] <- priceCand}


      margins   <- 1 - mc/priceCand
      revenues  <- calcShares(object,preMerger,revenue=TRUE)[subset]
      elasticities     <- t(elast(object,preMerger)[subset,subset])

      thisFOC <- revenues * diag(owner) + as.vector((elasticities * owner) %*% (margins * revenues))

      return(thisFOC)
    }

    ## Find price changes that set FOCs equal to 0
    minResult <- BBsolve(priceStart,FOC,quiet=TRUE,control=object@control.equ,...)

    if(minResult$convergence != 0){warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBsolve' reports: '",minResult$message,"'")}


    if(isMax){

      hess <- genD(FOC,minResult$par) #compute the numerical approximation of the FOC hessian at optimium
      hess <- hess$D[,1:hess$p]
      hess <- hess * (owner>0)   #0 terms not under the control of a common owner

      state <- ifelse(preMerger,"Pre-merger","Post-merger")

      if(any(eigen(hess)$values>0)){warning("Hessian of first-order conditions is not positive definite. ",state," price vector may not maximize profits. Consider rerunning 'calcPrices' using different starting values")}
    }


    priceEst[subset]        <- minResult$par
    names(priceEst) <- object@labels

    return(priceEst)

  }
)


#'@rdname Prices-Methods
#'@export
setMethod(
  f= "calcPrices",
  signature= "Auction2ndLogit",
  definition=function(object,preMerger=TRUE,exAnte=FALSE){

    nprods <- length(object@shares)

    if(preMerger){
      owner <- object@ownerPre
      mc <- object@mcPre
    }
    else{
      owner <- object@ownerPost
      mc <- object@mcPost}

    margins <- calcMargins(object,preMerger,exAnte=FALSE)

    prices <- margins + mc

    if(exAnte){


      prices <- prices * calcShares(object, preMerger=preMerger,revenue=FALSE)
    }


    names(prices) <- object@labels
    return(prices)
  }
)


#'@rdname Prices-Methods
#'@export
setMethod(
  f= "calcPrices",
  signature= "LogitCap",
  definition=function(object,preMerger=TRUE,isMax=FALSE,subset,...){




    if(preMerger){
      owner <- object@ownerPre
      mc    <- object@mcPre
      capacities <- object@capacitiesPre
    }
    else{
      owner <- object@ownerPost
      mc    <- object@mcPost
      capacities <- object@capacitiesPost
    }

    nprods <- length(object@shares)
    if(missing(subset)){
      subset <- rep(TRUE,nprods)
    }

    if(!is.logical(subset) || length(subset) != nprods ){stop("'subset' must be a logical vector the same length as 'shares'")}

    if(any(!subset)){
      owner <- owner[subset,subset]
      mc    <- mc[subset]
      priceStart <- priceStart[subset]
      capacities <- capacities[subset]
    }


    priceEst <- rep(NA,nprods)

    ##Define system of FOC as a function of prices
    FOC <- function(priceCand){

      if(preMerger){ object@pricePre[subset] <- priceCand}
      else{          object@pricePost[subset] <- priceCand}


      margins          <- 1 - mc/priceCand
      revenues         <- calcShares(object, preMerger = preMerger, revenue = TRUE)
      quantities       <- calcQuantities(object, preMerger = preMerger)
      revenues         <- revenues[subset]
      elasticities     <- elast(object,preMerger)[subset,subset]

      thisFOC <- revenues * diag(owner) + as.vector(t(elasticities * owner) %*% (margins * revenues))
      constraint <- ifelse(is.finite(capacities), (quantities - capacities) /object@insideSize, 0)



      measure <- ifelse( constraint != 0, thisFOC + constraint + sqrt(thisFOC^2 + constraint^2), thisFOC)

      return(measure)
    }


    ## Find price changes that set FOCs equal to 0
    minResult <- BBsolve(object@priceStart,FOC,quiet=TRUE,control=object@control.equ,...)

    if(minResult$convergence != 0){warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBsolve' reports: '",minResult$message,"'")}

    priceEst[subset]        <- minResult$par
    names(priceEst) <- object@labels

    if(isMax){

      hess <- genD(FOC,minResult$par) #compute the numerical approximation of the FOC hessian at optimium
      hess <- hess$D[,1:hess$p]
      hess <- hess * (owner>0)   #0 terms not under the control of a common owner

      state <- ifelse(preMerger,"Pre-merger","Post-merger")

      if(any(eigen(hess)$values>0)){warning("Hessian of first-order conditions is not positive definite. ",state," price vector may not maximize profits. Consider rerunning 'calcPrices' using different starting values")}
    }

    return(priceEst)

  }
)


#'@rdname Prices-Methods
#'@export
setMethod(
  f= "calcPrices",
  signature= "Linear",
  definition=function(object,preMerger=TRUE,subset,...){

    slopes    <- object@slopes
    intercept <- object@intercepts
    priceStart<- object@priceStart

    if(preMerger){
      owner <- object@ownerPre
      mc    <- object@mcPre
    }
    else{
      owner <- object@ownerPost
      mc    <- object@mcPost
    }


    nprods <- length(object@quantities)
    if(missing(subset)){
      subset <- rep(TRUE,nprods)
    }

    if(!is.logical(subset) || length(subset) != nprods ){stop("'subset' must be a logical vector the same length as 'quantities'")}



    ##first try the analytic solution
    #      prices <-
    #          solve((slopes*diag(owner)) + (t(slopes)*owner)) %*%
    #              ((t(slopes)*owner) %*% mc - (intercept*diag(owner)))
    #
    #       prices <- as.vector(prices)
    #       quantities <- as.vector(intercept + t(slopes) %*% prices)
    #
    #      ##use the numeric solution if analytic solution yields negative quantities
    #      if(any(subset) || any(quantities<0)){
    #
    #          if(any(quantities<0)){
    #           warning("Equilibrium prices yield negative equilibrium quantities. Recomputing equilibrium prices  under the restriction that equilbrium quantities must be non-negative")
    #                               }
    #          else if(any(subset)){
    #            warning("Elements of 'subset' are  FALSE. Computing equilbrium under the restriction that these products have 0 sales")
    #          }

    FOC <- function(priceCand){

      if(preMerger){ object@pricePre  <- priceCand}
      else{          object@pricePost <- priceCand}

      margins   <- priceCand - mc
      quantities  <- calcQuantities(object,preMerger)

      thisFOC <- quantities*diag(owner) + (t(slopes)*owner) %*% margins
      thisFOC[!subset] <- quantities[!subset] #set quantity equal to 0 for firms not in subset

      return(as.vector(crossprod(thisFOC)))

    }

    ##Find starting value that always meets boundary conditions
    ##startParm <- as.vector(solve(slopes) %*% (-intercept + 1))

    minResult <- constrOptim(object@priceStart,FOC,grad=NULL,ui=slopes,ci=-intercept,...)

    if(!isTRUE(all.equal(minResult$convergence,0,check.names=FALSE))){
      warning("'calcPrices' solver may not have successfully converged.'constrOptim' reports: '",minResult$message,"'")
    }

    prices <- minResult$par

    #}

    names(prices) <- object@labels

    return(prices)


  }
)

#'@rdname Prices-Methods
#'@export
setMethod(
  f= "calcPrices",
  signature= "LogLin",
  definition=function(object,preMerger=TRUE,subset,...){

    if(preMerger){
      owner <- object@ownerPre
      mc    <- object@mcPre
    }
    else{
      owner <- object@ownerPost
      mc    <- object@mcPost
    }

    nprods <- length(object@quantities)
    if(missing(subset)){
      subset <- rep(TRUE,nprods)
    }

    if(!is.logical(subset) || length(subset) != nprods ){stop("'subset' must be a logical vector the same length as 'quantities'")}


    FOC <- function(priceCand){

      if(preMerger){ object@pricePre <- priceCand}
      else{          object@pricePost <- priceCand}


      margins    <- 1 - mc/priceCand
      quantities <- calcQuantities(object,preMerger,revenue=TRUE)
      revenues   <- priceCand*quantities
      elasticities     <- t(elast(object,preMerger))

      thisFOC <- revenues * diag(owner) + as.vector((elasticities * owner) %*% (margins * revenues))
      thisFOC[!subset] <- revenues[!subset] #set quantity equal to 0 for firms not in subset


      return(thisFOC)
    }

    minResult <- BBsolve(object@priceStart,FOC,quiet=TRUE,control=object@control.equ,...)

    if(minResult$convergence != 0){warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBSolve' reports: '",minResult$message,"'")}


    priceEst        <- minResult$par
    names(priceEst) <- object@labels
    return(priceEst)

  }
)


#'@rdname Prices-Methods
#'@export
setMethod(
  f= "calcPrices",
  signature= "AIDS",
  definition=function(object,preMerger=TRUE,...){


    ##if(any(is.na(object@prices)){warning("'prices' contains missing values. AIDS can only predict price changes, not price levels")}

    if(preMerger){prices <- object@prices}
    else{ prices <- object@prices * (1 + object@priceDelta)}

    names(prices) <- object@labels
    return(prices)
  }
)

#'@rdname Prices-Methods
#'@export
setMethod(
  f= "calcPrices",
  signature= "VertBargBertLogit",
  definition=function(object,preMerger=TRUE,...){
    
    up <- object@up
    down <- object@down
    
    
    subsetDown <- down@subset
    subsetUp <- up@subset
    alpha <- down@slopes$alpha
    meanval <- down@slopes$meanval[subsetDown]
    bargparm <- up@bargpower[subsetDown]
    
    priceStartUp <- up@priceStart
    priceStartDown <- down@priceStart
    priceStart <- c(priceStartUp[subsetDown],priceStartDown[subsetDown])
    
    
    nprods <- length(meanval)
    
    FOC <-function(priceCand
                      #, preMerger=TRUE  
                      #subset=rep(TRUE,length(mkt$down$meanval))
                      #,useEst = FALSE, level = "all"
                      ){
      
       
        
      
      priceCandUp= rep(NA, 1,nprods)[subsetDown] 
      priceCandUp <- priceCand[1:length(priceCandUp)]
      priceCandDown <- priceCand[-(1:length(priceCandUp))]
      
      #v <- mkt$vertical
      
      if(preMerger){
        ownerUp <- up@ownerPre[subsetDown,subsetDown]
        ownerDown <- down@ownerPre[subsetDown,subsetDown]
        
        down@pricePre <- priceCandDown
        mcUp <- up@mcPre[subsetDown]
        mcDown <- down@mcPre[subsetDown]
        
      }
      
      else{
        
        ownerUp <- up@ownerPost[subsetDown,subsetDown]
        ownerDown <- down@ownerPost[subsetDown,subsetDown]
        
        down@pricePost <- priceCandDown
        mcUp <- up@mcPost[subsetDown]
        mcDown <- down@mcPost[subsetDown]
      
        }
      
       
        
      
      
      
      shareCandDown  <- calcShares(down, preMerger=preMerger,revenue=FALSE)
      #marginsDownCand <- calcMargins(down, preMerger=preMerger)*priceCandDown
      
     
      elast <-  -alpha*tcrossprod(shareCandDown)
      diag(elast) <- alpha*shareCandDown + diag(elast)
      elast.inv <- try(solve(ownerDown * elast),silent=TRUE)
      if(class(elast.inv) == "try-error"){elast.inv <- MASS::ginv(ownerDown * elast)}
      marginsDownCand <- -as.vector(elast.inv %*% shareCandDown)
      
      
      
      div <- tcrossprod(1/(1-shareCandDown),shareCandDown)
      diag(div) <- -1
      
      
      #if(length(v) > 0 && !preMerger){
      #  marginsDownCand <-  marginsDownCand - elast.inv %*% ( (v$ownerPost.down * elast) %*% (priceCandUp-mcUp) )
      #  #upFOC <- (ownerUp * div) %*% (priceCandUp-mcUp) - (v$ownerPostLambda.down * div) %*% marginsDownCand
      #  marginsUpCand <- as.vector(solve(ownerUp * div) %*% (v$ownerPostLambda.down * div) %*% marginsDownCand) 
      #}
      #else{
        
        marginsUpCand <-  (1-bargparm)/bargparm * as.vector(MASS::ginv(ownerUp * div) %*% (ownerDown * div) %*% marginsDownCand)
        #upFOC<- as.vector(bargparm *((ownerUp * div) %*% (priceCandUp-mcUp))  -  (1-bargparm) * ((ownerDown * div) %*% marginsDownCand))
        
        
      #}
      
      upFOC <-  priceCandUp-mcUp - marginsUpCand
      downFOC <- priceCandDown - priceCandUp - mcDown - marginsDownCand
      
      
      thisFOC= c(downFOC,upFOC)
      
      return(thisFOC)
    }
    
   
      minResult   <- nleqslv(priceStart,FOC)$x
      
      minResultUp <- rep(NA, length(priceStartUp))
      minResultDown <- rep(NA, length(priceStartDown))
      
        
        
      minResultUp[subsetDown] <- minResult[1:length(priceStartUp[subsetDown])] 
      minResultDown[subsetDown] <- minResult[-(1:length(priceStartUp[subsetDown]))]
    
   
      return(list(up=minResultUp,down=minResultDown)) 
  }
)
