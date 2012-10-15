setClass(
         Class = "LogLin",
         contains="Linear",
         prototype=prototype(
         symmetry=FALSE
         ),
          validity=function(object){

             ## Sanity Checks


             nprods <- length(object@prices)
             if(any(is.na(object@margins))){
                 stop("'margins' cannot contain NA values")
                 }
             if(nprods != length(object@priceStart)){
                 stop("'priceStart' must have the same length as 'prices'")}
 })


setMethod(
 f= "calcSlopes",
 signature= "LogLin",
 definition=function(object){



     margins <- object@margins
     quantities <- object@quantities
     prices <- object@prices

     revenues <- prices * quantities

     nprods <- length(margins)

     diversion <- object@diversion * tcrossprod(quantities,1/quantities)

     slopes <- matrix(margins * revenues,ncol=nprods, nrow=nprods,byrow=TRUE)
     slopes <- revenues / rowSums(slopes * diversion * object@ownerPre)
     slopes <- -t(slopes * diversion)

     dimnames(slopes) <- list(object@labels,object@labels)

     intercept <- as.vector(log(quantities) - slopes %*% log(prices))
     names(intercept) <-  object@labels

     object@slopes <- slopes
     object@intercepts <- intercept

     return(object)


 }
          )


setMethod(
 f= "calcPrices",
 signature= "LogLin",
 definition=function(object,preMerger=TRUE,...){

     if(preMerger){owner <- object@ownerPre}
     else{owner <- object@ownerPost}
     mc <- calcMC(object,preMerger)

 FOC <- function(priceCand){

     if(preMerger){ object@pricePre <- priceCand}
     else{          object@pricePost <- priceCand}


     margins   <- 1 - mc/priceCand
     revenues  <- calcShares(object,preMerger,revenue=TRUE)
     elasticities     <- elast(object,preMerger)

     thisFOC <- revenues + as.vector(t(elasticities * owner) %*% (margins * revenues))

     return(thisFOC)
 }

 minResult <- nleqslv(object@priceStart,FOC,...)

 if(minResult$termcd != 1){warning("'calcPrices' nonlinear solver may not have successfully converged. 'nleqslv' reports: '",minResult$message,"'")}

 priceEst        <- minResult$x
 names(priceEst) <- object@labels
 return(priceEst)

}
 )


setMethod(
 f= "calcQuantities",
 signature= "LogLin",
 definition=function(object,preMerger=TRUE,...){


     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

     slopes    <- object@slopes
     intercept <- object@intercepts
     quantities <- exp(intercept) * apply(prices^t(slopes),2,prod)
     names(quantities) <- object@labels

     return(quantities)

}
 )




setMethod(
 f= "elast",
 signature= "LogLin",
 definition=function(object,preMerger=TRUE,market=FALSE){

    if(market){

        quantities <-  calcQuantities(object,preMerger)
        prices     <-  calcPrices(object,preMerger)
        elast      <-  sum(t(t(object@slopes * quantities) * 1/prices)) / sum(quantities) * sum(quantities * prices / sum(quantities))


       }

    else{
        elast    <- object@slopes
        dimnames(elast) <- list(object@labels,object@labels)
    }

      return(elast)

}
 )

setMethod(
 f= "CV",
 signature= "LogLin",
 definition=function(object){
     stop("CV method not currently available for 'LogLin' Class")

 })


setMethod(
 f= "calcPriceDeltaHypoMon",
 signature= "LogLin",
 definition=function(object,prodIndex){


     mc <- calcMC(object,TRUE)[prodIndex]
     pricePre <- object@pricePre

     calcMonopolySurplus <- function(priceCand){

         pricePre[prodIndex] <- priceCand
         object@pricePre <- pricePre
         quantityCand <- calcQuantities(object,TRUE)


         surplus <- (priceCand-mc)*quantityCand[prodIndex]


         return(sum(surplus))
     }


     minResult <- optim(object@prices[prodIndex],calcMonopolySurplus,
                              method = "L-BFGS-B",lower = 0,
                              control=list(fnscale=-1))

     pricesHM <- minResult$par
     priceDelta <- pricesHM/pricePre[prodIndex] - 1

     return(priceDelta)

 })





loglinear <- function(prices,quantities,margins,diversions,
                     ownerPre,ownerPost,
                     mcDelta=rep(0,length(prices)),
                     priceStart=prices,
                     labels=paste("Prod",1:length(prices),sep=""),...
                     ){




    shares=quantities/sum(quantities)


    if(missing(diversions)){
        diversions <-  tcrossprod(1/(1-shares),shares)
        diag(diversions) <- -1
    }


    result <- new("LogLin",prices=prices, quantities=quantities,margins=margins,
                  shares=shares,mcDelta=mcDelta, priceStart=priceStart,
                  ownerPre=ownerPre,diversion=diversions,
                  ownerPost=ownerPost, labels=labels)


    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate Demand Slope Coefficients
    result <- calcSlopes(result)

    ## Calculate pre and post merger equilibrium prices
    result@pricePre  <- calcPrices(result,TRUE,...)
    result@pricePost <- calcPrices(result,FALSE,...)


   return(result)

}

