setClass(
         Class = "Linear",
         contains="Antitrust",
         representation=representation(
         prices           = "vector",
         quantities       = "vector",
         margins          = "vector",
         mc               = "vector",
         pricePre         = "vector",
         pricePost        = "vector",
         diversion        = "anyMatrix",
         symmetry         = "logical"
         ),
          prototype=prototype(

          pricePre      =  vector(),
          pricePost     =  vector(),
          mc            =  vector()


        ),
         validity=function(object){

             ## Sanity Checks

             nprods <- length(object@shares) # count the number of products


             if(nprods != length(object@quantities) ||
                nprods != length(object@margins) ||
                nprods != length(object@prices)){
                 stop("'prices', 'quantities', 'margins', and 'shares' must all be vectors with the same length")}

             if(any(object@prices<0,na.rm=TRUE))             stop("'prices' values must be positive")
             if(any(object@quantities<0,na.rm=TRUE))          stop("'quantities' values must be positive")
             if(any(object@margins<0 | object@margins>1,na.rm=TRUE)) stop("'margins' values must be between 0 and 1")

             if(!all(diag(object@diversion) == 1)){ stop("'diversion' diagonal elements must all equal 1")}
             if(any(abs(object@diversion) > 1)){ stop("'diversion' elements must be between -1 and 1")}

             if(nprods != nrow(object@diversion) ||
                nprods != ncol(object@diversion)){
                 stop("'diversions' must be a square matrix")
             }

             if(length(object@symmetry)!=1){stop("'symmetry' must equal 'TRUE' or 'FALSE'")}

         }
 )


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

     nprods <- length(margins)

     diag(diversion) <- -1


      if(!symmetry){


          slopes <- matrix(margins * prices,ncol=nprods, nrow=nprods,byrow=TRUE)
          slopes <- 1/rowSums(slopes * diversion * ownerPre) * quantities
          slopes <- -t(slopes * diversion)


      }

     else{
         existMargins <- which(!is.na(margins))
         nMargins <-  length(margins[existMargins])
         revenues <- prices*quantities
         k <- existMargins[1] ## choose a diagonal demand parameter corresponding to a provided margin

          minD <- function(s){
              beta <- -s*diversion[k,]/diversion[,k]
              slopesCand <- matrix(beta,ncol=nprods,nrow=nprods,byrow=TRUE)
              slopesCand <- slopesCand*t(diversion)
              elast <- slopesCand * tcrossprod(1/quantities,prices)

              marginsCand <- -1 * as.vector(solve(t(elast * ownerPre)) %*% revenues) / revenues

              measure <- sqrt(sum((margins - marginsCand)^2,na.rm=TRUE))/sqrt(nMargins)
              return(measure)
          }


         minSlope <- optimize(minD,c(-1e12,0))$minimum

         slopes <-  -minSlope*diversion[k,]/diversion[,k]
         slopes <- matrix(slopes,ncol=nprods,nrow=nprods,byrow=TRUE)
         slopes <- slopes*t(diversion)
     }




     dimnames(slopes) <- list(object@labels,object@labels)


     intercept <- as.vector(quantities - slopes %*% prices)

     if(!symmetry &&
        !isTRUE(all.equal(slopes[upper.tri(slopes)],slopes[lower.tri(slopes)]))){
         warning("Matrix of demand slopes coefficients is not symmetric. Demand parameters may not be consistent with utility maximization theory.")}

     if(any(intercept<0))   {warning(  "Some demand intercepts have been found to be negative")}
     if(any(diag(slopes)>0)){warning(  "Some own-slope coefficients have been found to be positive")}

     return(cbind(intercept,slopes))


 }
          )



## Create a function to recover marginal cost using
## demand parameters and supplied prices
setMethod(
          f= "calcMC",
          signature= "Linear",
          definition= function(object){

              object@pricePre <- object@prices


              marginPre <- calcMargins(object,TRUE)

              mc <- (1 - marginPre) * object@prices
              names(mc) <- object@labels

              return(mc)
          }
          )



setMethod(
 f= "calcPrices",
 signature= "Linear",
 definition=function(object,preMerger=TRUE,...){



     if(preMerger){
         mcDelta <- rep(0,length(object@mcDelta))
         owner <- object@ownerPre
     }

     else{
         mcDelta <- object@mcDelta
         owner <- object@ownerPost
     }

     slopes    <- object@slopes[,-1]
     intercept <- object@slopes[,1]
     mc <- object@mc

     prices <- solve( slopes +  t(slopes*owner))
     prices <- prices %*% (crossprod(slopes * owner,mc * (1+mcDelta)) -  intercept)
     prices <- as.vector(prices)

     quantities <- as.vector(intercept + slopes%*%prices)


     if(any(quantities<0)){
         warning("The linear demand model has predicted that some equilibrium levels of output are negative. Re-simulating with the constraint that all equilibrium levels of output are non-negative.")

         FOC <- function(priceCand){

             quantity <- intercept + slopes%*%priceCand
             thisFOC  <- quantity + t(slopes*owner) %*% (priceCand - mc * (1+mcDelta))

             return(as.vector(crossprod(thisFOC)))

         }

         ##Find starting value that always meets boundary conditions
         startParm <- as.vector(solve(slopes) %*% (-intercept + 1))
         minResult <- constrOptim(startParm,FOC,grad=NULL,ui=slopes,ci=-intercept,...)

          if(minResult$convergence != 0){
              warning("'calcPrices' solver may not have successfully converged. Returning unconstrained result. 'constrOptim' reports: '",minResult$message,"'")}

          else{prices <- minResult$par}

     }


     names(prices) <- object@labels

     return(prices)


 }
          )



setMethod(
 f= "calcQuantities",
 signature= "Linear",
 definition=function(object,preMerger=TRUE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

      slopes    <- object@slopes[,-1]
      intercept <- object@slopes[,1]


     quantities <- as.vector(intercept+slopes %*% prices)
     names(quantities) <- object@labels

     return(quantities)

}
 )

setMethod(
 f= "calcShares",
 signature= "Linear",
 definition=function(object,preMerger=TRUE){

     quantities <- calcQuantities(object,preMerger)

     return(quantities/sum(quantities))
 }
)


setMethod(
 f= "elast",
 signature= "Linear",
 definition=function(object,preMerger=TRUE){

       if(preMerger){ prices <- object@pricePre}
       else{          prices <- object@pricePost}

      slopes    <- object@slopes[,-1]


      quantities <-  calcQuantities(object,preMerger)

      elast <- slopes * tcrossprod(1/quantities,prices)
      dimnames(elast) <- list(object@labels,object@labels)

      return(elast)

}
 )

setMethod(
 f= "CV",
 signature= "Linear",
 definition=function(object){

     slopes    <- object@slopes[,-1]

      if(!isTRUE(all.equal(slopes[upper.tri(slopes)],slopes[upper.tri(t(slopes))]))){
                  stop("price coefficient matrix must be symmetric in order to calculate compensating variation. Suggest setting 'symmetry=TRUE'")
              }

     intercept <- object@slopes[,1]
     pricePre  <- object@pricePre
     pricePost <- object@pricePost

     result <- sum(intercept*(pricePost-pricePre)) + .5 * as.vector(pricePost%*%slopes%*%pricePost - pricePre%*%slopes%*%pricePre)

     return(result)
 })

setMethod(
 f= "calcPriceDeltaHypoMon",
 signature= "Linear",
 definition=function(object,prodIndex){

     nprods <- length(prodIndex)
     intercept <- object@slopes[,1]
     slopes <- object@slopes[,-1]
     mc <- object@mc[prodIndex]
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

         maxResult <- constrOptim(startParm,calcMonopolySurplus,
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

     priceDelta <- pricesHM/pricePre[prodIndex] - 1

     return(priceDelta)

 })



linear <- function(prices,quantities,margins, diversions=NULL, symmetry=TRUE,
                     ownerPre,ownerPost,
                     mcDelta=rep(0,length(prices)),
                     labels=paste("Prod",1:length(prices),sep=""),
                   ...
                     ){

    shares <- quantities/sum(quantities)

    if(is.null(diversions)){
        diversions <- tcrossprod(1/(1-shares),shares)
        diag(diversions) <- 1
    }


     result <- new("Linear",prices=prices, quantities=quantities,margins=margins,
                   shares=shares,mcDelta=mcDelta,
                   ownerPre=ownerPre,diversion=diversions, symmetry=symmetry,
                   ownerPost=ownerPost, labels=labels)


     ## Convert ownership vectors to ownership matrices
     result@ownerPre  <- ownerToMatrix(result,TRUE)
     result@ownerPost <- ownerToMatrix(result,FALSE)

     ## Calculate Demand Slope Coefficients
     result@slopes <- calcSlopes(result)

    ##Calculate constant marginal costs
     result@mc <- calcMC(result)

    ## Calculate Pre and Post merger equilibrium prices
    result@pricePre  <- calcPrices(result,TRUE)
    result@pricePost <- calcPrices(result,FALSE)


   return(result)

}


