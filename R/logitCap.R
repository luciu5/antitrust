
setClass(
         Class   = "LogitCap",
         contains="Logit",
         representation=representation(
         mktSize               = "numeric",
         capacities           = "numeric"

         ),


         validity=function(object){





             nprods <- length(object@shares)


             if(nprods != length(object@capacities)){
                 stop("'prices', 'capacities' must all be vectors with the same length")}

             if(any(is.na(object@capacities) |
                    !is.finite(object@capacities) |
                    object@capacities<0 ,na.rm=TRUE)){stop("'capacities' values must be positive, finite numbers")}

             if(length(object@mktSize)!=1 || isTRUE(object@mktSize<0)){stop("mktSize must be a positive number")}

             if(any(object@mktSize*object@shares > object@capacities)){stop("utilization is greater than capacity")}

             if(all(object@mktSize*object@shares == object@capacities)){stop("utilization cannot equal capacity for all products")}

             if(any(is.na(object@margins[object@mktSize*object@shares == object@capacities]))){
                 stop("'margins' cannot equal NA for capacity constrained products")
                 }

                return(TRUE)

         })


setMethod(
          f= "calcSlopes",
          signature= "LogitCap",
          definition=function(object){

              ## Uncover Demand Coefficents


              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              margins      <-  object@margins
              prices       <-  object@prices
              capacities  <-  object@capacities/object@mktSize
              idx          <-  object@normIndex

              if(is.na(idx)){
                  idxShare <- 1 - object@shareInside
                  idxPrice <- 0
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


                  FOC <- revenues + (elast * ownerPre * notBinds) %*% (margins * revenues)

                  ## omit the FOCs of single product, capacity constrained firms
                  measure <- sum(as.vector(FOC[!singleConstrained])^2,na.rm=TRUE)

                  return(measure)
              }

              minAlpha <- optimize(minD,c(-1e6,0))$minimum


              meanval <- log(shares) - log(idxShare) - minAlpha * (prices - idxPrice)

              names(meanval)   <- object@labels

              object@slopes    <- list(alpha=minAlpha,meanval=meanval)


              return(object)
          }
          )

setMethod(
 f= "calcQuantities",
 signature= "LogitCap",
 definition=function(object,preMerger=TRUE){

     quantities <- object@mktSize * calcShares(object,preMerger,revenue=FALSE)

     return(quantities)

 }
          )

## compute margins
setMethod(
 f= "calcMargins",
 signature= "LogitCap",
 definition=function(object,preMerger=TRUE){


     if( preMerger) {
         margins <- object@margins #capacity-constrained margins not identified -- use supplied margins
         constrained <- object@capacities == object@mktSize*object@shares

         owner  <- object@ownerPre
         revenue<- calcShares(object,preMerger,revenue=TRUE)[!constrained]
         elast <-  elast(object,preMerger)
         margins[!constrained] <-  -1 * as.vector(solve(t(elast*owner)[!constrained,!constrained]) %*% revenue) / revenue

     }

     else{
         prices <- object@pricePost
         mc     <- calcMC(object,preMerger)

         margins <- 1 - mc/prices
     }


     names(margins) <- object@labels

     return(as.vector(margins))
     }

 )



setMethod(
 f= "calcPrices",
 signature= "LogitCap",
 definition=function(object,preMerger=TRUE,isMax=FALSE,...){


     capacities <- object@capacities

     if(preMerger){owner <- object@ownerPre}
     else{owner <- object@ownerPost}

     mc <- calcMC(object,preMerger)


     ##Define system of FOC as a function of prices
     FOC <- function(priceCand){

         if(preMerger){ object@pricePre <- priceCand}
         else{          object@pricePost <- priceCand}


         margins          <- 1 - mc/priceCand
         quantities       <- calcQuantities(object,preMerger)
         revenues         <- quantities * priceCand
         elasticities     <- elast(object,preMerger)

         thisFOC <- revenues + as.vector(t(elasticities * owner) %*% (margins * revenues))
         constraint <- quantities - capacities

         measure <- thisFOC + constraint + sqrt(thisFOC^2 + constraint^2)

         return(measure)
     }

     ## Find price changes that set FOCs equal to 0
     minResult <- nleqslv(object@priceStart,FOC,...)

      if(minResult$termcd != 1){warning("'calcPrices' nonlinear solver may not have successfully converged. 'nleqslv' reports: '",minResult$message,"'")}

     priceEst        <- minResult$x
     names(priceEst) <- object@labels
  if(isMax){

         hess <- genD(FOC,priceEst) #compute the numerical approximation of the FOC hessian at optimium
         hess <- hess$D[,1:hess$p]
         hess <- hess * (owner>0)   #0 terms not under the control of a common owner

         state <- ifelse(preMerger,"Pre-merger","Post-merger")

         if(any(eigen(hess)$values>0)){warning("Hessian of first-order conditions is not positive definite. ",state," price vector may not maximize profits. Consider rerunning 'calcPrices' using different starting values")}
     }

     return(priceEst)

 }
 )



setMethod(
 f= "calcPricesHypoMon",
 signature= "LogitCap",
 definition=function(object,prodIndex,...){


     mc       <- calcMC(object,TRUE)[prodIndex]
     pricePre <- object@pricePre

      FOC <- function(priceCand){

          thisPrice <- pricePre
          thisPrice[prodIndex] <- priceCand

          if(preMerger){ object@pricePre <- thisPrice}
          else{          object@pricePost <- thisPrice}


          margins          <- 1 - mc/priceCand
          quantities       <- calcQuantities(object,preMerger)[prodIndex]
          revenues         <- quantities * priceCand
          elasticities     <- elast(object,preMerger)[prodIndex,prodIndex]

          thisFOC <- revenues + as.vector(t(elasticities) %*% (margins * revenues))
          constraint <- quantities - object@capacities[prodIndex]

          measure <- thisFOC + constraint + sqrt(thisFOC^2 + constraint^2)

         return(measure)
      }



      ## Find price changes that set FOCs equal to 0
     minResult <- nleqslv(object@priceStart[prodIndex],FOC,...)

     if(minResult$termcd != 1){warning("'calcPricesHypoMon' nonlinear solver may not have successfully converged. 'nleqslv' reports: '",minResult$message,"'")}


     pricesHM <- minResult$x
      #priceDelta <- pricesHM/pricePre[prodIndex] - 1
      #names(priceDelta) <- object@labels[prodIndex]
     names(priceHM) <- object@labels[prodIndex]

     return(priceHM)

 })



logit.cap <- function(prices,shares,margins,
                      ownerPre,ownerPost,
                      capacities,
                      mktSize=sum(capacities),
                      normIndex=ifelse(sum(shares)<1,NA,1),
                      mcDelta=rep(0,length(prices)),
                      priceStart = prices,
                      isMax=FALSE,
                      labels=paste("Prod",1:length(prices),sep=""),
                      ...
                      ){


    ## Create LogitCap  container to store relevant data
    result <- new("LogitCap",prices=prices, shares=shares,
                  margins=margins,capacities=capacities, mktSize=mktSize,
                  normIndex=normIndex,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  mcDelta=mcDelta,
                  priceStart=priceStart,shareInside=sum(shares),
                  labels=labels)

    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate Demand Slope Coefficients
    result <- calcSlopes(result)


    ## Solve Non-Linear System for Price Changes
    result@pricePre  <- calcPrices(result,preMerger=TRUE,isMax=isMax,...)
    result@pricePost <- calcPrices(result,preMerger=FALSE,isMax=isMax,...)


    return(result)

}

