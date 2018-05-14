
setClass(
         Class   = "LogitCap",
         contains="Logit",
         representation=representation(
         capacitiesPre           = "numeric",
         capacitiesPost          = "numeric"

         ),


         validity=function(object){





             nprods <- length(object@shares)


             if(nprods != length(object@capacitiesPre)){
                 stop("'shares', 'capacitiesPre' must all be vectors with the same length")}
             if(length(object@capacitiesPost) != length(object@capacitiesPre)){
               stop("'capacitiesPre', 'capacitiesPost', must be vectors with the same length")}
             

             if(any(is.na(object@capacitiesPre) |
                    #!is.finite(object@capacitiesPre) |
                    object@capacitiesPre<0 ,na.rm=TRUE)){stop("'capacitiesPre' values must be positive numbers")}


             if(any(is.na(object@capacitiesPost) |
                    #!is.finite(object@capacitiesPost) |
                    object@capacitiesPost<0 ,na.rm=TRUE)){stop("'capacitiesPost' values must be positive numbers")}
             
             if(is.na(object@insideSize)){stop("'insideSize' must eqal the total pre-merger units sold in the market")}
             
             if(any(object@insideSize*object@shares > object@capacitiesPre)){warning("utilization is greater than capacity")}

             if(identical(object@insideSize*object@shares,object@capacitiesPre)){warning("utilization equal capacity for all products")}

             if(any(is.na(object@margins[object@insideSize*object@shares == object@capacitiesPre]))){
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
              capacities  <-  object@capacitiesPre/object@insideSize
              idx          <-  object@normIndex

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

              minAlpha <- optimize(minD,c(-1e6,0),
                                   tol=object@control.slopes$reltol)$minimum


              meanval <- log(shares) - log(idxShare) - minAlpha * (prices - idxPrice)

              names(meanval)   <- object@labels

              object@slopes    <- list(alpha=minAlpha,meanval=meanval)
              object@mktSize <- object@insideSize / sum(shares)
              


              return(object)
          }
          )



## compute margins
setMethod(
 f= "calcMargins",
 signature= "LogitCap",
 definition=function(object,preMerger=TRUE){

   margins <- object@margins #capacity-constrained margins not identified -- use supplied margins

     if( preMerger) {
       capacities <- object@capacitiesPre
       
     }
    else{
      
      capacities <- object@capacitiesPost
    }
       
         
        
         quantities <- calcQuantities(object, preMerger=TRUE)
         constrained <-  abs(capacities - quantities) < 1e-5

         owner  <- object@ownerPre
         revenue<- calcShares(object,preMerger,revenue=TRUE)[!constrained]
         elast <-  elast(object,preMerger)
         margins[!constrained] <-  -1 * as.vector(MASS::ginv(t(elast*owner)[!constrained,!constrained]) %*% revenue) / revenue

     

     names(margins) <- object@labels

     return(as.vector(margins))
     }

 )



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

          margins          <- 1 - mc/priceCand
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



logit.cap <- function(prices,shares,margins,
                      ownerPre,ownerPost,
                      capacitiesPre = rep(Inf , length(shares)),
                      capacitiesPost = capacitiesPre,
                      insideSize,
                      normIndex=ifelse(sum(shares)<1,NA,1),
                      mcDelta=rep(0,length(prices)),
                      subset=rep(TRUE,length(prices)),
                      priceOutside=0,
                      priceStart = prices,
                      isMax=FALSE,
                      control.slopes,
                      control.equ,
                      labels=paste("Prod",1:length(prices),sep=""),
                      ...
                      ){


    ## Create LogitCap  container to store relevant data
    result <- new("LogitCap",prices=prices, shares=shares,
                  margins=margins,
                  capacitiesPre=capacitiesPre,
                  capacitiesPost=capacitiesPost, 
                  insideSize=insideSize,
                  normIndex=normIndex,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  mcDelta=mcDelta,
                  subset=subset,
                  priceOutside=priceOutside,
                  priceStart=priceStart,shareInside=sum(shares),
                  labels=labels)

    if(!missing(control.slopes)){
      result@control.slopes <- control.slopes
    }
    if(!missing(control.equ)){
      result@control.equ <- control.equ
    }
     
    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate Demand Slope Coefficients
    result <- calcSlopes(result)


    ## Calculate marginal cost
    result@mcPre <-  calcMC(result,TRUE)
    result@mcPost <- calcMC(result,FALSE)

    ## Solve Non-Linear System for Price Changes
    result@pricePre  <- calcPrices(result,preMerger=TRUE,isMax=isMax,...)
    result@pricePost <- calcPrices(result,preMerger=FALSE,isMax=isMax,subset=subset,...)


    return(result)

}

