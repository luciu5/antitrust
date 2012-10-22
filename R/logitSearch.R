
setClass(
         Class   = "LogitSearch",
         contains="Logit",

         representation=representation(
         parmsStart="numeric",
         searchSetsPre ="matrix",
         searchSetsPost="matrix",
         locations="numeric",
         distances="numeric"
         ),
          prototype=prototype(

          searchSetsPre   =  matrix(),
          searchSetsPost  =  matrix()

         ),
           validity=function(object){



               if(length(object@locations) != length(object@shares)){
                   stop("'locations' length must be less than or equal to the number of products")
               }
               if (length(unique(object@locations)) != length(object@distances)){
                   stop("'distances' length must equal the number of distinct locations")
               }

               if(object@shareInside==1){
                   stop("Sum of 'shares' must be less than 1")
                   }
           }

         )

setMethod(
          f= "calcSearchSets",
          signature= "LogitSearch",
          definition=function(object,preMerger=TRUE){

              locations <- unique(object@locations)


              numLocations <- length(locations)

              result <- do.call("cbind",
                            lapply(1:numLocations,
                                   function(x){
                                       c=combn(locations,x)
                                       r=matrix(0,ncol=ncol(c),nrow=numLocations)
                                       r[1:nrow(c),]=c
                                       return(r)}
                                   )
                                )

              return(result)

          }
          )

setMethod(
          f= "calcSlopes",
          signature= "LogitSearch",
          definition=function(object){

              ## Uncover Demand Coefficents


              ownerPre     <-  object@ownerPre
              shares       <-  object@shares

              nprods        <- length(shares)

              margins      <-  object@margins
              naMargins    <- !is.na(object@margins)

              prices       <-  object@prices

              object@pricePre <- prices

              minD <- function(theta){

                  alpha <- theta[1]
                  tau   <- theta[2]
                  meanval <- theta[-(1:2)]

                  object@slopes <- list(alpha=alpha,meanval=meanval,tau=tau)
                  sharesCand    <- calcShares(object,preMerger=TRUE)

                  revenuesCand <- sharesCand * prices

                  ## the following returns the elasticity TRANSPOSED
                  elast <- -alpha *  matrix(prices * sharesCand,ncol=nprods,nrow=nprods)
                  diag(elast) <- alpha*prices + diag(elast)


                  marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% revenuesCand) / revenuesCand

                  moments <- c(margins - marginsCand,shares - sharesCand)
                  measure <- sum(moments^2,na.rm=TRUE)

                  return(measure)
              }

               ## Constrain optimizer to look  alpha <0, tau < 0
              lowerB <- upperB <- rep(Inf,nprods + 2)
              lowerB <- -lowerB
              upperB[1:2] <- 0


              minTheta <- optim(object@parmsStart,minD,method="L-BFGS-B",lower= lowerB,upper=upperB)
              print(minTheta$value)
              meanval <- minTheta$par[-(1:2)]
              names(meanval) <- object@labels

              object@slopes    <- list(alpha=minTheta$par[1],meanval=meanval,tau=minTheta$par[2])


              return(object)
          }
          )

setMethod(
 f= "calcShares",
 signature= "LogitSearch",
 definition=function(object,preMerger=TRUE,revenue=FALSE){

      if(preMerger){
          prices     <- object@pricePre
          searchSets <- object@searchSetsPre
      }
      else{
          prices     <- object@pricePost
          searchSets <- object@searchSetsPost
      }

      alpha    <-  object@slopes$alpha
      meanval  <-  object@slopes$meanval
      tau      <-  object@slopes$tau
      distances <- object@distances
      locations <- object@locations
      isInside <- as.numeric(object@shareInside<1)


      IV <- searchSets

      searchCost <- tau * colSums((searchSets > 0) * distances)

      util <- exp(meanval + alpha*prices)

      ivNest <- tapply(util,locations,sum)
      ivNest <- ivNest[as.vector(searchSets)]



      IV[IV!=0] <- ivNest



      shareCond <-  sapply(1:ncol(searchSets),function(x){locations %in% searchSets[,x]}) #expensive!
      shareCond <- shareCond * util
      shareCond <- t(t(shareCond) / (isInside + colSums(IV)))

      IV        <- log(isInside + colSums(IV)) + searchCost
      shareSet  <- exp(IV)/sum(exp(IV))


      shares <- rowSums(t(t(shareCond)*shareSet))

      if(revenue){shares <- prices*shares/sum(prices*shares,1-sum(shares))}

      names(shares) <- object@labels

      return(shares)




}
 )








setMethod(
          f= "CV",
          signature= "LogitSearch",
          definition=function(object){


              names(result) <- NULL
              return(result)

 })



logit.search <- function(prices,shares,margins,
                         ownerPre,ownerPost,
                         locations,
                         distances,
                         normIndex=ifelse(sum(shares) < 1,NA,1),
                         mcDelta=rep(0,length(prices)),
                         priceStart = prices,
                         isMax=FALSE,
                         parmsStart,
                         labels=paste("Prod",1:length(prices),sep=""),
                        ...
                        ){


    if(missing(locations)){
        if(is.matrix(ownerPre)){
        locations <- unique(ownerPre)
        locations <- locations * 1: nrow(locations)
        locations <- apply(locations,2,max)

    }
        else{ locations <- as.numeric(as.factor(ownerPre))
          }
    }
    if(missing(distances)){
        distances <- rep(1,length(unique(locations)))
    }


    if(missing(parmsStart)){

        parmsStart <- -runif(length(shares)+2)
                   }

    ## Create LogitNests  container to store relevant data
    result <- new("LogitSearch",prices=prices, margins=margins,
                  shares=shares,mcDelta=mcDelta,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  locations=locations,
                  distances=distances,
                  normIndex=normIndex,
                  priceStart=priceStart,
                  shareInside=sum(shares),
                  parmsStart=parmsStart,
                  labels=labels)

    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    result@searchSetsPre  <- calcSearchSets(result,TRUE)
    result@searchSetsPost <- calcSearchSets(result,FALSE)
    ## Calculate Demand Slope Coefficients
    result <- calcSlopes(result)

    ## Solve Non-Linear System for Price Changes
    result@pricePre  <- calcPrices(result,preMerger=TRUE,isMax=isMax,...)
    result@pricePost <- calcPrices(result,preMerger=FALSE,isMax=isMax,...)

    return(result)

}

