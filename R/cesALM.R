setClass(
         Class   = "CESALM",
         contains="CES",
         representation=representation(
          parmsStart="numeric"
         ),
         prototype=prototype(
         normIndex         =  NA,
         control.slopes = list( 
           factr = 1e7 
         )
         ),

         validity=function(object){



             nMargins  <- length(object@margins[!is.na(object@margins)])

             if(nMargins<2){stop("At least 2 elements of 'margins' must not be NA in order to calibrate demand parameters")}

             if(object@shareInside!=1){
                 stop(" sum of 'shares' must equal 1")
             }

              if(length(object@parmsStart)!=2){
                 stop("'parmsStart' must a vector of length 2")
                 }
         }
         )

setMethod(
          f= "calcSlopes",
          signature= "CESALM",
          definition=function(object){

              ## Uncover Demand Coefficents

              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              margins      <-  object@margins
              prices       <-  object@prices

              nprods <- length(object@shares)

              ##identify which products have enough margin information
              ##  to impute Bertrand margins
              #isMargin    <- matrix(margins,nrow=nprods,ncol=nprods,byrow=TRUE)
              #isMargin[ownerPre==0]=0
              #isMargin    <- !is.na(rowSums(isMargin))

              minD <- function(theta){

                  gamma <- theta[1]
                  sOut  <- theta[2]

                  probs <- shares * (1 - sOut)
                  elasticity <- (gamma - 1 ) * matrix(probs,ncol=nprods,nrow=nprods)
                  diag(elasticity) <- -gamma + diag(elasticity)
                  
                  revenues <- probs * prices
                  marginsCand <- -1 * as.vector(MASS::ginv(elasticity * ownerPre) %*% (revenues * diag(ownerPre))) / revenues
                  measure <- sum((margins - marginsCand)^2,na.rm=TRUE)

                 
                  return(measure)
              }

               ## Constrain optimizer to look  alpha > 1,  0 < sOut < 1
              lowerB <- c(1,0)
              upperB <- c(Inf,.99999)

              minGamma <- optim(object@parmsStart,minD,
                                method="L-BFGS-B",
                                lower= lowerB,upper=upperB,
                                control=object@control.slopes)$par

              if(isTRUE(all.equal(minGamma[2],0,check.names=FALSE))){stop("ERROR: Estimated outside share is close to 0. Use `ces' function instead")}
              
              
              meanval <- log(shares * (1 - minGamma[2])) - log(minGamma[2]) + (minGamma[1] - 1) * (log(prices) - log(object@priceOutside))
              meanval <- exp(meanval)
              
              
              
              names(meanval)   <- object@labels


              object@slopes      <- list(alpha=1/(1 - minGamma[2]) - 1,gamma=minGamma[1],meanval=meanval)
              object@shareInside <- 1-minGamma[2]

              return(object)

          }

          )


ces.alm <- function(prices,shares,margins,
                      ownerPre,ownerPost,
                      mcDelta=rep(0,length(prices)),
                      subset=rep(TRUE,length(prices)),
                      priceOutside=1,
                      priceStart = prices,
                      isMax=FALSE,
                      parmsStart,
                      control.slopes,
                      control.equ,
                      labels=paste("Prod",1:length(prices),sep=""),
                      ...
                      ){


    if(missing(parmsStart)){
        parmsStart <- rep(.1,2)
        nm <- which(!is.na(margins))[1] 
        parmsStart[1] <- 1/(margins[nm]*(1-shares[nm])) - shares[nm]/(1-shares[nm]) #ballpark alpha for starting values
    }

   
  
    ## Create CES  container to store relevant data
    result <- new("CESALM",prices=prices, shares=shares,
                  margins=margins,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  mcDelta=mcDelta,
                  subset=subset,
                  priceOutside=priceOutside,
                  priceStart=priceStart,
                  shareInside=sum(shares),
                  parmsStart=parmsStart,
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


