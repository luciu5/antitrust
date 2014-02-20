setClass(
         Class = "AIDS",
         contains="Linear",
         representation=representation(
         priceStart  = "numeric",
         priceDelta       = "numeric",
         mktElast         = "numeric"
         ),
        prototype=prototype(
        priceDelta       =  numeric(),
        mktElast         =  numeric()
        ),

         validity=function(object){



             nprods <- length(object@shares)

             if(!isTRUE(all.equal(sum(object@shares),1))){
                 stop("The sum of 'shares' values must equal 1")}

             if(length(object@margins[!is.na(object@margins)])<2){
                 stop("'margins' must contain at least two non-missing margins in order to calibrate demand parameters")
             }

             if(!isTRUE(all.equal(rowSums(object@diversion,na.rm=TRUE),rep(0,nprods)))){ stop("'diversion' rows must sum to 0")}


             ## Need to write a check that tests if the margins for all the firm's products is present

         }

             )








setMethod(
 f= "calcSlopes",
 signature= "AIDS",
 definition=function(object){


     ## Uncover linear demand slopes
     shares     <- object@shares
     prices     <- object@prices
     margins    <- object@margins
     diversion  <- object@diversion
     labels     <- object@labels
     ownerPre   <- object@ownerPre

     existMargins <- which(!is.na(margins))

     k <- existMargins[1] # choose a diagonal demand parameter corresponding to a provided margin

     minD <- function(s){


         coefCand <- s[1]
         mktElast <- s[2]

         diagCand <- coefCand*diversion[k,]/diversion[,k]
         slopesCand <- -diversion*diagCand # slopes symmetric, no need to transpose


         elast <- t(slopesCand/shares) + shares * (mktElast + 1) #Caution: returns TRANSPOSED elasticity matrix
         diag(elast) <- diag(elast) - 1

         marginsCand <- -1 * as.vector(ginv(elast * ownerPre) %*% (shares * diag(ownerPre))) / shares


         measure <- sum((margins - marginsCand)^2,na.rm=TRUE)

         return(measure)
     }

     ##constrain mktElast, coef to be negative and mktElast to be less than the own price elasticity in absolute value
     ui <- matrix(c(-1/shares[k],0,-1,1 - shares[k],-1,0),ncol=2,nrow=3)
     ci <- c(shares[k] - 1,0,0)

     parmsStart <- c(-1 - shares[k]*(1-shares[k]),-2) # These starting values always satisfy constraints

     minParam <- constrOptim(parmsStart,minD,ui=ui,ci=ci,grad=NULL)$par

     diagB <- minParam[1]*diversion[k,]/diversion[,k]
     B       <- -diversion*diagB # slopes symmetric, no need to transpose


     dimnames(B) <- list(object@labels,object@labels)

     object@slopes <- B

     if(abs(minParam[2])>5){warning("'mktElast' estimate is large.")}
     object@mktElast <- minParam[2]
     object@intercepts <- as.vector(shares - B%*%log(prices))
     names(object@intercepts) <- object@labels


     return(object)

 }

      )



setMethod(
 f= "calcPriceDelta",
 signature= "AIDS",
 definition=function(object,isMax=FALSE,subset=TRUE,...){


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
     minResult <- BBsolve(object@priceStart,FOC,quiet=TRUE,...)


     if(minResult$convergence != 0){warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBsolve' reports: '",minResult$message,"'")}

     if(isMax){

         hess <- genD(FOC,minResult$par) #compute the numerical approximation of the FOC hessian at optimium
         hess <- hess$D[,1:hess$p]


         if(any(eigen(hess)$values>0)){warning("Hessian of first-order conditions is not positive definite. Price vector may not maximize profits. Consider rerunning 'calcPrices' using different starting values")}
     }

     deltaPrice <- exp(minResult$par)-1
     names(deltaPrice) <- object@labels

     return(deltaPrice)
 }
          )


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

setMethod(
 f= "calcShares",
 signature= "AIDS",
 definition=function(object,preMerger=TRUE,revenue=TRUE){

     if(!revenue &&
        any(is.na(object@prices))
        ){
         warning("'prices' contains missing values. Some results are missing")
        }

     prices <- calcPrices(object,preMerger)
     shares <- object@shares


     if(!preMerger){
         shares <-  shares + as.vector(object@slopes %*% log(object@priceDelta + 1))
     }

      if(!revenue){shares <- (shares/prices)/sum(shares/prices)}

     names(shares) <- object@labels
     return(shares)
 }
 )




setMethod(
          f= "elast",
          signature= "AIDS",
          definition=function(object,preMerger=TRUE,market=FALSE){


               if(market){

                   return(object@mktElast)

               }

               else{
                   shares <- calcShares(object,preMerger)

                   elast <- t(object@slopes/shares) + shares * (object@mktElast + 1) #Caution: returns TRANSPOSED elasticity matrix
                   diag(elast) <- diag(elast) - 1
                   dimnames(elast) <-  list(object@labels,object@labels)

                   return(t(elast))

               }
           }
          )




setMethod(
 f= "diversion",
 signature= "AIDS",
 definition=function(object,preMerger=TRUE,revenue=TRUE){

     if(revenue){
         diversion <- -t(object@slopes)/diag(object@slopes)
         dimnames(diversion) <-  list(object@labels,object@labels)
         return(diversion)
     }

     else{callNextMethod(object,preMerger,revenue)}

     }

          )

setMethod(
 f= "diversionHypoMon",
 signature= "AIDS",
 definition=function(object){

   return(diversion(object,revenue=FALSE))

          })

setMethod(
 f= "upp",
 signature= "AIDS",
 definition=function(object){

     if(any(is.na(object@prices))){stop("UPP cannot be calculated without supplying values to 'prices'")}

     else{return(callNextMethod(object))}
     })


setMethod(
 f= "cmcr",
 signature= "AIDS",
 definition=function(object){


     ownerPre  <- object@ownerPre
     ownerPost <- object@ownerPost

     isParty <- rowSums( abs(ownerPost - ownerPre) ) > 0

     sharesPre <- calcShares(object,TRUE)
     sharesPre <- tcrossprod(1/sharesPre,sharesPre)

     marginPre <- calcMargins(object,TRUE)


     elastPre  <- t(elast(object,TRUE))

     divPre    <- elastPre/diag(elastPre)


     Bpost      <- divPre * sharesPre * ownerPost
     marginPost <- -1 * as.vector(ginv(Bpost) %*% (diag(ownerPost)/diag(elastPre))
                                  )

     cmcr <- (marginPost - marginPre)/(1 - marginPre)
     names(cmcr) <- object@labels

     cmcr <- cmcr[isParty]

    return(cmcr * 100)
}
 )




setMethod(
          f= "CV",
          signature= "AIDS",
          definition=function(object,totalRevenue){

              ## computes compensating variation using the closed-form solution found in LaFrance 2004, equation 10

              if(any(is.na(object@prices))){stop("Compensating Variation cannot be calculated without supplying values to 'prices'")}


              slopes <- object@slopes
              intercepts <- object@intercepts
              pricePre <- log(object@pricePre)
              pricePost <- log(object@pricePost)

              result <- sum(intercepts*(pricePost-pricePre)) + .5 * as.vector(t(pricePost)%*%slopes%*%pricePost - t(pricePre)%*%slopes%*%pricePre)


              if(missing(totalRevenue)){
                  warning("'totalRevenue' is missing. Calculating CV as a percentage change in (aggregate) income")
                  return(result*100)}

               else{
                  return(totalRevenue*(exp(result)-1))
              }

          }
          )


setMethod(
 f= "calcMargins",
 signature= "AIDS",
 definition=function(object,preMerger=TRUE){

     priceDelta <- object@priceDelta
     ownerPre   <- object@ownerPre
     shares     <- calcShares(object,TRUE)

     elastPre <-  t(elast(object,TRUE))
     marginPre <-  -1 * as.vector(ginv(elastPre * ownerPre) %*% (shares * diag(ownerPre))) / shares

     if(preMerger){
         names(marginPre) <- object@labels
         return(marginPre)}

     else{

         marginPost <- 1 - ((1 + object@mcDelta) * (1 - marginPre) / (priceDelta + 1) )
         names(marginPost) <- object@labels
         return(marginPost)
     }

}
 )



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



setMethod(
 f= "calcPricesHypoMon",
 signature= "AIDS",
 definition=function(object,prodIndex,...){


     priceDeltaHM <- calcPriceDeltaHypoMon(object,prodIndex,...)

     prices <- object@prices[prodIndex] * (1 + priceDeltaHM)


     return(prices)

 })



setMethod(
 f= "show",
 signature= "AIDS",
 definition=function(object){

     print(object@priceDelta*100)

}
 )




setMethod(
 f= "summary",
 signature= "AIDS",
 definition=function(object,revenue=TRUE,parameters=FALSE,digits=2,...){


     curWidth <-  getOption("width")

     isParty <- as.numeric(rowSums( abs(object@ownerPost - object@ownerPre))>0)
     isParty <- factor(isParty,levels=0:1,labels=c(" ","*"))


     pricePre   <-  object@pricePre
     pricePost  <-  object@pricePost

     outPre  <-  calcShares(object,TRUE,revenue) * 100
     outPost <-  calcShares(object,FALSE,revenue) * 100

     mcDelta <- object@mcDelta

     outDelta <- (outPost/outPre - 1) * 100

     priceDelta   <-  object@priceDelta * 100



         results <- data.frame(priceDelta=priceDelta,sharesPre=outPre,
                               sharesPost=outPost,outputDelta=outDelta)



     if(!any(is.na(pricePre))){
         results <- cbind(pricePre=pricePre,pricePost=pricePost,results)
        }

     if(sum(abs(mcDelta))>0) results <- cbind(results,mcDelta=mcDelta)

     rownames(results) <- paste(isParty,object@labels)

     cat("\nMerger simulation results under '",class(object),"' demand:\n\n",sep="")


     options("width"=100) # this width ensures that everything gets printed on the same line
     print(round(results,digits),digits=digits)
     options("width"=curWidth) #restore to current width

     cat("\n\tNotes: '*' indicates merging parties' products. Deltas are percent changes.\n")
     if(revenue){cat("\tOutput is based on revenues.\n")}
     else{cat("\tOutput is based on units sold.\n")}
     results <- cbind(isParty, results)

     cat("\n\nShare-Weighted Price Change:",round(sum(outPost/100*priceDelta),digits),sep="\t")
     cat("\nShare-Weighted CMCR:",round(sum(cmcr(object)*outPost[isParty=="*"])/sum(outPost[isParty=="*"]),digits),sep="\t")
     if(!any(is.na(pricePre))){
         cat("\nCompensating Variation (CV):",round(CV(object,...),digits),sep="\t")
         }
     cat("\n\n")




     if(parameters){

         cat("\nAggregate Elasticity Estimate:",round(object@mktElast,digits),sep="\t")
         cat("\n\n")
         cat("\nDemand Parameter Estimates:\n\n")
         print(round(object@slopes,digits))
         cat("\n\n")


          if(.hasSlot(object,"intercepts") && all(!is.na(object@intercepts))){

              cat("\nIntercepts:\n\n")
              print(round(object@intercepts,digits))
              cat("\n\n")

             }
         if(hasMethod("getNestsParms",class(object))){
             cat("\nNesting Parameter Estimates:\n\n")
              print(round(getNestsParms(object),digits))

             cat("\n\n")
             }


     }

     rownames(results) <- object@labels
     return(invisible(results))


     })



aids <- function(shares,margins,prices,diversions,
                 ownerPre,ownerPost,
                 mcDelta=rep(0, length(shares)),
                 subset=rep(TRUE, length(shares)),
                 priceStart=runif(length(shares)),
                 isMax=FALSE,
                 labels=paste("Prod",1:length(shares),sep=""),
                 ...){



    if(missing(prices)){ prices <- rep(NA_real_,length(shares))}

    if(missing(diversions)){
        diversions <- tcrossprod(1/(1-shares),shares)
        diag(diversions) <- -1
    }



    ## Create AIDS container to store relevant data
    result <- new("AIDS",shares=shares,mcDelta=mcDelta,subset=subset,
                  margins=margins, prices=prices, quantities=shares,
                  ownerPre=ownerPre,ownerPost=ownerPost,
                  diversion=diversions,
                  priceStart=priceStart,labels=labels)

    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate Demand Slope Coefficients
    result <- calcSlopes(result)

    ## Solve Non-Linear System for Price Changes
    result@priceDelta <- calcPriceDelta(result,isMax=isMax,subset=subset,...)


    ## Calculate marginal cost
    result@mcPre <-  calcMC(result,TRUE)
    result@mcPost <- calcMC(result,FALSE)


    ## Calculate Pre and Post merger equilibrium prices
    result@pricePre  <- calcPrices(result,TRUE)
    result@pricePost <- calcPrices(result,FALSE)


    return(result)
}



