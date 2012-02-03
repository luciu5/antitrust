#setwd("h:/TaragIC/AntiTrustRPackage/antitrust/R")
#source("Classes.R")
#source("Methods.R")



setClass(
         Class = "AIDS",
         contains="Linear",
         representation=representation(
         priceDeltaStart  = "numeric",
         priceDelta       = "numeric",
         mktElast         = "numeric"
         ),
        prototype=prototype(
        priceDelta       =  numeric(),
        mktElast         =  numeric()
        ),

         validity=function(object){

             ## Sanity Checks

             nprods <- length(object@shares)
             if(nprods != length(object@priceDeltaStart)){
                 stop("'priceDeltaStart' must have the same length as 'shares'")}

             if(!isTRUE(all.equal(sum(object@shares),1))){
                 stop("The sum of 'shares' values must equal 1")}

             if(length(object@margins[!is.na(object@margins)])<2){
                 stop("'margins' must contain at least two non-missing margins in order to calibrate demand parameters")
             }

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

         marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% shares) / shares


         measure <- sum((margins - marginsCand)^2,na.rm=TRUE)

         return(measure)
     }

     ui <- matrix(c(-1/shares[k],0,1 - shares[k],-1),ncol=2,nrow=2)
     ci <- c(shares[k] - 1,1)

     parmsStart <- c(-1 - shares[k]*(1-shares[k]),-2) # These starting values always satisfy constraints

     minParam <- constrOptim(parmsStart,minD,ui=ui,ci=ci,grad=NULL)$par

     diagB <- minParam[1]*diversion[k,]/diversion[,k]
     B       <- -diversion*diagB # slopes symmetric, no need to transpose


     dimnames(B) <- list(object@labels,object@labels)

     object@slopes <- B

     if(abs(minParam[2])>5){warning("'mktElast' estimate is large.")}
     object@mktElast <- minParam[2]
     object@intercepts <- as.vector(shares - B%*%prices)
     names(object@intercepts) <- object@labels


     return(object)

 }

      )



setMethod(
 f= "calcPriceDelta",
 signature= "AIDS",
 definition=function(object,...){


     ## Calculate premerger margins
     marginPre <- calcMargins(object,TRUE)


     ##Define system of FOC as a function of priceDelta
     FOC <- function(priceDelta,object){

         sharePost <- object@shares +  as.vector(object@slopes %*% priceDelta)

         elastPost <- t(object@slopes/sharePost) + sharePost * (object@mktElast + 1) #Caution: returns TRANSPOSED elasticity matrix
         diag(elastPost) <- diag(elastPost) - 1

         marginPost <- 1 - ((1 + object@mcDelta) * (1 - marginPre) / exp(priceDelta))


         thisFOC <- sharePost + as.vector((elastPost*object@ownerPost) %*% (sharePost*marginPost))
         return(thisFOC)

     }




     ## Find price changes that set FOCs equal to 0
     minResult <- nleqslv(object@priceDeltaStart,FOC,object=object,...)

     if(minResult$termcd != 1){warning("'calcPriceDelta' nonlinear solver may not have successfully converged. 'nleqslv' reports: '",minResult$message,"'")}


     deltaPrice <- (exp(minResult$x)-1)
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
     else{ prices <- object@prices * (1 - object@priceDelta)}

     }
          )

setMethod(
 f= "calcShares",
 signature= "AIDS",
 definition=function(object,preMerger=TRUE){

     if(preMerger){return(object@shares)}

     else{

         shares <-  object@shares + as.vector(object@slopes %*% log(object@priceDelta + 1))
         return(shares)
     }

 }
 )




setMethod(
          f= "elast",
          signature= "AIDS",
          definition=function(object,preMerger=TRUE){


              shares <- calcShares(object,preMerger)

              elast <- t(object@slopes/shares) + shares * (object@mktElast + 1) #Caution: returns TRANSPOSED elasticity matrix
              diag(elast) <- diag(elast) - 1
              dimnames(elast) <-  list(object@labels,object@labels)

              return(t(elast))

          }
          )




setMethod(
 f= "diversion",
 signature= "AIDS",
 definition=function(object,preMerger=TRUE){

     diversion <- -t(object@slopes)/diag(object@slopes)
     dimnames(diversion) <-  list(object@labels,object@labels)


     return(diversion)

     }

          )

setMethod(
 f= "cmcr",
 signature= "AIDS",
 definition=function(object){

     sharesPre <- calcShares(object,TRUE)
     sharesPre <- tcrossprod(1/sharesPre,sharesPre)

     marginPre <- calcMargins(object,TRUE)


     elastPre  <- t(elast(object,TRUE))

     divPre    <- elastPre/diag(elastPre)


    Bpost      <- divPre * sharesPre * object@ownerPost
    marginPost <- -1 * as.vector(solve(Bpost) %*% (1/diag(elastPre)))

    cmcr <- (marginPost - marginPre)/(1 - marginPre)
    names(cmcr) <- object@labels

    return(cmcr * 100)
}
 )




setMethod(
          f= "CV",
          signature= "AIDS",
          definition=function(object){

              ## computes compensating variation using the closed-form solution found in LaFrance 2004, equation 10

              if(any(is.na(object@prices))){stop("Compensating Variation cannot be calculated without supplying values to 'prices'")}


              slopes <- object@slopes
              intercepts <- object@intercepts
              pricePre <- object@pricePre
              pricePost <- object@pricePost

              ePre <-   sum(intercepts * log(pricePre))    + .5 * as.vector(log(pricePre)    %*% slopes %*% log(pricePre) )
              ePost <-  sum(intercepts * log(pricePost)) + .5 * as.vector(log(pricePost) %*% slopes %*% log(pricePost) )

              return( exp(ePost) - exp(ePre)  )

          }
          )


setMethod(
 f= "calcMargins",
 signature= "AIDS",
 definition=function(object,preMerger=TRUE){

     priceDelta <- object@priceDelta
     shares     <- calcShares(object,TRUE)

     elastPre <-  t(elast(object,TRUE))
     marginPre <-  -1 * as.vector(solve(elastPre*object@ownerPre) %*% shares) / shares

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

     nprods <- length(prodIndex)
     priceDeltaOld <- object@priceDelta
     shares <- object@shares
     mktElast <- object@mktElast
     slopes <- object@slopes

      ## Calculate premerger margins
     marginPre <- calcMargins(object,TRUE)[prodIndex]


     ##Define system of FOC as a function of priceDelta
     FOC <- function(priceDelta){

         priceCand <- priceDeltaOld
         priceCand[prodIndex] <- priceDelta
         shareCand <- shares +  as.vector(slopes %*% priceCand)

         elastPost <- t(slopes/shareCand) + shareCand * (mktElast + 1) #Caution: returns TRANSPOSED elasticity matrix
         diag(elastPost) <- diag(elastPost) - 1
         elastPost <- elastPost[prodIndex,prodIndex]

         marginPost <- 1 - (1 - marginPre) / exp(priceDelta)

         shareCand <-  shareCand[prodIndex]
         thisFOC <- shareCand + as.vector(elastPost %*% (shareCand*marginPost))
         return(thisFOC)

     }




     ## Find price changes that set FOCs equal to 0
     minResult <- nleqslv(object@priceDeltaStart[prodIndex],FOC,...)

     if(minResult$termcd != 1){warning("'calcPriceDeltaHypoMon' nonlinear solver may not have successfully converged. 'nleqslv' reports: '",minResult$message,"'")}


     deltaPrice <- (exp(minResult$x)-1)

     return(deltaPrice)

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
 definition=function(object){


     curWidth <-  getOption("width")

     isParty <- as.numeric(colSums( object@ownerPost - object@ownerPre)>0)
     isParty <- factor(isParty,levels=0:1,labels=c(" ","*"))


     pricePre   <-  object@pricePre
     pricePost  <-  object@pricePost

     outPre  <-  calcShares(object,TRUE) * 100
     outPost <-  calcShares(object,FALSE) * 100

     mcDelta <- object@mcDelta

     outDelta <- (outPost/outPre - 1) * 100

     priceDelta   <-  object@priceDelta * 100



         results <- data.frame(priceDelta=priceDelta,outputPre=outPre,
                               outputPost=outPost,outputDelta=outDelta)



     if(!any(is.na(pricePre))){
         results <- cbind(pricePre=pricePre,pricePost=pricePost,results)
        }

     if(sum(abs(mcDelta))>0) results <- cbind(results,mcDelta=mcDelta)

     rownames(results) <- paste(isParty,object@labels)

      cat("\nMerger Simulation Results:\n\n")

     options("width"=100) # this width ensures that everything gets printed on the same line
     print(results,digits=2)
     options("width"=curWidth) #restore to current width

     cat("\n\tNotes: '*' indicates merging parties' products. Deltas are percent changes.\n")
     results <- cbind(isParty, results)

     cat("\n\nShare-Weighted Price Change:",round(sum(outPost/100*priceDelta),2),sep="\t")
     cat("\nShare-Weighted CMCR:",round(sum(cmcr(object)*outPost/100),2),sep="\t")


     cat("\n\nAggregate Elasticity Estimate:",round(object@mktElast,2),sep="\t")
     cat("\n\n")

     rownames(results) <- object@labels
     return(invisible(results))


     })





aids <- function(shares,margins,prices,diversions,
                 ownerPre,ownerPost,
                 mcDelta=rep(0, length(shares)),
                 priceDeltaStart=runif(length(shares)),
                 labels=paste("Prod",1:length(shares),sep=""),
                 ...){



    if(missing(prices)){ prices <- rep(NA,length(shares))}

    if(missing(diversions)){
        diversions <- tcrossprod(1/(1-shares),shares)
        diag(diversions) <- -1
    }



    ## Create AIDS container to store relevant data
    result <- new("AIDS",shares=shares,mcDelta=mcDelta
                  ,margins=margins, prices=prices, quantities=shares,
                  ownerPre=ownerPre,ownerPost=ownerPost,
                  diversion=diversions,
                  priceDeltaStart=priceDeltaStart,labels=labels)

    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate Demand Slope Coefficients
    result <- calcSlopes(result)

    ## Solve Non-Linear System for Price Changes
    result@priceDelta <- calcPriceDelta(result,...)


    ## Calculate Pre and Post merger equilibrium prices
    result@pricePre  <- calcPrices(result,TRUE)
    result@pricePost <- calcPrices(result,FALSE)

    ##Calculate constant marginal costs
    result@mc <- calcMC(result)



    return(result)
}



