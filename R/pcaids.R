#setwd("h:/TaragIC/AntiTrustRPackage/antitrust/R")
#source("Classes.R")
#source("Methods.R")



setClass(
         Class = "PCAIDS",
         contains="Antitrust",
         representation=representation(

         knownElast      = "numeric",
         knownElastIndex = "numeric",
         mktElast        = "numeric",
         priceDeltaStart = "vector",
         priceDelta      = "vector",
         diversion       = "anyMatrix"
         ),
        prototype=prototype(

        priceDelta      =  vector()

        ),

         validity=function(object){

             ## Sanity Checks


             nprods <- length(object@shares)

             if(!(object@knownElastIndex %in% seq(1,nprods)) ){
                 stop("'knownElastIndex' value must be between 1 and the length of 'shares'")}
             if(nprods != length(object@mcDelta)){
                 stop("'mcDelta' must have the same length as 'shares'")}
             if(nprods != length(object@priceDeltaStart)){
                 stop("'priceDeltaStart' must have the same length as 'shares'")}
             if(object@knownElast>0 || object@mktElast > 0 ){
                 stop("All elasticities must be negative")}
             if(abs(object@knownElast) <= abs(object@mktElast)){
                 stop("'mktElast' must be less than 'knownElast' in absolute value")}

             if(!all(diag(object@diversion) == 1)){ stop("'diversion' diagonal elements must all equal 1")}
             if(any(abs(object@diversion) > 1)){ stop("'diversion' elements must be between -1 and 1")}

             if(nprods != nrow(object@diversion) ||
                nprods != ncol(object@diversion)){
                 stop("'diversions' must be a square matrix")
             }

         }

             )








setMethod(
 f= "calcSlopes",
 signature= "PCAIDS",
 definition=function(object){


     ## Uncover linear demand slopes from shares, knownElast and mktElast
     ## Since demand is linear IN LOG PRICE, model assumes that slopes remain unchanged following merger

     shares    <- object@shares
     diversion <- object@diversion
     labels    <- object@labels

     diag(diversion) <- -1
     diversion <- -1*t(diversion)

     idx        <- object@knownElastIndex
     shareKnown <- shares[idx]

     b <- rep(1,length(shares))

     bKnown <- shareKnown * (object@knownElast + 1 - shareKnown * (object@mktElast + 1))

     b[-idx] <- qr.solve(diversion[,-idx], -1*diversion[,idx])
     b <- b * bKnown
     B <- t(t(diversion) * b)
     dimnames(B) <- list(labels,labels)

     return(B)
 }

      )




setMethod(
 f= "calcPriceDelta",
 signature= "PCAIDS",
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
 f= "calcShares",
 signature= "PCAIDS",
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
 signature= "PCAIDS",
 definition=function(object,preMerger=TRUE){


     shares <- calcShares(object,preMerger)

     elast <- t(object@slopes/shares) + shares * (object@mktElast + 1) #Caution: returns TRANSPOSED elasticity matrix
     diag(elast) <- diag(elast) - 1
     dimnames(elast) <-  list(object@labels,object@labels)

     return(t(elast))

}
 )






setMethod(
 f= "cmcr",
 signature= "PCAIDS",
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
          signature= "PCAIDS",
          definition=function(object,prices=NULL){

              ## computes compensating variation using the closed-form solution found in LaFrance 2004, equation 10

              if(is.null(prices) ||
                 length(prices) != length(object@shares) ||
                 any(prices<=0 || is.na(prices)) ){
                  stop("'prices' must be a vector of positive numbers whose length equals 'shares'")}


              slopes <- object@slopes

              ## The following test should be never be flagged (it is true by construction in PCAIDS)
               if(!isTRUE(all.equal(slopes[upper.tri(slopes)],slopes[upper.tri(t(slopes))]))){
                  stop("log-price coefficient matrix must be symmetric in order to calculate compensating variation")
              }


              shares <- object@shares
              priceDelta <- 1+object@priceDelta
              intercept <- as.vector(shares - slopes %*% log(prices))
              pricePost <-  prices*priceDelta

              ePre <-   sum(intercept * log(prices))    + .5 * as.vector(log(prices)    %*% slopes %*% log(prices) )
              ePost <-  sum(intercept * log(pricePost)) + .5 * as.vector(log(pricePost) %*% slopes %*% log(pricePost) )

              return( exp(ePost) - exp(ePre)  )

          }
          )


setMethod(
 f= "calcMargins",
 signature= "PCAIDS",
 definition=function(object,preMerger=TRUE){

     priceDelta <- object@priceDelta
     shares     <- object@shares

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
 f= "show",
 signature= "PCAIDS",
 definition=function(object){

     print(object@priceDelta*100)

}
 )


setMethod(
 f= "summary",
 signature= "PCAIDS",
 definition=function(object){



         outPre  <-  calcShares(object,TRUE) * 100
         outPost <-  calcShares(object,FALSE) * 100



         outDelta <- (outPost/outPre - 1) * 100

         priceDelta   <-  object@priceDelta * 100


         results <- data.frame(priceDelta=priceDelta,outputPre=outPre,
                               outputPost=outPost,outputDelta=outDelta)

         rownames(results) <- object@labels

         cat("\nMerger Simulation Results (Deltas are Percent Changes):\n\n")
         print(round(results,2))

     cat("\nMerger Simulation Results (Deltas are Percent Changes):\n\n")
     print(round(results,2))
     cat("\n\nShare-Weighted Price Change:",round(sum(outPost*priceDelta/100),2),sep="\t")
     cat("\nShare-Weighted CMCR:",round(sum(cmcr(object)*outPost/100),2),sep="\t")
     cat("\n\n")
         return(invisible(results))

     })


setMethod(
 f= "calcPriceDeltaHypoMon",
 signature= "PCAIDS",
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



pcaids <- function(shares,knownElast,mktElast=-1,diversions=NULL,
                   ownerPre,ownerPost,
                   knownElastIndex=1,
                   mcDelta=rep(0, length(shares)),
                   priceDeltaStart=runif(length(shares)),
                   labels=paste("Prod",1:length(shares),sep=""),
                   ...){



    if(is.null(diversions)){
        diversions <- tcrossprod(1/(1-shares),shares)
        diag(diversions) <- 1
    }

  ## Create PCAIDS container to store relevant data
    result <- new("PCAIDS",shares=shares,mcDelta=mcDelta
                  ,knownElast=knownElast,mktElast=mktElast,
                  ownerPre=ownerPre,ownerPost=ownerPost,
                  knownElastIndex=knownElastIndex,
                  diversion=diversions,
                  priceDeltaStart=priceDeltaStart,labels=labels)

    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate Demand Slope Coefficients
    result@slopes <- calcSlopes(result)

    ## Solve Non-Linear System for Price Changes
    result@priceDelta <- calcPriceDelta(result,...)

    return(result)
}



