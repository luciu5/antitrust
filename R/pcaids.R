#setwd("h:/TaragIC/AntiTrustRPackage/antitrust/R")
#source("Classes.R")
#source("Methods.R")


## Check if Consumer Surplus Calculation is correct


setClass(
         Class = "PCAIDS",
         contains="Antitrust",
         representation=representation(

         knownElast      = "numeric",
         knownElastIndex = "numeric",
         mktElast        = "numeric",
         priceDeltaStart = "vector",
         priceDelta      = "vector"

         ),
        prototype=prototype(

        priceDelta      =  vector()

        ),

         validity=function(object){

             ## Sanity Checks


             nprods <- length(object@shares)

             if(!(object@knownElastIndex %in% seq(1,nprods)) ){
                 stop("'knownElastIndex' values must be between 1 and the length of 'shares'")}
             if(nprods != length(object@mcDelta)){
                 stop("'mcDelta' must have the same length as 'shares'")}
             if(nprods != length(object@priceDeltaStart)){
                 stop("'priceDeltaStart' must have the same length as 'shares'")}
             if(object@knownElast>0 || object@mktElast > 0 ){
                 stop("All elasticities must be negative")}
             if(abs(object@knownElast) <= abs(object@mktElast)){
                 stop("'mktElast' must be less than 'knownElast' in absolute value")}


         }

             )








setMethod(
 f= "calcSlopes",
 signature= "PCAIDS",
 definition=function(object){


     ## Uncover linear demand slopes from shares, knownElast and mktElast
     ## Since demand is linear IN LOG PRICE, model assumes that slopes remain unchanged following merger



     shareKnown <- object@shares[object@knownElastIndex]


     bKnown <- shareKnown * (object@knownElast + 1 - shareKnown * (object@mktElast + 1))
     beta <- ((1-object@shares) / (1-shareKnown)) * (object@shares / shareKnown)
     b <- beta * bKnown
     B <- tcrossprod(-object@shares, 1/(1-object@shares))
     diag(B) <- 1
     B <- t(t(B) * b)
     dimnames(B) <- list(object@labels,object@labels)

     return(B)
 }

      )


setGeneric (
 name= "calcPriceDelta",
 def=function(object,...){standardGeneric("calcPriceDelta")}
 )

setMethod(
 f= "calcPriceDelta",
 signature= "PCAIDS",
 definition=function(object,...){

     require(BB) #needed to solve nonlinear system of firm FOC

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
     minResult <- BBsolve(object@priceDeltaStart,FOC,object=object,...)

     deltaPrice <- (exp(minResult$par)-1)
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
 f= "calcMargins",
 signature= "PCAIDS",
 definition=function(object,preMerger=TRUE){

     elastPre <-  t(elast(object,TRUE))
     marginPre <-  -1 * as.vector(solve(elastPre*object@ownerPre) %*% object@shares) / object@shares

     if(preMerger){return(marginPre)}

     else{

         marginPost <- 1 - ((1 + object@mcDelta) * (1 - marginPre) / (object@priceDelta + 1) )
         return(marginPost)
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
 f= "diversion",
 signature= "PCAIDS",
 definition=function(object,preMerger=TRUE){

    elasticity <- elast(object,preMerger)
    shares <-  calcShares(object,preMerger)

    diversion <- -1 * t(elasticity) / diag(elasticity)
    diag(diversion) <- 1
    diversion <- diversion * tcrossprod(1/shares,shares)

    return(diversion)
}
 )








#setMethod(
# f= "deltaCS",
# signature= "PCAIDS",
# definition=function(object){


#    deltaLogPrice <-  log(object@priceDelta + 1)
#    sharesPre <- calcShares(object,TRUE)
#    sharesPost <- calcShares(object,FALSE)
#    deltaShares <- sharesPost - sharesPre

#     deltaCS <- .5 * (deltaLogPrice * 100) * (deltaShares * 100)
#     names(deltaCS) <- object@labels

#     return(deltaCS)

#}
# )




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


     sharesPre   <-  calcShares(object,TRUE)
     sharesPost  <-  calcShares(object,FALSE)
     sharesDelta <- (sharesPost - sharesPre)/sharesPre * 100

     results <- data.frame(
                           priceDelta=object@priceDelta * 100,sharesPre=sharesPre,
                           sharesPost=sharesPost,sharesDelta=sharesDelta)

     rownames(results) <- object@labels

     cat("\nMerger Simulation Results (Deltas are Percent Changes):\n\n")
     print(round(results,2))

    ## cat("\n\nPercent Change in Consumer Surplus:\n\n")
    ## print(round(deltaCS(object),2))

}
 )









pcaids <- function(shares,knownElast,mktElast=-1,ownerPre,ownerPost,
                   knownElastIndex=1,
                   mcDelta=rep(0, length(shares)),
                   priceDeltaStart=runif(length(shares)),
                   labels=paste("Prod",1:length(shares),sep=""),
                   ...){



  ## Create PCAIDS container to store relevant data
    result <- new("PCAIDS",shares=shares,mcDelta=mcDelta
                  ,knownElast=knownElast,mktElast=mktElast,
                  ownerPre=ownerPre,ownerPost=ownerPost,knownElastIndex=knownElastIndex,
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



