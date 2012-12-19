

setClass(

         Class = "Bertrand",
         contains="Antitrust",
         representation=representation(
         shares       = "numeric",
         mcDelta      = "numeric",
         slopes       = "matrixOrList"
         ),
         prototype=prototype(

         slopes          = matrix()

         ),
         validity=function(object){



             nprods <- length(object@labels)

             if(nprods != length(object@shares)){
                 stop("'labels' must have the same length as 'shares'")}

             if(any(object@shares < 0 | object@shares > 1,na.rm=TRUE)){
                 stop("'shares' values must be between 0 and 1")}

             if(sum(object@shares) > 1){
                 stop("The sum of 'shares' values must be less than or equal to 1")}


             ## if(any(object@mcDelta < 0 | object@mcDelta > 1,na.rm=TRUE)){
             ##    stop("'mcDelta' values must be between 0 and 1")}

             if(nprods != length(object@mcDelta)){
                 stop("'mcDelta' must have the same length as 'shares'")}

             if(any(object@mcDelta>0,na.rm=TRUE)){
                 warning("positive values of 'mcDelta' imply an INCREASE in marginal costs")}
             if(any(abs(object@mcDelta)>1,na.rm=TRUE)){
                 warning("Values of 'mcDelta' greater than 1 in absolute value imply a marginal cost change greater than 100%")}

             return(TRUE)

         }

         )


##
## Bertrand Methods
##

## Generate a bunch of generic functions


setGeneric (
 name= "calcSlopes",
 def=function(object){standardGeneric("calcSlopes")}
 )

setGeneric (
 name= "calcQuantities",
 def=function(object,...){standardGeneric("calcQuantities")}
 )

setGeneric (
 name= "calcShares",
 def=function(object,...){standardGeneric("calcShares")}
 )

setGeneric (
name= "calcMC",
 def=function(object,...){standardGeneric("calcMC")}
 )

setGeneric (
name= "calcMargins",
 def=function(object,...){standardGeneric("calcMargins")}
 )


setGeneric (
 name= "calcPricesHypoMon",
 def=function(object,...){standardGeneric("calcPricesHypoMon")}
 )
setGeneric (
 name= "calcPriceDeltaHypoMon",
 def=function(object,...){standardGeneric("calcPriceDeltaHypoMon")}
 )
setGeneric (
 name= "HypoMonTest",
 def=function(object,...){standardGeneric("HypoMonTest")}
 )
#setGeneric (
# name= "calcSearchSets",
# def=function(object,...){standardGeneric("calcSearchSets")}
# )

setGeneric (
 name= "isMax",
 def=function(object,...){standardGeneric("isMax")}
 )

setGeneric (
 name= "elast",
 def=function(object,...){standardGeneric("elast")}
 )


setGeneric (
 name= "diversion",
 def=function(object,...){standardGeneric("diversion")}
 )

setGeneric (
 name= "diversionHypoMon",
 def=function(object,...){standardGeneric("diversionHypoMon")}
 )

setGeneric (
 name= "CV",
 def=function(object,...){standardGeneric("CV")}
 )

setGeneric (
 name= "hhi",
 def=function(object,...){standardGeneric("hhi")}
 )

setGeneric (
 name= "upp",
 def=function(object,...){standardGeneric("upp")}
 )


setGeneric (
 name= "cmcr",
 def=function(object,...){standardGeneric("cmcr")}
 )


setGeneric (
 name= "getNestsParms",
 def=function(object,...){standardGeneric("getNestsParms")}
 )

setGeneric (name= "summary")


## Create some methods for the Bertrand Class




## Method to compute HHI
setMethod(
 f= "hhi",
 signature= "Bertrand",
 definition=function(object,preMerger=TRUE,revenue=FALSE){

     if(preMerger){owner <- object@ownerPre}
     else{owner <- object@ownerPost}

     control <- owner>0              #assumes that a firm can set prices
                                     #on products over which it has partial ownership

     weights <- crossprod(control,owner)
     weights <- t(t(weights)/diag(weights)) # divide each element by its corresponding diagonal

     shares <- calcShares(object,preMerger,revenue) *100

     result <- as.vector(shares %*% weights %*% shares)




     return(result)



 }
)



## compute margins
setMethod(
 f= "calcMargins",
 signature= "Bertrand",
 definition=function(object,preMerger=TRUE){



     if( preMerger) {

         owner  <- object@ownerPre
         revenue<- calcShares(object,preMerger,revenue=TRUE)

         elast <-  elast(object,preMerger)
         margins <-  -1 * as.vector(solve(t(elast*owner)) %*% (revenue * diag(owner))) / revenue


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




## Create a method to recover marginal cost using
## demand parameters and supplied prices
setMethod(
          f= "calcMC",
          signature= "Bertrand",
          definition= function(object,preMerger=TRUE){

              object@pricePre <- object@prices


              marginPre <- calcMargins(object,TRUE)

              mc <- (1 - marginPre) * object@prices

              if(!preMerger){
                  mc <- mc*(1+object@mcDelta)
              }

              names(mc) <- object@labels

              return(as.vector(mc))
          }
          )




##summarize method

setMethod(
 f= "summary",
 signature= "Bertrand",
 definition=function(object,revenue=TRUE,shares=TRUE,parameters=FALSE,digits=2,...){

     curWidth <-  getOption("width")


     pricePre   <-  object@pricePre
     pricePost  <-  object@pricePost
     priceDelta <- (pricePost/pricePre - 1) * 100

     if(!shares && hasMethod("calcQuantities",class(object))){
         outPre  <-  calcQuantities(object,TRUE)
         outPost <-  calcQuantities(object,FALSE)

         if(revenue){
             outPre <- pricePre*outPre
             outPost <- pricePost*outPost
         }

         sumlabels=paste("quantity",c("Pre","Post"),sep="")
     }

     else{
         if(!shares){warning("'shares' equals FALSE but 'calcQuantities' not defined. Reporting shares instead of quantities")}

         outPre  <-  calcShares(object,TRUE,revenue) * 100
         outPost <-  calcShares(object,FALSE,revenue) * 100

         sumlabels=paste("shares",c("Pre","Post"),sep="")
     }

     mcDelta <- object@mcDelta

     outDelta <- (outPost/outPre - 1) * 100


     isParty <- as.numeric(colSums( object@ownerPost - object@ownerPre)>0)
     isParty <- factor(isParty,levels=0:1,labels=c(" ","*"))

     results <- data.frame(pricePre=pricePre,pricePost=pricePost,
                           priceDelta=priceDelta,outputPre=outPre,
                           outputPost=outPost,outputDelta=outDelta)

     colnames(results)[colnames(results) %in% c("outputPre","outputPost")] <- sumlabels

     if(sum(abs(mcDelta))>0) results <- cbind(results,mcDelta=mcDelta)


     rownames(results) <- paste(isParty,object@labels)

     sharesPost <- calcShares(object,FALSE,revenue)

     cat("\nMerger simulation results under '",class(object),"' demand:\n\n",sep="")

     options("width"=100) # this width ensures that everything gets printed on the same line
     print(round(results,digits),digits=digits)
     options("width"=curWidth) #restore to current width

     cat("\n\tNotes: '*' indicates merging parties' products. Deltas are percent changes.\n")
     if(revenue){cat("\tOutput is based on revenues.\n")}
     else{cat("\tOutput is based on units sold.\n")}

     results <- cbind(isParty, results)

     cat("\n\nShare-Weighted Price Change:",round(sum(sharesPost*priceDelta),digits),sep="\t")
     cat("\nShare-Weighted CMCR:",round(sum(cmcr(object)*sharesPost[isParty=="*"])/sum(sharesPost[isParty=="*"]),digits),sep="\t")

     ##Only compute upp if prices are supplied
     uppExists <- tryCatch(upp(object),error=function(e) FALSE)
     if(is.logical(uppExists)){
     cat("\nShare-Weighted Net UPP:",round(sum(upp(object)*sharesPost),digits),sep="\t")}

     ##Only compute CV if prices  are supplied
     cvExists <- tryCatch(CV(object,...),error=function(e) FALSE)
     if(cvExists){
     cat("\nCompensating Variation (CV):",round(CV(object,...),digits),sep="\t")}

     cat("\n\n")


     if(parameters){

         cat("\nDemand Parameter Estimates:\n\n")
         if(is.list(object@slopes)){
             print(lapply(object@slopes,round,digits=digits))
             }
         else{
         print(round(object@slopes,digits))
         }
         cat("\n\n")

         if(.hasSlot(object,"intercepts")){

             cat("\nIntercepts:\n\n")
             print(round(object@intercepts,digits))
             cat("\n\n")

             }

             if(.hasSlot(object,"constraint") && object@constraint){cat("\nNote: (non-singleton) nesting parameters are constrained to be equal")}
             cat("\n\n")

         }

     rownames(results) <- object@labels
     return(invisible(results))

 })



## Method to compute diversion
setMethod(
 f= "diversion",
 signature= "Bertrand",
 definition=function(object,preMerger=TRUE,revenue=FALSE){


     output  <-  calcShares(object,preMerger,revenue)

     elasticity <- elast(object,preMerger)



     if(revenue){
         diversion <- t(elasticity) / (diag(elasticity)-1)
         diversion <- -1 * diversion / diag(diversion)
     }

     else{
         diversion <- -1 * t(elasticity) / diag(elasticity)
     }

     diversion <- diversion * tcrossprod(1/output,output)
     dimnames(diversion) <-  list(object@labels,object@labels)

     return(diversion)
}
 )

##Method to compute Compensating Marginal Cost Reduction
setMethod(
 f= "cmcr",
 signature= "Bertrand",
 definition=function(object){

     isParty <- colSums( object@ownerPost - object@ownerPre) > 0

     ##Compute pre-merger margins
     marginPre  <- calcMargins(object,TRUE)


     ##compute post-merger margins evaluated at pre-merger prices
     object@ownerPre <- object@ownerPost
     marginPost <- calcMargins(object,TRUE)

     cmcr <- (marginPost - marginPre)/(1 - marginPre)
     names(cmcr) <- object@labels

     cmcr <- cmcr[isParty]
     return(cmcr * 100)
}
 )

##Method to compute upp
setMethod(
          f= "upp",
          signature= "Bertrand",
         definition=function(object){

             isParty     <- colSums( object@ownerPost - object@ownerPre) > 0

             ownerPre    <- object@ownerPre[isParty,isParty]
             ownerPost   <- object@ownerPost[isParty,isParty] #this will typically be a matrix of ones

             is1         <- ownerPre[1,] > 0 # identify which products belong to the first merging party listed

             upp         <- rep(NA,length(is1))
             result      <- rep(0,length(isParty))

             div         <- diversion(object,preMerger=TRUE)[isParty,isParty]
             price       <- object@pricePre[isParty]

             mcPre       <- calcMC(object,preMerger=TRUE)[isParty]
             mcPost      <- calcMC(object,preMerger=FALSE)[isParty]
             mcDelta     <- mcPost - mcPre

             margin      <- 1 - mcPre/price


             ## weight diversion ratios by price ratios and ownership matrices ##
             priceRatio = tcrossprod(1/price, price)
             Bpre  =  -1 * div * priceRatio * ownerPre
             Bpost =  -1 * div * priceRatio * ownerPost


             D1 <- (solve(Bpre[is1,is1,drop=FALSE])   %*%  Bpost[is1,!is1,drop=FALSE]) %*% margin[!is1] # owner 1 gross upp
             D2 <- (solve(Bpre[!is1,!is1,drop=FALSE]) %*%  Bpost[!is1,is1,drop=FALSE]) %*% margin[is1]  # owner 2 gross upp

             upp[is1]  <- -as.vector(D1)
             upp[!is1] <- -as.vector(D2)

             result[isParty] <- upp*price + mcDelta #net UPP

             names(result) <- object@labels

             return(result)

         }
          )




## Use the Hypothetical Monopolist Test to determine whether a candidate market satisfies a SSNIP.
setMethod(
 f= "calcPriceDeltaHypoMon",
 signature= "Bertrand",
          definition=function(object,prodIndex,...){


              pricesHM <-  calcPricesHypoMon(object,prodIndex,...)

              pricesDelta <- pricesHM/object@pricePre[prodIndex] - 1

              return(pricesDelta)

          })

setMethod(
 f= "HypoMonTest",
 signature= "Bertrand",
          definition=function(object,prodIndex,ssnip=.05,...){

              ownerPre <- object@ownerPre
              nprods   <- ncol(ownerPre)
              pricesDelta <- rep(0,nprods)

              if(missing(prodIndex) || any(prodIndex>nprods | prodIndex <1 ) ){
                  stop("'prodIndex' must be a vector of product indices between 1 and ",nprods)
              }

              if(length(ssnip)>1 || ssnip<0 | ssnip>1 ){stop("'ssnip' must be a number between 0 and 1")}

              isParty <- colSums( object@ownerPost - ownerPre)>0 #identify which products belong to the merging parties

              if(length(intersect(which(isParty),prodIndex))==0){
                  stop("'prodIndex' does not contain any of the merging parties' products. Add at least one of the following indices: ",
                       paste(which(isParty),collapse=","))
                  }



              pricesDelta[prodIndex] <-  calcPriceDeltaHypoMon(object,prodIndex,...)


              result <- max(pricesDelta[isParty]) > ssnip

              return( result)
          }

              )


setMethod(
 f= "diversionHypoMon",
 signature= "Bertrand",
          definition=function(object,prodIndex,...){

              object@pricePre[prodIndex] <- calcPricesHypoMon(object,prodIndex,...)

              return(diversion(object,preMerger=TRUE,revenue=TRUE))



              }
          )
