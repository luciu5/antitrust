## Generate a bunch of generic functions

setGeneric (
 name= "ownerToMatrix",
 def=function(object,...){standardGeneric("ownerToMatrix")}
 )

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
name= "calcMargins",
 def=function(object,...){standardGeneric("calcMargins")}
 )

setGeneric (
 name= "calcPrices",
 def=function(object,...){standardGeneric("calcPrices")}
 )

setGeneric (
 name= "calcPriceDelta",
 def=function(object,...){standardGeneric("calcPriceDelta")}
 )
setGeneric (
 name= "calcPriceDeltaHypoMon",
 def=function(object,...){standardGeneric("calcPriceDeltaHypoMon")}
 )
setGeneric (
 name= "elast",
 def=function(object,...){standardGeneric("elast")}
 )

setGeneric (
 name= "defineMarkets",
 def=function(object,...){standardGeneric("defineMarkets")}
 )


setGeneric (
 name= "diversion",
 def=function(object,...){standardGeneric("diversion")}
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
 name= "cmcr",
 def=function(object,...){standardGeneric("cmcr")}
 )


setGeneric (
 name= "getNestsParms",
 def=function(object,...){standardGeneric("getNestsParms")}
 )

setGeneric (name= "summary")


## Create some methods for the Antitrust Class

## Method to compute price changes
setMethod(
 f= "calcPriceDelta",
 signature= "Antitrust",
 definition=function(object){

     pricePre  <- object@pricePre
     pricePost <- object@pricePost

     priceDelta <- pricePost/pricePre - 1
     names(priceDelta) <- object@labels

     return(priceDelta)

 }
)


## Method to compute HHI
setMethod(
 f= "hhi",
 signature= "Antitrust",
 definition=function(object,preMerger=TRUE){

     shares=calcShares(object,preMerger) *100

     return(sum(shares^2))

 }
)

## print method
setMethod(
 f= "show",
 signature= "Antitrust",
 definition=function(object){

    print(calcPriceDelta(object)*100)

}
 )


## compute margins
setMethod(
 f= "calcMargins",
 signature= "Antitrust",
 definition=function(object,preMerger=TRUE){

     priceDelta <- calcPriceDelta(object)
     shares     <- object@shares

     elastPre <-  t(elast(object,TRUE))
     marginPre <-  -1 * as.vector(solve(elastPre*object@ownerPre) %*% shares) / shares


     if(preMerger){
         names(marginPre) <- object@labels
         return(marginPre)
     }

     else{

         marginPost <- 1 - ((1 + object@mcDelta) * (1 - marginPre) / (priceDelta + 1) )
         names(marginPost) <- object@labels
         return(marginPost)
     }

}
 )


##summarize method

setMethod(
 f= "summary",
 signature= "Antitrust",
 definition=function(object){

     curWidth <-  getOption("width")

     if(hasMethod("calcQuantities",class(object))){
         outPre  <-  calcQuantities(object,TRUE)
         outPost <-  calcQuantities(object,FALSE)
     }

     else{
         outPre  <-  calcShares(object,TRUE) * 100
         outPost <-  calcShares(object,FALSE) * 100
     }

     mcDelta <- object@mcDelta

     outDelta <- (outPost/outPre - 1) * 100

     pricePre   <-  object@pricePre
     pricePost  <-  object@pricePost
     priceDelta <- (pricePost/pricePre - 1) * 100

     isParty <- as.numeric(colSums( object@ownerPost - object@ownerPre)>0)
     isParty <- factor(isParty,levels=0:1,labels=c(" ","*"))

     results <- data.frame(pricePre=pricePre,pricePost=pricePost,
                           priceDelta=priceDelta,outputPre=outPre,
                           outputPost=outPost,outputDelta=outDelta)

     if(sum(abs(mcDelta))>0) results <- cbind(results,mcDelta=mcDelta)


     rownames(results) <- paste(isParty,object@labels)

     sharesPost <- calcShares(object,FALSE)

     cat("\nMerger Simulation Results:\n\n")

     options("width"=100) # this width ensures that everything gets printed on the same line
     print(results,digits=2)
     options("width"=curWidth) #restore to current width

     cat("\n\tNotes: '*' indicates merging parties' products. Deltas are percent changes.\n")
     results <- cbind(isParty, results)

     cat("\n\nShare-Weighted Price Change:",round(sum(sharesPost*priceDelta),2),sep="\t")
     cat("\nShare-Weighted CMCR:",round(sum(cmcr(object)*sharesPost),2),sep="\t")
     cat("\n\n")

     rownames(results) <- object@labels
     return(invisible(results))

 })


## create ownership matrix
setMethod(
 f= "ownerToMatrix",
signature= "Antitrust",
definition=function(object,preMerger=TRUE){


    ## transform ownerPre and ownerPost vectors into matrices, when applicable

    if(preMerger) thisOwner <- object@ownerPre
    else{         thisOwner <- object@ownerPost}



    if(is.vector(thisOwner) || is.factor(thisOwner)){

        nprod <- length(object@shares)
        owners <- as.numeric(factor(thisOwner))
        thisOwner <- matrix(0,ncol=nprod,nrow=nprod)


        for( o in unique(owners)){
            thisOwner [owners == o, owners == o] = 1
        }


    }


    return(thisOwner)

}
)


## Method to compute diversion
setMethod(
 f= "diversion",
 signature= "Antitrust",
 definition=function(object,preMerger=TRUE){

     if(hasMethod("calcQuantities",class(object))){
         output  <-  calcQuantities(object,preMerger)
     }

     else{
         output  <-  calcShares(object,preMerger)
     }

    elasticity <- elast(object,preMerger)

    diversion <- -1 * t(elasticity) / diag(elasticity)
    diag(diversion) <- 1
    diversion <- diversion * tcrossprod(1/output,output)

    return(diversion)
}
 )

##Method to compute Compensating Marginal Cost Reduction
setMethod(
 f= "cmcr",
 signature= "Antitrust",
 definition=function(object){

    pricePre <- object@pricePre
    pricePre <- tcrossprod(1/pricePre,pricePre)

    elastPre <- elast(object,TRUE)
    divPre <- diversion(object,TRUE)

    Bpre =  -1 * divPre * pricePre * object@ownerPre;  diag(Bpre) = 1
    Bpost = -1 * divPre * pricePre * object@ownerPost; diag(Bpost) = 1

    marginPre <- -1 * as.vector(solve(Bpre)  %*% (1/diag(elastPre)))
    marginPost <-     as.vector(solve(Bpost) %*%  Bpre %*% marginPre)

    cmcr <- (marginPost - marginPre)/(1 - marginPre)
    names(cmcr) <- object@labels

    return(cmcr * 100)
}
 )


## Use the Hypothetical Monopolist Test to create the set of candidate markets
## that satisfy a SSNIP

setMethod(
 f= "defineMarkets",
 signature= "Antitrust",
 definition=function(object,startProd=NULL,ssnip=.05,supplyResponse=FALSE,...){

     ownerPre <- object@ownerPre
     nprods   <- ncol(ownerPre)
     isParty <- colSums( object@ownerPost - object@ownerPre)>0 #identify which products belong to the merging parties

     if(length(ssnip)>1 || ssnip<0 | ssnip>1 ){stop("'ssnip' must be a number between 0 and 1")}


     if(is.null(startProd)){
         startProd <- which(isParty)[1]
         warning("'startProd' not specified. Setting equal to ",startProd)
     }

     else if(!startProd %in% which(isParty)){stop("'startProd' must equal the position number of a merging parties' product")}

     startDiversion <- diversion(object,TRUE)[startProd,]
     closestSubs <- order(startDiversion,decreasing = TRUE)

     result <- matrix(NA,ncol=nprods+1,nrow=nprods)


     for( i in 1:length(closestSubs)){
         candMonopolist <- closestSubs[1:i]

         if(supplyResponse){
             candOwner <- ownerPre

             ##create ownership structure for
             ## products owned by hypothetical monopolist

             candOwner[candMonopolist,] <- candOwner[,candMonopolist] <- 0
             candOwner[candMonopolist,candMonopolist] <- 1
             object@ownerPost <- candOwner

             if(hasMethod("calcPrices",class(object))){
                 object@pricePost <- rep(NA,nprods) #allows priceDelta to equal NA if next line fails
                 try(object@pricePost <- calcPrices(object,FALSE,...),silent=TRUE)
             }

             priceDelta <-  calcPriceDelta(object)[candMonopolist]

         }

         else{
             priceDelta <- rep(NA,length(candMonopolist))  #allows priceDelta to equal NA if next line fails
             try(priceDelta <-  calcPriceDeltaHypoMon(object,candMonopolist,...),silent=TRUE)}

         result[i,candMonopolist] <- 1
         result[i,nprods+1] <- max(priceDelta[which(isParty[candMonopolist])])


     }

     colnames(result) <- c(object@labels,"MaxPriceDelta")
     rownames(result) <- 1:nprods
     result[,"MaxPriceDelta"] <- result[,"MaxPriceDelta"]*100

     missingResults <- is.na(result[,nprods+1])

     if(any(missingResults)){
         warning("Price equilibria for some candidate markets could not be calculated.")
     }

     #cat("\nCandidate Markets:\n\n")

     result <- result[!missingResults & result[,nprods+1] >= ssnip*100,]
     #print(round(result))

     #cat(paste("\n\tNotes: 'MaxPriceDelta' is in percent changes.\n"))



     return(result)
 }
          )
