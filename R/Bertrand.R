

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
  name= "calcProducerSurplus",
  def=function(object,...){standardGeneric("calcProducerSurplus")}
)
setGeneric (
  
  name= "calcProducerSurplusGrimTrigger",
  def=function(object,...){standardGeneric("calcProducerSurplusGrimTrigger")}
)
setGeneric (
 name= "HypoMonTest",
 def=function(object,...){standardGeneric("HypoMonTest")}
 )
#setGeneric (
# name= "calcSearchSets",
# def=function(object,...){standardGeneric("calcSearchSets")}
# )

#setGeneric (
# name= "isMax",
# def=function(object,...){standardGeneric("isMax")}
# )


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
         margins <-  -1 * as.vector(solve(t(elast)*owner) %*% (revenue * diag(owner))) / revenue


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




## compute margins
setMethod(
  f= "calcProducerSurplus",
  signature= "Bertrand",
  definition=function(object,preMerger=TRUE){
    
    
    mc     <- calcMC(object,preMerger)
    if( preMerger) {prices <- object@pricePre}
    else{prices <- object@pricePost}
    
    if(hasMethod("calcQuantities",class(object))){
      output <- calcQuantities(object,preMerger)
    }
    else{
      warning("'calcQuantities' method not defined for class ",class(object),". Using 'calcShares' instead")
      output <- calcShares(object,preMerger,revenue=FALSE)
    }
    
    ps <- (prices - mc) * output
    names(ps) <- object@labels
    
    return(ps)
  }
  
)
## Coordinated effects using Grim Trigger strategies
## and optimal deffection

setMethod(
  f= "calcProducerSurplusGrimTrigger",
  signature= "Bertrand",
  definition= function(object,coalition,discount,preMerger=TRUE){
    nprod <- length(object@labels)
    if(!is.numeric(coalition) ||
         length(coalition) > nprod ||
         !coalition %in% 1:nprod){
      stop ("'coalition' must be a vector of product indices no greater than the number of products")
    }
    if(any(discount<0 | discount>=1) ||
       length(discount) != length(coalition)){
      stop ("'discount' must be a vector of values between 0 and 1 and whose length equals that of 'coalition'")
    }
    
    psPunish <- calcProducerSurplus(object,preMerger)
    
    if(preMerger){
      ownerPre <- object@ownerPre
      object@ownerPre[coalition,coalition] <- 1
      ownerCoalition <- object@ownerPre
      object@pricePre <- calcPrices(object,preMerger)
      psCoord  <- calcProducerSurplus(object,preMerger)
      
      psDefect <- psPunish
      
      #test to see if c defects
      for(c in coalition){
        object@ownerPre <- ownerCoalition
        object@ownerPre[c,] <- ownerPre[c,]
        object@pricePre <- calcPrices(object,preMerger)
        psDefect[c] <- calcProducerSurplus(object,preMerger)[c]
      }
      
    }
    else{
      
      ownerPost <- object@ownerPost
      object@ownerPost[coalition,coalition] <- 1
      ownerCoalition <- object@ownerPost
      object@pricePost <- calcPrices(object,preMerger)
      psCoord  <- calcProducerSurplus(object,preMerger)
      
      psDefect <- psPunish
      
      for(c in coalition){
        object@ownerPost <- ownerCoalition
        object@ownerPost[c,] <- ownerPost[c,]
        object@pricePost <- calcPrices(object,preMerger)
        psDefect[c] <- calcProducerSurplus(object,preMerger)[c]
      }
    }

    result<-cbind(Discount=discount,Coord=psCoord,Defect=psDefect,Punish=psPunish)
    rownames(result) <- object@labels
    
    result <- result[coalition,]
    
    result <- cbind(result,IC=result[,"Coord"] > result[,"Defect"]*(1-discount) + result[,"Punish"]*discount)
    return(result)
  }
)



##plot method

setMethod(
  f= "plot",
  signature= "Bertrand",
  definition=function(x,scale=.1){
    object=x
    nprods=length(object@labels)
    pricePre=object@pricePre
    pricePost=object@pricePost
    mcPre=calcMC(object,preMerger=TRUE)
    mcPost=calcMC(object,preMerger=FALSE)
    labels=object@labels
    isParty <- rowSums( abs(object@ownerPost - object@ownerPre))>0
    isParty <- ifelse(isParty,"*","")
    labels  <- paste(isParty,labels,sep="")
    
    if(hasMethod("calcQuantities",class(object))){
      outPre=calcQuantities(object,preMerger=TRUE)
      outPost=calcQuantities(object,preMerger=FALSE)
    }
    else{
      outPre=calcShares(object,preMerger=TRUE)
      outPost=calcShares(object,preMerger=FALSE)
    }
    
    prices<-quantPre<-quantPost<-prod<-NULL
    
    plotThis <- function(price,idx,preMerger){
      thisObj=object
      if(preMerger){thisObj@pricePre[idx]=price}
      else{thisObj@pricePost[idx]=price}
      
      if(hasMethod("calcQuantities",class(thisObj))){
        return(calcQuantities(thisObj,preMerger=preMerger)[idx])
      }
      else{return(calcShares(thisObj,preMerger=preMerger)[idx])}
              }
    
    for(i in 1:nprods){
      thesePrices=seq((1-scale)*min(mcPre[i],mcPost[i]),(1+scale)*max(pricePre[i],pricePost[i]),length.out=100)
      quantPre=c(quantPre,sapply(thesePrices,plotThis,idx=i,preMerger=TRUE))
      quantPost=c(quantPost,sapply(thesePrices,plotThis,idx=i,preMerger=FALSE))
      prices=c(prices,thesePrices)
      prod=c(prod,rep(labels[i],length(thesePrices)))
    }
    
    resultPre=data.frame(output=quantPre,price=prices,prod=prod,Demand="pre-merger")
    resultPost=data.frame(output=quantPost,price=prices,prod=prod,Demand="post-merger")
    result=rbind(resultPre,resultPost)
    result=result[result$output>0,]
    
    equilibria=data.frame(output=c(outPre,
                                   outPost),
                          price=c(pricePre,pricePost),
                          mc=c(mcPre,mcPost),
                          prod=labels,
                          Demand=rep(c("pre-merger","post-merger"),each=nprods))
    equilibria$Cost=equilibria$Demand

    
    thisPlot=ggplot(result,(aes(output,price,color=Demand,group=Demand))) + geom_line() + theme_bw() + theme(legend.position="bottom", legend.direction="horizontal",legend.title=element_blank())
    thisPlot=thisPlot + facet_wrap(~prod,scales="free_x") 
    
    thisPlot=thisPlot + geom_vline(aes(xintercept = output,group=Demand,colour=Demand),linetype=3,equilibria[,c("output","Demand","prod")] )
    thisPlot=thisPlot + geom_hline(aes(yintercept = price,group=Demand,colour=Demand),linetype=3,equilibria[,c("price","Demand","prod")] )
    thisPlot=thisPlot + geom_point(aes(y=price,x=output,color=Demand,group=Demand),equilibria)
    
    thisPlot=thisPlot + geom_hline(aes(yintercept = mc,group=Cost,color=Cost),data=equilibria[,c("mc","Cost","prod")],show_guide=TRUE)
   
    #if(!isTRUE(all.equal(mcPre,mcPost))){
    #  thisPlot=thisPlot + geom_hline(aes(yintercept = mc), color="orange",data=data.frame(mc=mcPost[mcPost!=mcPre],prod=labels[mcPost!=mcPre]),show_guide=TRUE)
     
    #}
    
    
    return(thisPlot)
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


     isParty <- as.numeric(rowSums( abs(object@ownerPost - object@ownerPre))>0)
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
     thisUPP <- tryCatch(upp(object),error=function(e) FALSE)
     if(!is.logical(thisUPP)){
     cat("\nShare-Weighted Net UPP:",round(sum(thisUPP*sharesPost[isParty=="*"])/sum(sharesPost[isParty=="*"]),digits),sep="\t")}

     ##Only compute CV if prices  are supplied
     thisCV <- tryCatch(CV(object,...),error=function(e) FALSE)
     if(!is.logical(thisCV)){
     cat("\nCompensating Variation (CV):",round(thisCV,digits),sep="\t")}

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

     isParty <- rowSums( abs(object@ownerPost - object@ownerPre) ) > 0

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

             isParty     <- rowSums( abs(object@ownerPost  - object@ownerPre) ) > 0

             ownerPre    <- object@ownerPre[isParty,isParty,drop=FALSE]
             ownerPost   <- object@ownerPost[isParty,isParty,drop=FALSE] #this will typically be a matrix of ones

             is1         <- unique(ceiling(ownerPre))[1,] > 0 # identify which products belong to the first merging party listed
             is2         <- unique(ceiling(ownerPre))[2,] > 0 # identify which products belong to the 2nd   merging party listed

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


             D1 <- (solve(Bpre[is1,is1,drop=FALSE])   %*%
                    ((diag(ownerPre)[is1]/diag(ownerPost)[is1]) * as.vector(Bpost[is1,is2,drop=FALSE] %*%
                                                                            margin[is2]))) # owner 1 gross upp
             D2 <- (solve(Bpre[is2,is2,drop=FALSE]) %*%
                    ((diag(ownerPre)[is2]/diag(ownerPost)[is2]) * as.vector(Bpost[is2,is1,drop=FALSE] %*%
                                                                              margin[is1])))  # owner 2 gross upp

             upp[is1]  <- -as.vector(D1)
             upp[is2] <- -as.vector(D2)

             result[isParty] <- upp*price + mcDelta #net UPP

             names(result) <- object@labels

             return(result[isParty])

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

              isParty <- rowSums( abs(object@ownerPost - ownerPre) )>0 #identify which products belong to the merging parties

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
