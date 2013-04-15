##Consider switching the order in which parameters are solved for: put distribution parms first


setClass(
         Class   = "Auction2ndCap",
         contains="Antitrust",
         representation=representation(
         capacities       = "numeric",
         margins          = "numeric",
         prices           = "numeric",
         reserve          = "numeric",
         shareInside      = "numeric",
         sellerCostCDF    = "character",
         sellerCostCDFLowerTail    = "logical",
         sellerCostPDF    = "function",
         sellerCostBounds = "numeric",
         sellerCostParms  = "numeric",
	       buyerValuation        = "numeric", #was buyerValuation
         reservePre       = "numeric",
         reservePost      = "numeric",
         mcDelta          = "numeric",
         parmsStart       = "numeric"

         ),
         prototype=prototype(
         reservePre      =  numeric(),
         reservePost     =  numeric(),
         buyerValuation       =  numeric(),
         sellerCostParms =  numeric(),
         sellerCostCDFLowerTail    = TRUE


          ),
         validity=function(object){

             nprods <- length(object@labels)
             #cdf    <- as.character(quote(object@sellerCostCDF))
             #print(cdf)
             cdf    <- object@sellerCostCDF
             if(is.na(object@reserve)){parmsStart <- object@parmsStart[-1]}
             else{parmsStart <- object@parmsStart}

             if(nprods != length(object@capacities)){
                 stop("'labels' must have the same length as 'capacities'")}


             if( cdf=="punif"){
                if(length(parmsStart)!=2 ||
                    parmsStart[1] >= parmsStart[2]){
                  if(is.na(object@reserve)){stop("The first element of parmsStart must be a starting value for 'reserve'")}
                  stop("For the Uniform distribution, 'parmsStart' must be a numeric vector of length ", 2 + is.na(object@reserve),
                        "whose final element must be greater than the next-to-last element")
                }
             }

             else if( cdf=="pexp"){
                if(length(parmsStart)!=1 ||
                    parmsStart[1] <=0){
                    if(is.na(object@reserve)){stop("The first element of parmsStart must be a starting value for 'reserve'")}
                    stop("For the Exponential distribution, 'parmsStart' must be a numeric vector of length ",
                         1 + is.na(object@reserve),"whose final element is greater than 0")
                 }
             }
             else if( cdf=="pweibull"){
                if(length(parmsStart)!=2 ||
                     parmsStart[1] <=0  ||
                     parmsStart[2] <=0 ){
                  if(is.na(object@reserve)){stop("The first element of parmsStart must be a starting value for 'reserve'")}
                    stop("For the Weibull distribution, 'parmsStart' must be a numeric vector of length ",
                         2 + is.na(object@reserve)," whose final 2 elements must be greater than 0")
                        }
                      }
             else if( cdf=="pgumbel"){
                if(length(parmsStart)!=2 ||
                     parmsStart[2] <=0){
                  if(is.na(object@reserve)){stop("The first element of parmsStart must be a starting value for 'reserve'")}
                    stop("For the Gumbel distribution, 'parmsStart' must be a numeric vector of length ",
                         2 + is.na(object@reserve)," whose final element must be greater than 0")
                    }
                  }
             else if( cdf=="pfrechet"){
                if(length(parmsStart)!=3 ||
                      parmsStart[2] <=0  ||
                      parmsStart[3] < 2  ){
                  if(is.na(object@reserve)){stop("The first element of parmsStart must be a starting value for 'reserve'")}
                      stop("For the Frechet distribution, 'parmsStart' must be a numeric vector of length ",
                           2 + is.na(object@reserve)," whose next-to-last element must be positive and whose final element must be at least 2")
                        }
                      }


             return(TRUE)
         }
         )

setGeneric (
 name= "calcOptimalReserve",
 def=function(object,...){standardGeneric("calcOptimalReserve")}
 )
setGeneric (
 name= "cdfG",
 def=function(object,...){standardGeneric("cdfG")}
 )
setGeneric (
  name= "calcExpectedPrice",
  def=function(object,...){standardGeneric("calcExpectedPrice")}
)
setGeneric (
 name= "calcExpectedLowestCost",
 def=function(object,...){standardGeneric("calcExpectedLowestCost")}
 )
setGeneric (
 name= "calcBuyerExpectedCost",
 def=function(object,...){standardGeneric("calcBuyerExpectedCost")}
 )
setGeneric (
 name= "calcExpectedSupplierProfits",
 def=function(object,...){standardGeneric("calcExpectedSupplierProfits")}
 )
setGeneric (
  name= "calcSellerCostParms",
  def=function(object,...){standardGeneric("calcSellerCostParms")}
)
setGeneric (
  name= "calcBuyerValuation",
  def=function(object,...){standardGeneric("calcBuyerValuation")}
)

setMethod(
    f= "calcSellerCostParms",
    signature= "Auction2ndCap",
    definition=function(object,...){

        ## method to calibrate seller cost distribution parameters
        sellerCostParms <- object@sellerCostParms
        reserve         <- object@reserve
        shareInside     <- object@shareInside
        margins         <- object@margins
        prices          <- object@prices
        parmsStart      <- object@parmsStart


        cdf <- object@sellerCostCDF

        minD <- function(parmsStart){

            if(is.na(reserve)){r <- parmsStart[1]; parmsStart <- parmsStart[-1]}
            else{r <- reserve}

            sellerCostParms <- parmsStart

            object@reservePre      <- r
            object@sellerCostParms <- sellerCostParms

            ## For uniform, frechet,  distribution bounds are a function
            ## of distribution parameters

            if(cdf == "punif"){
                object@sellerCostBounds <- parmsStart }
           else if(cdf == "pfrechet"){
               object@sellerCostBounds[1] <- parmsStart[1] }





            ##calculate each bidder's profit margin, conditional on totCap bidder winning
            thisInShare <- cdfG(object,preMerger=TRUE)
            thisMargin  <- calcExpectedSupplierProfits(object,preMerger=TRUE)/calcShares(object,preMerger=TRUE) # why is this not divided by probability of a firm having a cost draw below the reserve?
            thisPrice   <- calcPrices(object,preMerger=TRUE)


            measure   <- c(margins - thisMargin/thisPrice,
                           1 - thisPrice/prices,
                           shareInside - thisInShare)

            return(sum(measure^2,na.rm=TRUE))

        }


    if(cdf == "punif") {

         ui = diag(length(parmsStart))
         if(is.na(reserve)){ui[1,nrow(ui)-1]=-1} #constrain reserve to be greater than cLower
         ui[nrow(ui),nrow(ui)-1]=-1 #constrain cLower to be less than cUpper
         ci = rep(0,length(parmsStart))


         result <- constrOptim(parmsStart,minD,grad=NULL,ui=ui,ci=ci,...)

       }
        else{
          lb <- ub <- rep(Inf,length(parmsStart))

          if( cdf=="pexp"){lb[1] <-0}
          else if( cdf=="pweibull"){ lb[1:2] <- 0} #shape,scale must be positive
          else if( cdf=="pgumbel"){  lb[2] <- 0} #scale must be positive
          else if( cdf=="pfrechet"){ lb[2] <- 0; lb[3] <- 2}  #scale must be positive, shape must be > 2 for finite variance

          if(is.na(reserve)){lb <- c(0,lb); ub <- c(Inf,ub)}

          if(length(parmsStart)>1){method="L-BFGS-B"}
          else{method="Brent"; ub=1e12} #'Brent' is equivalent to using optimize for 1D problems

          result <- optim(parmsStart,minD,method=method,lower=lb,upper=ub,...)

        }

        result <- result$par
        if(is.na(reserve)){object@reserve <- result[1]; result <- result[-1]}
        object@sellerCostParms <- result

        if(cdf == "punif"){
          object@sellerCostBounds <- result }
        else if(cdf == "pfrechet"){
          object@sellerCostBounds[1] <- result[1] }


        return(object)
  }
)

setMethod(
  f= "calcBuyerValuation",
  signature= "Auction2ndCap",
  definition=function(object){

    ## Use FOC from buyers cost minimization problem
    ## to uncover buyer cost parameter
    capacities <- object@capacities
    totCap       <- sum(capacities)
    reserve    <- object@reserve

    object@reservePre <- reserve
    cdfF = match.fun(object@sellerCostCDF)
    pdfF = object@sellerCostPDF

    sellerCostParms <- c(list(reserve),as.list(object@sellerCostParms))
    fc = do.call(pdfF,sellerCostParms)

    sellerCostParms <- c(sellerCostParms,
                         lower.tail=as.list(object@sellerCostCDFLowerTail))
    Fc = do.call(cdfF,sellerCostParms)

    gc <- totCap*fc*(1-Fc)^(totCap-1)

    expectedPrice  <- calcExpectedPrice(object,preMerger=TRUE)
    partialSupplierProfits <- (1-Fc)^(totCap-capacities) - (1-Fc)^totCap
    partialSupplierProfits <- sum(partialSupplierProfits)/gc

    result <- reserve + partialSupplierProfits

    return(result)

  }
)
setMethod(
          f= "calcOptimalReserve",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE,lower,upper){

            if(missing(lower)){lower <- max(object@sellerCostBounds[1],0)}
            if(missing(upper)){upper <- object@buyerValuation}

            minD <- function(r){
              if(preMerger){object@reservePre <- r}
              else{object@reservePost <- r}
                  calcBuyerExpectedCost(object,preMerger=preMerger)
            }

            res <- optimize(
              f  = minD,
              lower = lower,
              upper = upper,
            )


              rStar <- res$minimum
              return(rStar)
          }
)


setMethod(
          f= "calcBuyerExpectedCost",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){


              shareInside <- cdfG(object,preMerger=preMerger)
              val  <- object@buyerValuation * (1-shareInside) + calcExpectedPrice(object,preMerger=preMerger)*shareInside

              return(val)
          })

setMethod(
          f= "calcPrices",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){

             # if(!preMerger){r      <- object@reservePost}
             # else{r <- object@reservePre}

            val <- calcExpectedSupplierProfits(object,preMerger=preMerger) + calcMC(object,preMerger=preMerger)

            #val <- val/cdfG(object,r,preMerger=preMerger)
            val <- val/calcShares(object,preMerger=preMerger)
            names(val) <- object@labels
            return(val)
          })

setMethod(
          f= "cdfG",
          signature= "Auction2ndCap",
          definition=function(object,c,preMerger=TRUE){

              if(missing(c)){
                  if(preMerger){
                      c <- object@reservePre
                      capacities <- sum(object@capacities)
                            }
                  else{
                      c <- object@reservePost
                      capacities <- sum(object@capacities*(1+object@mcDelta))
                   }


              }

              else{

                  if(preMerger){capacities <- object@capacities}
                  else{capacities <- tapply(object@capacities*(1+object@mcDelta),object@ownerPost,sum)}


              }


              cdfF = match.fun(object@sellerCostCDF)
              sellerCostParms <- c(list(c),as.list(object@sellerCostParms),
                                   lower.tail=as.list(object@sellerCostCDFLowerTail))

              Fc = do.call(cdfF,sellerCostParms)
              retval = 1-(1-Fc)^capacities

              if(!preMerger && length(capacities)>1){

                  temp <- rep(NA, length(object@ownerPre))
                  temp[object@ownerPre == object@ownerPost] <- retval
                  retval <- temp

              }
                  return(retval)

              }
          )

setMethod(
  f= "calcExpectedPrice",
  signature= "Auction2ndCap",
  definition=function(object,preMerger=TRUE){

    val <- calcExpectedLowestCost(object,preMerger=preMerger) + sum(calcExpectedSupplierProfits(object,preMerger=preMerger),na.rm=TRUE)/cdfG(object,preMerger=preMerger)
    return(val)
  }
)

setMethod(
  f= "calcMC",
  signature= "Auction2ndCap",
  definition=function(object,t,preMerger=TRUE){


    cdfF <- match.fun(object@sellerCostCDF)
    pdfF <- object@sellerCostPDF
    sellerCostBounds <-object@sellerCostBounds



    if(preMerger) {
      capacities <- object@capacities
      r    <- object@reservePre
    }
    else {
      capacities <- tapply(object@capacities*(1+object@mcDelta),object@ownerPost,sum)
      r    <- object@reservePost
    }

    totCap <- sum(capacities)

    if(missing(t)){t <- capacities}



    ## The expected production cost
    ecIntegrand = function(c,t){
      sellerCostParms <- c(list(c),as.list(object@sellerCostParms))

      fc = do.call(pdfF,sellerCostParms)

      sellerCostParms <- c(sellerCostParms,
                           lower.tail=as.list(object@sellerCostCDFLowerTail))
      Fc = do.call(cdfF,sellerCostParms)

      retval = t*c*fc*(1-Fc)^(totCap-1)
      retval = ifelse(is.finite(retval),retval,0)

      return(retval)
    }

    result <- sapply(
      t,
      function(t.i) {
        if( r < sellerCostBounds[2]) {
            retval = integrate(ecIntegrand,lower=sellerCostBounds[1],upper=r, stop.on.error = FALSE,t=t.i)$value
        }
        else {
          retval = integrate(ecIntegrand,lower=sellerCostBounds[1],upper=sellerCostBounds[2], stop.on.error = FALSE,t=t.i)$value
        }

        return(retval)
      })


    if(!preMerger && length(t)>1){

        temp <- rep(NA, length(object@ownerPre))
        temp[object@ownerPre == object@ownerPost] <- result
        result <- temp

        }

    return(result)

  }
)


setMethod(
          f= "calcExpectedLowestCost",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){


              if(preMerger){capacities <- object@capacities}
              else{capacities <- object@capacities*(1+object@mcDelta)}

              num    <- calcMC(object,t=sum(capacities),preMerger=preMerger)


              retval <- num/cdfG(object, preMerger=preMerger)
              return(retval)


          })


setMethod(
          f= "calcExpectedSupplierProfits",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){

             sellerCostBounds <-object@sellerCostBounds
              if(preMerger){r    <- object@reservePre}
              else{r    <- object@reservePost}


              if(preMerger) { capacities = object@capacities }
              else {          capacities = tapply(object@capacities*(1+object@mcDelta),object@ownerPost,sum) }

              totCap = sum(capacities)

              espIntegrand = function(c,t){
                  sellerCostParms <- c(list(c),as.list(object@sellerCostParms),
                                       lower.tail=as.list(object@sellerCostCDFLowerTail))
                  Fc <- do.call(match.fun(object@sellerCostCDF),sellerCostParms)
                  val <- (1-Fc)^(totCap-t)-(1-Fc)^totCap
              }


              retval <- sapply(
                                capacities,
                                function(t.i) {

                                  if( r < sellerCostBounds[2]) {
                                      retval <- integrate(espIntegrand,lower=sellerCostBounds[1],upper=r,
                                                          stop.on.error = FALSE,t=t.i)$value
                                  }
                                  else {
                                    retval <- integrate(espIntegrand,lower=sellerCostBounds[1],upper=sellerCostBounds[2],
                                                        stop.on.error = FALSE,t=t.i)$value
                                  }

                                    return(retval)
                                })

             if(!preMerger){

                 temp <- rep(NA, length(object@ownerPre))
                 temp[object@ownerPre == object@ownerPost] <- retval
                 retval <- temp

             }

              return(retval)
          })

setMethod(
          f= "calcShares",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){


              ownerPre  <- object@ownerPre
              ownerPost  <- object@ownerPost

              shareInside <- cdfG(object,preMerger=preMerger)

              if(preMerger){
                  return((shareInside*object@capacities)/sum(object@capacities))
              }
              else{

                  capacities <- object@capacities*(1+object@mcDelta)
                  result <- rep(NA, length(object@ownerPre))
                  result[object@ownerPre == ownerPost] <- (shareInside*tapply(capacities,ownerPost,sum))/sum(capacities)
                  return(result)
              }


          }
)


##summarize method

setMethod(
 f= "summary",
 signature= "Auction2ndCap",
 definition=function(object,parameters=FALSE,digits=2,...){

     curWidth <-  getOption("width")


     pricePre   <-  object@pricePre
     pricePost  <-  object@pricePost
     priceDelta <- (pricePost/pricePre - 1) * 100


         outPre  <-  calcShares(object,TRUE) * 100
         outPost <-  calcShares(object,FALSE) * 100

         sumlabels=paste("shares",c("Pre","Post"),sep="")


     #mcDelta <- object@mcDelta

     outDelta <- (outPost/outPre - 1) * 100


     isParty <- object@ownerPost != object@ownerPre
     isParty <- c(object@ownerPre[isParty],object@ownerPost[isParty])
     isParty <- factor(ifelse(object@ownerPre %in% isParty,1,0),levels=0:1,labels=c(" ","*"))

     results <- data.frame(pricePre=pricePre,pricePost=pricePost,
                           priceDelta=priceDelta,outputPre=outPre,
                           outputPost=outPost,outputDelta=outDelta)

     colnames(results)[colnames(results) %in% c("outputPre","outputPost")] <- sumlabels

     #if(sum(abs(mcDelta))>0) results <- cbind(results,mcDelta=mcDelta)


     rownames(results) <- paste(isParty,object@labels)

     sharesPost <- calcShares(object,FALSE)

     cat("\nMerger simulation results under '",class(object),"' demand:\n\n",sep="")

     options("width"=100) # this width ensures that everything gets printed on the same line
     print(round(results,digits),digits=digits)
     options("width"=curWidth) #restore to current width

     cat("\n\tNotes: '*' indicates merging parties' products. Deltas are percent changes.\n")

     cat("\tOutput is based on units sold.\n")

     results <- cbind(isParty, results)

     cat("\n\n% Change In Expected Price:",round((calcExpectedPrice(object,FALSE)-calcExpectedPrice(object,TRUE))/calcExpectedPrice(object,TRUE)*100,digits),sep="\t")


     cat("\n\n")


     if(parameters){

         cat("\nSupplier Cost Distribution Parameters:\n\n")

         print(round(object@sellerCostParms,digits))

         cat("\nBuyer Valuation:\n\n")
         print(round(object@buyerValuation,digits))

             cat("\n\n")

         }

     rownames(results) <- object@labels
     return(invisible(results))

 })




## Create Constructor Function

auction2nd.cap <- function(capacities, margins,prices,reserve=NA,shareInside=1,
                           sellerCostCDF=c("punif","pexp","pweibull","pgumbel","pfrechet"),
                           parmsStart,
                           ownerPre,ownerPost,
                           mcDelta=rep(0,length(capacities)), reserve.fixed=FALSE,
                           labels=paste("Firm",1:length(capacities),sep=""),...
                          ){


    sellerCostCDF <- match.arg(sellerCostCDF)
    lower.tail    <- TRUE
    sellerCostPDF <- match.fun(paste("d",substring(sellerCostCDF,2,),sep=""))

    if(is.na(reserve)){reserve <- NA_real_}

    if (missing(parmsStart)){
      avgPrice =  mean(prices,na.rm=TRUE)
      minPrice =  min(prices,na.rm=TRUE)
      if(sellerCostCDF=="punif"){parmsStart <- sellerCostBounds <- range(prices,na.rm=TRUE)} #uniform on price range
      else if(sellerCostCDF=="pexp"){parmsStart=1/avgPrice; }
      else if(sellerCostCDF=="pweibull"){ parmsStart=c(avgPrice,avgPrice)}
      else if(sellerCostCDF=="pgumbel"){parmsStart=c(avgPrice,avgPrice) }
      else if(sellerCostCDF=="pfrechet"){parmsStart=c(minPrice,minPrice,2) }


      if(is.na(reserve) || missing(reserve)){
        parmsStart<-c(max(prices,na.rm=TRUE),parmsStart)
      }
    }


    if(sellerCostCDF=="punif"){sellerCostBounds <- parmsStart[-1]}
    else if(sellerCostCDF=="pexp"){sellerCostBounds=c(0,Inf)}
    else if(sellerCostCDF=="pweibull"){ sellerCostBounds=c(0,Inf)}
    else if(sellerCostCDF=="pgumbel"){ sellerCostBounds=c(-Inf,Inf) ;lower.tail=FALSE}
    else if(sellerCostCDF=="pfrechet"){ sellerCostBounds=c(parmsStart[2],Inf) ; lower.tail=FALSE}



    result <- new("Auction2ndCap",capacities=capacities,
                  margins=margins,
                  prices=prices,
                  reserve=reserve,
                  shareInside=shareInside,
                  sellerCostCDF=sellerCostCDF,
                  sellerCostCDFLowerTail=lower.tail,
                  sellerCostBounds=sellerCostBounds,
                  sellerCostPDF=sellerCostPDF,
                  ownerPre=ownerPre,ownerPost=ownerPost,
                  mcDelta=mcDelta,
                  parmsStart=parmsStart,
                  labels=labels)


    ## Calibrate seller cost parameters
    result                 <- calcSellerCostParms(result,...)

    ## Calibrate buyer cost parameter
    result@buyerValuation       <- calcBuyerValuation(result)

    ## Compute pre- and post-merger reserves
    if(reserve.fixed){ result@reservePre  <-  result@reservePost <- result@reserve}
    else{
        result@reservePre      <- calcOptimalReserve(result,preMerger=TRUE) #Find Buyer Reserve pre-merger
        result@reservePost     <- calcOptimalReserve(result,preMerger=FALSE) #Find Buyer Reserve pre-merger
    }

    result@pricePre        <- calcPrices(result,preMerger=TRUE)
    result@pricePost        <- calcPrices(result,preMerger=FALSE)

    return(result)
    }



