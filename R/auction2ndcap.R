
setClass(
         Class   = "Auction2ndCap",
         contains="Antitrust",
         representation=representation(
         capacities       = "numeric",
         reserve          = "numeric",
         sellerCostCDF    = "function",
         sellerCostCDFLowerTail    = "logical",
         sellerCostPDF    = "function",
         sellerCostParms  = "numeric",
	       buyerCost        = "numeric", #was buyerValuation
         reservePre       = "numeric",
         reservePost      = "numeric",
         parmsStart       = "numeric"

         ),
         prototype=prototype(
         reservePre      =  numeric(),
         reservePost     =  numeric(),
         parmsStart      =  numeric(),
         buyerCost       =  numeric(),
         sellerCostCDFLowerTail    = TRUE


          ),
         validity=function(object){

             nprods <- length(object@labels)
             cdf    <- as.character(quote(object@sellerCostCDF))


             if(nprods != length(object@capacities)){
                 stop("'labels' must have the same length as 'capacities'")}

             
             if( (cdf=="punif"){
                if(length(object@parmsStart)!=2 ||
                    object@parmsStart[1] >= object@parmsStart[2]){
                  stop("For the Uniform distribution, 'sellerCostParms' must a numeric vector of length 2 whose first element is less than the second element")
                }
             }
               
             else if( (cdf=="pexp"){
                if(length(object@parmsStart)!=1 ||
                    object@parmsStart[1] <=0){
                    stop("For the Exponential distribution, 'sellerCostParms' must a numeric vector of length 1 whose first element is greater than 0")
                 }
             }
             else if( (cdf=="pweibull"){
                if(length(object@parmsStart)!=2 ||
                     object@parmsStart[1] <=0  ||
                     object@parmsStart[2] <=0 ){
                    stop("For the Weibull distribution, 'sellerCostParms' must a numeric vector of length 1 whose first element is greater than 0")
                        }
                      }
             else if( (cdf=="pgumbel"){
                if(length(object@parmsStart)!=2 ||
                     object@parmsStart[2] <=0){
                    stop("For the Gumbel distribution, 'sellerCostParms' must a numeric vector of length 1 whose first element is greater than 0")
                    }
                  }
             else if( (cdf=="pfrechet"){
                if(length(object@parmsStart)!=3 ||
                      object@parmsStart[2] <=0  ||
                      object@parmsStart[3] <=0  ){
                      stop("For the Frechet distribution, 'sellerCostParms' must a numeric vector of length 1 whose first element is greater than 0")
                        }
                      }       
            else if( (cdf=="pgev"){
                 if(length(object@parmsStart)!=3 ||
                    object@parmsStart[1] <=0){
                    stop("For the GEV distribution, 'sellerCostParms' must a numeric vector of length 1 whose first element is greater than 0")
                    }
                  }                 
             if(length(object@buyerCost)!=1){
                 stop("'buyerCost' must be a length-1 numeric vector")
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
  name= "calcBuyerCost",
  def=function(object,...){standardGeneric("calcBuyerCost")}
)

setMethod(
    f= "calcSellerCostParms",
    signature= "Auction2ndCap",
    definition=function(object){

        ## method to calibrate seller cost distribution parameters
        sellerCostParms <- object@sellerCostParms
        reserve         <- object@reserve
        shareInside     <- object@shareInside
        margins         <- object@margins
        prices          <- object@prices
        parmsStart      <- object@parmsStart


        minD <- function(parmsStart){

            if(is.na(reserve)){r <- parmsStart[1]; parmsStart <- parmsStart[-1]}
            else{r <- reserve}

            sInside <- parmsStart[1]
            sellectCostParms <- parmsStart[-1]

            object@reservePre      <- r
            object@sellerCostParms <- sellerCostParms


            ##calculate each bidder's profit margin, conditional on that bidder winning
            thisMargin  <- calcExpectedSupplierProfits(object)/calcShares(object)
            thisPrice   <- calcPrices(object)
            thisInShare <- cdfG(object)

            predicted   <- c(thisMargin/thisPrice,thisInShare)
            observed    <- c(margins,shareInside)

            return(sum((predicted-observed)^2,na.rm=TRUE))

        }

        result <- optim(parmsStart,minD)

  }
)

setMethod(
  f= "calcBuyerCost",
  signature= "Auction2ndCap",
  definition=function(object){
    
    ## Use FOC from buyers cost minimization problem
    ## to uncover buyer cost parameter
    capacities <- object@capacities
    tHat       <- sum(capacities)
    reserve    <- object@reserve
    
    object@reservePre <- reserve
    cdfF = object@sellerCostCDF
    pdfF = object@sellerCostPDF
    
    sellerCostParms <- c(reserve,as.list(object@sellerCostParms))
    fc = do.call(pdfF,sellerCostParms)
    
    sellerCostParms <- c(sellerCostParms,
                         lower.tail=object@sellerCostCDFLowerTail)
    Fc = do.call(cdfF,sellerCostParms)
    
    gc <- tHat*fc*(1-Fc)^(tHat-1)
    
    expectedPrice  <- calcPrices(object,preMerger=TRUE)
    partialSupplierProfits <- (1-Fc)^(tHat-capacities) - (1-Fc)^tHat
    partialSupplierProfits <- sum(partialSupplierProfits)/gc
    
    result <- reserve - 2*expectedPrice + partialSupplierProfits
    
    return(max(result,0))

  }
)
setMethod(
          f= "calcOptimalReserve",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE,lower=-1e12,upper=-lower){

            #sellerCostCDF = as.character(quote(object@sellerCostCDF))



            minD <- function(r){
              if(preMerger){object@reservePre <- r}
              else{object@reservePost <- r}
                  calcBuyerExpectedCost(object,preMerger)
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
            

              val  <- object@buyerCost * (1-cdfG(object,preMerger=preMerger)) + calcPrices(object,preMerger)*cdfG(object,preMerger=preMerger)

              return(val)
          })

setMethod(
          f= "calcPrices",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){

              val <- calcExpectedLowestCost(object) + sum(calcExpectedSupplierProfits(object,preMerger))/cdfG(object,preMerger=preMerger)

              return(val)
          })

setMethod(
          f= "cdfG",
          signature= "Auction2ndCap",
          definition=function(object,c,preMerger=TRUE){

              if(missing(c)){
                  if(preMerger){c <- object@reservePre}
                  else{c <- object@reservePost}

                  capacities <- sum(object@capacities)
              }
              else{capacities <- object@capacities}


              cdfF = object@sellerCostCDF
              sellerCostParms <- c(c,as.list(object@sellerCostParms),
                                   lower.tail=object@sellerCostCDFLowerTail)

              Fc = do.call(cdfF,sellerCostParms)
              retval = 1-(1-Fc)^capacities

                  return(retval)

              }
          )

setMethod(
          f= "calcExpectedLowestCost",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){

              cdfF <- object@sellerCostCDF
              pdfF <- object@sellerCostPDF

              if(preMerger) {
                  tVec <- object@capacities
                  r    <- object@reservePre
              }
              else {
                  tVec <- tapply(object@capacities,object@ownerPost,sum)
                  r    <- object@reservePost
              }

              tHat = sum(tVec)

               ## The expected lowest production cost
               elcIntegrand = function(c){
                   sellerCostParms <- c(c,as.list(object@sellerCostParms),
                                        lower.tail=object@sellerCostCDFLowerTail)

                  fc = do.call(pdfF,sellerCostParms)
                  Fc = do.call(cdfF,sellerCostParms)
                  retval = c*tHat*fc*(1-Fc)^(tHat-1)

                  return(retval)
               }



                  num = integrate(elcIntegrand,lower=-Inf,upper=r)$value

                  retval = num/cdfG(object,r)
                  return(retval)


          })


setMethod(
          f= "calcExpectedSupplierProfits",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){

              if(preMerger){r    <- object@reservePre}
              else{r    <- object@reservePost}


              if(preMerger) { tVec = object@capacities }
              else { tVec = tapply(object@capacities,object@ownerPost,sum) }

              tHat = sum(tVec)

              espIntegrand = function(c,t){
                  sellerCostParms <- c(c,as.list(object@sellerCostParms),
                                   lower.tail=object@sellerCostCDFLowerTail)
                  Fc = do.call(object@sellerCostCDF,sellerCostParms)
                val = (1-Fc)^(tHat-t)-(1-Fc)^tHat
              }


              retval <- sapply(
                                tVec,
                                function(t.i) {
                                        retval = integrate(espIntegrand,lower=-Inf,upper=r,t.i)$value
                                    return(retval)
                                })

              return(retval)
          })

setMethod(
          f= "calcShares",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){

              capacities <- object@capacities
              ownerPost  <- object@ownerPost
              cdfF = object@sellerCostCDF
              sellerCostParms <- as.list(object@sellerCostParms)
              tHat <- sum(capacities)



              if(preMerger){

                  sellerCostParms <- c(object@reservePre,sellerCostParms,
                                   lower.tail=object@sellerCostCDFLowerTail)
                  Fc = do.call(cdfF,sellerCostParms)
                  shareOutside <- 1 - (1-Fc^tHat)
                  return((shareOutside*capacities)/sum(capacities))
              }
              else{
                  sellerCostParms <- c(object@reservePost,sellerCostParms,
                                   lower.tail=object@sellerCostCDFLowerTail)
                  Fc = do.call(cdfF,sellerCostParms)
                  return((shareOutside*tapply(capacities,ownerPost,sum))/sum(capacities))
              }

          }
)





## Create Constructor Function

auction2nd.cap <- function(capacities, margins,prices, shareInside=1, reserve=10*max(prices),
                           sellerCostCDF=c("punif","pexp","pweibull","pgumbel","pfrechet","pgev"),
                           parmsStart,
                           ownerPre,ownerPost,
                           labels=paste("Firm",1:length(capacities),sep="")
                          ){

##Gumbel is the log of the weibull and as such may be added
##Frechet is 1/Weibull and may be added. GEV works for the same reason
##add package evd

    sellerCostCDF <- match.arg(sellerCostCDF)
    lower.tail    <- TRUE
    sellerCostPDF <- match.fun(paste("d",substring(sellerCostCDF,2,),sep=""))



    if (missing(parmsStart)){
      if(sellerCostCDF=="punif"){parmsStart=c(0,1)} #uniform on 0,1

      else if(sellerCostCDF=="pexp"){parmsStart=1}

      else if(sellerCostCDF=="pweibull"){parmsStart=c(1,5)} #symmetric, unimodal distribution
      else if(sellerCostCDF=="pgumbel"){parmsStart=c(1,5); lower.tail=FALSE} #symmetric, unimodal distribution
      else if(sellerCostCDF=="pfrechet"){parmsStart=c(0,1,5); lower.tail=FALSE} #symmetric, unimodal distribution located @ 0
      else if(sellerCostCDF=="pgev"){parmsStart=c(1,5,0); lower.tail=FALSE} #symmetric, unimodal distribution located @ 0

    }

    result <- new("Auction2ndCap",capacities=capacities,
                  reserve=reserve,
                  shareInside=shareInside,
                  sellerCostCDF=match.fun(sellerCostCDF),
                  sellerCostCDFLowerTail=lower.tail,
                  sellerCostPDF=sellerCostPDF,
                  sellerCostParms=sellerCostParms,
                  ownerPre=ownerPre,ownerPost=ownerPost,
                  parmsStart=parmsStart,
                  labels=labels)


    result@sellerCostParms <- calcSellerCostParms(result)
    result@buyerCost       <- calcBuyerCost(result)
    
    result@reservePre      <- calcOptimalReserve(result,preMerger=TRUE) #Find Buyer Reserve pre-merger
    result@reservePost     <- calcOptimalReserve(result,preMerger=FALSE) #Find Buyer Reserve pre-merger
    
    result@pricePre        <- calcPrices(result,preMerger=TRUE) # Calculate premerger expected price
    result@pricePost       <- calcPrices(result,preMerger=FALSE) # Calculate postmerger expected price

    return(result)
    }



