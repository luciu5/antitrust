
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
         sellerCostParms  = "numeric",
	       buyerCost        = "numeric", #was buyerValuation
         reservePre       = "numeric",
         reservePost      = "numeric",
         parmsStart       = "numeric"

         ),
         prototype=prototype(
         reservePre      =  numeric(),
         reservePost     =  numeric(),
         buyerCost       =  numeric(),
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
                      parmsStart[3] <=0  ){
                  if(is.na(object@reserve)){stop("The first element of parmsStart must be a starting value for 'reserve'")}
                      stop("For the Frechet distribution, 'parmsStart' must be a numeric vector of length ", 
                           2 + is.na(object@reserve)," whose final 2 elements must be greater than 0")
                        }
                      }       
            else if( cdf=="pgev"){
                 if(length(parmsStart)!=3 ||
                    parmsStart[1] <=0){
                   if(is.na(object@reserve)){stop("For the GEV distribution, 'parmsStart' must be a numeric vector of length 4. The first element of parmsStart must be a starting value for 'reserve'. The second element must be positive")}
                    stop("For the GEV distribution, 'parmsStart' must be a numeric vector of length  ", 
                         3," whose 1st element must be greater than 0")
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
  name= "calcBuyerCost",
  def=function(object,...){standardGeneric("calcBuyerCost")}
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


            ##calculate each bidder's profit margin, conditional on that bidder winning
            thisMargin  <- calcExpectedSupplierProfits(object)/calcShares(object)
            thisPrice   <- calcPrices(object)/calcShares(object)
            thisInShare <- cdfG(object)

            predicted   <- c(thisMargin,thisPrice,thisInShare)
            observed    <- c(margins*prices,prices,shareInside)

            return(sum((predicted-observed)^2,na.rm=TRUE))

        }

       if(cdf == "punif") {
         ui = diag(length(parmsStart) + is.na(reserve))
         ui[nrow(ui),nrow(ui)-1]=-1 #constrain cLower to be less than cUpper
         ci = rep(0,length(parmsStart) + is.na(reserve))
         
         result <- constrOptim(parmsStart,minD,grad=NULL,ui=ui,ci=ci,...)
       }
        else{
          lb <- ub <- rep(Inf,length(parmsStart))
          
          if( cdf=="pexp"){lb[1] <-0}
          else if( cdf=="pweibull"){ lb[1:2] <- 0}
          else if( cdf=="pgumbel"){  lb[2] <- 0}
          else if( cdf=="pfrechet"){ lb[2:3] <- 0}       
          else if( cdf=="pgev"){ lb[1] <- 0}                 
          
          if(is.na(reserve)){lb <- c(0,lb); ub <- c(Inf,ub)}
          
          if(length(parmsStart)>1){method="L-BFGS-B"}
          else{method="Brent"; ub=1e12} #'Brent' is equivalent to using optimize for 1D problems
          
          result <- optim(parmsStart,minD,method=method,lower=lb,upper=ub,...)
          
        }

        object@sellerCostParms <- 
        if(is.na(reserve)){object@reserve <- result$par[1]; object@sellerCostParms <-  result$par[-1]}
        else{object@sellerCostParms <- result$par}
        
        return(object)
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
    cdfF = match.fun(object@sellerCostCDF)
    pdfF = object@sellerCostPDF
    
    sellerCostParms <- c(list(reserve),as.list(object@sellerCostParms))
    fc = do.call(pdfF,sellerCostParms)
    
    sellerCostParms <- c(sellerCostParms,
                         lower.tail=as.list(object@sellerCostCDFLowerTail))
    Fc = do.call(cdfF,sellerCostParms)
    
    gc <- tHat*fc*(1-Fc)^(tHat-1)
    
    expectedPrice  <- calcExpectedPrice(object,preMerger=TRUE)
    partialSupplierProfits <- (1-Fc)^(tHat-capacities) - (1-Fc)^tHat
    partialSupplierProfits <- sum(partialSupplierProfits)/gc
    
    result <- reserve - 2*expectedPrice + partialSupplierProfits
    
    return(max(result,0,na.rm=TRUE))

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
            

              val  <- object@buyerCost * (1-cdfG(object,preMerger=preMerger)) + calcExpectedPrice(object,preMerger)*cdfG(object,preMerger=preMerger)

              return(val)
          })

setMethod(
          f= "calcPrices",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){

            ## is this correct? why are all other firm's profits a function of firm's price?
              #val <- calcExpectedLowestCost(object) + sum(calcExpectedSupplierProfits(object,preMerger))/cdfG(object,preMerger=preMerger)
            ## I think this is the ex ante price that firm i sets. Should I be using the lowest cost 
            ## across all firms, or just that firm's lowest cost. I think the latter-- check with Nate.
            val <- calcExpectedSupplierProfits(object,preMerger) + calcExpectedLowestCost(object,preMerger)*cdfG(object,preMerger=preMerger)
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


              cdfF = match.fun(object@sellerCostCDF)
              sellerCostParms <- c(list(c),as.list(object@sellerCostParms),
                                   lower.tail=as.list(object@sellerCostCDFLowerTail))

              Fc = do.call(cdfF,sellerCostParms)
              retval = 1-(1-Fc)^capacities

                  return(retval)

              }
          )

setMethod(
  f= "calcExpectedPrice",
  signature= "Auction2ndCap",
  definition=function(object,preMerger=TRUE){
  
    val <- calcExpectedLowestCost(object,preMerger) + sum(calcExpectedSupplierProfits(object,preMerger))/cdfG(object,preMerger=preMerger)
    return(val)  
  }
)

setMethod(
          f= "calcExpectedLowestCost",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){

              cdfF <- match.fun(object@sellerCostCDF)
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
                   sellerCostParms <- c(list(c),as.list(object@sellerCostParms))
                   fc = do.call(pdfF,sellerCostParms)
                   sellerCostParms <- c(sellerCostParms,
                                        lower.tail=as.list(object@sellerCostCDFLowerTail))
                  Fc = do.call(cdfF,sellerCostParms)
                  retval = c*tHat*fc*(1-Fc)^(tHat-1)

                  return(retval)
               }



              
                  num = integrate(elcIntegrand,lower=-Inf,upper=r)$value

                  retval = num/cdfG(object)
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
                  sellerCostParms <- c(list(c),as.list(object@sellerCostParms),
                                   lower.tail=as.list(object@sellerCostCDFLowerTail))
                  Fc = do.call(match.fun(object@sellerCostCDF),sellerCostParms)
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
              cdfF = match.fun(object@sellerCostCDF)
              sellerCostParms <- as.list(object@sellerCostParms)
              tHat <- sum(capacities)



              if(preMerger){

                  sellerCostParms <- c(list(object@reservePre),sellerCostParms,
                                   lower.tail=as.list(object@sellerCostCDFLowerTail))
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

auction2nd.cap <- function(capacities, margins,prices,reserve=10*max(prices),shareInside=1,
                           sellerCostCDF=c("punif","pexp","pweibull","pgumbel","pfrechet","pgev"),
                           parmsStart,
                           ownerPre,ownerPost,
                           labels=paste("Firm",1:length(capacities),sep=""),...
                          ){


    sellerCostCDF <- match.arg(sellerCostCDF)
    lower.tail    <- TRUE
    sellerCostPDF <- match.fun(paste("d",substring(sellerCostCDF,2,),sep=""))

    if(is.na(reserve)){reserve <- NA_real_} 

    if (missing(parmsStart)){
      
      if(sellerCostCDF=="punif"){parmsStart=range(prices)} #uniform on price range
      else if(sellerCostCDF=="pexp"){parmsStart=1}
      else if(sellerCostCDF=="pweibull"){parmsStart=c(1,5)} #symmetric, unimodal distribution
      else if(sellerCostCDF=="pgumbel"){parmsStart=c(1,5); lower.tail=FALSE} #symmetric, unimodal distribution
      else if(sellerCostCDF=="pfrechet"){parmsStart=c(0,1,5); lower.tail=FALSE} #symmetric, unimodal distribution located @ 0
      else if(sellerCostCDF=="pgev"){parmsStart=c(1,5,0); lower.tail=FALSE} #symmetric, unimodal distribution located @ 0

    
      if(is.na(reserve)){
        parmsStart<-c(10*max(prices),parmsStart)
      }
    }

    result <- new("Auction2ndCap",capacities=capacities,
                  margins=margins,
                  prices=prices,
                  reserve=reserve,
                  shareInside=shareInside,
                  sellerCostCDF=sellerCostCDF,
                  sellerCostCDFLowerTail=lower.tail,
                  sellerCostPDF=sellerCostPDF,
                  ownerPre=ownerPre,ownerPost=ownerPost,
                  parmsStart=parmsStart,
                  labels=labels)


    ## Calibrate seller cost parameters
    result                 <- calcSellerCostParms(result,...)
    
    ## Calibrate buyer cost parameter
    result@buyerCost       <- calcBuyerCost(result)
    
    ## Compute pre- and post-merger reserves
    result@reservePre      <- calcOptimalReserve(result,preMerger=TRUE) #Find Buyer Reserve pre-merger
    result@reservePost     <- calcOptimalReserve(result,preMerger=FALSE) #Find Buyer Reserve pre-merger
    
    result@pricePre        <- calcPrices(result,preMerger=TRUE)
    result@pricePost        <- calcPrices(result,preMerger=FALSE)
    
    return(result)
    }



