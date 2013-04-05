
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
           print(parmsStart)
            if(is.na(reserve)){r <- parmsStart[1]; parmsStart <- parmsStart[-1]}
            else{r <- reserve}

            sellerCostParms <- parmsStart

            object@reservePre      <- r
            object@sellerCostParms <- sellerCostParms
           
           ## For uniform. frechet, and gev, distribution bounds are a function 
           ## of distribution parameters
           
           if(cdf == "punif"){
                object@sellerCostBounds <- parmsStart }
           else if(cdf == "pfrechet"){
                object@sellerCostBounds[1] <- parmsStart[1] }
           else if (cdf=="pgev"){
              object@sellerCostBounds <-c(ifelse(parmsStart[3]>0,parmsStart[1]-parmsStart[2]/parmsStart[3],-Inf),
                                          ifelse(parmsStart[3]<0,parmsStart[1]-parmsStart[2]/parmsStart[3],Inf)
                                          ) 
            
            }

            
            thisShare <- calcShares(object)

            ##calculate each bidder's profit margin, conditional on totCap bidder winning
            thisMargin  <- calcExpectedSupplierProfits(object)/thisShare
            thisPrice   <- calcPrices(object)/thisShare
            thisInShare <- cdfG(object)

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
          else if( cdf=="pfrechet"){ lb[2:3] <- 0}  #shape,scale must be positive     
          else if( cdf=="pgev"){ lb[2] <- 0} #scale must be positive                 
          
          if(is.na(reserve)){lb <- c(0,lb); ub <- c(Inf,ub)}
          
          if(length(parmsStart)>1){method="L-BFGS-B"}
          else{method="Brent"; ub=1e12} #'Brent' is equivalent to using optimize for 1D problems
          
          result <- optim(parmsStart,minD,method=method,lower=lb,upper=ub,...)
          
        }
        print(result)
        result <- result$par
        if(is.na(reserve)){object@reserve <- result[1]; result <- result[-1]}
        object@sellerCostParms <- result
       
        if(cdf == "punif"){
          object@sellerCostBounds <- result }
        else if(cdf == "pfrechet"){
          object@sellerCostBounds[1] <- result[1] }
        else if (cdf=="pgev"){
          object@sellerCostBounds <-c(ifelse(result[3]>0,result[1]-result[2]/result[3],-Inf),
                                      ifelse(result[3]<0,result[1]-result[2]/result[3],Inf)
          ) 
          
        }
        
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
            if(missing(upper)){upper <- object@buyerCost}

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
              val  <- object@buyerCost * (1-shareInside) + calcExpectedPrice(object,preMerger=preMerger)*shareInside

              return(val)
          })

setMethod(
          f= "calcPrices",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){

            val <- calcExpectedSupplierProfits(object,preMerger=preMerger) + calcMC(object,preMerger=preMerger)
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
              else{
               
                capacities <- object@capacities
                #if(length(c)!=length(capacities)){stop("'c' must have the same length as 'capacities")}
              }


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
  
    val <- calcExpectedLowestCost(object,preMerger=preMerger) + sum(calcExpectedSupplierProfits(object,preMerger=preMerger))/cdfG(object,preMerger=preMerger)
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
      capacities <- tapply(object@capacities,object@ownerPost,sum)
      r    <- object@reservePost
    }
   
    if(missing(t)){t <- capacities}
    
    ## The expected lowest production cost
    elcIntegrand = function(c,t){
      sellerCostParms <- c(list(c),as.list(object@sellerCostParms))
     
      fc = do.call(pdfF,sellerCostParms)
     
      sellerCostParms <- c(sellerCostParms,
                           lower.tail=as.list(object@sellerCostCDFLowerTail))
      Fc = do.call(cdfF,sellerCostParms)
      retval = c*t*fc*(1-Fc)^(t-1)
      retval = ifelse(is.finite(retval),retval,0)
      
      return(retval)
    }
    
    result <- sapply(
      t,
      function(t.i) {
        if( r < sellerCostBounds[2]) {
            retval = integrate(elcIntegrand,lower=sellerCostBounds[1],upper=r, stop.on.error = FALSE,t=t.i)$value 
        } 
        else {
          retval = integrate(elcIntegrand,lower=sellerCostBounds[1],upper=sellerCostBounds[2], stop.on.error = FALSE,t=t.i)$value
        }
       
        return(retval)
      })
    #result = integrate(elcIntegrand,lower=sellerCostBounds[1],upper=r)$value
    
    return(result)
    
  }
)


setMethod(
          f= "calcExpectedLowestCost",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){
            
            
              num    <- calcMC(object,t=sum(object@capacities),preMerger=preMerger)

  
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
              else {          capacities = tapply(object@capacities,object@ownerPost,sum) }

              totCap = sum(capacities)

              espIntegrand = function(c,t){
                  sellerCostParms <- c(list(c),as.list(object@sellerCostParms),
                                   lower.tail=as.list(object@sellerCostCDFLowerTail))
                  Fc = do.call(match.fun(object@sellerCostCDF),sellerCostParms)
                 val = (1-Fc)^(totCap-t)-(1-Fc)^totCap
              }


              retval <- sapply(
                                capacities,
                                function(t.i) {
                                 
                                  if( r < sellerCostBounds[2]) {
                                    retval = integrate(espIntegrand,lower=sellerCostBounds[1],upper=r, stop.on.error = FALSE,t=t.i)$value 
                                  } 
                                  else {
                                    retval = integrate(espIntegrand,lower=sellerCostBounds[1],upper=sellerCostBounds[2], stop.on.error = FALSE,t=t.i)$value
                                  }
                                        
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
             
              shareOutside <- cdfG(object,preMerger=preMerger)

              if(preMerger){ 
                  return((shareOutside*capacities)/sum(capacities))
              }
              else{
                  
                  return((shareOutside*tapply(capacities,ownerPost,sum))/sum(capacities))
              }

          }
)





## Create Constructor Function

auction2nd.cap <- function(capacities, margins,prices,reserve=NA,shareInside=1,
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
      
      if(sellerCostCDF=="punif"){parmsStart <- sellerCostBounds <- range(prices,na.rm=TRUE)} #uniform on price range
      else if(sellerCostCDF=="pexp"){parmsStart=1; }
      else if(sellerCostCDF=="pweibull"){ parmsStart=c(1,5)} #symmetric, unimodal distribution
      else if(sellerCostCDF=="pgumbel"){parmsStart=c(1,5) } #symmetric, unimodal distribution
      else if(sellerCostCDF=="pfrechet"){parmsStart=c(0,1,5) } #symmetric, unimodal distribution located @ 0
      else if(sellerCostCDF=="pgev"){parmsStart=c(1,5,0)} #symmetric, unimodal distribution located @ 0

    
      if(is.na(reserve) || missing(reserve)){
        parmsStart<-c(max(prices,na.rm=TRUE),parmsStart)
      }
    }
    
    
    if(sellerCostCDF=="punif"){sellerCostBounds <- parmsStart[-1]} 
    else if(sellerCostCDF=="pexp"){sellerCostBounds=c(0,Inf)}
    else if(sellerCostCDF=="pweibull"){ sellerCostBounds=c(0,Inf)} 
    else if(sellerCostCDF=="pgumbel"){ sellerCostBounds=c(-Inf,Inf) ;lower.tail=FALSE} 
    else if(sellerCostCDF=="pfrechet"){ sellerCostBounds=c(parmsStart[2],Inf) ; lower.tail=FALSE} 
    else if(sellerCostCDF=="pgev"){sellerCostBounds=c(ifelse(parmsStart[4]>0,parmsStart[2]-parmsStart[3]/parmsStart[4],-Inf),
                                                      ifelse(parmsStart[4]<0,parmsStart[2]-parmsStart[3]/parmsStart[4],Inf)); lower.tail=FALSE} 
    
    

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
                  parmsStart=parmsStart,
                  labels=labels)


    ## Calibrate seller cost parameters
    #result                 <- calcSellerCostParms(result,...)
    
    ## Calibrate buyer cost parameter
    #result@buyerCost       <- calcBuyerCost(result)
    
    ## Compute pre- and post-merger reserves
    #result@reservePre      <- calcOptimalReserve(result,preMerger=TRUE) #Find Buyer Reserve pre-merger
    #result@reservePost     <- calcOptimalReserve(result,preMerger=FALSE) #Find Buyer Reserve pre-merger
    
    #result@pricePre        <- calcPrices(result,preMerger=TRUE)
    #result@pricePost        <- calcPrices(result,preMerger=FALSE)
    
    return(result)
    }



