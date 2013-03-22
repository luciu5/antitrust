
setClass(
         Class   = "Auction2ndCap",
         contains="Antitrust",
         representation=representation(
         capacities       = "numeric",
         sellerCostCDF    = "function",
         sellerCostPDF    = "function",
         sellerCostParms    = "numeric",
	       buyerCost        = "numeric", #was buyerValuation
         buyerReserve     = "numeric"
         ),
         prototype=prototype(
          buyerReserve     =  numeric()

          ),
         validity=function(object){

             nprods <- length(object@labels)
            

             if(nprods != length(object@capacities)){
                 stop("'labels' must have the same length as 'capacities'")}

             if(length(object@sellerCostParms)!=2 ||
                object@sellerCostParms[1] >= object@sellerCostParms[2]){
                 stop("'sellerCostParms' must a numeric vector of length 2 whose first element is less than the second element")
             }

             if(length(object@buyerCost)!=1 ||
                !object@buyerCost %in% object@sellerCostParms){
                 stop("'buyerCost' must be a length-1 numeric vector within the range of 'sellerCostParms'")
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
 name= "calcSaleProb",
 def=function(object,...){standardGeneric("calcSaleProb")}
 )




setMethod(
          f= "calcOptimalReserve",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){


              c.LB = object@sellerCostParms[1]
              c.UB = object@sellerCostParms[2]


              minD <- function(r){
                  object@buyerReserve <- r
                  calcBuyerExpectedCost(object,preMerger)
              }

              res <- optimize(
                              f  = minD,
                              lower = c.LB,
                              upper = c.UB,
                              )

              rStar <- res$minimum
              return(rStar)
          }
)


setMethod(
          f= "calcBuyerExpectedCost",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){

              val  <- object@buyerCost * (1-cdfG(object)) + calcPrices(object,preMerger)*cdfG(object)

              return(val)
          })

setMethod(
          f= "calcPrices",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){

              val <- calcExpectedLowestCost(object) + sum(calcExpectedSupplierProfits(object,preMerger))/cdfG(object)

              return(val)
          })

setMethod(
          f= "cdfG",
          signature= "Auction2ndCap",
          definition=function(object,c){

              if(missing(c)){ c <- object@buyerReserve}

              tHat <- sum(object@capacities)

              cdfF = object@sellerCostCDF
              c.LB = object@sellerCostParms[1]
              c.UB = object@sellerCostParms[2]

               if(c < c.LB) { retval = 0 }
               else if(c > c.UB) { retval = 1 }
               else {
                   Fc = do.call(cdfF,list(c))
                   retval = 1-(1-Fc)^tHat
               }
                  return(retval)

              }
          )

setMethod(
          f= "calcExpectedLowestCost",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){

              r    <- object@buyerReserve
              cdfF <- object@sellerCostCDF
              pdfF <- object@sellerCostPDF
              c.LB <- object@sellerCostParms[1]
              c.UB <- object@sellerCostParms[2]

              if(preMerger) { tVec = object@capacities }
              else { tVec = tapply(object@capacities,object@ownerPost,sum) }

              tHat = sum(tVec)

               ## The expected lowest production cost
               elcIntegrand = function(c){
                   fc = do.call(pdfF,list(c))
                   Fc = do.call(cdfF,list(c))
                   retval = c*tHat*fc*(1-Fc)^(tHat-1)
                   return(retval)
               }


                  if(r < c.UB) {
                      num = integrate(elcIntegrand,lower=c.LB,upper=r)$value
                  } else {
                      num = integrate(elcIntegrand,lower=c.LB,upper=c.UB)$value
                  }
                  retval = num/cdfG(object,r)
                  return(retval)


          })


setMethod(
          f= "calcExpectedSupplierProfits",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){

              r    <- object@buyerReserve
              c.LB <- object@sellerCostParms[1]
              c.UB <- object@sellerCostParms[2]

              if(preMerger) { tVec = object@capacities }
              else { tVec = tapply(object@capacities,object@ownerPost,sum) }

              tHat = sum(tVec)

              espIntegrand = function(c,t){
                  Fc = do.call(object@sellerCostCDF,list(c))
                  val = (1-Fc)^(tHat-t)-(1-Fc)^tHat
              }


              retval <- sapply(
                                tVec,
                                function(t.i) {
                                    if( r < c.UB ) {
                                        retval = integrate(espIntegrand,lower=c.LB,upper=r,t.i)$value
                                    } else {
                                        retval = integrate(espIntegrand,lower=c.LB,upper=c.UB,t.i)$value
                                    }
                                    return(retval)
                                })

              return(retval)
          })

setMethod(
          f= "calcSaleProb",
          signature= "Auction2ndCap",
          definition=function(object){


          }
)





## Create Constructor Function

auction2nd.cap <- function(capacities,buyerCost=0,
                           sellerCostCDF=c(punif,pexp,pweibull)
                           sellerCostParms,
                           ownerPre,ownerPost,
                           labels=paste("Firm",1:length(capacities),sep="")
                          ){

  
    sellerCostCDF <- match.arg(sellerCostCDF)
    sellerCostPDF <- match.fun(paste("d",substring(sellerCostCDF,2,),sep=""))
    
    if(missing(sellerCostParms)){
      
      
    }

   
    result <- new("Auction2ndCap",capacities=capacities,
                  buyerCost=buyerCost,
                  sellerCostCDF=sellerCostCDF,
                  sellerCostPDF=sellerCostPDF,
                  sellerCostParms=sellerCostParms,
                  ownerPre=ownerPre,ownerPost=ownerPost,
                  labels=labels)


    result@buyerReserve  <- calcOptimalReserve(result,preMerger=TRUE) #Find Buyer Reserve
    result@pricePre      <- calcPrices(result,preMerger=TRUE) # Calculate premerger expected price
    result@pricePost     <- calcPrices(result,preMerger=FALSE) # Calculate postmerger expected price

    return(result)
    }



