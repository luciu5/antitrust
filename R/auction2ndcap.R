
setClass(
         Class   = "Auction2ndCap",
         contains="Antitrust",
         representation=representation(
         capacities       = "numeric",
         sellerCostCDF    = "function",
         sellerCostPDF    = "function",
         rangeCostDist    = "numeric",
	 buyerCost        = "numeric", #was buyerValuation
         buyerReserve     = "numeric"
         ),
         prototype=prototype(
          buyerReserve     =  numeric()

          ),
         validity=function(object){

             nprods <- length(object@labels)
             ##cdfName <- paste(quote(object@sellerCostCDF))

             ##if(substring(cdfName,1,1) != "p"){
             ##    stop("'sellerCostCDF' function must be a density beginning  with the letter 'p'")
             ##}

             if(nprods != length(object@capacities)){
                 stop("'labels' must have the same length as 'capacities'")}

             if(length(object@rangeCostDist)!=2 ||
                object@rangeCostDist[1] >= object@rangeCostDist[2]){
                 stop("'rangeCostDist' must a numeric vector of length 2 whose first element is less than the second element")
             }

             if(length(object@buyerCost)!=1 ||
                !object@buyerCost %in% object@rangeCostDist){
                 stop("'buyerCost' must be a length-1 numeric vector within the range of 'rangeCostDist'")
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


              c.LB = object@rangeCostDist[1]
              c.UB = object@rangeCostDist[2]


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
              c.LB = object@rangeCostDist[1]
              c.UB = object@rangeCostDist[2]

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
              c.LB <- object@rangeCostDist[1]
              c.UB <- object@rangeCostDist[2]

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
              c.LB <- object@rangeCostDist[1]
              c.UB <- object@rangeCostDist[2]

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
                           sellerCostCDF=punif,
                           sellerCostPDF=dunif,
                         rangeCostDist=c(0,1e12),
                         ownerPre,ownerPost,
                         labels=paste("Firm",1:length(capacities),sep="")
                         ){

    ##pdfName <- paste("d",substring(sellerCostCDF,2,),sep="")

    ##if(length(getAnywhere(pdfName)$objs) > 0){
    ##    sellerCostPDF <- getAnywhere(pdfName)$objs[[1]]
    ##}
    ##else{stop("density function for CDF ",cdfName," not found")}

    ## Create Logit  container to store relevant data
    result <- new("Auction2ndCap",capacities=capacities,
                  buyerCost=buyerCost,
                  sellerCostCDF=sellerCostCDF,
                  sellerCostPDF=sellerCostPDF,
                  rangeCostDist=rangeCostDist,
                  ownerPre=ownerPre,ownerPost=ownerPost,
                  labels=labels)


    result@buyerReserve  <- calcOptimalReserve(result,preMerger=TRUE) #Find Buyer Reserve
    result@pricePre      <- calcPrices(result,preMerger=TRUE) # Calculate premerger expected price
    result@pricePost     <- calcPrices(result,preMerger=FALSE) # Calculate postmerger expected price

    return(result)
    }


#test <- auction2nd.cap(capacities=c(.5,.5),
#         buyerCost=1,
#         sellerCostCDF=punif,
#         rangeCostDist=0:1,
#         ownerPre=factor(0:1),ownerPost=factor(rep(0,2)),
#         labels=as.character(0:1)
#         )

