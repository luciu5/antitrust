#' @title Merger Simulation With Debt
#' @name Debt-Functions
#' @aliases debt logit.debt
#' @description Calibrates consumer demand using Logit
#' and then simulates the price effect of a merger between two firms
#' under the assumption that all firms in the market are playing a
#' differentiated products Bertrand pricing game with Debt.
#' @description Let f denote the number of products. Let n denote the number of markets
#' @param shares A n x f  matrix of product (quantity) shares. Values must be
#'   between 0 and 1 and rows should sum to 1.
#' @param prices A n x f  matrix of product prices.
#' @param margins A n x f  matrix of product margins in proportions (bounded between 0 and 1), some of which may
#'   equal NA.
#' @param ownerPre A vector of length f whose values
#'   indicate which firm produced a product pre-merger.
#' @param ownerPost A vector of length f whose values
#'   indicate which firm produced a product post-merger.
#' @param mcDelta An n x f matrix of where each element equals the
#'   proportional change in a downstream firm's product-level marginal costs due to
#'     the merger. Default is 0, which assumes that the merger does not
#'     affect any products' marginal cost.
#' @param insideSize An integer equal to total pre-merger units sold.
#'   If shares sum to one, this also equals the size of the market.
#' @param subset A n x f  matrix  where each element equals TRUE if
#' the product indexed by that element should be included in the
#' post-merger simulation and FALSE if it should be excluded.Default is a
#' matrix where each element is TRUE.
#' @param priceOutside A length 1 vector indicating the price of the
#' outside good. Default is 0.
#' @param priceStart A n x f  matrix of starting values used to solve for
#' downstream equilibrium prices. Default is the \sQuote{prices} matrix.
#' @param control.slopes A list of  \code{\link{optim}}
#' control parameters passed to the calibration routine
#' optimizer (typically the \code{calcSlopes} method).
#' @param control.equ A list of  \code{\link[BB]{BBsolve}} control parameters
#' passed to the non-linear equation solver (typically the \code{calcPrices} method).
#' @param labels A k-length vector of labels. Default is "Prod#", where
#' \sQuote{#} is a number between 1 and the length of \sQuote{prices}.
#' @param ... Additional options to feed to the \code{\link[BB]{BBsolve}}
#' optimizer used to solve for equilibrium prices.
#'
#' @details Using product prices, quantity shares and all of the
#' product margins from at least one firm, \code{logit} is able to
#' recover the price coefficient and product mean valuations in a
#' Logit demand model. \code{logit} then uses these
#' calibrated parameters to simulate a merger between two firms.
#'
#'
#' @return When 'supplyDown' equals "bertand", \code{vertical.barg} returns an instance of class
#' \code{\linkS4class{VertBargBertLogit}}. When 'supplyDown' equals "2nd", \code{vertical.barg} returns an instance of class
#' \code{\linkS4class{VertBarg2ndLogit}} 
#' @author Charles Taragin \email{ctaragin+antitrustr@gmail.com}
#' @references 
#' @examples 
#' ## 3 single-product firms, each selling in 2 markets. 
#' ## Firms 1 and 2 merge
#' nprods <- 3
#' nmarkets <- 2
#' marketId <- rep(1:nmarkets,each=nprods)
#' shares <- c(0.0734, 0.6339, 0.2927, 0.5557, 0.3224, 0.1219)
#' prices <- c(107.3837, 156.9111, 131.0415, 118.6787, 108.4132, 99.1364)
#' margins <-c(23.2785, 66.7227, 21.971, 57.1583, 26.6393, 52.3467)/prices
#' ownerPre <- ownerPost <- as.character(rep(c(1, 2, 3),each=nmarkets))
#' ownerPost[ownerPost=='2'] <- '1'
#' debt <- c(14.6943, 30.3841, 16.3995)
#' insideSize <- 1
#'
#'

#' simres_debt <- logit.debt(prices,shares,margins,
#'                              ownerPre,ownerPost,debtPre=debt,marketID=marketId,
#'                              insideSize=insideSize)
#' 
#' 
#' print(simres_debt)
#' summary(simres_debt)
#' 

#' @include DebtClasses.R
NULL

#' @rdname Debt-Functions
#' @export

logit.debt <- function(prices,shares,margins,
                       ownerPre,ownerPost,
                       debtPre = rep(0, length(unique(ownerPre))),
                       debtPost = debtPre, 
                       marketID=as.character(rep(1,length(prices))),
                       productID,
                       insideSize,
                       density=c("beta","uniform"),
                       normIndex=1,
                       focal=1,
                       mcDelta=rep(0, length(prices)),
                       subset=rep(TRUE, length(prices)),
                       priceOutside=0,
                       priceStart = prices,
                       parmsStart=c(rep(-1, length(unique(marketID))),c(1,1)),
                       shareOutParm=c(1,1),
                       isMax=FALSE,
                       control.slopes,
                       control.equ,
                       labels,
                       ...
){
  
  density=match.arg(density)
  if(missing(productID)){
    nFirmProds <- tapply(prices,ownerPre,length)
    productID <- unlist(sapply(nFirmProds,function(x){return(1:x)}),use.names = FALSE)
    productID <- paste0(ownerPre,productID)
  }
  
  mktData <- data.frame(marketID,productID,prices,shares,margins,ownerPre,ownerPost,mcDelta,subset)
  
  ## fill in all missing market/owner/product combinations 
  ## necessary for reshaping 
  allPerms <- expand.grid(marketID=unique(marketID),
                          ownerPre=unique(ownerPre),
                          productID=unique(productID)
  )
  
  mktData <- merge(allPerms,mktData)
  
  mktData <- mktData[with(mktData,order(marketID,ownerPre,productID)),]

  if(missing(labels)){labels <- with(mktData,paste(marketID,ownerPre,productID,sep=":"))}
  nMarkets= length(unique(marketID))
  
  hFocalIntegral <- function(u,s1=shareOutParm[1],s2=shareOutParm[2]){
    beta(s1+1,s2+1)/beta(s1,s2)*
      pbeta(u,shape1=s1+1,shape2=s2+1,lower.tail = TRUE)
  }
  
  gFocalIntegral <- function(u,s1=shareOutParm[1],s2=shareOutParm[2]){
    beta(s1,s2+1)/beta(s1,s2)*
      pbeta(u,shape1=s1,shape2=s2+1,lower.tail = TRUE)
  }
  tOtherIntegral <- function(u){NULL}
  rOtherIntegral <- function(u){NULL}
  
  if(nMarkets >1 && density=="uniform"){
    
    mktCnt=2
    hFormula <- "u^3 -u^4/2~u"
    gFormula <- "u^2 -u^3/3~u"
    tFormula <- "u~u"
    
    while (mktCnt<=nMarkets) {
      
      h_anti_deriv <- deparse1(body(mosaicCalc::antiD(as.formula(hFormula))))
      h_anti_deriv <- gsub("\\s*\\+\\s*C","",h_anti_deriv,perl = TRUE)
      hFormula <- paste0(h_anti_deriv,"~u")
      
      g_anti_deriv <- deparse1(body(mosaicCalc::antiD(as.formula(gFormula))))
      g_anti_deriv <- gsub("\\s*\\+\\s*C","",g_anti_deriv,perl = TRUE)
      gFormula <- paste0(g_anti_deriv,"~u")
      
      t_anti_deriv <- deparse1(body(mosaicCalc::antiD(as.formula(tFormula))))
      t_anti_deriv <- gsub("\\s*\\+\\s*C","",t_anti_deriv,perl = TRUE)
      tFormula <- paste0(t_anti_deriv,"~u")
      
      mktCnt <- mktCnt+1
    }
    
    
    hFocalIntegral <- eval(parse(text=paste0("function(u){","beta(2,2)*",h_anti_deriv,"}")))
    gFocalIntegral <- eval(parse(text=paste0("function(u){","beta(1,2)*",g_anti_deriv,"}")))
    tOtherIntegral <- eval(parse(text=paste0("function(u){","beta(2,2)*",t_anti_deriv,"}")))
    rOtherIntegral <- eval(parse(text=paste0("function(u){","beta(1,2)*",t_anti_deriv,"}")))
  }
  
  
  
  ## Create LogitDebt  container to store relevant data
  result <- new("LogitDebt",prices=mktData$prices, shares=mktData$shares,
                margins=mktData$margins,
                debtPre=debtPre,
                debtPost=debtPost,
                insideSize=insideSize,
                normIndex=normIndex,
                focal=focal,
                density=density,
                shareOutParm=shareOutParm,
                ownerPre=mktData$ownerPre,
                ownerPost=mktData$ownerPost,
                mcDelta=mktData$mcDelta,
                subset=mktData$subset,
                marketID=as.character(mktData$marketID),
                productID=as.character(mktData$productID),
                #priceOutside=priceOutside,
                priceStart=priceStart,
                parmsStart=parmsStart,
                hFocalIntegral=hFocalIntegral,
                gFocalIntegral=gFocalIntegral,
                tOtherIntegral=tOtherIntegral,
                rOtherIntegral=rOtherIntegral,
                labels=labels)
  
  if(!missing(control.slopes)){
    result@control.slopes <- control.slopes
  }
  if(!missing(control.equ)){
    result@control.equ <- control.equ
  }
  
  ## Convert ownership vectors to ownership matrices
  #result@ownerPre  <- ownerToMatrix(result,TRUE)
  #result@ownerPost <- ownerToMatrix(result,FALSE)
  return(result)
  ## Calculate Demand Slope Coefficients
  result <- calcSlopes(result)
  
  
  ## Calculate marginal cost
  result@mcPre <-  calcMC(result,TRUE)
  result@mcPost <- calcMC(result,FALSE)
  
  ## Solve Non-Linear System for Price Changes
  result@pricePre  <- calcPrices(result,preMerger=TRUE,isMax=isMax,...)
  result@pricePost <- calcPrices(result,preMerger=FALSE,isMax=isMax,subset=subset,...)
  
  
  return(result)
  
}
