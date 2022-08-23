#' @title Bertrand Calibration and Merger Simulation With Logit, CES and AIDS Demand
#' @name Bertrand-Functions
#' @aliases Bertrand
#' bertrand
#' bertrand.alm
#'
#' @description Calibrates consumer demand using either a
#' Logit, CES, or AIDS demand system and then simulates the
#' prices effect of a merger between two firms
#' under the assumption that all firms in the market
#' are playing a Nash-Bertrand price setting game.
#' @description Let k denote the number of products produced by all firms below.
#'
#' @param demand A character vector indicating which demand system to use.
#' Currently allows logit (default), ces, or aids.
#' @param prices A length k vector product prices. Default is missing, in
#' which case demand intercepts are not calibrated.
#' @param quantities A length k vector of product quantities.
#' @param margins A length k vector of product margins. All margins must
#' be either be between 0 and 1, or NA.
#' @param diversions A k x k matrix of diversion ratios with diagonal
#' elements equal to -1. Default is missing, in which case diversion
#' according to revenue share is assumed.
#' @param ownerPre EITHER a vector of length k whose values
#' indicate which firm produced a product before the merger OR a
#' k x k matrix of pre-merger ownership shares.
#' @param ownerPost EITHER a vector of length k whose values
#' indicate which firm produced a product after the merger OR
#' a k x k matrix of post-merger ownership shares.
#' @param mktElast A negative number equal to the industry pre-merger
#' price elasticity. Default is NA.
#' @param insideSize Size of all units included in the market. For logit,
#' this defaults to total quantity, while for aids and ces
#' this defaults to total revenues.
#' @param mcDelta A vector of length k where each element equals the
#' proportional change in a product's marginal costs due to
#' the merger. Default is 0, which assumes that the merger does not
#' affect any products' marginal cost.
#' @param subset A vector of length k where each element equals TRUE if
#' the product indexed by that element should be included in the
#' post-merger simulation and FALSE if it should be excluded. Default is a
#' length k vector of TRUE.
#' @param parmStart \code{aids} only. A vector of length 2 who elements equal to an
#' initial guess for "known" element of the diagonal of the demand matrix and the market elasticity.
#' @param priceOutside A postive real number equal to the price of the outside good.
#' Default either equals 1 for Logit demand or 0 for CES demand.
#' @param priceStart A vector of length k who elements equal to an
#' initial guess of the proportional change in price caused by the
#' merger.  For aids, the default is to draw k random elements from
#' a [0,1] uniform distribution. For ces and logit, the default is prices.
#' @param isMax If TRUE, checks to see whether computed price equilibrium
#' locally maximizes firm profits and returns a warning if not. Default is FALSE.
#' @param control.slopes A list of  \code{\link{optim}}  control parameters passed
#' to the calibration routine optimizer (typically the \code{calcSlopes} method).
#' @param control.equ A list of  \code{\link[BB]{BBsolve}} control parameters passed
#' to the non-linear equation solver (typically the \code{calcPrices} method).
#' @param labels A k-length vector of labels.
#' @param ... Additional options to feed to the \code{\link[BB]{BBsolve}}
#' optimizer used to solve for equilibrium prices.
#'
#' @details
#' The main purpose of this function is to provide a more convenient front-end
#' for the \code{aids}, \code{logit.alm} and \code{ces} functions.
#'
#' Using price, and quantity, information for all products
#' in each market, as well as margin information for at least
#' one products in each market, \code{bertrand.alm} is able to
#' recover the slopes and intercepts of either a Logit, CES, or AIDS demand
#' system. These parameters are then used to simulate the price
#' effects of a merger between
#' two firms under the assumption that the firms are playing a
#' simultaneous price setting game.
#'
#'
#' \sQuote{ownerPre} and \sQuote{ownerPost} values will typically be equal to either 0
#' (element [i,j] is not commonly owned) or 1 (element [i,j] is commonly
#'                                             owned), though these matrices may take on any value between 0 and 1 to
#' account for partial ownership.
#'
#' @return \code{bertrand.alm} returns an instance of class \code{\linkS4class{LogitALM}},
#' \code{\linkS4class{CESALM}}, or \code{\linkS4class{AIDS}},
#' depending upon the value of the ``demand'' argument.
#' @author Charles Taragin \email{ctaragin+antitrustr@gmail.com}
#' @include Auction2ndLogitFunctions.R
NULL

#'@rdname Bertrand-Functions
#'@export
bertrand.alm <- function(
  demand = c("logit","ces","aids"),
  prices,quantities,margins,
  ownerPre,ownerPost,
  mktElast = NA_real_,
  insideSize = ifelse(demand == "logit",sum(quantities,na.rm=TRUE), sum(prices*quantities,na.rm=TRUE)),
  diversions,
  mcDelta=rep(0,length(prices)),
  subset=rep(TRUE,length(prices)),
  priceOutside=ifelse(demand== "logit",0, 1),
  priceStart=prices,
  isMax=FALSE,
  parmStart,
  control.slopes,
  control.equ,
  labels=paste("Prod",1:length(prices),sep=""),
  ...){


  demand <- match.arg(demand)


  shares_revenue <- shares_quantity <- quantities/sum(quantities)



  if(all(!is.na(prices))) shares_revenue <- prices*shares_quantity/sum(prices*shares_quantity)

  if(demand == "aids"){

    if(missing(prices)){ prices <- rep(NA_real_,length(shares_revenue))}

    if(missing(parmStart)) parmStart <- rep(NA_real_,2)

    if(missing(priceStart)) priceStart <- runif(length(shares_revenue))

    if(missing(diversions)){
      diversions <- tcrossprod(1/(1-shares_revenue),shares_revenue)
      diag(diversions) <- -1


    }

  }

  else if (demand %in% c("logit","ces")){

    if(missing(parmStart)){
      parmStart <- rep(.1,2)
      nm <- which(!is.na(margins))[1]
      if(demand == "logit"){
        parmStart[1] <- -1/(margins[nm]*prices[nm]*(1-shares_quantity[nm])) #ballpark alpha for starting values
      }
      else{parmStart[1] <- 1/(margins[nm]*(1-shares_revenue[nm])) - shares_revenue[nm]/(1-shares_revenue[nm])} #ballpark gamma for starting values
    }
    if(missing(priceStart)) priceStart <- prices
  }






  result <-   switch(demand,
                     aids=new("AIDS",shares=shares_revenue,mcDelta=mcDelta,subset=subset,
                              margins=margins, prices=prices, quantities=shares_revenue,  mktElast = mktElast,
                              insideSize = insideSize,
                              ownerPre=ownerPre,ownerPost=ownerPost, parmStart=parmStart,
                              diversion=diversions,
                              priceStart=priceStart,labels=labels),

                     logit=  new("LogitALM",prices=prices, shares=shares_quantity,
                                 margins=margins,
                                 ownerPre=ownerPre,
                                 ownerPost=ownerPost,
                                 mktElast = mktElast,
                                 mcDelta=mcDelta,
                                 subset=subset,
                                 priceOutside=priceOutside,
                                 priceStart=priceStart,
                                 shareInside= sum(shares_quantity),
                                 parmsStart=parmStart,
                                 insideSize = insideSize,
                                 labels=labels),

                     ces = new("CESALM",prices=prices, shares=shares_revenue,
                               margins=margins,
                               ownerPre=ownerPre,
                               ownerPost=ownerPost,
                               mktElast = mktElast,
                               mcDelta=mcDelta,
                               subset=subset,
                               priceOutside=priceOutside,
                               priceStart=priceStart,
                               shareInside=sum(shares_revenue),
                               parmsStart=parmStart,
                               insideSize =insideSize,
                               labels=labels),

                     linear=new("Linear",prices=prices, quantities=shares_quantity,margins=margins,
                                shares=shares_quantity,mcDelta=mcDelta, subset=subset,
                                ownerPre=ownerPre,diversion=diversions, symmetry=TRUE,
                                ownerPost=ownerPost, priceStart=priceStart,labels=labels)
  )


  if(!missing(control.slopes)){
    result@control.slopes <- control.slopes
  }
  if(!missing(control.equ)){
    result@control.equ <- control.equ
  }

  ## Convert ownership vectors to ownership matrices
  result@ownerPre  <- ownerToMatrix(result,TRUE)
  result@ownerPost <- ownerToMatrix(result,FALSE)

  ## Calculate Demand Slope Coefficients
  result <- calcSlopes(result)

  ## Solve Non-Linear System for Price Changes (AIDS only)
  if (demand == "aids"){
    result@priceDelta <- calcPriceDelta(result,isMax=isMax,subset=subset,...)
  }


  ## Calculate marginal cost
  result@mcPre <-  calcMC(result,TRUE)
  result@mcPost <- calcMC(result,FALSE)



  ## Solve Non-Linear System for Price Changes
  result@pricePre  <- calcPrices(result,preMerger=TRUE,isMax=isMax,...)
  result@pricePost <- calcPrices(result,preMerger=FALSE,isMax=isMax,subset=subset,...)

  return(result)

}
