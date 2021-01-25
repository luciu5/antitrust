#' @title Merger Simulation With User-Supplied Demand Parameters
#' @name Sim-Functions
#' @aliases sim
#' @description Simulates the price effects of a merger between two firms
#' with user-supplied demand parameters under the
#' assumption that all firms in the market are playing a
#' differentiated products Bertrand pricing game.
#' @description Let k denote the number of products produced by all firms below.
#'
#' @param prices A length k vector of product prices.
#' @param demand A character string indicating the type of demand system
#'   to be used in the merger simulation. Supported demand systems are
#'   linear (\sQuote{Linear}), log-linear(\sQuote{LogLin}), logit (\sQuote{Logit}), nested logit
#'   (\sQuote{LogitNests}), ces (\sQuote{CES}), nested CES (\sQuote{CESNests}) and capacity
#'   constrained Logit (\sQuote{LogitCap}).
#' @param demand.param  See Below.
#' @param ownerPre EITHER a vector of length k whose values
#'   indicate which firm produced a product pre-merger OR
#'   a k x k matrix of pre-merger ownership shares.
#' @param ownerPost EITHER a vector of length k whose values
#'   indicate which firm produced a product after the merger OR
#'   a k x k matrix of post-merger ownership shares.
#' @param nests A length k vector identifying the nest that each
#'   product belongs to. Must be supplied when \sQuote{demand} equals \sQuote{CESNests} and
#'   \sQuote{LogitNests}.
#' @param capacities A length k vector of product capacities. Must be
#'   supplied when \sQuote{demand} equals \sQuote{LogitCap}.
#' @param mcDelta A vector of length k where each element equals the
#'   proportional change in a product's marginal costs due to
#'   the merger. Default is 0, which assumes that the merger does not
#'   affect any products' marginal cost.
#' @param subset A vector of length k where each element equals TRUE if
#'   the product indexed by that element should be included in the
#'   post-merger simulation and FALSE if it should be excluded.Default is a
#'   length k vector of TRUE.
#' @param priceOutside A length 1 vector indicating the price of the
#'   outside good. This option only applies to the \sQuote{Logit} class and its child classes
#'   Default for \sQuote{Logit},\sQuote{LogitNests}, and \sQuote{LogitCap} is 0,
#'   and for \sQuote{CES} and \sQuote{CesNests} is 1.
#' @param priceStart A length k vector of starting values used to solve for
#'   equilibrium price. Default is the \sQuote{prices} vector for all values of
#'   demand except for \sQuote{AIDS}, which is set equal to a vector of 0s.
#' @param labels A k-length vector of labels. Default is \dQuote{Prod#}, where
#'   \sQuote{#} is a number between 1 and the length of \sQuote{prices}.
#' @param ... Additional options to feed to the
#'       optimizer used to solve for equilibrium prices.
#'
#' @details Using user-supplied demand parameters,
#' \code{sim} simulates the effects of a merger in a market where
#' firms are playing a differentiated products pricing game.
#'
#' If \sQuote{demand} equals \sQuote{Linear}, \sQuote{LogLin}, or
#' \sQuote{AIDS}, then \sQuote{demand.param} must be a
#' list containing \sQuote{slopes}, a k x k matrix of slope coefficients, and
#' \sQuote{intercepts}, a length-k vector of intercepts. Additionally, if
#' \sQuote{demand} equals \sQuote{AIDS}, \sQuote{demand.param} must contain \sQuote{mktElast}, an
#' estimate of aggregate market elasticity.  For \sQuote{Linear}
#' demand models, \code{sim} returns an error if any intercepts are
#' negative, and for both \sQuote{Linear}, \sQuote{LogLin}, and \sQuote{AIDS} models,  \code{sim}
#' returns an error if not all diagonal elements of the slopes matrix are
#' negative.
#'
#' If \sQuote{demand} equals \sQuote{Logit} or \sQuote{LogitNests}, then
#' \sQuote{demand.param} must equal a list containing
#' \itemize{
#'   \item{alpha}{The price coefficient.}
#'   \item{meanval}{A length-k vector of mean valuations \sQuote{meanval}. If
#'     none of the values of \sQuote{meanval} are zero, an outside good is assumed
#'     to exist.}
#' }
#' If demand equals \sQuote{CES} or \sQuote{CESNests}, then
#' \sQuote{demand.param} must equal a list containing
#'
#' \itemize{
#'   \item{gamma}{ The price coefficient,}
#'   \item{alpha}{The coefficient on the numeraire good. May instead be
#'     calibrated using \sQuote{shareInside},}
#'   \item{meanval}{A length-k vector of mean valuations \sQuote{meanval}. If
#'     none of the values of \sQuote{meanval} are zero, an outside good is assumed
#'     to exist,}
#'   \item{shareInside}{ The budget share of all products in the
#'     market. Default is 1, meaning that all consumer wealth is spent on
#'     products in the market. May instead be specified using \sQuote{alpha}.}
#'
#' }
#'
#' @return \code{sim} returns an instance of the class specified by the
#' \sQuote{demand} argument.
#' @seealso The S4 class documentation for: \code{\linkS4class{Linear}},
#' \code{\linkS4class{AIDS}}, \code{\linkS4class{LogLin}}, \code{\linkS4class{Logit}},
#' \code{\linkS4class{LogitNests}}, \code{\linkS4class{CES}}, \code{\linkS4class{CESNests}}
#' @author Charles Taragin \email{ctaragin@ftc.gov}
#'
#' @examples ## Calibration and simulation results from a merger between Budweiser and
#' ## Old Style. Note that the in the following model there is no outside
#' ## good; BUD's mean value has been normalized to zero.
#'
#' ## Source: Epstein/Rubenfeld 2004, pg 80
#'
#'
#' prodNames <- c("BUD","OLD STYLE","MILLER","MILLER-LITE","OTHER-LITE","OTHER-REG")
#' ownerPre <-c("BUD","OLD STYLE","MILLER","MILLER","OTHER-LITE","OTHER-REG")
#' ownerPost <-c("BUD","BUD","MILLER","MILLER","OTHER-LITE","OTHER-REG")
#' nests <- c("Reg","Reg","Reg","Light","Light","Reg")
#'
#' price    <- c(.0441,.0328,.0409,.0396,.0387,.0497)
#'
#' demand.param=list(alpha=-48.0457,
#'                   meanval=c(0,0.4149233,1.1899885,0.8252482,0.1460183,1.4865730)
#' )
#'
#' sim.logit <- sim(price,demand="Logit",demand.param,ownerPre=ownerPre,ownerPost=ownerPost)
#'
#'
#'
#' print(sim.logit)           # return predicted price change
#' summary(sim.logit)         # summarize merger simulation
#'
#' elast(sim.logit,TRUE)      # returns premerger elasticities
#' elast(sim.logit,FALSE)     # returns postmerger elasticities
#'
#' diversion(sim.logit,TRUE)  # return premerger diversion ratios
#' diversion(sim.logit,FALSE) # return postmerger diversion ratios
#'
#'
#' cmcr(sim.logit)            #calculate compensating marginal cost reduction
#' upp(sim.logit)             #calculate Upwards Pricing Pressure Index
#'
#' CV(sim.logit)              #calculate representative agent compensating variation
#'
#' @include LogitFunctions.R
NULL


#'@rdname Sim-Functions
#'@export
sim <- function(prices,demand=c("Linear","AIDS","LogLin","Logit","CES","LogitNests","CESNests","LogitCap"),demand.param,
                ownerPre,ownerPost,nests, capacities,
                mcDelta=rep(0,length(prices)),
                subset=rep(TRUE,length(prices)),
                priceOutside,
                priceStart,
                labels=paste("Prod",1:length(prices),sep=""),...){

  demand <- match.arg(demand)
  nprods <- length(prices)

  if(missing(priceStart)){
    if(demand=="AIDS"){priceStart <- runif(nprods)}
    else{              priceStart <- prices}
  }

  ## Create placeholders values to fill required Class slots

  shares <- margins <- rep(1/nprods,nprods)


  if(!missing(nests)){nests <- factor(nests,levels=unique(nests))}


  ## general checks
  if(!is.list(demand.param)){stop("'demand.param' must be a list.")}

  ## Checks for discrete choice models
  if(demand %in% c("CESNests","LogitNests","CES","Logit","LogitCap")){

    if(!("meanval" %in% names(demand.param))){
      stop("'demand.param' does not contain 'meanval'.")
    }
    if(length(demand.param$meanval) != nprods || any(is.na(demand.param$meanval))){
      stop("'meanval' must be a length-k vector of product mean valuations. NAs not allowed.")
    }

    if(demand %in% c("LogitNests","Logit","LogitCap")){

      ## An outside option is assumed to exist if all mean valuations are non-zero
      if(all(demand.param$meanval!=0)){
        normIndex <- NA
        shares <- rep(1/(nprods+1),nprods)
      }
      else{
        normIndex <- which(demand.param$meanval==0)

        if(length(normIndex)>1){
          warning("multiple values of meanval are equal to zero. Normalizing with respect to the first product with zero mean value.")
          normIndex <- normIndex[1]
        }

      }

      if(!("alpha" %in% names(demand.param))   ||
         length(demand.param$alpha) != 1     ||
         isTRUE(demand.param$alpha>0)){
        stop("'demand.param' does not contain 'alpha' or 'alpha' is not a negative number.")
      }

      shareInside <- sum(shares)
      if(missing(priceOutside)){priceOutside <- 0}

    }


    else  if(demand %in% c("CESNests","CES")){
      if(!("gamma" %in% names(demand.param))   ||
         length(demand.param$gamma) != 1     ||
         isTRUE(demand.param$gamma<0)){
        stop("'demand.param' does not contain 'gamma' or 'gamma' is not a positive number.")
      }


      ## uncover Numeraire Coefficients
      if(!("alpha" %in% names(demand.param)) &&
         !("shareInside" %in% names(demand.param))){
        warning("'demand.param' does not contain either 'alpha' or 'shareInside'. Setting shareInside=1 and alpha=NULL.")
        shareInside=1
        demand.param$alpha=NULL
      }

      else if("shareInside" %in% names(demand.param)){
        shareInside=demand.param$shareInside
        demand.param$shareInside <- NULL

        if(shareInside<1) {demand.param$alpha <- 1/shareInside -1}
        else{ demand.param$alpha <- NULL}


      }


      ## An outside option is assumed to exist if all mean valuations are non-zero
      if(all(demand.param$meanval!=1)){
        normIndex <- NA
        shares <- rep(1/(nprods+1),nprods)
      }
      else{
        normIndex <- which(demand.param$meanval==1)

        if(length(normIndex)>1){
          warning("multiple values of meanval are equal to one. Normalizing with respect to the first product with  mean value equal to 1.")
          normIndex <- normIndex[1]
        }

      }


      if(missing(priceOutside)){priceOutside <- 1}
    }

    if(demand %in% c("CESNests","LogitNests")){

      if(!("sigma" %in% names(demand.param))){
        stop("'demand.param' does not contain 'sigma'.")
      }

      if(missing(nests) ||
         length(nests)!= nprods ){stop("When 'demand' equals 'CESNests' or 'LogitNests', 'nests' must equal a vector whose length equals the number of products.")}

      if(length(demand.param$sigma)==1){
        constraint=TRUE
        demand.param$sigma <- rep(demand.param$sigma,nlevels(nests))
      }
      else{constraint=FALSE}


      if(nlevels(nests) != length(demand.param$sigma)){
        stop("The number of nests in 'nests' must either equal the number of nesting parameters in 'demand.param$sigma'.")}

    }


  }


  ## Checks for Linear-demand style models
  if(demand %in% c("Linear","LogLin","AIDS")){

    if(!("slopes" %in% names(demand.param))){stop("'demand.param' does not contain 'slopes'")}
    if(!("intercepts" %in% names(demand.param))){stop("'demand.param' does not contain 'intercepts'")}

    if(!(is.matrix(demand.param$slopes))   ||
       ncol(demand.param$slopes)!=nprods   ||
       nrow(demand.param$slopes)!=nprods   ||
       any(diag(demand.param$slopes)>0)){
      stop("'slopes' must be a k x k matrix of slope coeficients whose diagonal elements must all be negative.")}
    if(!is.vector(demand.param$intercepts)     ||
       length(demand.param$intercepts)!=nprods ||
       isTRUE(any(demand.param$intercepts<0,na.rm=TRUE))){
      stop("'intercepts' must be a length-k vector whose elements are all non-negative")
    }

    if (demand == "AIDS" &&
        !("mktElast" %in% names(demand.param))){
      warning("'demand.param' does not contain 'mktElast'. Setting 'mktElast' equal to -1")
      demand.param$mktElast=-1

    }

  }






  ## Create constructors for each demand system specified in the 'demand' parameter

  if(demand == "CESNests"){

    result <- new(demand,prices=prices, shares=shares,margins=margins,
                  mcDelta=mcDelta,
                  subset=subset,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  nests=nests,
                  normIndex=normIndex,
                  parmsStart=c(demand.param$gamma,demand.param$sigma),
                  priceStart=priceStart,
                  constraint=constraint,
                  shareInside=shareInside,labels=labels)

  }

  else if(demand == "LogitNests"){

    result <- new(demand,prices=prices, shares=shares,margins=margins,
                  mcDelta=mcDelta,
                  subset=subset,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  nests=nests,
                  normIndex=normIndex,
                  parmsStart=c(demand.param$alpha,demand.param$sigma),
                  priceStart=priceStart,
                  constraint=constraint,
                  shareInside=shareInside,labels=labels)

  }


  else if(demand %in% c("Logit","CES")){


    result <- new(demand,prices=prices, shares=shares,
                  margins=margins,
                  normIndex=normIndex,
                  mcDelta=mcDelta,
                  subset=subset,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  priceStart=priceStart,shareInside=shareInside,
                  labels=labels)

  }


  else if(demand == "LogitCap"){

    if(!("mktSize" %in% names(demand.param))){
      if(!missing(capacities) ){
        warning("'demand.param' does not contain 'mktSize'. Setting 'mktSize' equal to the sum of 'capacities'.")
        mktSize <- sum(capacities)
      }
      else{stop("'demand.param' does not contain 'mktSize'")}
    }
    else{mktSize <- demand.param$mktSize}


    shares <- capacities/mktSize
    shares <- shares/sum(shares)

    result <- new(demand, prices=prices, shares=shares,
                  margins=margins,capacities=capacities, mktSize=mktSize,
                  normIndex=normIndex,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  mcDelta=mcDelta,
                  subset=subset,
                  priceStart=priceStart,shareInside=shareInside,
                  labels=labels)
  }


  else if(demand == "Linear"){



    result <- new(demand,prices=prices, quantities=shares,margins=margins,
                  shares=shares,mcDelta=mcDelta,  subset=subset,
                  ownerPre=ownerPre,diversion=-diag(nprods),
                  symmetry=identical(demand.param$slopes,t(demand.param$slopes)),
                  ownerPost=ownerPost, priceStart=priceStart,labels=labels)

  }

  else if(demand == "AIDS"){

    ## find the market elasticity that best explains user-supplied intercepts and prices

    aidsShares    <- as.vector(demand.param$intercepts + demand.param$slopes %*% log(prices)) # AIDS needs actual shares for prediction
    aidsDiv       <- tcrossprod(1/(1-aidsShares),aidsShares)
    diag(aidsDiv) <- -1

    result <- new(demand,prices=prices, quantities=shares,margins=margins,
                  shares=aidsShares,
                  mcDelta=mcDelta,  subset=subset,mktElast=demand.param$mktElast,
                  ownerPre=ownerPre,diversion=aidsDiv,
                  priceStart=priceStart,
                  ownerPost=ownerPost, labels=labels)

  }



  else if(demand == "LogLin"){


    result <- new(demand,prices=prices, quantities=shares,margins=margins,
                  shares=shares,mcDelta=mcDelta, subset=subset, priceStart=priceStart,
                  ownerPre=ownerPre,diversion=-diag(nprods),
                  ownerPost=ownerPost, labels=labels)

  }


  if(demand %in% c("Linear","LogLin","AIDS")){
    result@slopes <- demand.param$slopes
    result@intercepts <- demand.param$intercepts
  }
  else{result@slopes=demand.param}


  ## Convert ownership vectors to ownership matrices
  result@ownerPre  <- ownerToMatrix(result,TRUE)
  result@ownerPost <- ownerToMatrix(result,FALSE)

  ## Calculate marginal cost
  result@mcPre     <-  calcMC(result,TRUE)
  result@mcPost    <-  calcMC(result,FALSE)

  if(demand == "AIDS"){
    ## Solve Non-Linear System for Price Changes
    result@priceDelta <- calcPriceDelta(result,...)
  }


  ## Solve Non-Linear System for Price Changes
  result@pricePre  <- calcPrices(result,TRUE,...)
  result@pricePost <- calcPrices(result,FALSE,subset=subset,...)


  return(result)
}
