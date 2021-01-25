#' @title Multi-product Cournot/Stackelberg Calibration and Merger Simulation With Linear or Log-Linear Demand
#' @name Cournot-Functions
#' @aliases cournot
#' Cournot
#' stackelberg
#' Stackelberg
#' @description Calibrates consumer demand for multiple products using either a
#' linear or log-linear demand system and then simulates the
#' prices effect of a merger between two multi-plant firms
#' under the assumption that all firms in the market are
#' playing either a Cournot or Stackelberg quantity setting game.
#' @description Let k denote the number of products and n denote the number of plants below.
#'
#' @param prices A length k vector product prices.
#' @param quantities An n x k matrix of product quantities.
#' All quantities must either be positive, or if the product is not produced by a plant, NA.
#' @param margins An n x k  matrix of product margins. All margins must
#' be either be between 0 and 1, or NA.
#' @param demand A length k character vector equal to "linear" if a product's
#' demand curve is assumed to be linear or "log" if a product's demand curve
#' is assumed to be log-linear.
#' @param cost A length n character vector equal to "linear" if a plant's
#' marginal cost curve is assumed to be linear or "constant" if a plant's marginal curve
#' is assumed to be constant. Returns an error if a multi-plant firm with constant marginal costs
#' does not have capacity constraints.
#' @param isLeaderPre An n x k logical matrix equal to TRUE if a firm is
#' a "leader" pre-merger for a particular product and FALSE otherwise.
#' Default is FALSE, which is equivalent to \code{cournot}.
#' @param isLeaderPost An n x k logical matrix equal to TRUE if a firm
#' is a "leader" pre-merger for a particular product and FALSE otherwise.
#' Default is FALSE, which is equivalent to \code{cournot}.
#' @param mcfunPre a length n list of functions that calculate a plant's
#' pre-merger marginal cost. If empty (the default), assumes  quadratic costs.
#' @param mcfunPost a length n list of functions that calculate a plant's
#' post-merger marginal cost. If empty (the default), equals \sQuote{mcfunPre}
#' @param vcfunPre a length n list of functions that calculate a plant's
#' pre-merger variable cost. If empty (the default), assumes quadratic variable costs.
#' @param vcfunPost a length n list of functions that calculate a plant's
#' post-merger variable cost. If empty (the default), equals \sQuote{vcfunPre}
#' @param dmcfunPre a length n list of functions that calculate the derivative
#' of a plant's pre-merger marginal cost. If empty (the default),
#' assumes  quadratic variable costs.
#' @param dmcfunPost a length n list of functions that calculate the derivative
#' of a plant's post-merger marginal cost. If empty (the default), equals \sQuote{mcfunPre}
#' @param capacitiesPre A length n numeric vector of pre-merger
#' plant capacities. Default is Inf.
#' @param capacitiesPost A length n numeric vector of post-merger
#' plant capacities. Default \sQuote{capacitiesPre}.
#' @param productsPre An n x k  matrix that equals TRUE if pre-merger,
#' a plant produces a product. Default is TRUE if 'quantities' is not NA.
#' @param productsPost An n x k  matrix that equals TRUE if post-merger,
#' a plant produces a product. Default equals \sQuote{productsPre}.
#' @param ownerPre EITHER a vector of length n whose values
#' indicate which plants are commonly owned pre-merger OR
#' an n x n matrix of pre-merger ownership shares.
#' @param ownerPost EITHER a vector of length n whose values
#' indicate which plants will be commonly owned after the merger OR
#' an n x n matrix of post-merger ownership shares.
#' @param mktElast A length k vector of product elasticities. Default is a length k vector of NAs
#' @param mcDelta A length n vector where each element equals the
#' proportional change in a firm's marginal costs due to
#' the merger. Default is 0, which assumes that the merger does not
#' affect any products' marginal cost.
#' @param quantityStart A length k vector of quantities used as the initial guess
#' in the nonlinear equation solver. Default is \sQuote{quantities}.
#' @param control.slopes A list of  \code{\link{optim}}  control parameters
#' passed to the calibration routine optimizer
#' (typically the \code{calcSlopes} method).
#' @param control.equ A list of  \code{\link[BB]{BBsolve}} control parameters
#' passed to the non-linear equation solver
#' (typically the \code{calcQuantities} method).
#' @param labels A list with 2 elements. The first element is a
#' vector of firm names, while the second element is a vector of products names.
#' Default is \sQuote{O1:On}, and \sQuote{P1:Pk}.
#' @param ... Additional options to feed to the solver. See below.
#'
#' @details Using price, and quantity, information for all products
#' in each market, as well as margin information for at least
#' one products in each market, \code{cournot} is able to
#' recover the slopes and intercepts of either a Linear or Log-linear demand
#' system as well as the cost parameters (see below for further details). 
#' These parameters are then used to simulate the price
#' effects of a merger between
#' two firms under the assumption that the firms are playing a
#' homogeneous products simultaneous quantity setting game.
#'
#' \code{stackelberg}, is similar to \code{cournot}, except that for a given product,
#' firms are either "leaders" or "followers". leaders gain a first mover advantage
#' over followers, which allows the leaders to anticipate how changes to their output
#' will effect the follower's output decisions. Firms can be the leader for
#' some products but the follower in others.
#'
#'
#' \sQuote{mcfunPre} and \sQuote{mcfunPost} are length n lists whose elements are \sQuote{R} functions 
#' that return a firm's marginal cost. The first argument of each function should be total firm quantities.
#'  By default, each firm is assumed to have quadratic costs with a firm-specific parameter calibrated from a firm's margin.
#'    \sQuote{vcfunPre} and \sQuote{vcfunPost} are similarly defined. 
#'    \sQuote{dmcfunPre} and \sQuote{dmcfunPost} are the changes in marginal cost and are only required for \code{stackelberg}.
#'
#'   \sQuote{ownerPre} and \sQuote{ownerPost} values will typically be equal to either 0
#'   (element [i,j] is not commonly owned) or 1 (element [i,j] is commonly
#'   owned), though these matrices may take on any value between 0 and 1 to
#'   account for partial ownership.
#'
#'
#'   Under linear demand and linear marginal costs, an analytic solution to the Cournot quantity game
#'   exists. However, this solution can at times produce negative
#'   equilibrium quantities. To accommodate this issue, \code{cournot}
#'   uses \code{\link[BB]{BBsolve}}  to
#'   find equilibrium quantities subject to a non-negativity constraint. \code{...} may
#'   be used to change the default options for \code{\link[BB]{BBsolve}}.
#'
#' @return \code{cournot} returns an instance of class \code{\linkS4class{Cournot}}.
#' \code{stackelberg} returns an instance of class \code{\linkS4class{Stackelberg}}.
#' @author Charles Taragin \email{ctaragin@ftc.gov}
#'
#' @examples ## Simulate a Cournot merger between two single-plant firms
#' ## producing a single product in a
#' ## 5-firm market with linear demand and quadratic costs
#'
#'
#'
#' n <- 5 #number of firms in market pre-merger
#' cap <- rnorm(n,mean = .5, sd = .1)
#' int <- 10
#' slope <- -.25
#'
#' B.pre.c = matrix(slope,nrow=n,ncol=n)
#' diag(B.pre.c) = 2* diag(B.pre.c) - 1/cap
#' quantity.pre.c = rowSums(solve(B.pre.c) * -int)
#' price.pre.c = int + slope * sum(quantity.pre.c)
#' mc.pre.c = quantity.pre.c/cap
#' vc.pre.c = quantity.pre.c^2/(2*cap)
#' margin.pre.c = 1 - mc.pre.c/price.pre.c
#' ps.pre.c = price.pre.c*quantity.pre.c - vc.pre.c
#'
#' mktQuant.pre.c = sum(quantity.pre.c)
#'
#' ## suppose firm 1 acquires firm 2
#' ## This model has a closed form solution
#' B.post.c = B.pre.c
#' B.post.c[1,2] = 2*B.post.c[1,2]
#' B.post.c[2,1] = 2*B.post.c[2,1]
#'
#' quantity.post.c = rowSums(solve(B.post.c) * -int)
#' price.post.c = int + slope * sum(quantity.post.c)
#' mc.post.c = quantity.post.c/cap
#' vc.post.c = quantity.post.c^2/(2*cap)
#' margin.post.c = 1 - mc.post.c/price.post.c
#' ps.post.c = price.post.c*quantity.post.c - vc.post.c
#'
#' mktQuant.post.c = sum(quantity.post.c, na.rm=TRUE)
#'
#' #check if merger is profitable for merging parties
#' isprofitable.c = ps.post.c - ps.pre.c
#' isprofitable.c= sum(isprofitable.c[1:2]) > 0
#'
#'
#' #prep inputs for Cournot
#' owner.pre <- diag(n)
#' owner.post <- owner.pre
#' owner.post[1,2] <- owner.post[2,1] <- 1
#'
#'
#'
#' result.c <- cournot(prices = price.pre.c,quantities = as.matrix(quantity.pre.c),
#'                     margins=as.matrix(margin.pre.c),
#'                     ownerPre=owner.pre,ownerPost=owner.post)
#'
#' print(result.c)           # return predicted price change
#' summary(result.c)         # summarize merger simulation
#'
#' ## check if 'cournot' yields the same result as closed-form solution
#' #print(all.equal(sum(result.c@quantityPre) , mktQuant.pre.c))
#' #print(all.equal(sum(result.c@quantityPost) , mktQuant.post.c))
#'
#'
#'
#'
#' ## Simulate a Stackelberg merger between two single-plant firms
#' ## producing a single product in a
#' ## 5-firm market with linear demand and quadratic costs.
#' ## Allow both merging parties to be followers pre-merger,
#' ## but assume that they become leaders post-merger.
#' ## Finally, assume that pre-merger, there is a single leader who ## remains a leader post-merger
#' ## Note: This example uses setup from the above Cournot example
#'
#' isLeader.pre = matrix(rep(FALSE,n), ncol=1)
#' isLeader.pre[n,] = TRUE
#' isLeader.post = isLeader.pre
#' isLeader.post[1:2,] = TRUE
#'
#' passthru.pre = matrix(-slope^2/(2*slope - 1/cap))
#' passthru.post = passthru.pre
#' passthru.pre[isLeader.pre] = 0
#' passthru.post[isLeader.post] = 0
#'
#' B.pre.s = matrix(slope,nrow=n,ncol=n)
#' diag(B.pre.s) = 2* diag(B.pre.s) - 1/cap
#' diag(B.pre.s)[n] = diag(B.pre.s)[n] + sum(passthru.pre)
#'
#' quantity.pre.s = rowSums(solve(B.pre.s) * ( -int))
#' price.pre.s = int + slope * sum(quantity.pre.s)
#' mc.pre.s = quantity.pre.s/cap
#' vc.pre.s = quantity.pre.s^2/(2*cap)
#' margin.pre.s = 1 - mc.pre.s/price.pre.s
#' ps.pre.s = price.pre.s*quantity.pre.s - vc.pre.s
#'
#' mktQuant.pre.s = sum(quantity.pre.s)
#'
#' ## suppose firm 1 acquires firm 2
#' ## This model has a closed form solution
#' B.post.s = matrix(slope,nrow=n,ncol=n)
#' diag(B.post.s) = 2* diag(B.post.s) - 1/cap
#' B.post.s[1,2] = 2*B.post.s[1,2]
#' B.post.s[1,1:2] = B.post.s[1,1:2]  + sum(passthru.post)
#' B.post.s[2,1] = 2*B.post.s[2,1]
#' B.post.s[2,1:2] = B.post.s[2,1:2]  + sum(passthru.post)
#' diag(B.post.s)[n] = diag(B.post.s)[n]  +  sum(passthru.post)
#'
#' quantity.post.s = rowSums(solve(B.post.s) * as.vector( -int ) )
#' price.post.s = int + slope * sum(quantity.post.s)
#' mc.post.s = quantity.post.s/cap
#' vc.post.s = quantity.post.s^2/(2*cap)
#' margin.post.s = 1 - mc.post.s/price.post.s
#' ps.post.s = price.post.s*quantity.post.s - vc.post.s
#'
#' mktQuant.post.s = sum(quantity.post.s, na.rm=TRUE)
#'
#' #check if merger is profitable for merging parties
#' isprofitable.s = ps.post.s - ps.pre.s
#' isprofitable.s = sum(isprofitable.s[1:2]) > 0
#'
#'
#' #prep inputs for Stackelberg
#' owner.pre <- diag(n)
#' owner.post <- owner.pre
#' owner.post[1,2] <- owner.post[2,1] <- 1
#'
#'
#'
#' result.s <- stackelberg(prices = price.pre.s,quantities = as.matrix(quantity.pre.s),
#'                         margins=as.matrix(margin.pre.s),ownerPre=owner.pre,
#'                         ownerPost=owner.post,
#'                         isLeaderPre = isLeader.pre, isLeaderPost = isLeader.post)
#'
#' print(result.s)           # return predicted price change
#' summary(result.s)         # summarize merger simulation
#'
#' ## check if 'stackelberg' yields the same result as closed-form solution
#' #print(all.equal(sum(result.s@quantityPre) , mktQuant.pre.s))
#' #print(all.equal(sum(result.s@quantityPost) , mktQuant.post.s))
#'
#' @include CMCRCournotFunctions.R

#'@rdname Cournot-Functions
#'@export
cournot <- function(prices,quantities,
                    margins = matrix(NA_real_ , nrow(quantities),ncol(quantities)),
                    demand = rep("linear",length(prices)),
                    cost   =   rep("linear",nrow(quantities)),
                    mcfunPre=list(),
                    mcfunPost=mcfunPre,
                    vcfunPre=list(),
                    vcfunPost=vcfunPre,
                    capacitiesPre = rep(Inf,nrow(quantities)),
                    capacitiesPost = capacitiesPre,
                    productsPre=!is.na(quantities),
                    productsPost=productsPre,
                    ownerPre,ownerPost,
                    mktElast = rep(NA_real_, length(prices)),
                    mcDelta =rep(0,nrow(quantities)),
                    quantityStart=as.vector(quantities),
                    control.slopes,
                    control.equ,
                    labels,
                    ...
){

  shares <- as.vector(quantities/sum(quantities))



  if(missing(labels)){
    if(is.null(dimnames(quantities))){
      rname <- paste0("O",1:nrow(quantities))
      cname <- paste0("P",1:ncol(quantities))
    }
    else{rname <- rownames(quantities)
    cname <- colnames(quantities)
    }
    labels <- list(rname,cname)

  }

  result <- new("Cournot",prices=prices, quantities=quantities,margins=margins,
                shares=shares,mcDelta=mcDelta, subset= rep(TRUE,length(shares)), demand = demand, cost=cost,
                mcfunPre=mcfunPre, mcfunPost=mcfunPost,vcfunPre=vcfunPre, vcfunPost=vcfunPost,
                capacitiesPre=capacitiesPre,capacitiesPost=capacitiesPost,
                ownerPre=ownerPre, mktElast = mktElast,productsPre=productsPre,productsPost=productsPost,
                ownerPost=ownerPost, quantityStart=quantityStart,labels=labels)


  if(!missing(control.slopes)){
    result@control.slopes <- control.slopes
  }

  
  if(!missing(control.equ)){
    result@control.equ <- control.equ
  }
  
  ## Convert ownership vectors to ownership matrices
  result@ownerPre  <- ownerToMatrix(result,TRUE)
  result@ownerPost <- ownerToMatrix(result,FALSE)

  ## Calculate Demand Slope Coefficients and Intercepts
  result <- calcSlopes(result)


  result@quantityPre  <- calcQuantities(result, preMerger = TRUE,...)
  result@quantityPost <- calcQuantities(result,preMerger = FALSE,...)

  result@pricePre  <- calcPrices(result, preMerger = TRUE)
  result@pricePost <- calcPrices(result,preMerger = FALSE)

  return(result)

}

#'@rdname Cournot-Functions
#'@export
stackelberg <- function(prices,quantities,margins,
                        demand = rep("linear",length(prices)),
                        cost   =   rep("linear",nrow(quantities)),
                        isLeaderPre = matrix(FALSE,ncol = ncol(quantities), nrow= nrow(quantities)),
                        isLeaderPost= isLeaderPre,
                        mcfunPre=list(),
                        mcfunPost=mcfunPre,
                        vcfunPre=list(),
                        vcfunPost=vcfunPre,
                        dmcfunPre=list(),
                        dmcfunPost=dmcfunPre,
                        capacitiesPre = rep(Inf,nrow(quantities)),
                        capacitiesPost = capacitiesPre,
                        productsPre=!is.na(quantities),
                        productsPost=productsPre,
                        ownerPre,ownerPost,
                        mcDelta =rep(0,nrow(quantities)),
                        quantityStart=as.vector(quantities),
                        control.slopes,
                        control.equ,
                        labels,
                        ...
){

  shares <- as.vector(quantities/sum(quantities))



  if(missing(labels)){
    if(is.null(dimnames(quantities))){
      rname <- paste0("O",1:nrow(quantities))
      cname <- paste0("P",1:ncol(quantities))
    }
    else{rname <- rownames(quantities)
    cname <- colnames(quantities)
    }
    labels <- list(rname,cname)

  }

  result <- new("Stackelberg",prices=prices, quantities=quantities,margins=margins,
                shares=shares,mcDelta=mcDelta, subset= rep(TRUE,length(shares)), demand = demand, cost = cost,
                mcfunPre=mcfunPre, mcfunPost=mcfunPost,vcfunPre=vcfunPre, vcfunPost=vcfunPost,
                dmcfunPre=dmcfunPre, dmcfunPost=dmcfunPost, isLeaderPre = isLeaderPre, isLeaderPost = isLeaderPost,
                ownerPre=ownerPre,productsPre=productsPre,productsPost=productsPost,
                capacitiesPre=capacitiesPre,capacitiesPost=capacitiesPost,
                ownerPost=ownerPost, quantityStart=quantityStart,labels=labels)


  if(!missing(control.slopes)){
    result@control.slopes <- control.slopes
  }

  
  if(!missing(control.equ)){
    result@control.equ <- control.equ
  }

  ## Convert ownership vectors to ownership matrices
  result@ownerPre  <- ownerToMatrix(result,TRUE)
  result@ownerPost <- ownerToMatrix(result,FALSE)

  ## Calculate Demand Slope Coefficients and Intercepts
  result <- calcSlopes(result)


  result@quantityPre  <- calcQuantities(result, preMerger = TRUE,...)
  result@quantityPost <- calcQuantities(result,preMerger = FALSE,...)

  result@pricePre  <- calcPrices(result, preMerger = TRUE)
  result@pricePost <- calcPrices(result,preMerger = FALSE)

  return(result)

}

