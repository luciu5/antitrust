\name{Cournot}
\alias{cournot}
\alias{Cournot}

\title{Multi-product Cournot Calibration and Merger Simulation With Linear or Log-Linear Demand}

\description{Calibrates consumer demand for multiple products using either a
  linear or log-linear demand system and then simulates the prices effect of a merger between two firms
  under the assumption that all firms in the market are playing a Cournot quantity setting game.}

\usage{


       cournot(prices,quantities,margins, 
                    demand = rep("linear",length(prices)),
                    mcfunPre=list(),
                    mcfunPost=mcfunPre,
                    productsPre=!is.na(quantities), 
                    productsPost=productsPre, 
                    ownerPre,ownerPost,
                    mcDelta =rep(0,nrow(quantities)),
                    quantityStart=as.vector(quantities),
                    control.slopes,
                    control.eq,
                    labels,
                   ...
                     )

}

\arguments{
  \item{}{Let k denote the number of products and n denote the number of firms.}
  \item{prices}{A length k vector product prices.}
  \item{quantities}{An n x k vector of product quantities. All quantities must            either be positive, or  if the product is not produced by a firm, NA. }
  \item{margins}{An n x k  matrix of product margins. All margins must
    be either be between 0 and 1, or NA.}
  \item{demand}{A length k character vector equal to "linear" if a product's              demand curve is assumed to be linear or "log" if a product's demand curve         is assumed to be log-linear.}
   \item{mcfunPre}{a length n list of functions that calculate a firm's pre-merger marginal cost. If empty (the default), assumes  quadratic costs. }
   \item{mcfunPost}{a length n list of functions that calculate a firm's post-merger marginal cost. If empty (the default), equals \sQuote{mcfunPre} )}
   \item{productsPre}{An n x k  matrix that equals TRUE if pre-merger, a firms producrs a product. Default is TRUE if 'quantities' is not NA. }
    \item{productsPre}{An n x k  matrix that equals TRUE if pre-merger, a firms producrs a product. Default equals \sQuote{productsPre}. }
  \item{ownerPre}{EITHER a vector of length n whose values
    indicate which firms produced a product pre-merger OR
    an n x n matrix of pre-merger ownership shares.}
  \item{ownerPost}{EITHER a vector of length n whose values
    indicate which firms will be commonly owned after the merger OR
    an n x n matrix of post-merger ownership shares.}
  \item{mcDelta}{A length n vector where each element equals the
    proportional change in a firm's marginal costs due to
    the merger. Default is 0, which assumes that the merger does not
    affect any products' marginal cost.}
   \item{quantityStart}{A length k vector of quantities used as the initial guess
     in the nonlinear equation solver. Default is \sQuote{quantities}.}
      \item{control.slopes}{A list of  \code{\link{optim}}  control parameters passed to the calibration routine optimizer (typically the \code{calcSlopes} method).}
    \item{control.equ}{A list of  \code{\link[BB]{BBsolve}} control parameters passed to the non-linear equation solver (typically the \code{calcQuantities} method).}
   \item{labels}{A list with 2 elements. The first element is a vector of firm names, while the second element is a vector of products names. Default is \sQuote{O1:On}, and \sQuote{P1:Pk}.}
   \item{...}{Additional options to feed to the solver. See below.}
 }



\details{Using price, and quantity, information for all products
  in each market, as well as margin information for at least
  one products in each market, \code{cournot} is able to
  recover the slopes and intercepts of either a Linear or Log-linear demand
  system. These parameters are then used to simulate the price
  effects of a merger between
  two firms under the assumption that the firms are playing a
  homogeneous products quantity setting game.



\sQuote{mcfunPre} and \sQuote{mcfunPost} are length n lists whose elements are \sQuote{R} functions that return a firm's marginal cost. The first argument of each function should be total firm quantities. By default, each firm is assumed to have quadratic costs with a firm-specific parameter calibrated from a firm's margin.  
  
  \sQuote{ownerPre} and \sQuote{ownerPost} values will typically be equal to either 0
  (element [i,j] is not commonly owned) or 1 (element [i,j] is commonly
  owned), though these matrices may take on any value between 0 and 1 to
  account for partial ownership.


  Under linear demand and linear marginal costs, an analytic solution to the Cournot quantity game
  exists. However, this solution can at times produce negative
  equilibrium quantities. To accommodate this issue, \code{cournot}
  uses \code{\link[BB]{BBsolve}}  to
  find equilibrium quantities subject to a non-negativity constraint. \code{...} may
  be used to change the default options for \code{\link[BB]{BBsolve}}.

  

}

\value{\code{cournot} returns an instance of class \code{\linkS4class{Cournot}}.}






\seealso{
  }

\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\references{
von Haefen, Roger (2002).
\dQuote{A Complete Characterization Of The Linear, Log-Linear, And Semi-Log
  Incomplete Demand System Models.}
\emph{Journal of Agricultural and Resource Economics}, \bold{27}(02).
\url{http://ideas.repec.org/a/ags/jlaare/31118.html}.
}


\examples{

## (Mcafee and Williams 1992)
## Simulate a merger between two single-plant firms 
## producing a single product in a
## 5-firm market with linear demand 



n <- 4 #number of firms in market pre-merger
cap <- rnorm(n,mean = 2, sd = .5)
int <- 10
slope <- .25 # absolute value of demand slope

betas.pre <- slope * cap/( 1 + slope * cap)
B.pre <- sum(betas.pre)
price.pre <- int/(1 + B.pre)

quantity.pre <- int*(betas.pre)/(slope * (1 + B.pre))
mc.pre = quantity.pre/cap
margin.pre = 1 - mc.pre/price.pre

mktQuant.pre = sum(quantity.pre)

## suppose firm 1 acquires firm 2
## This model has a closed form solution 
betas.post = betas.pre 
betas.post[1] = (betas.post[1] + betas.post[2] - 2*betas.post[1]*betas.post[2])/(1 - betas.post[1]*betas.post[2])
betas.post[2] = NA
B.post = sum(betas.post, na.rm=TRUE)
price.post <- int/(1 + B.post)
quantity.post <- int*(betas.post)/(slope * (1 + B.post))

mktQuant.post = sum(quantity.post, na.rm=TRUE)



#prep inputs for Cournot
owner.pre <- diag(n)
owner.post <- owner.pre
owner.post[1,2] <- owner.post[2,1] <- 1



result.c <- cournot(price.pre,as.matrix(quantity.pre),as.matrix(margin.pre),ownerPre=owner.pre,ownerPost=owner.post)

print(result.c)           # return predicted price change
summary(result.c)         # summarize merger simulation




}
