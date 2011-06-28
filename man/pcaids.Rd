\name{pcaids}
\alias{pcaids}
\alias{pcaids.nests}
\title{(Nested) PCAIDS Calibration and Merger Simulation}

\description{Calibrates consumer demand using (nested) PCAIDS and then simulates the price effect of a merger between two firms
  under the assumption that all firms in the market are playing a differentiated products Bertrand game.}

\usage{pcaids(shares,knownElast,mktElast=-1,
              diversions=NULL,
              ownerPre,ownerPost,
              knownElastIndex=1,
              mcDelta=rep(0, length(shares)),
              priceDeltaStart=runif(length(shares)),
              labels=paste("Prod",1:length(shares),sep=""),
              ...)

pcaids.nests(shares,margins,knownElast,mktElast=-1,
              ownerPre,ownerPost,
              nests=rep(1,length(shares)),
              knownElastIndex=1,
              mcDelta=rep(0, length(shares)),
              priceDeltaStart=runif(length(shares)),
              nestsParmStart = NULL,
              labels=paste("Prod",1:length(shares),sep=""),
              ...)

}

\arguments{
  Let k denote the number of products produced by all firms.
  \item{shares}{A length k vector of product revenue shares. All shares must
  be between 0 and 1.}
  \item{knownElast}{A negative number equal to the pre-merger own-price
    elasticity for any of the k products.}
  \item{mktElast}{A negative number equal to the industry pre-merger
    price elasticity. Default is -1.}
   \item{diversions}{A k x k matrix of diversion ratios with diagonal
    elements equal to 1. Default is NULL, in which case diversion
    according to revenue share is assumed.}
  \item{ownerPre}{EITHER a vector of length k whose values
      indicate which firm produced a product before the merger OR a
      k x k matrix of pre-merger ownership shares.}
  \item{ownerPost}{EITHER a vector of length k whose values
    indicate which firm produced a product after the merger OR
    a k x k matrix of post-merger ownership shares.}
   \item{margins}{A length k vector of product margins. NA values permitted.}
  \item{knownElastIndex}{An integer equal to the position of the
    'knownElast' product in the 'shares' vector. Default is 1, which
    assumes that the own-price elasticity of the first product is
    known. }
  \item{nests}{A length k vector identifying which nest a product
  belongs to. Default is that all products belong to a single nest.}
  \item{mcDelta}{A vector of length k where each element equals the
    proportional change in a product's marginal costs due to
    the merger. Default is 0, which assumes that the merger does not
  affect any products' marginal cost.}
  \item{priceDeltaStart}{A vector of length k who elements equal to an
  initial guess of the proportional change in price caused by the
  merger. The default is to draw k random elements from a [0,1] uniform
  distribution.}
\item{nestsParmStart}{A vector of starting values used to solve for
  price coefficient and nest parameters. If 'NULL' then the random
  draws with the appropriate restrictions are employed.}
  \item{labels}{A k-length vector of labels.}
  \item{...}{Additional options to feed to the \code{\link[nleqslv]{nleqslv}}
  optimizer used solve the model.}

}



\details{Using product market shares, the own-price elasticity of any single
  product, and the industry-wide own-price elasticity, \code{pcaids} is able to
  recover the slopes in a Proportionally Calibrated Almost Ideal Demand
  System (PCAIDS) and then use these slopes to simulate a merger between
  two firms under the assumption that all firms in the market are
  playing a differentiated Bertrand pricing game.

  By default, \code{pcaids} assumes that diversion between the products
  in the market occurs according to revenue share. This assumption may be relaxed
  by setting 'diversions' equal to a k x k matrix of diversion
  ratios. The diagonal of this matrix must equal 1.

  \code{pcaids.nests} extends \code{pcaids} by allowing products
  to be grouped into nests. Although products within the same nest still
  have the independence of irrelevant alternatives (IIA) property,
  products in different nests do not. Note that the 'diversions'
  argument is absent from \code{pcaids.nests}.

  \code{pcaids.nests}  assumes that the share diversion between nests is symmetric
 (i.e for 2 nests A and B, the diversion from A to B is the same as B to
 A). Therefore, if there are w nests, \eqn{2\le w \le k}{2 <= w <=k}, then the model
 must estimate \eqn{w(w-1)/2} distinct nesting parameters. To accomplish
 this, \code{pcaids.nests} uses margin information to produce estimates of
 the nesting parameters. It is important to note that the number of
 supplied margins must be at least as great as the number of nesting
 parameters in order for PCAIDS to work.

 The nesting parameters are constrained to be between 0 and
 1. Therefore, one way to test the validity of the nesting structure is
 to check whether the nesting parameters are between 0 and 1. The value
 of the nesting parameters may be obtained from calling either the 'summary' or
 'getNestsParms' functions.



}

\value{\code{pcaids} returns an instance of class \code{"\linkS4class{PCAIDS}"}.
  \code{pcaids.nests} returns an instance of \code{"\linkS4class{PCAIDSNests}"}, a
  child class of \code{"\linkS4class{PCAIDS}".}
}




\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\references{
  Epstein, Roy and  Rubinfeld, Daniel (2004).
``Merger Simulation with Brand-Level Margin Data: Extending PCAIDS
with Nests.''
\emph{The B.E. Journal of Economic Analysis \& Policy}, \bold{advances.4}(1), pp. 2.

  Epstein, Roy and Rubinfeld, Daniel (2004).
  ``Effects of Mergers Involving Differentiated Products.''
}


\examples{

## Simulate a merger between two single-product firms A and B in a
## three-firm market (A, B, C). This example assumes that the merger is between
## the firms A and B and that A's own-price elasticity is
## known.
## Source: Epstein and Rubinfeld (2004), pg 9, Table 2.

shares= c(.2,.3,.5)
knownElast=-3
mktElast=-1

## Define ownership using a vector of firm identities
ownerPre <- c("A","B","C")
ownerPost <- c("A","A","C")

## Alternatively, ownership could be defined using matrices
#ownerPre=diag(1,length(shares))
#ownerPost=ownerPre
#ownerPost[1,2] <- ownerPost[2,1] <- 1



result <- pcaids(shares,knownElast,mktElast,ownerPre=ownerPre,ownerPost=ownerPost,labels=ownerPre)

print(result)           # return predicted price change
summary(result)         # summarize merger simulation

elast(result,TRUE)      # returns premerger elasticities
elast(result,FALSE)     # returns postmerger elasticities

diversion(result,TRUE)  # return premerger diversion ratios
diversion(result,FALSE) # return postmerger diversion ratios


cmcr(result)            #calculate compensating marginal cost reduction

##Suppose that you also know pre-merger prices
price    <- c(2.9,3.4,2.2)

CV(result,price)        #calculate compensating variation






}