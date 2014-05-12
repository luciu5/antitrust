\name{aids}
\alias{aids}
\alias{pcaids}
\alias{pcaids.nests}
\title{(Nested) AIDS Calibration and Merger Simulation}

\description{Calibrates consumer demand using (nested) AIDS and then simulates the price effect of a merger between two firms
  under the assumption that all firms in the market are playing a differentiated products Bertrand game.}

\usage{


aids(shares,margins,prices,diversions,
                 ownerPre,ownerPost,
                 mcDelta=rep(0, length(shares)),
                 subset=rep(TRUE, length(shares)),
                 priceStart=runif(length(shares)),
                 isMax=FALSE,
                 labels=paste("Prod",1:length(shares),sep=""),
                 ...)

pcaids(shares,knownElast,mktElast=-1,
              prices,diversions,
              ownerPre,ownerPost,
              knownElastIndex=1,
              mcDelta=rep(0, length(shares)),
              subset=rep(TRUE, length(shares)),
              priceStart=runif(length(shares)),
              isMax=FALSE,
              labels=paste("Prod",1:length(shares),sep=""),
              ...)

pcaids.nests(shares,margins,knownElast,mktElast=-1,
              prices,ownerPre,ownerPost,
              nests=rep(1,length(shares)),
              knownElastIndex=1,
              mcDelta=rep(0, length(shares)),
              subset=rep(TRUE, length(shares)),
              priceStart=runif(length(shares)),
              isMax=FALSE,
              nestsParmStart,
              labels=paste("Prod",1:length(shares),sep=""),
              ...)

}

\arguments{
  \item{}{Let k denote the number of products produced by all firms.}
  \item{shares}{A length k vector of product revenue shares. All shares must
    be between 0 and 1.}
  \item{margins}{A length k vector of product margins. All margins must
    be either be between 0 and 1, or NA.}
  \item{prices}{A length k vector product prices. Default is missing, in
  which case demand intercepts are not calibrated.}
  \item{knownElast}{A negative number equal to the pre-merger own-price
    elasticity for any of the k products.}
  \item{mktElast}{A negative number equal to the industry pre-merger
    price elasticity. Default is -1.}
  \item{diversions}{A k x k matrix of diversion ratios with diagonal
    elements equal to -1. Default is missing, in which case diversion
    according to revenue share is assumed.}
  \item{ownerPre}{EITHER a vector of length k whose values
      indicate which firm produced a product before the merger OR a
      k x k matrix of pre-merger ownership shares. }
  \item{ownerPost}{EITHER a vector of length k whose values
    indicate which firm produced a product after the merger OR
    a k x k matrix of post-merger ownership shares.}
  \item{knownElastIndex}{An integer equal to the position of the
    \sQuote{knownElast} product in the \sQuote{shares} vector. Default is 1, which
    assumes that the own-price elasticity of the first product is
    known. }
  \item{nests}{A length k vector identifying which nest a product
  belongs to. Default is that all products belong to a single nest.}
  \item{mcDelta}{A vector of length k where each element equals the
    proportional change in a product's marginal costs due to
    the merger. Default is 0, which assumes that the merger does not
    affect any products' marginal cost.}
  \item{subset}{A vector of length k where each element equals TRUE if
  the product indexed by that element should be included in the
  post-merger simulation and FALSE if it should be excluded.Default is a
  length k vector of TRUE.}
  \item{priceStart}{A vector of length k who elements equal to an
  initial guess of the proportional change in price caused by the
  merger. The default is to draw k random elements from a [0,1] uniform
  distribution.}
\item{isMax}{If TRUE, checks to see whether computed price equilibrium
  locally maximizes firm profits and returns a warning if not. Default is FALSE.}
\item{nestsParmStart}{A vector of starting values used to solve for
  price coefficient and nest parameters. If missing then the random
  draws with the appropriate restrictions are employed.}
  \item{labels}{A k-length vector of labels.}
  \item{...}{Additional options to feed to the \code{\link[BB]{BBsolve}}
  optimizer used to solve for equilibrium prices.}

}



\details{Using product market revenue shares and all of the
  product product margins from at least two firms, \code{aids} is able to
  recover the slopes in a proportionally calibrated Almost Ideal Demand System (AIDS)
  without income effects. \code{aids} then uses these slopes to simulate the
  price effects of a merger between
  two firms under the assumption that all firms in the market are
  playing a differentiated Bertrand pricing game.

  If prices are also supplied, \code{aids} is able to recover the
  intercepts from the AIDS demand system. Intercepts are helpful because
  they can be used to simulate pre- and post-merger price \emph{levels} as well as price
  \emph{changes}. Whatsmore, the intercepts are necessary in order to
  calculate compensating variation.


  \code{aids} assumes that diversion between the products
  in the market occurs according to revenue share. This assumption may be relaxed
  by setting \sQuote{diversions} equal to a k x k matrix of diversion
  ratios. The diagonal of this matrix must equal -1, the off-diagonal
  elements must be between 0 and 1, and the rows must sum to 1.

   \code{pcaids} is almost identical to \code{aids}, but instead of
   assuming that at least two margins are known, \code{pcaids} assumes
  that the own-price elasticity of any single
  product, and the industry-wide own-price elasticity, are
  known. Demand intercepts cannot be recovered using \code{pcaids}.



  \code{pcaids.nests} extends \code{pcaids} by allowing products
  to be grouped into nests. Although products within the same nest still
  have the independence of irrelevant alternatives (IIA) property,
  products in different nests do not. Note that the \sQuote{diversions}
  argument is absent from \code{pcaids.nests}.

  \code{pcaids.nests}  assumes that the share diversion between nests is symmetric
 (i.e for 2 nests A and B, the diversion from A to B is the same as B to
 A). Therefore, if there are \eqn{w} nests, \eqn{2\le w \le k}{2 <= w <=k}, then the model
 must estimate \eqn{w(w-1)/2} distinct nesting parameters. To accomplish
 this, \code{pcaids.nests} uses margin information to produce estimates of
 the nesting parameters. It is important to note that the number of
 supplied margins must be at least as great as the number of nesting
 parameters in order for PCAIDS to work.

 The nesting parameters are constrained to be between 0 and
 1. Therefore, one way to test the validity of the nesting structure is
 to check whether the nesting parameters are between 0 and 1. The value
 of the nesting parameters may be obtained from calling either the \sQuote{summary} or
 \sQuote{getNestsParms} functions.



}

\value{\code{aids} returns an instance of class
  \code{\linkS4class{AIDS}}, a child class of \code{\linkS4class{Linear}}.
  \code{pcaids} returns an instance of class
  \code{\linkS4class{PCAIDS}}, while \code{pcaids.nests} returns an
  instance of \code{\linkS4class{PCAIDSNests}}. Both are
  children of the \code{\linkS4class{AIDS}} class.
}



\seealso{\code{\link{linear}} for a demand system based on quantities
  rather than revenue shares.}

\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\references{
  Epstein, Roy and  Rubinfeld, Daniel (2004).
\dQuote{Merger Simulation with Brand-Level Margin Data: Extending PCAIDS
  with Nests.}
\emph{The B.E. Journal of Economic Analysis \& Policy}, \bold{advances.4}(1), pp. 2.

  Epstein, Roy and Rubinfeld, Daniel (2004).
  \dQuote{Effects of Mergers Involving Differentiated Products.}
}


\examples{

## Simulate a merger between two single-product firms A and B in a
## three-firm market (A, B, C). This example assumes that the merger is between
## the firms A and B and that A's own-price elasticity is
## known.
## Source: Epstein and Rubinfeld (2004), pg 9, Table 2.


prices    <- c(2.9,3.4,2.2) ## optional for aids, unnecessary for pcaids
shares    <- c(.2,.3,.5)


## The following are used by aids but not pcaids
## only two of the margins are required to calibrate the demand parameters
margins <- c(0.33, 0.36, 0.44)

## The following are used by pcaids, but not aids
knownElast<- -3
mktElast  <- -1


## Define ownership using a vector of firm identities
ownerPre <- c("A","B","C")
ownerPost <- c("A","A","C")

## Alternatively, ownership could be defined using matrices
#ownerPre=diag(1,length(shares))
#ownerPost=ownerPre
#ownerPost[1,2] <- ownerPost[2,1] <- 1


## AIDS: the following assumes both prices and margins are known.
##       Prices are not needed to estimate price changes



result.aids <- aids(shares,margins,prices,ownerPre=ownerPre,ownerPost=ownerPost,labels=ownerPre)



print(result.aids)           # return predicted price change
summary(result.aids)         # summarize merger simulation

elast(result.aids,TRUE)      # returns premerger elasticities
elast(result.aids,FALSE)     # returns postmerger elasticities

diversion(result.aids,TRUE)  # return premerger diversion ratios
diversion(result.aids,FALSE) # return postmerger diversion ratios


cmcr(result.aids)            #calculate compensating marginal cost reduction
upp(result.aids)             #calculate Upwards Pricing Pressure Index


## Implement the Hypothetical Monopolist Test
## for products A and B using a 5\% SSNIP

HypoMonTest(result.aids,prodIndex=1:2)


CV(result.aids)        #calculate compensating variation as a percent of
                       #representative consumer income
                       #CV can only be calculated if prices are supplied

CV(result.aids,14.5e12)     #calculate compensating variation in dollars
                            #14.5e12 is an estimate of total US GDP


## Get a detailed description of the 'AIDS' class slots
showClass("AIDS")

## Show all methods attached to the 'AIDS' Class
showMethods(classes="AIDS")

## Show which class have their own 'elast' method
showMethods("elast")

## Show the method definition for 'elast' and Class 'AIDS'
getMethod("elast","AIDS")



## PCAIDS: the following assumes that only one product's elasticity is
##       known as well as the market elasticity.



result.pcaids <- pcaids(shares,knownElast,mktElast,
                        ownerPre=ownerPre,ownerPost=ownerPost,
                        labels=ownerPre)

print(result.pcaids)           # return predicted price change
summary(result.pcaids)         # summarize merger simulation

elast(result.pcaids,TRUE)      # returns premerger elasticities
elast(result.pcaids,FALSE)     # returns postmerger elasticities

diversion(result.pcaids,TRUE)  # return premerger diversion ratios
diversion(result.pcaids,FALSE) # return postmerger diversion ratios


cmcr(result.pcaids)            #calculate compensating marginal cost reduction


## Implement the Hypothetical Monopolist Test
## for products A and B using a 5\% SSNIP

HypoMonTest(result.pcaids,prodIndex=1:2)



## Nested PCAIDS: in addition to the PCAIDS information requirements,
##                users must supply the nesting structure as well as margin infromation.

nests <- c('H','L','L') # product A assigned to nest H, products B and C assigned to nest L


result.pcaids.nests <- pcaids.nests(shares,knownElast,mktElast,margins=margins,
                                    nests=nests,ownerPre=ownerPre,
                                    ownerPost=ownerPost,labels=ownerPre)




}
