\\name{pcaids.nests}
\alias{pcaids.nests}
\title{Nested PCAIDS Calibration and Merger Simulation}

\description{Calibrates consumer demand using Nested PCAIDS and then simulates the price effect of a merger between two firms
  under the assumption that all firms in the market are playing a differentiated products Bertrand game.}

\usage{pcaids(shares,margins,knownElast,mktElast=-1,
              ownerPre,ownerPost,
              nests=rep(1,length(shares)),
              knownElastIndex=1,
              mcDelta=rep(0, length(shares)),
              priceDeltaStart=runif(length(shares)),
              labels=paste("Prod",1:length(shares),sep=""),
              ...)}

\arguments{
  Let k denote the number of products produced by all firms.
  \item{shares}{A length k vector of product market shares.}
  \item{margins}{A length k vector of product margins. NA values permitted.}
  \item{knownElast}{A negative number equal to the pre-merger own-price
    elasticity for any of the k products.}
  \item{mktElast}{A negative number equal to the industry pre-merger
    price elasticity. Default is -1.}
  \item{ownerPre}{{EITHER a vector of length k whose values
      indicate which firm produced a product before the merger OR a
      k x k matrix of pre-merger ownership shares.}
  \item{ownerPost}{EITHER a vector of length k whose values
    indicate which firm produced a product after the merger OR
    a k x k matrix of post-merger ownership shares}
  \item{knownElastIndex}{An integer equal to the position number of the
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
  merger. The default is to draw k random elements from a [0,1] uniform distribution.}
  \item{labels}{A k-length vector of labels}
  \item{...}{Additional options to feed to the \code{\link[BB]{BBsolve}}
  optimizer used solve the model}

}



\details{Nested PCAIDS extends \code{\link{aids}} by allowing products
  to be grouped into nests. Although products within the same nest still
  have the independence of irrelevant alternatives (IIA) property,
  products in different nests do not.

 The model assumes that the share diversion between nests is symmetric
 (i.e for 2 nests A and B, the diversion from A to B is the same as B to
 A). Therefore, if there are $w$ nests, $2\le w \le k$, then the model
 must estimate $w(w-1)/2$ distinct nesting parameters. To accomplish
 this, Nested PCAIDS uses margin information to produce estimates of
 the nesting parameters. It is important to note that the number of
 supplied margins must be at least as great as the number of nesting
 parameters in order for PCAIDS to work.

 The nesting parameters are constrained to be between 0 and
 1. Therefore, one way to test the validity of the nesting structure is
 to check whether the nesting parameters are between 0 and 1. The value
 of the nesting parameters may be obtained from calling either summary or
 getNestsParms.
}

\value{An object of class "PCAIDSNests"}






\seealso{}

\author{Charles Taragin <email: charles.Taragin@usdoj.gov>}

\references{}


\examples{

## Simulate a merger between two single-product firms A and B in a
## three-firm market (A, B, C). This example assumes that the merger is between
## the firms A and B and that A's own-price elasticity is
## known.
## Source: Epstein and Rubinfeld (2004), pg 9, Table 2.
## Note: Epstein and Rubinfeld assume that the nesting parameter is .5,
## while this example calculates the nesting parameter from the provided
## margins

shares= c(.2,.3,.5)
nests = c (1,2,1)
margins=c(.333,.481,NA)
knownElast=-3
mktElast=-1

## Define ownership using a vector of firm identities
ownerPre <- c("A","B","C")
ownerPost <- c("A","A","C")

## Alternatively, ownership could be defined using matrices
#ownerPre=diag(1,length(shares))
#ownerPost=ownerPre
#ownerPost[1,2] <- ownerPost[2,1] <- 1


result <- pcaids.nests(shares,margins,knownElast,mktElast,ownerPre,ownerPost,nests)

print(result)           # return predicted price change
summary(result)         # summarize merger simulation

elast(result,TRUE)      # returns premerger elasticities
elast(result,FALSE)     # returns postmerger elasticities

diversion(result,TRUE)  # return premerger diversion ratios
diversion(result,FALSE) # return postmerger diversion ratios

deltaCS(result)         # returns the percent change in consumer surplus
                        # from the merger

cmcr(result)            #calculate compensating marginal cost reduction



## Simulate a merger between two firms in a
## 5-firm market.  Firm A produces 2 products, while the remaining firms
## each produce a single product. This example assumes that the merger is between
## Firm A and B and that the firm B's own-price elasticity may
## be uncovered from its margin
Source: Epstein and Rubinfeld (2004), pg 22, Table 8.



shares= c(.1,.075,.125,.15,.25,.3)
nests = c ("Popular","Prestige","Budget","Popular","Budget","Popular")
margins=c(.4,.55,.45,NA,NA,NA)
knownElast=-1/.45 # firm B's own elasticity
knownElastIndex = 2
mktElast=-1
ownerPre= c("A","A","B","C","D","E")
ownerPost=ownerPre
ownerPost[3] = "A"


result <- pcaids.nests(shares,margins,knownElast,mktElast,ownerPre,ownerPost,nests)



}