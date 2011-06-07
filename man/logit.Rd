\name{logit}
\alias{logit}
\alias{logit.nests}
\title{(Nested) Logit Demand Calibration and Merger Simulation}

\description{Calibrates consumer demand using (Nested) Logit and then simulates the price effect of a merger between two firms
  under the assumption that all firms in the market are playing a differentiated products Bertrand game.}

\usage{logit(prices,shares,margins,
                ownerPre,ownerPost,
                normIndex=ifelse(sum(shares)<1,0,1),
                mcDelta=rep(0,length(prices)),
                priceStart = prices,
                labels=paste("Prod",1:length(prices),sep=""),
                ...
                )

logit.nests(prices,quantities,margins,
                ownerPre,ownerPost,
                nests=rep(1,length(shares)),
                shareInside = 1,
                normIndex=ifelse(shareInside < 1,0,1),
                mcDelta=rep(0,length(prices)),
                priceStart = prices,
                parmsStart=NULL,
                labels=paste("Prod",1:length(prices),sep=""),
                ...
                )

             }

\arguments{
  \item{prices}{A length k vector of product prices.}
  \item{quantities}{A length k vector of product quantities.}
  \item{margins}{A length k vector of product margins, some of which may
    equal NA.}
  \item{nests}{A length k vector identifying the nest that each
    product belongs to.}
   \item{ownerPre}{{EITHER a vector of length k whose values
      indicate which firm produced a product pre-merger OR
      a k x k matrix of pre-merger ownership shares}
  \item{ownerPost}{EITHER a vector of length k whose values
    indicate which firm produced a product after the merger OR
    a k x k matrix of post-merger ownership shares}
  \item{shareInside}{The share of customers that purchase any of the
  products included in the 'prices' vector. Default of 1 assumes
  that no other products are purchased.}
  \item{normIndex}{An integer specifying the product index against which
     the mean values of all other products are normalized. Default is 0
  if 'shareInside' is less than 1 and 1 if 'shareInside' equals 1.}
  \item{mcDelta}{A vector of length k where each element equals the
    proportional change in a product's marginal costs due to
    the merger. Default is 0, which assumes that the merger does not
    affect any products' marginal cost.}
  \item{priceStart}{A length k vector of starting values used to solve for
    equilibrium price. Default is the 'prices' vector.}
  \item{labels}{A k-length vector of labels. Default is "Prod#", where
    '#' is a number between 1 and the length of 'prices'.}

}



\details{Using product prices, shares and (at least) all of the
  product margins from any firm, 'logit' is able to
  recover the price coefficient and product mean valuations in a
  Logit demand model. 'logit' then uses these
  calibrated parameters to simulate a merger between two firms.

  'Logit.nests' is identical to 'ces' except that it includes the 'nests'
  argument which may be used to assign products to different
  nests. Nests are useful becuase they allow for richer substitution
  patterns between products. Products within the same nest are assumed
  to be closer substitutes than products in different nests. The degree
  of substitutability between products located in different nests is
  controlled by the value of the nesting parameter sigma. The vector of
  sigmas is calibrated from the prices, quantities, and margins supplied
  by the user.


}

\value{An object of class "Logit" or "LogitNests". "LogitNests" is a child
  class of "CES", which is a child class of "Antitrust".}






\seealso{}

\author{Charles Taragin <email: charles.Taragin@usdoj.gov>}

\references{}


\examples{

## Simulate a merger between two single-product firms A and B in a
## three-firm market (A, B, C). This example assumes that the merger is between
## the firms A and B.


price    <- c(2.9,3.4,2.2)
quantity <- c(650,998,1801)
margin   <- c(1/2.3,NA,NA)

## Define ownership using a vector of firm identities
ownerPre <- c("A","B","C")
ownerPost <- c("A","A","C")

## Alternatively, ownership could be defined using matrices
#ownerPre=diag(1,length(shares))
#ownerPost=ownerPre
#ownerPost[1,2] <- ownerPost[2,1] <- 1



result <-  ces(price,quantity,margin,ownerPre=owner.pre,ownerPost=owner.post)

print(result)           # return predcited price change
summary(result)         # summarize merger simulation

elast(result,TRUE)      # returns premerger elasticities
elast(result,FALSE)     # returns postmerger elasticities

diversion(result,TRUE)  # return premerger diversion ratios
diversion(result,FALSE) # return postmerger diversion ratios


cmcr(result)            #calculate compensating marginal cost reduction








}