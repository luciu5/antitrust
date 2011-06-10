\name{logit}
\alias{logit}
\alias{logit.nests}
\title{(Nested) Logit Demand Calibration and Merger Simulation}

\description{Calibrates consumer demand using (Nested) Logit and then simulates the price effect of a merger between two firms
  under the assumption that all firms in the market are playing a
  differentiated products Bertrand pricing game.}

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
  \item{shares}{A length k vector of product shares. Values must be
  between 0 and 1.}
  \item{quantities}{A length k vector of product quantities.}
  \item{margins}{A length k vector of product margins, some of which may
    equal NA.}
  \item{nests}{A length k vector identifying the nest that each
    product belongs to.}
   \item{ownerPre}{{EITHER a vector of length k whose values
      indicate which firm produced a product pre-merger OR
      a k x k matrix of pre-merger ownership shares.}
  \item{ownerPost}{EITHER a vector of length k whose values
    indicate which firm produced a product after the merger OR
    a k x k matrix of post-merger ownership shares.}
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
  \item{...}{Additional options to feed to the \code{\link[nleqslv]{nleqslv}}
  optimizer used solve the model}

}



\details{Using product prices, shares and (at least) all of the
  product margins from any firm, 'logit' is able to
  recover the price coefficient and product mean valuations in a
  Logit demand model. 'logit' then uses these
  calibrated parameters to simulate a merger between two firms.

  'logit.nests' is identical to 'logit' except that it includes the 'nests'
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






\seealso{\code{\link{ces}}}

\author{Charles Taragin <email: charles.taragin@usdoj.gov>}

\references{}


\examples{

## Calibration and simulation results from a merger between Budweiser and
## Old Style.
## Source: Epstein/Rubenfeld 2004, pg 80


prodNames <- c("BUD","OLD STYLE","MILLER","MILLER-LITE","OTHER-LITE","OTHER-REG")
ownerPre <-c("BUD","OLD STYLE","MILLER","MILLER","OTHER-LITE","OTHER-REG")
ownerPost <-c("BUD","BUD","MILLER","MILLER","OTHER-LITE","OTHER-REG")
nests <- c("Reg","Reg","Reg","Light","Light","Reg")

price    <- c(.0441,.0328,.0409,.0396,.0387,.0497)
shares   <- c(.066,.172,.253,.187,.099,.223)
margins <- c(.3830,.5515,.5421,.5557,.4453,.3769)

names(price) <-
    names(shares) <-
    names(margins) <-
    prodNames

result.logit <- logit(price,shares,margins,ownerPre=ownerPre,ownerPost=ownerPost)



print(result.logit)           # return predicted price change
summary(result.logit)         # summarize merger simulation

elast(result.logit,TRUE)      # returns premerger elasticities
elast(result.logit,FALSE)     # returns postmerger elasticities

diversion(result.logit,TRUE)  # return premerger diversion ratios
diversion(result.logit,FALSE) # return postmerger diversion ratios


cmcr(result.logit)            #calculate compensating marginal cost reduction

CV(result.logit)              #calculate average compensating variation


}