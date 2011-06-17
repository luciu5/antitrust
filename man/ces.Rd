\name{ces}
\alias{ces}
\alias{ces.nests}
\title{(Nested) Constant Elasticity of Substitution Demand Calibration and Merger Simulation}

\description{Calibrates consumer demand using (Nested) Constant Elasticity of
  Substitution (CES) and then simulates the price effect of a merger between two firms
  under the assumption that all firms in the market are playing a
  differentiated products Bertrand pricing game.}

\usage{ces(prices,shares,margins,
                ownerPre,ownerPost,
                shareInside = 1,
                normIndex=1,
                mcDelta=rep(0,length(prices)),
                priceStart = prices,
                labels=paste("Prod",1:length(prices),sep=""),
                ...
                )

ces.nests(prices,shares,margins,
                ownerPre,ownerPost,
                nests=rep(1,length(shares)),
                shareInside = 1,
                normIndex=1,
                mcDelta=rep(0,length(prices)),
                priceStart = prices,
                parmsStart=NULL,
                labels=paste("Prod",1:length(prices),sep=""),
                ...
                )

             }

\arguments{
  \item{prices}{A length k vector of product prices.}
  \item{shares}{A length k vector of product revenue shares.}
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
  \item{shareInside}{The proportion that a typical consumer spends on
     all products included in the 'prices' vector. Only needed to calculate
     compensating variation.
     Default is 1, meaning that all of a consumer's
     income is spent on products within the market.}
  \item{normIndex}{An integer specifying the product index against which
     the mean values of all other products are normalized. Default is 1.}
  \item{mcDelta}{A vector of length k where each element equals the
    proportional change in a product's marginal costs due to
    the merger. Default is 0, which assumes that the merger does not
    affect any products' marginal cost.}
  \item{priceStart}{A length k vector of starting values used to solve for
    equilibrium price. Default is the 'prices' vector.}
  \item{labels}{A k-length vector of labels. Default is "Prod#", where
    '#' is a number between 1 and the length of 'prices'.}

}



\details{Using product prices, quantities and (at least) all of the
  product margins from any firm, 'ces' is able to
  recover the price coefficient and product mean valuations in a
  Constant Elasticity of Substitution demand model. 'ces' then uses these
  calibrated parameters to simulate a merger between two firms under the
  assumption that that all firms in the market are playing a
  differentiated products Bertrand pricing game.

  'ces.nests' is identical to 'ces' except that it includes the 'nests'
  argument which may be used to assign products to different
  nests. Nests are useful becuase they allow for richer substitution
  patterns between products. Products within the same nest are assumed
  to be closer substitutes than products in different nests. The degree
  of substitutability between products located in different nests is
  controlled by the value of the nesting parameter sigma. The vector of
  sigmas is calibrated from the prices, quantities, and margins supplied
  by the user.


}

\value{'ces' returns an instance  of class "CES", a child class of
  "Logit".  'ces.nests' returns an instance of class "CESNests", a child
  class of "CES".


  For all classes, the 'print' function produces the predicted percentage price change
  due to the merger. The 'summary' function provides a more detailed
  analysis of the merger simulation.

  An object of class 'CES' is a class instance containing the
  following slots:

  \item{priceStart}{See above.}
  \item{totExp}{See above.}
  \item{normIndex}{See above.}
  \item{prices}{See above.}
  \item{quantities}{See above.}
  \item{margins}{See above.}
  \item{pricesPre}{See above.}
  \item{pricesPost}{See above.}
  \item{shareInside}{See above.}
  \item{shares}{See above.}
  \item{mcDelta}{See above}
  \item{slopes} {a list containing the coefficient on the numeraire
  ('alpha'),  the coefficient on price ('gamma'), and the vector of mean
  valuations ('meanval')}
  \item{ownerPre}{See above.}
  \item{ownerPost}{See above.}
  \item{labels}{See above.}

  An object of class 'PCAIDSNests' contains all of the above slots as
  well as

  \item{nests}{See above.}
  \item{parmsStart}{See above.}

}






\seealso{}

\author{Charles Taragin <email: charles.Taragin@usdoj.gov>}

\references{}


\examples{


## Calibration and simulation results from a merger between Budweiser and
## Old Style.
## Source: Epstein/Rubenfeld 2004, pg 80

prodNames <- c("BUD","OLD STYLE","MILLER","MILLER-LITE","OTHER-LITE","OTHER-REG")
ownerPre <-c("BUD","OLD STYLE","MILLER","MILLER","OTHER-LITE","OTHER-REG")
ownerPost <-c("BUD","BUD","MILLER","MILLER","OTHER-LITE","OTHER-REG")
nests <- c("R","R","R","L","L","R")

price    <- c(.0441,.0328,.0409,.0396,.0387,.0497)
shares   <- c(.071,.137,.251,.179,.093,.269)
margins  <- c(.3830,.5515,.5421,.5557,.4453,.3769)

names(price) <-
    names(shares) <-
    names(margins) <-
    prodNames

result.ces <- ces(price,shares.quantity,margins.logit,ownerPre=ownerPre,ownerPost=ownerPost)

print(result.ces)           # return predicted price change
summary(result.ces)         # summarize merger simulation

elast(result.ces,TRUE)      # returns premerger elasticities
elast(result.ces,FALSE)     # returns postmerger elasticities

diversion(result.ces,TRUE)  # return premerger diversion ratios
diversion(result.ces,FALSE) # return postmerger diversion ratios


cmcr(result.ces)            #calculate compensating marginal cost reduction

CV(result.ces)              #calculate representative agent compensating variation






}