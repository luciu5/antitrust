\name{linear}
\alias{linear}
\alias{loglog}

\title{Linear and Log-Log Demand Calibration and Merger Simulation}

\description{Calibrates consumer demand using either a
  linear or log-log demand system and then simulates the prices effect of a merger between two firms
  under the assumption that all firms in the market are playing a differentiated products Bertrand game.}

\usage{


       linear(prices,quantities,margins,
                     diversion=NULL,
                     ownerPre,ownerPost,
                     mcDelta=rep(0,length(prices)),
                     labels=paste("Prod",1:length(prices),sep="")
                     )


       loglog(prices,quantities,margins,
                     diversion=NULL,
                     ownerPre,ownerPost,
                     mcDelta=rep(0,length(prices)),
                     priceStart=prices,
                     labels=paste("Prod",1:length(prices),sep="")
}

\arguments{
  Let k denote the number of products produced by all firms.
  \item{prices}{A length k vector of product prices}
  \item{quantities}{A length k vector of product quantities}
  \item{margins}{A length k vector of product margins}
  \item{diversion}{A k x k matrix of diversion ratios with diagonal
  elements equal to 1. Default is NULL, in which case diversion
  according to quantity share is assumed.}
  \item{ownerPre}{{EITHER a vector of length k whose values
      indicate which firm produced a product pre-merger OR
      a k x k matrix of pre-merger ownership shares}
  \item{ownerPost}{EITHER a vector of length k whose values
    indicate which firm produced a product after the merger OR
    a k x k matrix of post-merger ownership shares}
  \item{mcDelta}{A vector of length k where each element equals the
    proportional change in a product's marginal costs due to
    the merger. Default is 0, which assumes that the merger does not
  affect any products' marginal cost.}
 \item{priceStart}{A length k vector of prices used as the initial guess
  in the nonlinear equation solver. Default is 'prices'.}
  \item{labels}{A k-length vector of labels. Default is "Prod#", where
  '#' is a number between 1 and the length of 'prices'}
}



\details{Using price, quantity, margin, and diversion information for all products
  in a market, 'linear' is able to
  recover the slopes and intercepts in a Linear Demand
  System and then use these demand parameters to simulate a merger between
  two firms. Likewise, 'loglog' uses the same information to uncover the
  slopes and intercepts in a Log-Log demand system, and then uses these
  demand parameters to simulate a merger.

  The elements of 'prices' and 'quantities' must be positive, and the
  elements of 'margins', and 'mcDelta'  must be between 0 and 1.

  'diversion' either equal NULL or a square matrix whose elements are be between
  -1 and 1. If 'diversion' equals NULL, then diversion according to quantity share is assumed. If
  a square matrix is supplied, the off-diagonal elements [i,j] of 'diversion'
  must equal an estimate of the diversion ratio from product i to product j
  (i.e. the estimated fraction of i's sales that go to j due to a small
  increase in i's price). Also, 'diversion'
  elements are positive if i and j are substitutes and negative if i and j are
  complements. The Diagonal elements must all equal 1.



  'ownerPre' and 'ownerPost' values will typically be equal to either 0
  (element [i,j] is not commonly owned) or 1 (element [i,j] is commonly
  owned), though these matrices may take on any value between 0 and 1 to
  account for partial ownership.




}

\value{'linear' returns an object  of class 'Linear' a child class of
  'Antitrust'. 'loglog' returns an object of
  class 'LogLog', which is a child class of 'Linear'

  For all classes, the 'print' function produces the predicted percentage price change
  due to the merger. The 'summary' function provides a more detailed
  analysis of the merger simulation.

  An object of class 'Linear' is a class instance containing the
  following slots:

  mc:        a length-k vector of  estimated marginal cost
  mcDelta:   see above
  slopes:    a k x k matrix of linear demand slope coefficients
  intercept: a length-k vector of linear demand intercepts
  ownerPre : see above
  ownerPost: see above
  labels   : see above
  diversion: see above

 An object of class 'LogLog' contains all of the above slots as well as

  priceStart: see above
}






\seealso{
  'elast' for calculating a k x k matrix of price elasticities
  'deltaCS' for calculating a length-k vector of changes in consumer
  surplus
  'cmcr' for calculating the compensating marginal cost reduction
  }

\author{Charles Taragin <email: charles.Taragin@usdoj.gov>}

\references{}


\examples{

## Simulate a merger between two single-product firms in a
## three-firm market with linear demand with diversions
## that are proportional to shares.
## This example assumes that the merger is between
## the first two firms



n <- 3 #number of firms in market
price    <- c(2.9,3.4,2.2)
quantity <- c(650,998,1801)


slopes <- matrix(
c(-2.3,	  0.18,	0.28,
   0.11, -2.4,	0.1,
   0.13,   .16,-2.7),ncol=n)


margin <- -1/diag(slopes)


#simulate merger between firms 1 and 2
owner.pre <- diag(n)
owner.post <- owner.pre
owner.post[1,2] <- owner.post[2,1] <- 1



result <- pclinear(price,quantity,margin,ownerPre=owner.pre,ownerPost=owner.post)

print(result)           # return predicted price change
summary(result)         # summarize merger simulation

elast(result,TRUE)      # returns premerger elasticities
elast(result,FALSE)     # returns postmerger elasticities

diversion(result,TRUE)  # returns premerger diversion ratios
diversion(result,FALSE) # returns postmeger diversion ratios

cmcr(result)            # returns the compensating marginal cost reduction




## Same as above, but with log-log instead of linear demand


result <- pcloglog(price,quantity,margin,ownerPre=owner.pre,ownerPost=owner.post)
print(result)





}