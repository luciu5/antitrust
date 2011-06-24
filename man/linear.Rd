\name{linear}
\alias{linear}
\alias{loglog}

\title{Linear and Log-Linear Demand Calibration and Merger Simulation}

\description{Calibrates consumer demand using either a
  linear or log-log demand system and then simulates the prices effect of a merger between two firms
  under the assumption that all firms in the market are playing a differentiated products Bertrand game.}

\usage{


       linear(prices,quantities,margins,
                     diversions=NULL,
                     ownerPre,ownerPost,
                     mcDelta=rep(0,length(prices)),
                     labels=paste("Prod",1:length(prices),sep="")
                     )


       loglog(prices,quantities,margins,
                     diversions=NULL,
                     ownerPre,ownerPost,
                     mcDelta=rep(0,length(prices)),
                     priceStart=prices,
                     labels=paste("Prod",1:length(prices),sep=""),
                     ...
             )
}

\arguments{
  Let k denote the number of products produced by all firms.
  \item{prices}{A length k vector product prices.}
  \item{quantities}{A length k vector of product quantities.}
  \item{margins}{A length k vector of product margins. All margins must
  be between 0 and 1.}
  \item{diversions}{A k x k matrix of diversion ratios with diagonal
  elements equal to 1. Default is NULL, in which case diversion
  according to quantity share is assumed.}
  \item{ownerPre}{EITHER a vector of length k whose values
      indicate which firm produced a product pre-merger OR
      a k x k matrix of pre-merger ownership shares.}
  \item{ownerPost}{EITHER a vector of length k whose values
    indicate which firm produced a product after the merger OR
    a k x k matrix of post-merger ownership shares.}
  \item{mcDelta}{A length k vector where each element equals the
    proportional change in a product's marginal costs due to
    the merger. Default is 0, which assumes that the merger does not
  affect any products' marginal cost.}
 \item{priceStart}{A length k vector of prices used as the initial guess
  in the nonlinear equation solver. Default is 'prices'.}
  \item{labels}{A k-length vector of labels. Default is "Prod#", where
  '#' is a number between 1 and the length of 'prices'.}
 \item{...}{Additional options to feed to the \code{\link[nleqslv]{nleqslv}}
  optimizer used solve the model.}
}



\details{Using price, quantity, margin, and diversion information for all products
  in a market, \code{linear} is able to
  recover the slopes and intercepts in a Linear Demand
  System and then use these demand parameters to simulate a merger between
  two firms under the assumption that the firms are playing a
  differentiated Bertrand pricing game.

  Similarly,  \code{loglog} uses the same information to uncover the
  slopes and intercepts in a Log-Log demand system, and then uses these
  demand parameters to simulate a merger two firms under the
  assumption that the firms are playing a
  differentiated Bertrand pricing game.


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

\value{\code{linear} returns an instance of class \code{"\linkS4class{Linear}"}.
  \code{loglog} returns an instance of \code{"\linkS4class{LogLog}"}, a
  child class of \code{"\linkS4class{Linear}"}.}






\seealso{
\code{\link{pcaids}} for a demand system that employs shares rather than quantities.
  }

\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\references{
von Haefen, Roger (2002).
``A Complete Characterization Of The Linear, Log-Linear, And Semi-Log
Incomplete Demand System Models.''
\emph{Journal of Agricultural and Resource Economics}, \bold{27}(02).
\url{http://ideas.repec.org/a/ags/jlaare/31118.html}.
}


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



result.linear <- linear(price,quantity,margin,ownerPre=owner.pre,ownerPost=owner.post)

print(result.linear)           # return predicted price change
summary(result.linear)         # summarize merger simulation

elast(result.linear,TRUE)      # returns premerger elasticities
elast(result.linear,FALSE)     # returns postmerger elasticities

diversion(result.linear,TRUE)  # returns premerger diversion ratios
diversion(result.linear,FALSE) # returns postmeger diversion ratios

cmcr(result.linear)            # returns the compensating marginal cost reduction




## Same as above, but with log-log instead of linear demand


result.log <- loglog(price,quantity,margin,ownerPre=owner.pre,ownerPost=owner.post)

print(result.log)
summary(result.log)         # summarize merger simulation

elast(result.log,TRUE)      # returns premerger elasticities
elast(result.log,FALSE)     # returns postmerger elasticities

diversion(result.log,TRUE)  # returns premerger diversion ratios
diversion(result.log,FALSE) # returns postmeger diversion ratios

cmcr(result.log)            # returns the compensating marginal cost reduction





}