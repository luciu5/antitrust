\name{auction2nd}
\alias{auction2nd}
\alias{auction2nd.cap}
\
\title{(Capacity Constrained) 2nd Price Auction Model}

\description{}

\usage{

auction2nd.cap(capacities,
               buyerCost=0,
               sellerCostCDF=punif,
               sellerCostPDF=dunif,
               rangeCostDist=c(0,1e12),
               ownerPre,ownerPost,
               labels=paste("Firm",1:length(capacities),sep="")
               )

             }
\arguments{
  Let k denote the number of firms bidding in the auction.
  \item{capacities}{A length k vector of product capacities.}
  \item{ownerPre}{EITHER a vector of length k whose values
    indicate which firm produced a product pre-merger OR
    a k x k matrix of pre-merger ownership shares.}
  \item{ownerPost}{EITHER a vector of length k whose values
    indicate which firm produced a product after the merger OR
    a k x k matrix of post-merger ownership shares.}
  \item{labels}{A k-length vector of labels. Default is "Firm", where
    `#' is a number between 1 and the length of `capacities'.}
}



\details{
}

\value{\code{auction2nd.cap} returns an instance of class
  \code{\linkS4class{Auction2ndCap}}.
}






\seealso{}

\author{Michael Sandfort \email{michael.sandfort@usdoj.gov}}

\references{
  Keith Waehrer and  Perry, Martin (2003).
  ``The Effects of Mergers in Open Auction Markets'',
  \emph{Rand Journal of Economics}, \bold{34(2)}, pp. 287-304.
}




\examples{

## Both bidders in a 2-bid auction merge
## Source: Waehrer and Perry (2003) pg. 297


cap = c(.5,.5)
bc = 1
costRange = 0:1
owner.pre = 0:1
owner.post = rep(0,2)

result <- auction2nd.cap(capacities=cap,
         buyerCost=bc,
         sellerCostCDF=punif,
         sellerCostPDF=dunif,
         rangeCostDist=costRange,
         ownerPre=owner.pre,ownerPost=owner.post
         )

print(result)           # return percent predicted price change in



}