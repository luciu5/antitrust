\name{upp}
\alias{upp}
\title{Upward Pricing Pressure (UPP) Index}

\description{Calculate the upward pricing presure index for each pair of
  products.}

\usage{upp(price, margin, diversion, ownerPre=NULL, efficiency=rep(0,length(price)), labels=paste("Prod",1:length(price),sep=""))}


\arguments{
  Let k_1 denote the number of products produced by merging firm 1 and
  k_2 denote the number of products produced by merging firm 2.
  \item{price}{A vector of product prices with length k_1 + k_2}
  \item{margin}{A vector of product margins with length k_1 + k_2}
  \item{diversion}{A (k_1 + k_2) x (k_1 + k_2) matrix of diversion ratios
    with diagonal elements equal to 1}
  \item{efficiency}{A vector with length k_1 + k_2  equal to the
  proportion of marginal cost efficiencies that a product will
  achieve. Default assumes no efficiencies}
  \item{ownerPre}{EITHER a vector of length k_1 + k_2 whose values
    indicate which of the merging parties produced a product pre-merger OR
    a (k_1 + k_2) x (k_1 + k_2) matrix whose elements (i,j) are NA if
    products i and j are owned by the same firm and 1 otherwise. NULL
  creates a matrix with diagonal elements equal to NA and off-diagonal
  elements equal to 1}
  \item{labels}{Value vector labels}
}

\details{All 'price' elements must be positive and all 'margin',
  'efficiency' and 'diversion' elements must be between 0 and 1. In
  addition, diagonal elements of 'diversion' must equal 1 and
  off-diagonal elements (i,j) of 'diversion'
  must equal the diversion ratio to product j from product i. Also, 'diversion'
  elements are positive if i and j are substitutes and negative if i and j are
  complements.

  'ownerPre' will typically be a vector set equal to 1 if a product is produced
  by firm 1 and 0 otherwise.  Optionally, 'ownerPre' may be set equal
  to a matrix  whose elements (i,j) equal NA if
  products i and j are owned by the same firm and 1 otherwise.
}

\value{A matrix whose values equal to the upwards pricing pressure
  index.}

\references{Farrell2010a}


\author{Charles Taragin}


\seealso{mcreduce}

\examples{

     ## Let k_1 = 1 and and k_2 = 2 ##

    p1 = 50;      margin1 = .3
    p2 = c(45,70); margin2 = c(.4,.6)
    isOne=c(1,0,0)
    diversions = matrix(c(1,.5,.01,.6,1,.1,.02,.2,1),ncol=3)

    upp(c(p1,p2), c(margin1,margin2), diversions, isOne)



}