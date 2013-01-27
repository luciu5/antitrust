\name{cmcr.cournot}
\alias{cmcr.cournot}
\alias{upp.cournot}
\title{Compensating Marginal Cost Reductions and Upwards Pricing
  Pressure (Cournot)}

\description{Calculate the average marginal cost reduction necessary to restore
  pre-merger prices, or the net Upwards Pricing Pressure in a two-product merger involving firms playing a
  homogeneous product Cournot pricing game.}

\usage{
cmcr.cournot(shares,mktElast)

upp.cournot(prices, margins, ownerPre,
             ownerPost=matrix(1,ncol=length(prices), nrow=length(prices)),
             mcDelta=rep(0,length(prices)),
             labels=paste("Prod",1:length(prices),sep=""))

}


\arguments{

  \item{shares}{A length-2 vector containing merging party quantity shares.}
  \item{mktElast}{A length-1 containing the industry elasticity.}
   \item{prices}{A length-2 vector of product prices.}
  \item{margins}{A length-2 vector of product margins.}
  \item{ownerPre}{EITHER a vector of length 2 whose values
  indicate which of the merging parties produced a product pre-merger OR
    a 2 x 2 matrix of pre-merger ownership shares.}
  \item{ownerPost}{A 2 x 2 matrix of post-merger
    ownership shares. Default is a 2 x 2 matrix of 1s.}
  \item{mcDelta}{A vector of length 2 where each element equals the
    proportional change in a product's marginal costs due to
    the merger. Default is 0, which assumes that the merger does not
    affect any products' marginal cost.}
  \item{labels}{A length-2 vector of product labels.}

}

\details{The `shares' vector must have 2 elements, and all `shares'
  elements must be between 0 and 1. The `mktElast' vector must
  have 1 non-negative  element.}
\value{A vector with 1 element whose value equals the percentage change
  in the products' average marginal costs that the merged firms must achieve in order to offset a
  price increase.}

\references{
 Froeb, Luke and Werden, Gregory (1998).
``A robust test for consumer welfare enhancing mergers among sellers
of a homogeneous product.''
\emph{Economics Letters}, \bold{58}(3), pp. 367 - 369.
}


\author{Charles Taragin}


\seealso{\code{\link{cmcr.bertrand}} for a differentiated products Bertrand version of this measure.}

\examples{

    shares=c(.05,.65)
    industryElast = 1.9


    cmcr.cournot(shares,industryElast)


     ## Calculate the necessary percentage cost reductions for various shares and
     ##  industry elasticities in a two-product merger where both firm
     ##  products have identical share (see Froeb and
     ##  Werden, 1998, pg. 369, Table 1)



    deltaHHI = c(100, 500, 1000, 2500, 5000) #start with change in HHI
    shares = sqrt(deltaHHI/(2*100^2)) #recover shares from change in HHI
    industryElast = 1:3

    result = matrix(nrow=length(deltaHHI),ncol=length(industryElast),
                    dimnames=list(deltaHHI,industryElast))

    for(s in 1:length(shares)){
        for(e in 1:length(industryElast)){


           result[s,e] = cmcr.cournot(rep(shares[s],2),industryElast[e])[1]

    }}

    print(round(result,1))
}