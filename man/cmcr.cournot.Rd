\name{cmcr.cournot}
\alias{cmcr.cournot}
\title{Compensating Marginal Cost Reductions (Cournot)}

\description{Calculate the average marginal cost reduction necessary to restore
  premerger prices in a two-product merger involving firms playing a
  homogeneous product Cournot pricing game.}

\usage{cmcr.cournot(share,industryElasticity)}


\arguments{

  \item{share}{A vector of length 1 containing merging party quantity shares}
  \item{IndustryElasticity}{A vector of length 1 containing industry
  elasticity.}

}

\details{The 'share' vector must have 2 elements, and all 'share'
  elements must be between 0 and 1. The 'industryElasticity' vector must
  have 1 non-negative  element.}
\value{A vector with 1 element whose value equals the percentage change
  in the products' average marginal costs that the merged firms must achieve in order to offset a
  price increase.}

\references{Froeb1998}


\author{Charles Taragin}


\seealso{cmcr.bertrand for a differentiated products Bertrand version of this measure}

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

    result = matrix(nrow=length(deltaHHI),ncol=length(industryElast),dimnames=list(deltaHHI,industryElast))

    for(s in 1:length(shares)){
        for(e in 1:length(industryElast)){


           result[s,e] = cmcr.cournot(rep(shares[s],2),industryElast[e])[1]

    }}

    print(round(result,1))
}