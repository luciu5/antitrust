\name{cmcr.bertrand}
\alias{cmcr.bertrand}
\title{Compensating Marginal Cost Reductions (Bertrand)}

\description{Calculate the marginal cost reductions necessary to restore
  premerger prices in a merger involving firms playing a differentiated
  product Bertrand pricing game.}

\usage{cmcr.bertrand(price, margin,  diversion, ownerPre, ownerPost=matrix(1,ncol=length(price), nrow=length(price)), labels=paste("Prod",1:length(price),sep=""))}


\arguments{
  Let k_1 denote the number of products produced by merging firm 1 and
  k_2 denote the number of products produced by merging firm 2.
  \item{price}{A vector of product prices with length k_1 + k_2}
  \item{margin}{A vector of product margins with length k_1 + k_2}
  \item{diversion}{A (k_1 + k_2) x (k_1 + k_2) matrix of diversion ratios
    with diagonal elements equal to 1}
  \item{ownerPre}{EITHER a vector of length k_1 + k_2 whose values
  indicate which of the merging parties produced a product pre-merger OR
    a (k_1 + k_2) x (k_1 + k_2) matrix of pre-merger ownership shares}
  \item{ownerPost}{A (k_1 + k_2) x (k_1 + k_2) matrix of post-merger
  ownership shares. Default is a (k_1 + k_2) x (k_1 + k_2) matrix of 1s}
  \item{labels}{Value vector labels}
}

\details{All 'price' elements must be positive, all 'margin'
  elements must be between 0 and 1, and all 'diversion' elements must be
  between 0 and 1 in absolute value. In
  addition, off-diagonal elements (i,j) of 'diversion'
  must equal an estimate of the diversion ratio from product i to product j
  (i.e. the estimated fraction of i's sales that go to j due to a small
  increase in i's price). Also, 'diversion'
  elements are positive if i and j are substitutes and negative if i and j are
  complements.

  'ownerPre' will typically be a vector whose values equal 1 if a product is produced
  by firm 1 and 0 otherwise, though other values including firm name are
  acceptable.  Optionally, 'ownerPre' may be set equal
  to a matrix of the merging firms pre-merger ownership
  shares. These ownership shares must be between 0 and 1.

  'ownerPost'  is an optional argument that should only be specified if
  one party to the acquisition is assuming partial control of the
  other party's assets.   'ownerPost' elements must be between 0 and
  1.
}

\value{A vector whose value equals the percentage change in each products'
  marginal costs that the merged firms must achieve in order to offset a
  price increase.}

\references{Werden1996}


\author{Charles Taragin}


\seealso{upp}

\examples{

     ## Let k_1 = 1 and and k_2 = 2 ##

    p1 = 50;      margin1 = .3
    p2 = c(45,70); margin2 = c(.4,.6)
    isOne=c(1,0,0)
    diversions = matrix(c(1,.5,.01,.6,1,.1,.02,.2,1),ncol=3)

    cmcr.bertrand(c(p1,p2), c(margin1,margin2), diversions, isOne)

    ## Result:      Prod1     Prod2     Prod3
    ##          67.225502 73.097206  8.424692



     ## Calculate the necessary percentage cost reductions for various margins and
     ## diversion ratios in a two-product merger where both products have
     ## equal prices and diversions (see Werden 1996, pg. 412, Table 1)


    margins = seq(.4,.7,.1)
    diversions = seq(.05,.25,.05)
    prices = rep(1,2) #assuming prices are equal, we can set product prices to 1
    isOne = c(1,0)
    result = matrix(ncol=length(margins),nrow=length(diversions),dimnames=list(diversions,margins))

    for(m in 1:length(margins)){
        for(d in 1:length(diversions)){

           dMatrix = diag(2)
           dMatrix[2,1] <- dMatrix[1,2] <- diversions[d]

           firmMargins = rep(margins[m],2)

           result[d,m] = cmcr.bertrand(prices, firmMargins, dMatrix, isOne)[1]

    }}

    print(round(result,1))
}