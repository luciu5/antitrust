\name{HHI}
\alias{HHI}
\title{Herfindahl-Hirschman Index}

\description{Calculate the Herfindahl-Hirschman Index with arbitrary
  ownership and control.}

\usage{
HHI(shares,
    owner=diag(length(shares)),
    control)
}


\arguments{
  Let k denote the number of products produced by the merging parties.
  \item{shares}{A length-k vector of product quantity shares.}
  \item{owner}{EITHER a vector  of length k whose values
  indicate which of the merging parties produced a product OR
    a k x k matrix of ownership shares. Default is a diagonal matrix,
  which assumes that each product is owned by a separate firm.}
\item{control}{EITHER a vector  of length k whose values
  indicate which of the merging parties have the ability to make pricing
  or output decisions OR
  a k x k matrix of
  control shares. Default is a k x k matrix equal to 1 if `owner' > 0
  and 0 otherwise.}
}

\details{All `shares' must be between 0 and 1. When `owner' is a matrix,
  the i,jth element of `owner' should equal the percentage of product j's profits earned by the owner
  of product i. When `owner' is a vector, \code{HHI} generates a k x k
  matrix of whose i,jth element equals 1 if products i and j are
  commonly owned and 0 otherwise. `control' works in a fashion similar
  to `owner'.
}

\value{
 \code{HHI} returns a number between 0 and 10,000
}

\references{
  Salop, Steven and O'Brien, Daniel (2000)
  ``Competitive Effects of Partial Ownership: Financial Interest and Corporate Control''
  67 Antitrust L.J. 559, pp. 559-614.
}


\author{Charles Taragin \email{charles.taragin@usdoj.gov}}


\seealso{\code{\link{other-methods}} for computing HHI following merger
  simulation.}

\examples{

## Consider a market with 5 products labeled 1-5. 1,2 are produced
## by Firm A, 2,3 are produced by Firm B, 3 is produced by Firm C.
## The pre-merger product market shares are

shares = c(.15,.2,.25,.35,.05)
owner  = c("A","A","B","B","C")
nprod  = length(shares)

HHI(shares,owner)

## Suppose that Firm A acquires a 75\% ownership stake in product 3, and
## Firm B get a 10\% ownership stake in product 1. Assume that neither
## firm cedes control of the product to the other.

owner <- diag(nprod)

owner[1,2] <- owner[2,1] <- owner[3,4] <- owner[4,3] <- 1
control <- owner
owner[1,1] <- owner[2,1] <- .9
owner[3,1] <- owner[4,1] <- .1
owner[1,3] <- owner[2,3] <- .75
owner[3,3] <- owner[4,3] <- .25

HHI(shares,owner,control)

## Suppose now that in addition to the ownership stakes described
## earlier, B receives 30\% of the control of product 1
control[1,1] <- control[2,1] <- .7
control[3,1] <- control[4,1] <- .3

HHI(shares,owner,control)


}