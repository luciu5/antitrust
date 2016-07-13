\name{auction2nd.logit}
\alias{auction2nd.logit}
\title{2nd Score Procurement Auction Model}
\description{Calibrates consumer demand using Logit and then simulates the price effect of a merger between two firms
  under the assumption that all firms in the market are playing a
  differentiated products 2nd score auction game.}

\usage{
auction2nd.logit(prices,shares,margins,
                  ownerPre,ownerPost,
                  normIndex=ifelse(sum(shares)<1,NA,1),
                  mcDelta=rep(0,length(prices)),
                  subset=rep(TRUE,length(prices)),
                  mcOutside = 0,
                  labels=paste("Prod",1:length(prices),sep="")
                  )


             }
\arguments{
  \item{}{Let k denote the number of products produced by all firms playing the
    Bertrand pricing game.}
  \item{prices}{A length k vector of product prices.}
  \item{shares}{A length k vector of product (quantity) shares. Values must be
    between 0 and 1.}
  \item{margins}{A length k vector of product margins, some of which may
    equal NA.}
  \item{normIndex}{An integer equalling the index (position) of the
    inside product whose mean valuation will be normalized to 1. Default
    is 1, unless \sQuote{shares} sum to less than 1, in which case the default is
    NA and an outside good is assumed to exist.}
  \item{ownerPre}{EITHER a vector of length k whose values
    indicate which firm produced a product pre-merger OR
    a k x k matrix of pre-merger ownership shares.}
  \item{ownerPost}{EITHER a vector of length k whose values
    indicate which firm produced a product after the merger OR
    a k x k matrix of post-merger ownership shares.}
  \item{mcDelta}{A vector of length k where each element equals the
    proportional change in a product's marginal costs due to
    the merger. Default is 0, which assumes that the merger does not
    affect any products' marginal cost.}
  \item{subset}{A vector of length k where each element equals TRUE if
     the product indexed by that element should be included in the
     post-merger simulation and FALSE if it should be excluded.Default is a
     length k vector of TRUE.}
  \item{mcOutside}{A length 1 vector indicating the marginal cost of the
    outside good. Default is 0.}
  \item{labels}{A k-length vector of labels. Default is "Prod#", where
    \sQuote{#} is a number between 1 and the length of \sQuote{prices}.}
}



\details{Using product prices, quantity shares and all of the
  product margins from at least one firm, \code{auction2nd.logit} is able to
  recover the price coefficient and product mean valuations in a
  Logit demand model. \code{auction2nd.logit} then uses these
  calibrated parameters to simulate a merger between two firms, under the assumption that firms are particpating in a 2nd score procurement auction.

  


}

\value{
   \code{auction2nd.logit} returns an instance of \code{\linkS4class{Auction2ndLogit}}, a child class of \code{\linkS4class{Logit}.}.
  }






\seealso{\code{\link{logit}}}

\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\references{
Miller, Nathan (2014). \dQuote{Modeling the effects of mergers in procurement} \emph{International Journal of Industrial Organization} , \bold{37}, pp. 201-208.

}




\examples{

## Calibration and simulation results from a merger between firms 2 and 3 
## of a 4-firm market
## Source: Miller 2014 backup materials http://www.nathanhmiller.org/research

share = c(0.29,0.40,0.28,0.03)

price = c(35.53,  154, 84.08, 53.16)
cost  =  c(NA, 101, NA, NA)

ownerPre <- ownerPost <- diag(length(share))

#Suppose products 2 and 3 merge
ownerPost[2,3] <- ownerPost[3,2] <- 1 

margin = 1 - cost/price

result <- auction2nd.logit(price,share,margin,
ownerPre=ownerPre,ownerPost=ownerPost,normIndex=2)


print(result)
summary(result,revenue=FALSE)

##re-run without any price information except Firm 2

price[-2] <- NA

result.noprice <- auction2nd.logit(price,share,margin,
ownerPre=ownerPre,ownerPost=ownerPost,normIndex=2)

print(result.noprice)
summary(result.noprice,revenue=FALSE)

## Get a detailed description of the 'Auction2ndLogit' class slots
showClass("Auction2ndLogit")

## Show all methods attached to the 'Auction2ndLogit' Class
showMethods(classes="Auction2ndLogit")

}
