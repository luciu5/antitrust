\name{collusion-methods}
\docType{methods}
\alias{calcProducerSurplusGrimTrigger,Bertrand-method}
\alias{calcProducerSurplusGrimTrigger}

\title{Methods For Evaluating Collusion}

\description{
This page describes methods that may be used to explore how a merger
affects firms' incentives to collude.
}

\usage{
\S4method{calcProducerSurplusGrimTrigger}{Bertrand}(object,coalition,
          discount,preMerger=TRUE,isCollusion=FALSE,...)
}

\arguments{
  \item{}{Let k denote the number of products in the market, and let c denote the number of firms in a coalition}
  \item{object}{An instance of one of the classes listed above.}
  \item{coalition}{A length c vector of integers indicating the index of
    the products participating in the coalition.}
  \item{discount}{A length k vector of values between 0 and 1 that
    represent the product-specific discount rate for all products
  produced by firms particiapting in the coalition. NAs are allowed}
  \item{preMerger}{TRUE returns pre-merger result, FALSE returns
    post-merger results. Default is TRUE.}
  \item{isCollusion}{TRUE recalculates demand and cost parameters under
    the assumption that the coalition specified in \sQuote{coalition} is
    operating pre-merger. FALSE (the default) uses demand
    and cost parameters calculated from the \sQuote{ownerPre} matrix.}
  \item{...}{Additional argument to pass to \code{calcPrices}}
}

\details{
  \code{calcProducerSurplusGrimTrigger} calculates \sQuote{preMerger} product
  producer surplus (as well as other statistics -- see below), under the
  assumption that firms are playing
  an N-player iterated Prisoner's Dilemma where in each period a
  coalition of firms decides whether to \emph{cooperate} with
  one another by setting the joint surplus maximizing price on some
  \sQuote{coalition} of their products, or \emph{defect} from the coalition by setting all
  of their products'
  prices to optimally undercut the prices of the coalition's
  products. Moreover, firms are assumed
  to play Grim Trigger strategies where each firm cooperates in the
  current period so long as \emph{every} firm in the coalition cooperated last period and
  defects otherwise. product level \sQuote{discount} rates are then employed to determine
  whether a firm's discounted surplus from remaining in the coalition are greater than
  its surplus from optimally undercutting the coalition prices' for one
  period plus its discounted surplus when all firms set Nash-Bertrand prices in
  all subsequent periods.

}

\value{
  \code{calcProducerSurplusGrimTrigger} returns a data frame with rows
  equal to the number of products produced by any firm participating in
  the coalition and the following 5 columns

  \itemize{
    \item{Discount:}{The user-supplied discount rate}
    \item{Coord:}{Single period producer surplus from coordinating}
    \item{Defect:}{Single period producer surplus from defecting}
    \item{Punish:}{Single period producer surplus from punishing using Bertrand price}
    \item{IC:}{TRUE if the discounted producer surplus from coordinating
    across all firm products are
      greater than the surplus from defecting across all firm products for one period and receiving
      discounted Bertrand surplus for all subsequent periods under Grim Trigger.}
  }

}



\keyword{methods}

