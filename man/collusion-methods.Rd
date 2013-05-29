\name{collusion-methods}
\docType{methods}
\alias{calcProducerSurplusGrimTrigger,Bertrand-method}
\alias{calcProducerSurplusGrimTrigger}

\title{Methods For Evaluating Collusion}

\description{
This page describes methods that may be used to determine how a merger
affects firms' incentives to collude.
}

\usage{
\S4method{calcProducerSurplusGrimTrigger}{Bertrand}(object,coalition,discount,preMerger=TRUE,...)
}

\arguments{
  \item{}{Let k denote the number of products in the market, and let c denote the number of firms in a coalition}
  \item{object}{An instance of one of the classes listed above.}
  \item{coalition}{A length c vector of integers indicating the index of
    the products participating in the coalition.}
  \item{discount}{A length c vector of values between 0 and 1 that
    represent the product specific discount rate.}
  \item{preMerger}{TRUE return pre-merger result, FALSE returns
    post-merger results. Default is TRUE.}
   \item{...}{Additional argument to pass to \code{calcPrices}}
}

\details{
  \code{calcProducerSurplusGrimTrigger} assumes that firms are playing
  an N-player iterated Prisoner's Dilemma where in each period  a
  user-specified 'coalition' of firms decides whether to \emph{cooperate} with
  one another by setting the joint profit maximizing price on some of
  their products, or \emph{defect} from the coalition by setting all
  of their products'
  prices to optimally undercut the prices of the coalition's
  products. Moreover, firms are assumed
  to play Grim Trigger strategies where each firm cooperates in the
  current period so long as \emph{every} firm cooperated last period and
  defects otherwise. 'discount' rates are then employed to determine
  whether a firm's discounted profits from remaining in the coalition are greater than
  its profits from optimally undercutting the coalition prices' for one
  period plus its discounted profits when all firms set Nash-Bertrand prices in
  all subsequent periods.

}

\value{
  \code{calcProducerSurplusGrimTrigger} returns a data frame with c rows
  and the following 5 columns

  \itemize{
    \item{Discount}{The user-supplied discount rate}
    \item{Coord}{Single period profits from coordinating}
    \item{Defect}{Single period profits from defecting}
    \item{Punish}{Single period profits from punishing using Bertrand price}
    \item{IC}{TRUE if the discounted profits from coordinating
    across all firm products are
      greater than the profits from defecting across all firm products for one period and receiving
      discounted Bertrand profits for all subsequent periods.}
  }

}



\keyword{methods}

