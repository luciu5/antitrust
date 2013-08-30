\name{other-methods}
\docType{methods}
\alias{other-methods}
\alias{calcShares}
\alias{calcShares,ANY-method}
\alias{calcShares,AIDS-method}
\alias{calcShares,CES-method}
\alias{calcShares,CESNests-method}
\alias{calcShares,Linear-method}
\alias{calcShares,Logit-method}
\alias{calcShares,LogitNests-method}
\alias{calcQuantities}
\alias{calcQuantities,ANY-method}
\alias{calcQuantities,Linear-method}
\alias{calcQuantities,LogLin-method}
\alias{calcQuantities,LogitCap-method}
\alias{calcPrices}
\alias{calcPrices,ANY-method}
\alias{calcPrices,AIDS-method}
\alias{calcPrices,Linear-method}
\alias{calcPrices,LogLin-method}
\alias{calcPrices,Logit-method}
\alias{calcPrices,LogitCap-method}
\alias{calcPriceDelta}
\alias{calcPriceDelta,ANY-method}
\alias{calcPriceDelta,Antitrust-method}
\alias{calcPriceDelta,AIDS-method}
\alias{calcProducerSurplus}
\alias{calcProducerSurplus,ANY-method}
\alias{calcProducerSurplus,Bertrand-method}
\alias{calcMargins}
\alias{calcMargins,ANY-method}
\alias{calcMargins,AIDS-method}
\alias{calcMargins,Bertrand-method}
\alias{calcMargins,LogitCap-method}
\alias{calcMC}
\alias{calcMC,ANY-method}
\alias{calcMC,Bertrand-method}
\alias{calcSlopes}
\alias{calcSlopes,ANY-method}
\alias{calcSlopes,AIDS-method}
\alias{calcSlopes,CES-method}
\alias{calcSlopes,CESNests-method}
\alias{calcSlopes,Linear-method}
\alias{calcSlopes,LogLin-method}
\alias{calcSlopes,Logit-method}
\alias{calcSlopes,LogitALM-method}
\alias{calcSlopes,LogitCap-method}
\alias{calcSlopes,LogitNests-method}
\alias{calcSlopes,LogitNestsALM-method}
\alias{calcSlopes,PCAIDS-method}
\alias{calcSlopes,PCAIDSNests-method}
\alias{hhi}
\alias{hhi,ANY-method}
\alias{hhi,Bertrand-method}
\alias{ownerToVec}
\alias{ownerToVec,ANY-method}
\alias{ownerToVec,Antitrust-method}
\alias{ownerToMatrix}
\alias{ownerToMatrix,ANY-method}
\alias{ownerToMatrix,Antitrust-method}
\alias{plot,Bertrand-method}
\alias{summary}
\alias{summary,ANY-method}
\alias{summary,AIDS-method}
\alias{summary,Bertrand-method}
\title{Other Useful \sQuote{Bertrand} Methods}
\description{
 Methods defined for the \sQuote{Bertrand} class and its child classes.
}
\usage{
\S4method{calcShares}{ANY}(object,preMerger=TRUE,revenue=FALSE)
\S4method{calcQuantities}{ANY}(object,preMerger=TRUE)
\S4method{calcPrices}{ANY}(object,preMerger=TRUE,subset,...)
\S4method{calcPriceDelta}{Antitrust}(object)
\S4method{calcPriceDelta}{AIDS}(object,isMax=FALSE,subset,...)
\S4method{calcProducerSurplus}{Bertrand}(object,preMerger=TRUE)
\S4method{calcMargins}{ANY}(object,preMerger=TRUE)
\S4method{calcMC}{Bertrand}(object,preMerger=TRUE)
\S4method{calcSlopes}{ANY}(object,preMerger=TRUE)
\S4method{hhi}{Bertrand}(object,preMerger=TRUE,revenue=FALSE)
\S4method{ownerToMatrix}{Antitrust}(object,preMerger=TRUE)
\S4method{ownerToVec}{Antitrust}(object,preMerger=TRUE)
\S4method{plot}{Bertrand}(x,scale=.1)
\S4method{summary}{Bertrand}(object,revenue=TRUE,shares=TRUE,parameters=FALSE,digits=2,...)


}

\arguments{
\item{}{Let k denote the number of products}
  \item{object}{An instance of one of the classes listed above.}
  \item{x}{Used only in plot method. Should always be set equal to object}
  \item{preMerger}{If TRUE, returns pre-merger outcome. If
    FALSE, returns post-merger outcome.  Default is TRUE.}
  \item{isMax}{If TRUE, uses numerical derivatives to determine if
  equilibrium price vector is a local maximum. Default is FALSE.}
  \item{revenue}{If TRUE, returns revenues. If FALSE,
    returns quantities. Default is TRUE}
   \item{subset}{A vector of length k where each element equals TRUE if
     the product indexed by that element should be included in the
     post-merger simulation and FALSE if it should be excluded.Default is a
     length k vector of TRUE.}
  \item{shares}{If TRUE, returns shares. If FALSE,
    returns levels. Default is TRUE}
  \item{parameters}{If TRUE, reports demand and cost parameters. Default
    is FALSE}
    \item{digits}{The number of significant digits to round
  printed results. Default is 2}
  \item{scale}{The proportion below marginal cost and above equilbrium price that should be plotted. Default is .1 }
  \item{...}{Arguments to be passed to non-linear solver, OR for
    \code{summary} to \code{CV}.}
}
\section{Methods}{
  \describe{

  \item{calcShares}{\code{signature(object=
	c(Linear,AIDS,Logit,LogitNests,CES,CESNests),
	preMerger=TRUE, revenue=FALSE)}}{

    Computes equilibrium product shares assuming that firms are playing a
    Nash-Bertrand pricing game. \sQuote{revenue} takes
    on a value of TRUE or FALSE, where TRUE calculates revenue shares,
    while FALSE calculates quantity shares.
  }
  \item{calcQuantities}{\code{signature(object=c(Linear,LogLin,LogitCap),preMerger=TRUE)}}{

    Computes equilibrium product quantities assuming that firms are playing a
    Nash-Bertrand pricing game.
  }
  \item{calcPrices}{\code{signature(object=c(Linear,LogLin,AIDS,Logit,LogitNests,LogitCap,CES,CESNests),preMerger=TRUE, ...)}}{
    Computes equilibrium product price levels assuming that firms are playing a
    Nash-Bertrand pricing game. \sQuote{...} may be used to feed additional
    options to the optimizer responsible for computing equilibrium
    prices. Typically, \code{\link[BB]{BBsolve}} is used, but see
    the appropriate document for further details.
  }
   \item{calcPriceDelta}{\code{signature(object=Antitrust)}}{
    Computes equilibrium price changes due to a merger assuming that firms are playing a
    Nash-Bertrand pricing game. This is a wrapper method for computing
    the difference between pre- and post-merger equilbrium prices
  }
   \item{calcPriceDelta}{\code{signature(object=AIDS,isMax=FALSE,...)}}{
    Computes equilibrium price changes due to a merger assuming that firms are playing a
    Nash-Bertrand pricing game. This method calls a non-linear
    equations solver to find a sequence of price changes that satisfy
    the Bertrand FOCs.
  }
   \item{calcProducerSurplus}{\code{signature(object=Bertrand,preMerger=TRUE)}}{
    Computes equilibrium producer surplus.
  }
  \item{calcMargins}{\code{signature(object=c(Bertrand,LogitCap),preMerger=TRUE)}}{
    Computes equilibrium product margins assuming that firms are playing a
    Nash-Bertrand pricing game. For "LogitCap", assumes firms are
    playing a Nash-Bertrand pricing game with capacity constraints.
  }
   \item{calcMC}{\code{signature(object=Bertrand,preMerger=TRUE)}}{
    Computes either pre- or post-merger marginal costs. Marginal costs
    are assumed to be constant. Post-merger marginal costs are equal to
    pre-merger marginal costs multiplied by 1+\sQuote{mcDelta}, a length-k
    vector of marginal cost changes. \sQuote{mcDelta} will typically be between 0
    and 1.
  }
   \item{calcSlopes}{\code{signature(object=c(Linear,LogLin,AIDS,PCAIDSNests,Logit,LogitNests,LogitCap,CES,CESNests),preMerger=TRUE)}}{

    Computes demand parameters assuming that firms are playing a
    Nash-Bertrand pricing game.
  }
   \item{hhi}{\code{signature(object=Bertrand,preMerger=TRUE,revenue=FALSE)}}{

    Computes the  Herfindahl-Hirschman Index (HHI) using simulated market
    shares and either pre- or post-merger ownership information.
  }
  \item{ownerToMatrix}{\code{signature(object=Antitrust,preMerger=TRUE)}}{

    converts a length-k ownership vector into a k x k ownership matrix
    where element i,j equals 1 if products i and j are commonly owned, and
    0 otherwise.
  }
   \item{ownerToVec}{\code{signature(object=Antitrust,preMerger=TRUE)}}{

    converts a k x k ownership matrix into a length-k ownership vector

  }
   \item{plot}{\code{signature(x,scale=.1)}}{
Use \code{\link{ggplot}} to plot pre- and post-merger demand, marginal cost and equilibria. \sQuote{scale} controls the amount above marginal cost and below equilbrium price that is plotted.
  }
  \item{show}{\code{signature(object=Antitrust)}}{
    Displays the percentage change in prices due to the merger.
  }
  \item{summary}{\code{signature(object=c(Bertrand,AIDS),revenue=TRUE,shares=TRUE,parameters=FALSE,digits=2,...)}}{
    Summarizes the effect of the merger, including price and revenue
    changes. Setting \sQuote{revenue} equal to FALSE reports quantities rather
    than revenues. Setting \sQuote{shares} to FALSE reports quantities
    rather than than shares (when possible). Setting \sQuote{parameters} equal to TRUE reports
    all demand parameters. \sQuote{digits} controls the number of significant
    digits reported in output. \sQuote{...} allows other arguments to be passed
    to a \code{CV} method.
  }
  \item{upp}{\code{signature(object)}}{
    Calculate the Upwards Pricing Pressure (upp) index.
  }
}
}

\keyword{methods}

