\name{other-methods}
\docType{methods}
\alias{other-methods}
\alias{calcMC}
\alias{calcMargins}
\alias{calcPrices}
\alias{calcPriceDelta}
\alias{calcQuantities}
\alias{calcShares}
\alias{calcSlopes}
\alias{getNestsParms}
\alias{ownerToMatrix}
\alias{ownerToVec}
\alias{calcMC-methods}
\alias{calcMargins-methods}
\alias{calcPrices-methods}
\alias{calcPriceDelta-methods}
\alias{calcQuantities-methods}
\alias{calcShares-methods}
\alias{calcSlopes-methods}
\alias{getNestsParms-methods}
\alias{ownerToMatrix-methods}
\alias{ownerToVec-methods}
\alias{getNestsParms,PCAIDSNests-method}
\alias{hhi,Bertrand-method}
\alias{summary,Bertrand-method}
\alias{summary,AIDS-method}
\alias{ownerToMatrix,Antitrust-method}
\alias{ownerToVec,Antitrust-method}
\alias{calcMC,Bertrand-method}
\alias{calcMargins,Bertrand-method}
\alias{calcMargins,AIDS-method}
\alias{calcMargins,LogitCap-method}
\alias{CalcPriceDelta,Antitrust-method}
\alias{CalcPriceDelta,AIDS-method}
\alias{calcQuantities,Linear-method}
\alias{calcQuantities,LogLin-method}
\alias{calcQuantities,LogitALM-method}
\alias{calcShares,Linear-method}
\alias{calcShares,AIDS-method}
\alias{calcShares,Logit-method}
\alias{calcShares,LogitNests-method}
\alias{calcShares,CES-method}
\alias{calcShares,CESNests-method}
\alias{calcSlopes,Linear-method}
\alias{calcSlopes,LogLin-method}
\alias{calcSlopes,AIDS-method}
\alias{calcSlopes,PCAIDS-method}
\alias{calcSlopes,PCAIDSNests-method}
\alias{calcSlopes,Logit-method}
\alias{calcSlopes,LogitALM-method}
\alias{calcSlopes,LogitNests-method}
\alias{calcSlopes,LogitCap-method}
\alias{calcSlopes,CES-method}
\alias{calcSlopes,CESNests-method}
\alias{calcPrices,Linear-method}
\alias{calcPrices,LogLin-method}
\alias{calcPrices,Logit-method}
\alias{calcPrices,LogitALM-method}
\alias{calcPrices,LogitNests-method}
\alias{calcPrices,LogitCap-method}
\alias{calcPrices,CES-method}
\alias{calcPrices,CESNests-method}
\alias{calcPrices,AIDS-method}

\title{Other Useful Methods}
\description{
 Below is a description of methods that users may find helpful.
}
\usage{
\S4method{calcShares}{Linear,AIDS,Logit,LogitNests,CES,CESNests}(object,preMerger=TRUE,
revenue=FALSE)

\S4method{calcQuantities}{Linear,LogLin}(object,preMerger=TRUE)

\S4method{calcPrices}{Linear,LogLin,AIDS,Logit,LogitNests,LogitCap,CES,CESNests}(object,preMerger=TRUE,
...)

\S4method{calcPriceDelta}{Bertrand,AIDS}(object,preMerger=TRUE,
...)

\S4method{calcMargins}{Bertrand,LogitCap}(object,preMerger=TRUE)

\S4method{calcMC}{Bertrand}(object,preMerger=TRUE)

\S4method{calcSlopes}{Linear,LogLin,AIDS,PCAIDSNests,Logit,LogitNests,LogitCap,CES,CESNests}(object,preMerger=TRUE)

\S4method{getNestsParms}{PCAIDSNests}(object)

\S4method{hhi}{Bertrand}(object,preMerger=TRUE,revenue=FALSE)

\S4method{ownerToMatrix}{Antitrust}(object,preMerger=TRUE)
\S4method{ownerToVec}{Antitrust}(object,preMerger=TRUE)
}

\arguments{
  \item{object}{An instance of one of the classes listed above.}
  \item{preMerger}{If TRUE, returns pre-merger outcome. If
    FALSE, returns post-merger outcome.  Default is TRUE.}
  \item{revenue}{If TRUE, returns revenue shares. If FALSE,
    returns quantity shares.}
%  \item{mcDelta}{A vector of length k where each element equals the
%    proportional change in a product's marginal costs due to
%    the merger. Default is 0, which assumes that the merger does not
%    affect any products' marginal cost.}
  \item{...}{Arguments to be passed to non-linear solver.}
}
\section{Methods}{
  \describe{

  \item{calcShares}{\code{signature(object=
	c(Linear,AIDS,Logit,LogitNests,CES,CESNests),
	preMerger=TRUE, revenue=FALSE)}}{

    Computes equilibrium product shares assuming that firms are playing a
    Nash-Bertrand pricing game. `revenue' takes
    on a value of TRUE or FALSE, where TRUE calculates revenue shares,
    while FALSE calculates quantity shares.
  }
  \item{calcQuantities}{\code{signature(object=
	c(Linear,LogLin),
	preMerger=TRUE)}}{

    Computes equilibrium product quantities assuming that firms are playing a
    Nash-Bertrand pricing game.
  }
  \item{calcPrices}{\code{signature(object=c(Linear,LogLin,AIDS,Logit,LogitNests,LogitCap,CES,CESNests),preMerger=TRUE, ...)}}{
    Computes equilibrium product price levels assuming that firms are playing a
    Nash-Bertrand pricing game. `...' may be used to feed additional
    options to the optimizer responsible for computing equilibrium
    prices. Typically, \code{\link[nleqslv]{nleqslv}} is used, but see
    the appropriate document for further details.
  }
   \item{calcPriceDelta}{\code{signature(object=c(Bertrand,AIDS),preMerger=TRUE)}}{
    Computes equilibrium price changes due to a merger assuming that firms are playing a
    Nash-Bertrand pricing game.
  }
  \item{calcMargins}{\code{signature(object=c(Bertrand,LogitCap),preMerger=TRUE)}}{
    Computes equilibrium product margins assuming that firms are playing a
    Nash-Bertrand pricing game. For "LogitCap", assumes firms are
    playing a Nash-Bertrand pricing game with capacity constraints.
  }
   \item{calcMC}{\code{signature(object=
	Bertrand,
	preMerger=TRUE)}}{
    Computes either pre- or post-merger marginal costs. Marginal costs
    are assumed to be constant. Post-merger marginal costs are equal to
    pre-merger marginal costs multiplied by 1+`mcDelta', a length-k
    vector of marginal cost changes. `mcDelta' will typically be between 0
    and 1.
  }
   \item{calcSlopes}{\code{signature(object=c(Linear,LogLin,AIDS,PCAIDSNests,Logit,LogitNests,LogitCap,CES,CESNests),preMerger=TRUE)}}{

    Computes demand parameters assuming that firms are playing a
    Nash-Bertrand pricing game.
  }
  \item{getNestsParms}{\code{signature(object=PCAIDSNests)}}{

    Returns a vector a matrix  of
    calibrated nesting parameters.
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
  \item{show}{\code{signature(object)}}{
    Displays the percentage change in prices due to the merger.
  }
  \item{summary}{\code{signature(object,revenue=TRUE,shares=TRUE,parameters=FALSE,digits=2)}}{
    Summarizes the effect of the merger, including price and revenue
    changes. Setting `revenue' equal to FALSE reports quantities rather
    than revenues. Setting `shares' to FALSE reports quantities
    rather than than shares (when possible). Setting `parameters' equal to TRUE reports
    all demand parameters. `digits' controls the number of significant
    digits reported in output.
  }
  \item{upp}{\code{signature(object)}}{
    Calculate the Upwards Pricing Pressure (upp) index.
  }
}
}

\keyword{methods}

