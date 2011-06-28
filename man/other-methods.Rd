\name{other-methods}
\docType{methods}
\alias{other-methods}
\alias{calcMargins}
\alias{calcShares}
\alias{calcQuantities}
\alias{calcPrices}
\alias{calcPriceDelta}
\alias{calcSlopes}
\alias{getNestsParms}

%\alias{calcMargins,Antitrust-method}
%\alias{calcMargins,PCAIDS-method}
%\alias{CalcPriceDelta,Antitrust-method}
%\alias{CalcPriceDelta,PCAIDS-method}
%\alias{CalcPriceDelta,PCAIDSNests-method}
%\alias{calcShares,Linear-method}
%\alias{calcShares,PCAIDS-method}
%\alias{calcShares,Logit-method}
%\alias{calcShares,LogitNests-method}
%\alias{calcShares,CES-method}
%\alias{calcShares,CESNests-method}
%\alias{calcSlopes,Linear-method}
%\alias{calcSlopes,LogLog-method}
%\alias{calcSlopes,PCAIDS-method}
%\alias{calcSlopes,Logit-method}
%\alias{calcSlopes,LogitNests-method}
%\alias{calcSlopes,CES-method}
%\alias{calcSlopes,CESNests-method}
%\alias{calcPrices,Linear-method}
%\alias{calcPrices,LogLog-method}
%\alias{calcPrices,Logit-method}
%\alias{calcPrices,LogitNests-method}
%\alias{calcPrices,CES-method}
%\alias{calcPrices,CESNests-method}
%\alias{calcPrices,CES-method}
%\alias{getNestsParms,PCAIDSNests-method}
%\alias{getNestsParms,LogitNests-method}
%\alias{getNestsParms,CESNests-method}

\title{Other Useful Methods}
\description{
 Below is a description of methods that users may find helpful.
}
\section{Methods}{
  \describe{

  For all methods containing a 'preMerger' argument, setting 'preMerger'
  equal to TRUE does the calculation using pre-merger equilibrium
  prices, while setting 'preMerger' equal to FALSE does the calculation
  using post-merger prices.


  \item{\code{calcShares(object =
      c(("Linear","PCAIDS","Logit","LogitNests","CES","CESNests")),
      preMerger=TRUE)}}{

    Computes equilibrum product shares assuming that firms are playing a
    Nash-Bertrand pricing game.
  }

  \item{\code{calcQuantities(object =
      c(("Linear","LogLog")),
      preMerger=TRUE)}}{

    Computes equilibrum product quantities assuming that firms are playing a
    Nash-Bertrand pricing game.
  }

  \item{\code{calcPrices(object =
      c("Linear","LogLog","Logit","LogitNests","CES","CESNests"),
      preMerger=TRUE)}}{

    Computes equilibrum product prices assuming that firms are playing a
    Nash-Bertrand pricing game.
  }

  \item{\code{calcMargins(object = c("Antitrust","PCAIDS"),
      preMerger=TRUE)}}{

    Computes equilibrum product margins assuming that firms are playing a
    Nash-Bertrand pricing game.
  }


  \item{\code{calcPriceDelta(object = c("Antitrust", "PCAIDS", "PCAIDSNests"))}}{

    Computes equilbrium price changes due to a merger assuming that firms are playing a
    Nash-Bertrand pricing game.
  }


  \item{\code{calcSlopes(object =
      c("Antitrust","Linear","LogLog","PCAIDS","PCAIDSNests","Logit","LogitNests","CES","CESNests"),
      preMerger=TRUE)}}{

    Computes demand parameters assuming that firms are playing a
    Nash-Bertrand pricing game.
  }

  \item{\code{getNestsParms(object = c("PCAIDSNests","LogitNests","CESNests"),
      preMerger=TRUE)}}{

    Returns a vector (or in the case of PCAIDSNests, a matrix) of
    calibrated nesting parameters.
  }


}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
