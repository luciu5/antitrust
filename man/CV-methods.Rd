\name{CV-methods}
\docType{methods}
\alias{CV-methods}
\alias{CV}
\alias{CV,ANY-method}
\alias{CV,AIDS-method}
\alias{CV,CES-method}
\alias{CV,CESNests-method}
\alias{CV,Linear-method}
\alias{CV,LogLin-method}
\alias{CV,Logit-method}
\alias{CV,LogitNests-method}

\title{Methods For Calculating Compensating Variation (CV)}
\description{
 Calculate the amount of money a consumer would need to
 be paid to be just as well off as they were before the merger.
}
\section{Methods}{
\describe{

  \item{\code{signature(object = c(Logit,LogitNests))}}{
    All the information needed to
    compute CV is already available within the Logit and Nested Logit classes.
  }


  \item{\code{signature(object = c(CES, CESNests), revenueInside)}}{
    The CV method for the "CES" and nested "CES" classes has an additional
    parameter, `revenueInside', which must be set equal to
    the total amount that consumers have spent on
    products  inside the market in order for CV to be calculated.
  }


  \item{\code{signature(object = AIDS , totalRevenue)}}{
    The CV method for "AIDS" has an additional parameter, `totalRevenue',
    which should aggregate income (e.g. GDP). If supplied computes CV in terms of dollars. If missing, CV is
    calculated as a percentage change in aggregate in income. must be set equal to the vector of pre-merger prices for all
    products in the market in order for CV to be calculated.
  }

  \item{\code{signature(object = c(Linear,LogLin))}}{
    Although no additional information is needed to calculate CV for
    either the "Linear" or "LogLin" classes, The CV method will fail if
    the appropriate restrictions on the demand parameters have not been imposed.
  }
}}
\keyword{methods}

