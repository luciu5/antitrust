\name{CV-methods}
\docType{methods}
\alias{CV-methods}
\alias{CV}
\alias{CV,CES-method}
\alias{CV,CESNests-method}
\alias{CV,Linear-method}
\alias{CV,LogLog-method}
\alias{CV,Logit-method}
\alias{CV,LogitNests-method}
\alias{CV,PCAIDS-method}
\title{Methods For Calculating Compensating Variation (CV)}
\description{
 Calculate the amount of money a consumer would need to
 be paid to be just as well off as they were before the merger.
}
\section{Methods}{
\describe{

  \item{\code{signature(object = c("Logit","LogitNests"))}}{
    All the information needed to
    compute CV is already available within the Logit and Nested Logit classes.
  }


  \item{\code{signature(object = c("CES", "CESNests"), revenueInside)}}{
    The CV method for the CES and nested CES classes has an additional
    parameter, 'revenueInside', which must be set equal to
    the total amount that consumers have spent on
    products  inside the market in order for CV to be calculated.
  }


  \item{\code{signature(object = "PCAIDS" , prices)}}{
    The CV method for "PCAIDS" has an additional parameter, 'prices',
    which must be set equal to the vector of pre-merger prices for all
    products in the market in order for CV to be calculated.
  }

  \item{\code{signature(object = c("Linear","LogLog"))}}{
    Although no additional information is needed to calculate CV for
    either the "Linear" or "LogLog" classes, The CV method will fail if
    the appropriate restrictions on the demand parameters have not been imposed.
  }
}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
