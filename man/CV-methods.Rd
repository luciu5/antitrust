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
\alias{CV,Auction2ndLogit-method}
\alias{CV,Cournot-method}

\title{Methods For Calculating Compensating Variation (CV)}
\description{
 Calculate the amount of money a consumer would need to
 be paid to be just as well off as they were before the merger.
}
\section{Methods}{
\describe{

  \item{\code{signature(object = c(Logit,LogitNests, CES, CESNests))}}{
    All the information needed to
    compute CV is already available within the Logit, Nested Logit CES and Nested CES classes. 
    In CES and Nested CES, CV cannot be calculated if the outside share cannot be inferred.
  }


  \item{\code{signature(object = AIDS)}}{
    If the  \sQuote{insideSize}  slot to the \dQuote{AIDS} class equals NA, CV is calculated as a percentage of 
    total expenditure (revenues) on products included in the simulation. Otherwise CV is calculated in terms of dollars. 
    Pre-merger prices for all products in the market must be supplied in order for CV to be calculated.
  }

  \item{\code{signature(object = c(Linear,LogLin))}}{
    Although no additional information is needed to calculate CV for
    either the \dQuote{Linear} or \dQuote{LogLin} classes, The CV method will fail if
    the appropriate restrictions on the demand parameters have not been imposed.
  }
}}
\keyword{methods}

