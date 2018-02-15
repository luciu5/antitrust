\name{antitrust_shiny}
\alias{antitrust_shiny}

\title{A Shiny Interface to the Antitrust Package}

\description{A Shiny Interface to the Antitrust Package}

\usage{
antitrust_shiny()
      
             }


\details{
  \code{antitrust_shiny} launches a shiny interface for the antitrust package. The shiny interface provides users with the ability to calibrate and simulate horizontal mergers using  many of the supply and demand models included in the antitrust package. 

Unfortunately, the interface does not report convergence issues that may occur during calibration or simulation. As such, users are encouraged to use the shiny interface as an exploratory tool, but to check all results in R.


}


\author{Charles Taragin \email{charles.taragin@usdoj.gov}}



\examples{
## launch shiny interface to antitrust package
\dontrun{
antitrust_shiny()
}
}
