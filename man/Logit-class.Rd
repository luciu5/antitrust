\name{Logit-class}
\Rdversion{1.1}
\docType{class}
\alias{Logit-class}
%\alias{calcPrices,Logit-method}
%\alias{calcShares,Logit-method}
%\alias{calcSlopes,Logit-method}
%\alias{CV,Logit-method}
%\alias{elast,Logit-method}

\title{Class \dQuote{Logit}}
\description{
   The \dQuote{Logit} class contains all the information needed to calibrate a Logit
  demand system and perform a merger simulation analysis under the assumption that
  firms are playing a differentiated products Bertrand pricing game.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by using the constructor function \code{\link{logit}}.
%%  ~~ describe objects here ~~
}



      
\section{Slots}{
  Let k denote the number of products produced by all firms.
  \describe{
    \item{\code{prices}:}{A length k vector of product prices.}
    \item{\code{margins}:}{A length k vector of product margins, some of which may
      equal NA.}
    \item{\code{normIndex}:}{An integer specifying the product index against which
     the mean values of all other products are normalized.}
    \item{\code{shareInside}:}{The share of customers that purchase any of the
      products included in the `prices' vector.}
    \item{\code{priceOutside}:}{The price of the outside good. Default
      is 0.}
     \item{\code{slopes}:}{A list containing the coefficient on price (\sQuote{alpha}) and the vector of mean
      valuations (\sQuote{meanval})}
       \item{\code{mktElast}:}{A length 1 vector of market elasticities. }
 \item{\code{priceStart}:}{A length-k vector of starting prices for the non-linear solver}
  \item{\code{insideSize}:}{A positive number equal to total pre-merger quantities (revenues for CES) 
                            for all products
                            included in the simulation.}
   \item{\code{mktSize}:}{A positive number equal to total quantities (revenues for CES) 
   pre-merger for all products in the simulations
                          as well as the outside good.}                          
                            
  
  }
}
\section{Extends}{
Class \code{\linkS4class{Bertrand}}, directly.
Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 2.
}
\section{Methods}{
  For all of methods containing the \sQuote{preMerger} argument, \sQuote{preMerger} takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger ownership structure, while FALSE invokes the method using
  the post-merger ownership structure.
  \describe{
    \item{\code{calcPrices}}{\code{signature(object = Logit, preMerger
	= TRUE,isMax=FALSE,...)}}{ Compute either pre-merger or post-merger equilibrium
      prices under the assumptions that consumer demand is Logit and firms play a differentiated product
      Bertrand Nash pricing game. When isMax equals TRUE, a check is run
      to determine if the calculated equilibrium price vector locally maximizes
      profits. \sQuote{...} may be used to change the
      default values of \code{\link[BB]{BBsolve}}, the non-linear equation solver. }
    \item{\code{\link{calcPriceDeltaHypoMon}}}{\code{signature(object = Logit,prodIndex,...)}}{ Calculates the price changes that a
   Hypothetical Monopolist would impose on its products
   relative to pre-merger prices. }
 \item{\code{calcShares}}{\code{signature(object = Logit, preMerger
	= TRUE,revenue = FALSE)}}{ Compute either pre-merger or post-merger equilibrium
      shares under the assumptions that consumer demand is Logit and firms play a differentiated product
      Bertrand Nash pricing game. \sQuote{revenue} takes
      on a value of TRUE or FALSE, where TRUE calculates revenue shares,
      while FALSE calculates quantity shares.}
    \item{\code{calcSlopes}}{\code{signature(object = Logit)}}{ Uncover Logit demand
      parameters. Assumes that firms are currently at equilibrium  in a
      differentiated product Bertrand Nash pricing game.}
    \item{\code{\link{CV}}}{\code{signature(object =
	Logit)}}{ Calculate the amount of money a representative consumer would need to
      be paid to be just as well off as they were before the merger.}
    \item{\code{\link{elast}}}{\code{signature(object = Logit, preMerger
	= TRUE)}}{ Computes a k x k matrix of own and cross-price elasticities.}
  }

}


\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\examples{
showClass("Logit")           # get a detailed description of the class
showMethods(classes="Logit") # show all methods defined for the class
}
\keyword{classes}
