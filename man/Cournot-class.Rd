\name{Cournot-class}
\Rdversion{1.0}
\docType{class}
\alias{Cournot-class}


\title{Class \dQuote{Cournot}}
\description{
  The \dQuote{Cournot} class is a building block used to create other classes
  in this package. As such, it is most likely to be useful for developers
  who wish to code their own merger calibration/simulation routines.


}
\section{Objects from the Class}{

Objects can be created by calls of the form \code{new("Cournot", ...)}.

}

\section{Slots}{
  Let k denote the number of products and n denote the number of plants.
  \describe{
    \item{\code{intercepts}:}{A length k vector containing the calibrated demand intercept.}
    \item{\code{mcfunPre}:}{A length n list whose elements equal a function that calculates a plant's pre-merger marginal cost. }
    \item{\code{mcfunPost}:}{A length n list whose elements equal a function that calculates a plant's post-merger marginal cost. }
     \item{\code{vcfunPre}:}{A length n list whose elements equal a function that calculates a plant's pre-merger variable cost. }
    \item{\code{vcfunPost}:}{A length n list whose elements equal a function that calculates a plant's post-merger variable cost. }
    \item{\code{prices}:}{A length k vector of product prices. }
    \item{\code{quantities}:}{An n x k matrix of plant quantities produced for each product.}
    \item{\code{margins}:}{An n x k matrix of plant product margins.}
    \item{\code{quantityPre}:}{An n x k matrix of predicted pre-merger quantities.}
    \item{\code{quantityPost}:}{An n x k matrix of predicted post-merger quantities.}
    \item{\code{quantityStart}:}{A length n x k vector of starting quantities for the non-linear solver.}
    \item{\code{productsPre}:}{An n x k logical matrix qhose elements are TRUE if a plant produces a product pre-merger and FALSE otherwise.}
    \item{\code{productsPost}:}{An n x k logical matrix qhose elements are TRUE if a plant produces a product post-merger and FALSE otherwise.}
    item{\code{capacitiesPre}:}{A length-n logical vector whose elements equal to pre-merger plant capacities. Infinite values are allowed.}
    \item{\code{capacitiesPost}:}{A length-n logical vector whose elements equal to post-merger plant capacities. Infinite values are allowed.}
  \item{\code{demand}:}{A length k character vector specifying whether product demand is linear ("linear") or log-linear ("log").}
}}
 \section{Extends}{
   Class \code{\linkS4class{Bertrand}}, directly.
   Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 2.
 }
\section{Methods}{
  Many of the methods described below contain a \sQuote{preMerger} and
  \sQuote{revenue} argument. The \sQuote{preMerger} takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger values, while FALSE invokes the method using
  the post-merger ownership structure. The \sQuote{revenue} argument also takes
  on a value of TRUE or FALSE, where TRUE invokes the method using
  revenues, while FALSE invokes the method using quantities
  \describe{
    \item{\code{calcMC}}{\code{signature(object,preMerger=TRUE)}}{ Calculates marginal cost for each
      product. }
    \item{\code{calcVC}}{\code{signature(object,preMerger=TRUE)}}{ Calculates variable cost for each
      product.}
    \item{\code{calcMargins}}{\code{signature(object, preMerger
	= TRUE)}}{ Compute either pre-merger or post-merger equilibrium
      margins under the assumption that firms play a multi-product
      Cournot Nash pricing game.}
    \item{\code{\link{cmcr}}}{\code{signature(object)}}{ Calculates compensated marginal cost reduction,
      the percentage decrease in the marginal costs of the merging
      parties' products needed to offset a post-merger price increase.}
    \item{\code{\link{HypoMonTest}}}{\code{signature(object,prodIndex,ssnip=.05,...)}}{
      \code{HypoMonTest} implements the Hypothetical
      Monopolist Test for a given \sQuote{ssnip}.}
    \item{\code{\link{calcPriceDeltaHypoMon}}}{\code{signature(object,prodIndex,...)}}{
      Compute the proportional difference in product prices between the
      prices of products in \sQuote{prodIndex} (i.e. prices set by the
      Hypothetical Monopolist) and prices set in the pre-merger Cournot
      equilibrium. \sQuote{...} may be used to pass arguments to the
      optimizer.}
    \item{\code{hhi}}{\code{signature(object, preMerger= TRUE,revenue=FALSE)}}{ Compute either the pre-merger or
      post-merger Herfindahl-Hirschman Index (HHI) under the assumption
      that firms play a multi-product Cournot Nash pricing game.}
    \item{\code{\link{diversion}}}{\code{signature(object, preMerger
	= TRUE)}}{ Computes a k x k matrix of diversion ratios.}
    \item{\code{summary}}{\code{signature(object,revenue=TRUE,shares=TRUE,parameters=FALSE,digits=2)}}{ Summarizes
      the effect of the merger, including price and revenue
      changes. Setting \sQuote{revenue} equal to FALSE reports quantities rather
      than revenues. Setting \sQuote{shares} to FALSE reports quantities
      rather than than shares (when possible). Setting \sQuote{parameters} equal to TRUE reports
      all demand parameters. \sQuote{digits} controls the number of significant
      digits reported in output.}
    \item{\code{cmcr}}{\code{signature(object)}}{
      Calculates Compensating Marginal Cost Reduction (CMCR).
    }
     }}

\author{Charles Taragin \email{charles.taragin@usdoj.gov}}


\examples{
showClass("Cournot")           # get a detailed description of the class
showMethods(classes="Cournot") # show all methods defined for the class
}
\keyword{classes}
