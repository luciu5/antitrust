\name{Bertrand-class}
\Rdversion{1.1}
\docType{class}
\alias{Bertrand-class}
%\alias{cmcr,Bertrand-method}
%\alias{diversion,Bertrand-method}
%\alias{hhi,Bertrand-method}
\alias{show,Bertrand-method}
%\alias{summary,Bertrand-method}


\title{Class "Bertrand"}
\description{
  The "Bertrand" class is a building block used to create other classes
  in this package. As such, it is most likely to be useful for developers
  who wish to code their own merger calibration/simulation routines.


}
\section{Objects from the Class}{

Objects can be created by calls of the form \code{new("Bertrand", ...)}.

}
\section{Slots}{
  Let k denote the number of products produced by all firms.
  \describe{
    \item{\code{shares}:}{A length k vector containing observed
      output. Depending upon the model, output will be measured in
      units sold, quantity shares, or revenue shares.}
    \item{\code{mcDelta}:}{A length k vector where each element equals the
    proportional change in a product's marginal costs due to
    the merger.}
    \item{\code{slopes}:}{A k x (k+1) matrix of linear demand intercepts and slope coefficients}
  }
}
\section{Methods}{
  Many of the methods described below contain a `preMerger' and
  `revenue' argument. The `preMerger' takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger values, while FALSE invokes the method using
  the post-merger ownership structure. The `revenue' argument also takes
  on a value of TRUE or FALSE, where TRUE invokes the method using
  revenues, while FALSE invokes the method using quantities
  \describe{
    \item{\code{calcMC}}{\code{signature(object,preMerger=TRUE)}}{ Calculates (constant) marginal cost for each
      product. For those classes that do not require prices, returns a
      length-k vector of NAs when prices are not supplied.}
    \item{\code{calcMargins}}{\code{signature(object, preMerger
	= TRUE)}}{ Compute either pre-merger or post-merger equilibrium
      margins under the assumption that firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{\code{\link{cmcr}}}{\code{signature(object)}}{ Calculates compensated marginal cost reduction,
      the percentage decrease in the marginal costs of the merging
      parties' products needed to offset a post-merger price increase.}
    \item{\code{\link{isMarket}}}{\code{signature(object,prodIndex,ssnip=.05,
	supplyResponse=FALSE,...)}}{
      \code{isMarket} uses the Hypothetical
      Monopolist Test to determine whether the product indexed by
      `prodIndex' constitute an antitrust market.  Setting `supplyResponse' equal to FALSE requires
      firms who own products not included in `prodIndex'
      to fix their prices at (predicted) pre-merger levels, while TRUE allows these
      firms to set their product prices optimally.}
    \item{\code{\link{smallestMarket}}}{\code{signature(object,startProd,prodOrder,ssnip=.05,supplyResponse=FALSE,...)}}{
      Beginning with the product
      indexed by `startProd', \code{smallestMarket} uses the Hypothetical
      Monopolist Test described in the Merger Guidelines to construct
      the smallest market that satisfies a
      `ssnip'. Products are placed under the Hypothetical Monopolist's
      control in the order dictated by `prodOrder'. If `prodOrder' is missing, products are placed under the
      Hypothetical Monopolist's control in order from greatest to least
      diversion from `startProd'. Setting `supplyResponse' equal to FALSE requires
      firms who own products not included in the Hypothetical Monopolist
      to fix their prices at (predicted) pre-merger levels, while TRUE allows these
      firms to set their product prices optimally.}
    \item{\code{hhi}}{\code{signature(object, preMerger
	= TRUE,revenue=FALSE)}}{ Compute either the pre-merger or post-merger
      Herfindahl–Hirschman Index (HHI) under the assumption that firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{\code{\link{diversion}}}{\code{signature(object, preMerger
	= TRUE)}}{ Computes a k x k matrix of diversion ratios.}
    \item{\code{summary}}{\code{signature(object,revenue=TRUE,shares=TRUE,parameters=FALSE,digits=2)}}{ Summarizes
      the effect of the merger, including price and revenue
      changes. Setting `revenue' equal to FALSE reports quantities rather
      than revenues. Setting `shares' to FALSE reports quantities
      rather than than shares (when possible). Setting `parameters' equal to TRUE reports
      all demand parameters. `digits' controls the number of significant
      digits reported in output.}
    \item{\code{upp}}{\code{signature(object)}}{
      Calculate the Upwards Pricing Pressure (upp) index.
    }
     }}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("Bertrand")
}
\keyword{classes}
