\name{setRetention}
\alias{setRetention}
\alias{getRetention}
\title{Product revenue retention}
\usage{
setRetention(object, retentionPre = NULL, retentionPost = NULL)

getRetention(object, preMerger = TRUE)
}
\arguments{
\item{object}{An antitrust model or \code{AntitrustFit}.}
\item{retentionPre, retentionPost}{Positive finite product-level revenue retention vectors. Scalars are recycled and named vectors are matched to product labels.}
\item{preMerger}{Return pre-merger retention when \code{TRUE}.}
}
\value{\code{setRetention()} returns the updated object; \code{getRetention()} returns its product-level retention vector.}
\description{Revenue retention is the share of product revenue received by the seller. Values below one represent ad valorem charges and values above one represent subsidies. The default is one for every product. An omitted argument preserves existing retention; for a new model, post-merger retention inherits pre-merger retention. Economic ownership remains separate. \code{setRetention()} changes metadata only: it does not recover costs, recalibrate, or solve prices. Prefer \code{revenueRetentionPre} in \code{specify()} or \code{calibrate()} and \code{revenueRetentionPost} in \code{simulate()} for fitted models. When changing retention on an already fitted legacy model, refresh its stored \code{mcPre} and \code{mcPost} from \code{calcMC()} before calling \code{calcPrices()}; price solvers use those stored effective costs.}
