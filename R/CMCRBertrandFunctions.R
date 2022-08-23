#' @title Compensating Marginal Cost Reductions and Upwards Pricing Pressure (Bertrand)
#' @name CMCRBertrand-Functions
#' @aliases cmcr.bertrand
#' upp.bertrand
#' upp
#' cmcr
#'
#' @description Calculate the marginal cost reductions necessary to restore
#' premerger prices (CMCR), or the net Upwards Pricing Pressure (UPP) in a
#' merger involving firms playing a differentiated products Bertrand pricing game.
#' @description Let k denote the number of products produced by the merging parties below.
#'
#' @param prices A length-k vector of product prices.
#' @param margins A length-k vector of product margins.
#' @param diversions A k x k matrix of diversion ratios
#' with diagonal elements equal to -1.
#' @param ownerPre EITHER a vector of length k whose values
#' indicate which of the merging parties produced a product pre-merger OR
#' a k x k matrix of pre-merger ownership shares.
#' @param ownerPost A k x k matrix of post-merger
#' ownership shares. Default is a k x k matrix of 1s.
#' @param mcDelta A vector of length k where each element equals the
#' proportional change in a product's marginal costs due to
#' the merger. Default is 0, which assumes that the merger does not
#' affect any products' marginal cost.
#' @param rel A length 1 character vector indicating whether CMCR should be calculated
#' relative to pre-merger cost (``cost'') or pre-merger price (``price''), 
#' Default is ``cost''.
#' @param labels A length-k vector of product labels.
#'
#' @details All \sQuote{prices} elements must be positive, all \sQuote{margins}
#' elements must be between 0 and 1, and all \sQuote{diversions} elements must be
#' between 0 and 1 in absolute value. In
#' addition, off-diagonal elements (i,j) of \sQuote{diversions}
#' must equal an estimate of the diversion ratio from product i to product j
#' (i.e. the estimated fraction of i's sales that go to j due to a small
#' increase in i's price). Also, \sQuote{diversions}
#' elements are positive if i and j are substitutes and negative if i and j are
#' complements.
#'
#' \sQuote{ownerPre} will typically be a vector whose values equal 1 if a product is produced
#' by firm 1 and 0 otherwise, though other values including firm name are
#' acceptable.  Optionally, \sQuote{ownerPre} may be set equal
#' to a matrix of the merging firms pre-merger ownership
#' shares. These ownership shares must be between 0 and 1.
#'
#' \sQuote{ownerPost}  is an optional argument that should only be specified if
#' one party to the acquisition is assuming partial control of the
#' other party's assets.   \sQuote{ownerPost} elements must be between 0 and
#' 1.
#'
#' @return \code{cmcr.bertrand} returns a length-k vector whose values
#' equal the percentage change in each products'
#' marginal costs that the merged firms must achieve in order to offset a
#' price increase.
#'
#' \code{upp.bertrand} returns a length-k vector whose values equal
#' the generalized pricing pressure (GePP) for each of the merging's parties' products,
#' net any efficiency claims. GePP is a generalization of Upwards Pricing Pressure (UPP)
#' that accomodates multi-product firms.
#'
#' @references Farrell, Joseph and Shapiro, Carl (2010).
#' \dQuote{Antitrust Evaluation of Horizontal Mergers: An Economic Alternative to
#' Market Definition.}
#' \emph{The B.E. Journal of Theoretical Economics}, \bold{10}(1), pp. 1-39.
#'
#' Jaffe, Sonia and Weyl Eric (2012).
#' \dQuote{The First-Order Approach to Merger Analysis.}
#' \emph{SSRN eLibrary}
#'
#' Werden, Gregory (1996).
#' \dQuote{A Robust Test for Consumer Welfare Enhancing Mergers Among Sellers
#' of Differentiated Products.}
#' \emph{The Journal of Industrial Economics}, \bold{44}(4), pp. 409-413.
#'
#' @author Charles Taragin \email{ctaragin+antitrustr@gmail.com}
#' @seealso \code{\link{cmcr.cournot}} for a homogeneous products Cournot
#' version of CMCR, and \code{\link{cmcr-methods}} for calculating
#' CMCR and UPP after calibrating demand system parameters and simulating a merger.
#'
#' @examples
#' ## Let k_1 = 1 and and k_2 = 2 ##
#'
#' p1 = 50;      margin1 = .3
#' p2 = c(45,70); margin2 = c(.4,.6)
#' isOne=c(1,0,0)
#' diversions = matrix(c(-1,.5,.01,.6,-1,.1,.02,.2,-1),ncol=3)
#'
#' cmcr.bertrand(c(p1,p2), c(margin1,margin2), diversions, isOne)
#' upp.bertrand(c(p1,p2), c(margin1,margin2), diversions, isOne)
#'
#'
#' ## Calculate the necessary percentage cost reductions for various margins and
#' ## diversion ratios in a two-product merger where both products have
#' ## equal prices and diversions (see Werden 1996, pg. 412, Table 1)
#'
#'
#' margins = seq(.4,.7,.1)
#' diversions = seq(.05,.25,.05)
#' prices = rep(1,2) #assuming prices are equal, we can set product prices to 1
#' isOne = c(1,0)
#' result = matrix(ncol=length(margins),nrow=length(diversions),dimnames=list(diversions,margins))
#'
#' for(m in 1:length(margins)){
#'   for(d in 1:length(diversions)){
#'
#'     dMatrix = -diag(2)
#'     dMatrix[2,1] <- dMatrix[1,2] <- diversions[d]
#'
#'     firmMargins = rep(margins[m],2)
#'
#'     result[d,m] = cmcr.bertrand(prices, firmMargins, dMatrix, isOne)[1]
#'
#'   }}
#'
#' print(round(result,1))
#'
#' @include CESFunctions.R
NULL

#'@rdname CMCRBertrand-Functions
#'@export
cmcr.bertrand <- function(prices, margins, diversions, ownerPre, 
                          ownerPost=matrix(1,ncol=length(prices), nrow=length(prices)),
                          rel=c("cost","price"),
                          labels=names(prices))
{



  rel=match.arg(rel)
  
  
  if(!(is.vector(prices) & is.vector(margins))){
    stop("'prices' and 'margins' must be vectors")}

  nprod = length(prices)


  if(nprod != length(margins)){
    stop("'prices'and 'margins' vectors must be the same length")}

  if(any(prices < 0,na.rm=TRUE)){ stop("'prices'  must be non-negative")}
  if(any(margins < 0 | margins > 1,na.rm=TRUE) ){ stop("'margins' vector elements  must be between 0 and 1")}

  if(!is.matrix(diversions)){ stop("'diversions'  must be a matrix")}
  if(!all(diag(diversions) == -1,na.rm=TRUE)){ stop("'diversions' diagonal elements must all equal -1")}
  if(any( abs(diversions) > 1,na.rm=TRUE)){ stop("'diversions' elements must be between -1 and 1")}
  if(ncol(diversions)!=nrow(diversions) ||
     ncol(diversions)!= nprod){
    stop("'diversions' must be a square matrix whose dimension equals the length of 'prices'")
  }

  if(!is.matrix(ownerPost)){ stop("'ownerPost'  must be a matrix")}
  if(any(ownerPost < 0 | ownerPost > 1,na.rm=TRUE)){ stop("'ownerPost' elements must be between 0 and 1")}
  if(
    !isTRUE(all.equal(colSums(unique(ownerPost)),rep(1,nprod),check.names=FALSE))){
    stop("The columns of the matrix formed from the unique rows of 'ownerPost' must sum to 1")
  }

  ## transform pre-merger ownership vector into matrix##
  if(is.vector(ownerPre) ||
     is.factor(ownerPre)){

    if(nprod != length(ownerPre)){
      stop("'prices'and 'ownerPre' vectors must be the same length")}

    owners <- as.numeric(factor(ownerPre))
    ownerPre <- matrix(0,ncol=nprod,nrow=nprod)

    for( o in unique(owners)){
      ownerPre[owners == o, owners == o] = 1
    }

    rm(owners)
  }




  else if(!is.matrix(ownerPre) ||
          ncol(ownerPre) != nrow(ownerPre) ||
          any(ownerPre < 0 | ownerPre > 1,na.rm=TRUE) ||
          ncol(ownerPre) != nprod ||
          !isTRUE(all.equal(colSums(unique(ownerPre)),rep(1,nprod),check.names=FALSE))
  ){
    stop("'ownerPre' must be a square matrix whose dimension equals the length of 'prices' and whose elements are between 0 and 1. Also, the columns of the matrix formed from the unique rows of 'ownerPre' must sum to 1")}


  if(is.null(labels)){labels <- paste("Prod",1:length(prices),sep="")}

  ## weight diversion ratios by price ratios and ownership matrices ##
  priceRatio = tcrossprod(1/prices, prices)
  Bpre =  -1 * diversions * priceRatio * ownerPre;
  Bpost = -1 * diversions * priceRatio * ownerPost;



  ##calculate post-merger margin ##
  marginPost = as.vector( solve(Bpost) %*%
                            ((diag(ownerPost)/diag(ownerPre)) * as.vector(Bpre %*% margins))
  )

  ## calculate changes in marginal cost relative to pre-merger prices
  mcDelta= (marginPost - margins)
  
  ## calculate changes in marginal cost relative to pre-merger costs
  if(rel=="cost"){mcDelta <- mcDelta/(1 - margins)}

  names(mcDelta) <- labels

  return(mcDelta*100)
}


#'@rdname CMCRBertrand-Functions
#'@export
upp.bertrand <- function(prices, margins, diversions, ownerPre,
                         ownerPost=matrix(1,ncol=length(prices), nrow=length(prices)),
                         mcDelta=rep(0,length(prices)),
                         labels=paste("Prod",1:length(prices),sep=""))
{

  if(!(is.vector(prices) & is.vector(margins))){
    stop("'prices' and 'margins'  must be vectors")}

  nprod = length(prices)


  if(nprod != length(margins)){
    stop("'prices'and 'margins' vectors must be the same length")}

  if(any(prices < 0,na.rm=TRUE)){ stop("'prices' must be non-negative")}
  if(any(margins < 0 | margins > 1,na.rm=TRUE) ){ stop("'margins' vector elements  must be between 0 and 1")}

  if(!is.matrix(diversions)){ stop("'diversions' must be a matrix")}
  if(!all(diag(diversions) == -1,na.rm=TRUE)){ stop("'diversions' diagonal elements must all equal -1")}
  if(any( abs(diversions) > 1,na.rm=TRUE)){ stop("'diversions' elements must be between -1 and 1")}
  if(ncol(diversions)!=nrow(diversions) ||
     ncol(diversions)!= nprod){
    stop("'diversions' must be a square matrix whose dimension equals the length of 'prices'")
  }


  if(!is.matrix(ownerPost)){ stop("'ownerPost' must be a matrix")}
  if(any(ownerPost < 0 | ownerPost > 1,na.rm=TRUE)){ stop("'ownerPost' elements must be between 0 and 1")}
  if(
    !isTRUE(all.equal(colSums(unique(ownerPost)),rep(1,nprod),check.names=FALSE))){
    stop("The columns of the matrix formed from the unique rows of 'ownerPost' must sum to 1")
  }
  if(any(mcDelta<0,na.rm=TRUE)){stop("'mcDelta' must be positive")}

  ## transform pre-merger ownership vector into matrix##
  if(is.vector(ownerPre)  ||
     is.factor(ownerPre)){

    if(nprod != length(ownerPre)){
      stop("'prices'and 'ownerPre' vectors must be the same length")}

    owners <- as.numeric(factor(ownerPre))
    ownerPre <- matrix(0,ncol=nprod,nrow=nprod)

    for( o in unique(owners)){
      ownerPre[owners == o, owners == o] = 1
    }

    rm(owners)
  }




  else if(!is.matrix(ownerPre) ||
          ncol(ownerPre) != nrow(ownerPre) ||
          any(ownerPre < 0 | ownerPre > 1,na.rm=TRUE) ||
          ncol(ownerPre) != nprod ||
          !isTRUE(all.equal(colSums(unique(ownerPre)),rep(1,nprod),check.names=FALSE))
  ){
    stop("'ownerPre' must be a square matrix whose dimension equals then length of 'prices' and whose elements are between 0 and 1 Also, the columns of the matrix formed from the unique rows of 'ownerPre' must sum to 1")
  }



  mcPre  <- prices*(1-margins)
  mcPost <- mcPre*(1+mcDelta)

  marginsPre <- margins
  marginsPost <- 1 - mcPost/prices

  ## weight diversion ratios by price ratios and ownership matrices ##
  priceRatio = tcrossprod(1/prices, prices)
  Bpre =   diversions * priceRatio * ownerPre
  Bpost =  diversions * priceRatio * ownerPost



  result <- as.vector((Bpost %*% marginsPost)/diag(ownerPost) - (Bpre %*% marginsPre)/diag(ownerPre))

  names(result) <- labels

  return(result*100) #net UPP

}
