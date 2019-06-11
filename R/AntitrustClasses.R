#' @import stats
#' @import methods
#' @import BB
#' @import numDeriv
#' @import rhandsontable
#'@title \dQuote{Antitrust} Classes
#'@name Antitrust-Class
#'@aliases Antitrust-class
#'matrixOrList-class
#'matrixOrVector-class
#'characterOrList-class
#'@description The \dQuote{Antitrust} class is a building block used to create other classes
#'in this package. As such, it is most likely to be useful for developers
#'who wish to code their own calibration/simulation routines.
#'@description Let k denote the number of products produced by all firms below.

#'@section Objects from the Class:
#'Objects can be created by calls of the form \code{new("Antitrust", ...)}.
#'@slot pricePre A length k vector of simulated pre-merger prices.
#'@slot pricePost A length k vector of simulated post-merger prices.
#'@slot ownerPre  A k x k matrix of pre-merger ownership shares.
#'@slot ownerPost  A k x k matrix of post-merger ownership shares.
#'@slot labels A length k vector of labels.
#'@slot control.slopes A list of \code{\link{optim}} control parameters passed to the calibration routine optimizer (typically the \code{calcSlopes} method).
#'@slot control.equ A list of \code{\link[BB]{BBsolve}} control parameters passed to the non-linear equation solver (typically the \code{calcPrices} method).

#'@section The \dQuote{matrixOrList}, \dQuote{matrixOrVector} and \dQuote{characterOrList} Classes:
#'The \dQuote{matrixOrList},\dQuote{matrixOrVector} and \dQuote{characterOrList} classes are
#'virtual classes used for validity checking in the \sQuote{ownerPre} and
#'\sQuote{ownerPost} slots of \dQuote{Antitrust} and the \sQuote{slopes} slot in
#'\dQuote{Bertrand}.
#'@author Charles Taragin \email{ctaragin@ftc.gov}
#'@examples showClass("Antitrust")           # get a detailed description of the class
#'@keywords classes
NULL

setClassUnion("matrixOrVector", c("matrix", "numeric","character","factor"))
setClassUnion("matrixOrList", c("numeric","matrix", "list"))
setClassUnion("characterOrList", c("character", "list"))

setClass(

  Class = "Antitrust",
  representation=representation(
    ownerPre     = "matrixOrVector",
    ownerPost    = "matrixOrVector",
    pricePre     = "numeric",
    pricePost    = "numeric",
    mcPre        = "numeric",
    mcPost       = "numeric",
    labels       = "characterOrList",
    cls        = "character",
    control.slopes = "list",
    control.equ = "list"
  ),
  prototype(
    pricePre  = numeric(),
    pricePost = numeric(),
    mcPre     = numeric(),
    mcPost    = numeric(),
    cls  =    character(),
    ##copied from 'optim' definition:
    control.slopes = list(
      reltol = sqrt(.Machine$double.eps)
    ),
    ## copied from 'BBsolve' definition:
    ## changed 'maxit' to 2000 from 1500
    ## changed tol to sqrt(.Machine$double.eps) < 1e-07 (default)
    control.equ = list(maxit = 2000, M = c(50, 10),
                       tol = sqrt(.Machine$double.eps), trace = FALSE,
                       triter = 10,
                       noimp = 100, NM = c(TRUE, FALSE))
  ),
  validity=function(object){



    if(is.list(object@labels)){ nprods <- length(object@labels[[1]])}
    else{nprods <- length(object@labels)}


    if(is.matrix(object@ownerPre)){

      if(nprods != ncol(object@ownerPre)){
        stop("The number of rows and columns in 'ownerPre' must equal the number of products (plants in Cournot)")}
      if(nrow(object@ownerPre) != ncol(object@ownerPre)){
        stop("'ownerPre' must be a square matrix ")}

      if(
        any(colSums(unique(object@ownerPre),na.rm=TRUE)>1)
      ){
        stop("The columns of the matrix formed from the unique rows of 'ownerPre' must sum to no more than 1")
      }
    }

    else if (nprods != length(object@ownerPre)) stop("'ownerPre' and 'labels' must be vectors of the same length")
    if(is.matrix(object@ownerPost)){
      if(nprods != ncol(object@ownerPost)){
        stop("The number of rows and columns in 'ownerPost' must equal the number of products")}
      if(nrow(object@ownerPost) != ncol(object@ownerPost)){
        stop("'ownerPost' must be a square matrix")}
      if(
        any(colSums(unique(object@ownerPost))>1,na.rm=TRUE)
      ){
        stop("The columns of the matrix formed from the unique rows of 'ownerPost' must sum to no more than 1")
      }
    }

    else if (nprods != length(object@ownerPost)) stop("'ownerPost' and 'labels' must be vectors of the same length")

    if(identical(object@ownerPre,object@ownerPost)) warning("'ownerPost' and 'ownerPre' are the same")
    return(TRUE)
  }

)
