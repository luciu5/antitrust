#' @title Herfindahl-Hirschman Index
#' @name HHI-Functions
#' @aliases HHI
#' @description Calculate the Herfindahl-Hirschman Index with arbitrary
#' ownership and control.
#' @description Let k denote the number of products produced by the merging parties below.
#' @param shares A length-k vector of product quantity shares.
#' @param owner EITHER a vector  of length k whose values
#' indicate which of the merging parties produced a product OR
#' a k x k matrix of ownership shares. Default is a diagonal matrix,
#' which assumes that each product is owned by a separate firm.
#' @param control EITHER a vector  of length k whose values
#' indicate which of the merging parties have the ability to make pricing
#' or output decisions OR a k x k matrix of
#' control shares. Default is a k x k matrix equal to 1 if \sQuote{owner} > 0
#' and 0 otherwise.
#' @details All \sQuote{shares} must be between 0 and 1. When \sQuote{owner} is a matrix,
#' the i,jth element of \sQuote{owner} should equal the percentage of
#' product j's profits earned by the owner
#' of product i. When \sQuote{owner} is a vector, \code{HHI} generates a k x k
#' matrix of whose i,jth element equals 1 if products i and j are
#' commonly owned and 0 otherwise. \sQuote{control} works in a fashion similar
#' to \sQuote{owner}.
#'
#' @return \code{HHI} returns a number between 0 and 10,000
#' @references Salop, Steven and O'Brien, Daniel (2000)
#' \dQuote{Competitive Effects of Partial Ownership: Financial Interest and Corporate Control}
#' 67 Antitrust L.J. 559, pp. 559-614.
#' @author Charles Taragin \email{ctaragin+antitrustr@gmail.com}
#' @seealso \code{\link{HHI-Methods}} for computing HHI following merger simulation.
#'
#' @examples ## Consider a market with 5 products labeled 1-5. 1,2 are produced
#' ## by Firm A, 2,3 are produced by Firm B, 3 is produced by Firm C.
#' ## The pre-merger product market shares are
#'
#' shares = c(.15,.2,.25,.35,.05)
#' owner  = c("A","A","B","B","C")
#' nprod  = length(shares)
#'
#' HHI(shares,owner)
#'
#' ## Suppose that Firm A acquires a 75\% ownership stake in product 3, and
#' ## Firm B get a 10\% ownership stake in product 1. Assume that neither
#' ## firm cedes control of the product to the other.
#'
#' owner <- diag(nprod)
#'
#' owner[1,2] <- owner[2,1] <- owner[3,4] <- owner[4,3] <- 1
#' control <- owner
#' owner[1,1] <- owner[2,1] <- .9
#' owner[3,1] <- owner[4,1] <- .1
#' owner[1,3] <- owner[2,3] <- .75
#' owner[3,3] <- owner[4,3] <- .25
#'
#' HHI(shares,owner,control)
#'
#' ## Suppose now that in addition to the ownership stakes described
#' ## earlier, B receives 30\% of the control of product 1
#' control[1,1] <- control[2,1] <- .7
#' control[3,1] <- control[4,1] <- .3
#'
#' HHI(shares,owner,control)
#'
#' @include CournotFunctions.R
NULL

#'@rdname HHI-Functions
#'@export
HHI <- function(shares,
                owner=diag(length(shares)),
                control){

  nprod <- length(shares)

  if(!(is.vector(shares))){
    stop("'shares' must be a vector")}

  if(any(shares < 0,na.rm=TRUE) ||
     any(shares >1,na.rm=TRUE) ) {
    stop("'shares'  must be between 0 and 1")}


  ## transform pre-merger ownership vector into matrix##
  if(is.vector(owner) ||
     is.factor(owner)){

    if(nprod != length(owner)){
      stop("'shares' and 'owner' vectors must be the same length")}

    owners <- as.numeric(factor(owner))
    owner <- matrix(0,ncol=nprod,nrow=nprod)

    for( o in unique(owners)){
      owner[owners == o, owners == o] = 1
    }

    rm(owners)

  }

  else if(!is.matrix(owner) ||
          ncol(owner) != nrow(owner) ||
          ncol(owner) != nprod ||
          any(owner < 0 | owner > 1,na.rm=TRUE)){

    stop("'owner' must be a square matrix whose dimensions equal the length of 'shares' and whose elements are between 0 and 1")
  }


  if(missing(control)){
    control <- owner>0
  }

  else if(is.vector(control) ||
          is.factor(control)){

    if(nprod != length(control)){
      stop("'shares' and 'control' vectors must be the same length")}

    controls <- as.numeric(factor(control))
    control <- matrix(0,ncol=nprod,nrow=nprod)

    for( c in unique(controls)){
      control[controls == c, controls == c] = 1
    }

    rm(controls)
  }

  else if (!is.matrix(control) ||
           ncol(control) != nrow(control) ||
           ncol(control) != nprod ||
           any(control < 0 | control > 1,na.rm=TRUE)
  ){
    stop("'control' must be a square matrix whose dimensions equal the length of 'shares' and whose elements are between 0 and 1")

  }

  weights <- crossprod(control,owner)
  weights <- t(t(weights)/diag(weights)) # divide each element by its corresponding diagonal

  shares <- shares*100
  result <- as.vector(shares %*% weights %*% shares)



  return(result)




}

