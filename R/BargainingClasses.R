#'@title \dQuote{Bargaining} Classes
#'@name Bargaining-Class
#'@include BertrandRUMClasses.R
#'@description The \dQuote{Bargaining} class is a building block used to create other classes
#'in this package. As such, it is most likely to be useful for developers
#'who wish to code their own calibration/simulation routines.
#'@description Let k denote the number of products produced by all firms below.

#'@section Objects from the Class:
#'Objects can be created by calls of the form \code{new("Bargaining", ...)}.
#'@slot bargpowerPre A length k vector of pre-merger bargaining power parameters.
#'@slot bargpowerPre A length k vector of post-merger bargaining power parameters.
#'@slot prices A length k vector of of observes prices.
#'@slot margins  A length k vector of of observes margins.


setClass(
  
  Class = "Bargaining",
  contains="Bertrand",
  representation=representation(
    bargpowerPre       = "numeric",
    bargpowerPost       = "numeric",
    prices           = "numeric",
    margins          = "numeric",
    priceStart       = "numeric"
  ),
  prototype=prototype(
    
    bargpowerPre          = numeric(),
    priceStart         = numeric()
    
  ),
  validity=function(object){
    
    nprods <- length(object@prices)
    
    if(length(object@margins) != nprods){stop("'margins' and 'prices' must have the same length")}
    
    if(
      !(all(object@margins >0) &&
        all(object@margins <=1))
    ){
      stop("elements of vector 'margins' must be between 0 and 1")
    }
    
    if(
      !(all(object@bargpowerPre >=0) &&
        all(object@bargpowerPre <=1))
    ){
      stop("elements of vector 'bargpowerPre' must be between 0 and 1")
    }
   
    if(
      !(all(object@bargpowerPost >=0) &&
        all(object@bargpowerPost <=1))
    ){
      stop("elements of vector 'bargpowerPost' must be between 0 and 1")
    }
    
    if(nprods != length(object@priceStart)){
      stop("'priceStart' must have the same length as 'prices'")}
    
    
  }
)

setClass(
  
  Class = "BargainingLogit",
  contains="Logit",
  representation=representation(
    bargpowerPre       = "numeric",
    bargpowerPost       = "numeric"
  ),
  prototype=prototype(
    
    bargpowerPre          = numeric()
  ),
  validity=function(object){
    
    if(
      !(all(object@bargpowerPre >=0) &&
        all(object@bargpowerPre <=1))
    ){
      stop("elements of vector 'bargpowerPre' must be between 0 and 1")
    }
    
    if(
      !(all(object@bargpowerPost >=0) &&
        all(object@bargpowerPost <=1))
    ){
      stop("elements of vector 'bargpowerPost' must be between 0 and 1")
    }
    
  }
    
)
