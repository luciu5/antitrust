#'@title \dQuote{Bargaining} Classes
#'@name Bargaining-Class
#'@include BertrandClasses.R
#'@description The \dQuote{Bargaining} class is a building block used to create other classes
#'in this package. As such, it is most likely to be useful for developers
#'who wish to code their own calibration/simulation routines.
#'@description Let k denote the number of products produced by all firms below.

#'@section Objects from the Class:
#'Objects can be created by calls of the form \code{new("Bargaining", ...)}.
#'@slot bargpower A length k vector of calibrated bargaining power parameters.
#'@slot prices A length k vector of of observes prices.
#'@slot margins  A length k vector of of observes margins.


setClass(
  
  Class = "Bargaining",
  contains="Bertrand",
  representation=representation(
    bargpower       = "numeric",
    prices           = "numeric",
    margins          = "numeric",
    priceStart       = "numeric"
  ),
  prototype=prototype(
    
    bargpower          = numeric(),
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
      !(all(object@bargpower >=0) &&
        all(object@bargpower <=1))
    ){
      stop("elements of vector 'bargpower' must be between 0 and 1")
    }
   
    if(nprods != length(object@priceStart)){
      stop("'priceStart' must have the same length as 'prices'")}
    
    
  }
)
