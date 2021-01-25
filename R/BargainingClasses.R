#'@title \dQuote{Bargaining} Classes
#'@name Bargaining-Classes
#'@aliases Bargaining-class
#'@include BertrandRUMClasses.R
#'
#'@description Each class contains all the information needed to calibrate a specific type of demand system and
#'perform a merger simulation analysis under the assumption that firms are playing a differentiated products Nash Bargaining  game.
#'
#'@description The \dQuote{Bargaining} class is a building block used to create other classes
#'in this package. As such, it is most likely to be useful for developers
#'who wish to code their own calibration/simulation routines.
#'@description The \dQuote{BargainingLogit} class has the information for a Nash Bargaining game with Logit demand.
#'@description Let k denote the number of products produced by all firms below.
#'@section Objects from the Class:
#'Objects can be created by calls of the form \code{new("Bargaining", ...)}.
#'@slot bargpowerPre A length k vector of pre-merger bargaining power parameters.
#'@slot bargpowerPre A length k vector of post-merger bargaining power parameters.
#'@slot prices A length k vector of of observes prices.
#'@slot margins  A length k vector of of observes margins.


#'@rdname Bargaining-Classes
#'@export
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
      !(all(object@margins >0,na.rm=TRUE) &&
        all(object@margins <=1,na.rm=TRUE))
    ){
      stop("elements of vector 'margins' must be between 0 and 1")
    }
    
    if(
      !(all(object@bargpowerPre >=0,na.rm=TRUE) &&
        all(object@bargpowerPre <=1,na.rm=TRUE))
    ){
      stop("elements of vector 'bargpowerPre' must be between 0 and 1")
    }
   
    if(
      !(all(object@bargpowerPost >=0,na.rm=TRUE) &&
        all(object@bargpowerPost <=1,na.rm=TRUE))
    ){
      stop("elements of vector 'bargpowerPost' must be between 0 and 1")
    }
    
    if(nprods != length(object@priceStart)){
      stop("'priceStart' must have the same length as 'prices'")}
    
    
  }
)


#'@rdname Bargaining-Classes
#'@export
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
      !(all(object@bargpowerPre >=0,na.rm = TRUE) &&
        all(object@bargpowerPre <=1,na.rm=TRUE))
    ){
      stop("elements of vector 'bargpowerPre' must be between 0 and 1")
    }
    
    if(
      !(all(object@bargpowerPost >=0,na.rm=TRUE) &&
        all(object@bargpowerPost <=1,na.rm=TRUE))
    ){
      stop("elements of vector 'bargpowerPost' must be between 0 and 1")
    }
    
  }
    
)
