#'@title \dQuote{Vertical} Classes
#'@name Vertical-Classes
#'@include BertrandRUMClasses.R
#'@description The \dQuote{Vertical} class is a building block used to create other classes
#'in this package. As such, it is most likely to be useful for developers
#'who wish to code their own calibration/simulation routines.
#'@description Let k denote the number of products produced by all firms below.

#'@section Objects from the Class:
#'Objects can be created by calls of the form \code{new("Vertical", ...)}.
#'@slot bargpower A length k vector of calibrated Vertical power parameters.
#'@slot prices A length k vector of of observes prices.
#'@slot margins  A length k vector of of observes margins.

#'@rdname Vertical-Classes
#'@export

setClass(
  
  Class = "VertBargBertLogit",
  representation=representation(
     up = "Bargaining",
     down = "Logit"
  )
 
)

setClass(
  
  Class = "VertBarg2ndLogit",
  representation=representation(
    up = "Bargaining",
    down = "Auction2ndLogit"
  )
  
)
