#'@title \dQuote{Vertical} Classes
#'@name Vertical-Classes
#'@aliases VertBargBertLogit VertBarg2ndLogit 

#'@description The \dQuote{Vertical} classes arebuilding blocks used to create other classes
#'in this package. As such, it is most likely to be useful for developers
#'who wish to code their own calibration/simulation routines.

#'@description The \dQuote{VertBargBertLogit} class has the information for a Vertical Supply Chain with Logit demand and a downstream Nash-Bertrand Pricing game.

#'@description The \dQuote{VertBarg2ndLogit} class has the information fora Vertical Supply Chain with Logit demand and a downstream 2nd Score Auction.

#'@description Let k denote the number of products produced by all firms below.
#'@slot bargpower A length k vector of calibrated Vertical power parameters.
#'@slot prices A length k vector of of observes prices.
#'@slot margins  A length k vector of of observes margins.

#'@author Charles Taragin \email{ctaragin@ftc.gov}
#'@include BertrandRUMClasses.R AuctionClasses.R BargainingClasses.R
#'@keywords classes
#'
NULL

#'@rdname Vertical-Classes
#@export
setClass(
  
  Class = "VertBargBertLogit",
  representation=representation(
     up = "Bargaining",
     down = "Logit",
     constrain="character"
  )
 
)

#'@rdname Vertical-Classes
#@export
setClass(
  
  Class = "VertBarg2ndLogit",
  contains = "VertBargBertLogit"
  ,
  representation=representation(
    up = "Bargaining",
    down = "Auction2ndLogit"
  )
  
)
