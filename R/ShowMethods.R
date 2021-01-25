#' @title Show Method
#' @name Show-Methods
#' @docType methods
#'
#' @aliases show,Antitrust-method
#' @param object An instance of the Antitrust class.
#' @description Displays the percentage change in prices due to the merger.
#'
#' @include SummaryMethods.R
#' @keywords methods
NULL

## print method
#'@rdname Show-Methods
#'@export
setMethod(
  f= "show",
  signature= "Antitrust",
  definition=function(object){

    res <- summary(object,market=TRUE)

    return(NULL)
  }
)
#'@rdname Show-Methods
#'@export
setMethod(
  f= "show",
  signature= "VertBargBertLogit",
  definition=function(object){
    
    res <- summary(object,market=TRUE)
    
    return(NULL)
  }
)
