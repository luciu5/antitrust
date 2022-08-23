#'@title \dQuote{Cournot} Classes
#'@name Cournot-classes
#'@aliases Cournot-class
#'Stackelberg-class
#'@description The \dQuote{Cournot} and \dQuote{Stackelberg} classes are building blocks used to create other classes
#'in this package. As such, they are most likely to be useful for developers
#'who wish to code their own merger calibration/simulation routines.
#'@description Note below that k is the number of products and n is the number of plants.
#'@section Objects from the Class:
#'For Cournot, objects can be created by calls of the form \code{new("Cournot", ...)}.
#'
#'For Stackelberg, objects can be created by calls of the form \code{new("Stackelberg", ...)}.

#'@slot intercepts A length k vector containing the calibrated demand intercept.
#'@slot mcfunPre A length n list whose elements equal a function that calculates a plant's pre-merger marginal cost.
#'@slot mcfunPost A length n list whose elements equal a function that calculates a plant's post-merger marginal cost.
#'@slot vcfunPre A length n list whose elements equal a function that calculates a plant's pre-merger variable cost.
#'@slot vcfunPost A length n list whose elements equal a function that calculates a plant's post-merger variable cost.
#'@slot prices A length k vector of product prices.
#'@slot quantities An n x k matrix of plant quantities produced for each product.
#'@slot margins An n x k matrix of plant product margins.
#'@slot quantityPre An n x k matrix of predicted pre-merger quantities.
#'@slot quantityPost An n x k matrix of predicted post-merger quantities.
#'@slot quantityStart A length n x k vector of starting quantities for the non-linear solver.
#'@slot productsPre An n x k logical matrix qhose elements are TRUE if a plant produces a product pre-merger and FALSE otherwise.
#'@slot productsPost An n x k logical matrix qhose elements are TRUE if a plant produces a product post-merger and FALSE otherwise.
#'@slot capacitiesPre A length-n logical vector whose elements equal to pre-merger plant capacities. Infinite values are allowed.
#'@slot capacitiesPost A length-n logical vector whose elements equal to post-merger plant capacities. Infinite values are allowed.
#'@slot demand A length k character vector specifying whether product demand is linear ("linear") or log-linear ("log").
#'@slot cost A length k character vector equal to "linear" if a plant's marginal cost curve is assumed to be linear
#'or "constant" if a plant's marginal curve is assumed to be constant.
#'Returns an error if a multi-plant firm with constant marginal costs does not have capacity constraints.
#'@slot mktElast A length k vector of market elasticities.
#'@slot dmcfunPre A length n list whose elements equal a function that calculates the derivative of a
#'plant's pre-merger marginal cost with respect to that plant's output. (Stackelberg only)
#'@slot dmcfunPost A length n list whose elements equal a function that calculates the derivative of a
#'plant's post-merger marginal cost with respect to that plant's output. (Stackelberg only)
#'@slot isLeaderPre An n x k logical matrix qhose elements are TRUE if a plant produces a product pre-merger and FALSE otherwise. (Stackelberg only)
#'@slot isLeaderPost An n x k logical matrix qhose elements are TRUE if a plant produces a product post-merger and FALSE otherwise. #'@slot dmcfunPre A length n list whose elements equal a function that calculates the derivative of a
#'plant's pre-merger marginal cost with respect to that plant's output. (Stackelberg only)
#'@section Extends:
#'Cournot:
#'Class \code{\linkS4class{Bertrand}}, directly.
#'Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 2.
#'
#'Stackelberg:
#'Class \code{\linkS4class{Cournot}}, directly.
#'Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Cournot}}, distance 2.
#'Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 3.
#'@author Charles Taragin \email{ctaragin+antitrustr@gmail.com}
#'@examples
#'showClass("Cournot")           # get a detailed description of the class
#'showClass("Stackelberg")           # get a detailed description of the class
#'@keywords classes
#'@include AuctionClasses.R
NULL

#'@rdname Cournot-classes
#' @export
setClass(

  Class = "Cournot",
  contains="Bertrand",
  representation=representation(

    intercepts       = "numeric",
    mcfunPre           = "list",
    mcfunPost           = "list",
    vcfunPre           = "list",
    vcfunPost           = "list",
    prices           = "vector",
    quantities       = "matrix",
    margins          = "matrix",
    quantityPre      = "matrix",
    quantityPost     = "matrix",
    capacitiesPre           = "numeric",
    capacitiesPost           = "numeric",
    quantityStart     = "numeric",
    productsPre      = "matrix",
    productsPost     = "matrix",
    demand           = "character",
    cost             = "character",
    mktElast        =  "numeric"
  ),
  prototype(
    intercepts    =  numeric(),
    mcfunPre =list(),
    mcfunPost=list(),
    vcfunPre =list(),
    vcfunPost=list(),
    mktElast = NA_real_
  ),
  validity=function(object){

    nplants <- nrow(object@quantities) # count the number of plants
    nprods  <- ncol(object@quantities)     # count the number of products


    if(nplants != length(object@capacitiesPre)){
      stop("capacitiesPre' must be a vector whose length equals the number of rows in 'quantities'")}

    if( any(object@capacitiesPre<0 ,na.rm=TRUE)){stop("'capacitiesPre' values must be positive, and not NA")}



    if(nplants != length(object@capacitiesPost)){
      stop("capacitiesPre' must be a vector whose length equals the number of rows in 'quantities'")}

    if(any(is.na(object@capacitiesPost) |
           object@capacitiesPost<0 ,na.rm=TRUE)){stop("'capacitiesPost' values must be positive, and not NA")}

    if(
      !length(object@mcfunPre) %in% c(nplants,0)) {stop("'mcfunPre' must be a list of functions whose length equals the number of plants")}
    if(length(object@mcfunPre) >0 && any(sapply(object@mcfunPre,class) != "function"))
    {stop("'mcfunPre' must be a list of functions")}

    if(!identical(length(object@mcfunPre),length(object@vcfunPre) )){stop("'mcfunPre' and 'vcfunPre' should be lists of the same length")}


    if(
      !length(object@mcfunPost) %in% c(nplants,0)) {stop("'mcfunPost' must be a list of functions whose length equals the number of plants")}
    if(length(object@mcfunPost) >0 && any(sapply(object@mcfunPost,class) != "function"))
    {stop("'mcfunPost' must be a list of functions")}

    if(
      !length(object@vcfunPre) %in% c(nplants,0)) {stop("'vcfunPre' must be a list of functions whose length equals the number of plants")}
    if(length(object@vcfunPre) >0 && any(sapply(object@vcfunPre,class) != "function"))
    {stop("'vcfunPre' must be a list of functions")}

    if(
      !length(object@vcfunPost) %in% c(nplants,0)) {stop("'vcfunPost' must be a list of functions whose length equals the number of plants")}
    if(length(object@vcfunPost) >0 && any(sapply(object@vcfunPost,class) != "function"))
    {stop("'vcfunPost' must be a list of functions")}

    if(!is.logical(object@productsPre)) stop("'productsPre' must be a logical matrix")
    if(!is.logical(object@productsPost)) stop("'productsPost' must be a logical matrix")

    if (!identical(dim(object@quantities), dim(object@margins))) stop("'margins' and 'quantities' must be matrices of the same dimension")
    if (!identical(dim(object@quantities), dim(object@productsPre))) stop("'productsPre' and 'quantities' must be matrices of the same dimension")
    if (!identical(dim(object@quantities), dim(object@productsPost))) stop("'productsPost' and 'quantities' must be matrices of the same dimension")


    if(!is.list(object@labels)) stop("'labels' must be a list")
    if (isTRUE(nplants != length(object@labels[[1]]))) stop("'labels' length must be a list whose first element is a vector whose length equals the number of plants")

    if (isTRUE(nprods != length(object@labels[[2]]))) stop("'labels' length must be a list whose 2nd element is a vector whose length equals the number of products")


    if(any(object@prices<=0,na.rm=TRUE))             stop("'prices' values must be positive")
    if(any(object@quantities<=0,na.rm=TRUE))          stop("'quantities' values must be positive")
    if(any(object@margins<0 | object@margins>1,na.rm=TRUE)) stop("'margins' values must be between 0 and 1")
    if(any(colSums(object@margins,na.rm=TRUE) == 0  & is.na(object@mktElast))) stop("at least one plant margin or a market elasticity must be supplied for each product")

    if(all(!(object@demand %in% c("linear","log")))){stop("'demand' must equal 'linear' or 'log'")}
    if(length(object@demand) != nprods) stop("the length of 'demand' must equal the number of products")

    if(all(!(object@cost %in% c("linear","constant")))){stop("'cost' must equal 'linear' or 'constant'")}
    if(length(object@cost) != nplants) stop("the length of 'cost' must equal the number of products")

    #if(any(rowSums(object@ownerPre) > 0 & !(is.finite(object@capacitiesPre) & is.finite(object@capacitiesPost))  & object@cost == "constant")){ stop("multi-plant firms with 'cost' equal to 'constant' must have finite pre- and post-merger capacity constraints")}

    if(length(object@prices) != nprods) stop("the length of 'prices' must equal the number of products")

    if(ncol(object@quantities) != nprods) stop("the number of columns in 'quantities' must equal the number of products")

    if(any(rowSums(object@quantities,na.rm=TRUE) > object@capacitiesPre)){
      stop("pre-merger plant output must be less than pre-merger capacity constraint")
    }



  }

)

#'@rdname Cournot-classes
#' @export
setClass(

  Class = "Stackelberg",
  contains="Cournot",
  representation=representation(
    isLeaderPre = "matrix",
    isLeaderPost = "matrix",
    dmcfunPre           = "list",
    dmcfunPost          = "list"
  ),
  prototype(
    isLeaderPre = matrix(),
    isLeaderPost = matrix(),
    dmcfunPre = list(),
    dmcfunPost=list()
  ),
  validity=function(object){

    nplants <- nrow(object@quantities) # count the number of plants
    nprods  <- length(object@prices)     # count the number of products

    if(!is.logical(object@isLeaderPre)) stop("'leaderPre' must be a logical matrix")
    if(!is.logical(object@isLeaderPost)) stop("'leaderPost' must be a logical matrix")

    if(!identical(dim(object@quantities), dim(object@isLeaderPre))){stop("'isLeaderPre' must be a logical matrix whose dimensions must equal 'quantities' ")}
    if(!identical(dim(object@quantities), dim(object@isLeaderPost))){stop("'isLeaderPost' must be a logical matrix whose dimensions must equal 'quantities'")}

    if(
      !length(object@dmcfunPre) %in% c(nplants,0)) {stop("'dmcfunPre' must be a list of functions whose length equals the number of plants")}
    if(length(object@dmcfunPre) >0 && any(sapply(object@dmcfunPre,class) != "function"))
    {stop("'dmcfunPre' must be a list of functions")}

    if(
      !length(object@dmcfunPost) %in% c(nplants,0)) {stop("'dmcfunPost' must be a list of functions whose length equals the number of plants")}
    if(length(object@dmcfunPost) >0 && any(sapply(object@dmfunPost,class) != "function"))
    {stop("'dmcfunPost' must be a list of functions")}


  }
)
