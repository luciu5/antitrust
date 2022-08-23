#'@title \dQuote{Bertrand} Classes
#'@name BertrandOther-Classes
#'@description The \dQuote{Bertrand} class is a building block used to create other classes
#'in this package. As such, it is most likely to be useful for developers
#'who wish to code their own merger calibration/simulation routines.

#'@description Each class below contains all the information needed to calibrate a specific type of demand system and
#'perform a merger simulation analysis under the assumption that firms are playing a differentiated products Bertrand pricing game.

#'@description The \dQuote{Linear} class has the information for a Linear demand system.
#'@description The \dQuote{LogLin} class has the information for a Log-Linear demand system.
#'@description The \dQuote{AIDS} class has the information for a AIDS demand system.
#'@description The \dQuote{PCAIDS} class has the information for a PCAIDS demand system
#'@description The \dQuote{PCAIDSNests} class has the information for a nested PCAIDS demand system
#'@description Below, let k denote the number of products produced by all firms.

#'@section Objects from the Class:
#'For Bertrand, objects can be created by calls of the form \code{new("Bertrand", ...)}.
#'
#'For Linear, objects can be created by using the constructor function \code{\link{linear}}.
#'
#'For LogLin, objects can be created by using the constructor function \code{\link{loglin}}.
#'
#'For AIDS, objects can be created by using the constructor function \code{\link{aids}}.
#'
#'For PCAIDS, objects can be created by using the constructor \code{\link{pcaids}}.
#'
#'For nested PCAIDS, objects can be created by using the constructor \code{\link{pcaids.nests}}.

#'@slot shares A length k vector containing observed output. Depending upon the model, output will be measured in
#'units sold, quantity shares, or revenue shares.
#'@slot mcDelta A length k vector where each element equals the proportional change in a product's marginal costs due to the merger.
#'@slot slopes A k x (k+1) matrix of linear demand intercepts and slope coefficients
#'@slot subset A vector of length k where each element equals TRUE if the product indexed by that element should be included in the
#'post-merger simulation and FALSE if it should be excluded.
#'@slot intercepts A length k vector of demand intercepts. (Linear only)
#'@slot prices A length k vector product prices. (Linear only)
#'@slot quantities A length k vector of product quantities. (Linear only)
#'@slot margins A length k vector of product margins. All margins must be between 0 and 1. (Linear only)
#'@slot diversion A k x k matrix of diversion ratios with diagonal elements equal to -1.
#'@slot priceStart A length k vector of prices used as the initial guess in the nonlinear equation solver. (Linear and AIDS only)
#'@slot symmetry If TRUE, requires the matrix of demand slope coefficients to be consistent with utility maximization theory.
#'Default is false. (Linear and LogLin only)
#'@slot insideSize A positive number equal to total pre-merger revenues for all products included in the simulation. (AIDS only)
#'@slot mktElast A negative number equal to the industry pre-merger price elasticity. (AIDS only)
#'@slot parmStart A length 2 vector who elements equal to an initial of a single diagonal element of the matrix of slope coefficients,
#'as well as the market elasticity. (AIDS only)
#'@slot priceDelta A length k vector containing the simulated price effects from the merger. (AIDS only)
#'@slot knownElast A negative number equal to the pre-merger own-price elasticity for any of the k products. (PCAIDS only)
#'@slot knownElastIndex An integer equal to the position of the `knownElast' product in the \sQuote{shares} vector. (PCAIDS only)
#'@slot nests A length k vector identifying which nest a product belongs to. (Nested PCAIDS only)
#'@slot nestsParms A length k vector containing nesting parameters. (Nested PCAIDS only)

#'@section Extends:
#'Bertrand: Class \code{\linkS4class{Antitrust}}, directly.
#'
#'Linear: Class \code{\linkS4class{Bertrand}}, directly.
#'Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 2.
#'
#'LogLin: Class \code{\linkS4class{Linear}}, directly.
#'Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Linear}}, distance 2.
#'Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 3.
#'
#'AIDS: Class \code{\linkS4class{Linear}}, directly. Class \code{\linkS4class{Bertrand}}, by class \dQuote{Linear}, distance 2.
#'
#'PCAIDS: Class \code{\linkS4class{AIDS}}, directly. Class \code{\linkS4class{Linear}}, by class \code{\linkS4class{AIDS}}, distance 2.
#'Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Linear}}, distance 3. Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 4.
#'
#'Nested PCAIDS: Class \code{\linkS4class{PCAIDS}}, directly. Class \code{\linkS4class{AIDS}}, by class \code{\linkS4class{PCAIDS}}, distance 2.
#'Class \code{\linkS4class{Linear}}, by class \code{\linkS4class{AIDS}}, distance 3. Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Linear}}, distance 4.
#'Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 5.

#'@author Charles Taragin \email{ctaragin+antitrustr@gmail.com}
#'@examples
#'showClass("Bertrand")           # get a detailed description of the class
#'showClass("Linear")           # get a detailed description of the class
#'showClass("LogLin")           # get a detailed description of the class
#'showClass("AIDS")           # get a detailed description of the class
#'showClass("PCAIDS")           # get a detailed description of the class
#'showClass("PCAIDSNests")           # get a detailed description of the class
#'@keywords classes
#'
#'@include AntitrustClasses.R
NULL

#'@rdname BertrandOther-Classes
#'@export
setClass(

  Class = "Bertrand",
  contains="Antitrust",
  representation=representation(
    shares       = "numeric",
    mcDelta      = "numeric",
    slopes       = "matrixOrList",
    subset       = "logical",
    diversion        = "matrix"
  ),
  prototype=prototype(

    slopes          = matrix(),
    mcDelta         = numeric(),
    subset          =  logical(),
    diversion       = matrix()

  ),
  validity=function(object){

    if(is.list(object@labels)){ nprods <- length(object@labels[[1]])}
    else{nprods <- length(object@labels)}


    if(!is.list(object@labels) &&
       (nprods != length(object@shares) ||
        nprods != length(object@subset))){
      stop("'labels', 'shares', and 'subset' must all have the same length")}

    if(any(object@shares < 0 | object@shares > 1,na.rm=TRUE)){
      stop("'shares' values must be between 0 and 1")}

    if(!(sum(object@shares,na.rm=TRUE) < 1 ||
         isTRUE(all.equal(sum(object@shares),1,check.names=FALSE, tolerance = 1e-6)))){
      stop("The sum of 'shares' values must be less than or equal to 1")}


    if(nprods != length(object@mcDelta) ||
       any(is.na(object@mcDelta))){
      stop("'mcDelta' must be a numeric vector with the same length as 'shares' and no element of 'mcDelta' can equal NA")}

    if(any(object@mcDelta>0,na.rm=TRUE)){
      warning("positive values of 'mcDelta' imply an INCREASE in marginal costs")}
    
    diversion <- object@diversion
    
    if(!all(is.na(diversion))){
      
    if(!isTRUE(all.equal(diag(diversion),rep(-1,nprods), check.names=FALSE))){ stop("'diversions' diagonal elements must all equal -1")}
    
    allhavezeros <- all(apply(diversion,1,function(x){any(x==0)}))
    if(allhavezeros){stop("every row of 'diversions' contains zeros. Cannot calibrate demand parameters!")}
    
    diag(diversion)=1
    if(any(diversion > 1 | diversion<0)){
      stop("'diversions' off-diagonal elements must be between 0 and 1")}
    
    if (!isTRUE(all.equal(rowSums(object@diversion,na.rm=TRUE),rep(0,nprods),check.names=FALSE,tolerance=1e-3)) &&
        any(rowSums(object@diversion,na.rm=TRUE)>0,na.rm=TRUE)){ stop("'diversions' rows cannot sum to greater than 0")}
    
    if(nprods != nrow(object@diversion) ||
       nprods != ncol(object@diversion)){
      stop("'diversions' must be a square matrix")
    }
    }

    
    return(TRUE)

  }

)


#'@rdname BertrandOther-Classes
#'@export
setClass(
  Class = "Linear",
  contains="Bertrand",
  representation=representation(
    intercepts       = "vector",
    prices           = "vector",
    quantities       = "numeric",
    margins          = "numeric",
    priceStart       = "numeric",
    symmetry         = "logical"
  ),
  prototype=prototype(
    intercepts    =  numeric(),
    symmetry      =  TRUE
  ),
  validity=function(object){



    nprods <- length(object@shares) # count the number of products
   

    if(nprods != length(object@quantities) ||
       nprods != length(object@margins) ||
       nprods != length(object@prices)){
      stop("'prices', 'quantities', 'margins', and 'shares' must all be vectors with the same length")}

    if(any(object@prices<0,na.rm=TRUE))             stop("'prices' values must be positive")
    if(any(object@quantities<0,na.rm=TRUE))          stop("'quantities' values must be positive")
    if(any(object@margins<0 | object@margins>1,na.rm=TRUE)) stop("'margins' values must be between 0 and 1")


    if(any(is.na(object@diversion))){stop("'diversions' matrix cannot contain NA")}
    
    if(nprods != length(object@priceStart)){
      stop("'priceStart' must have the same length as 'shares'")}


    if(!is.logical(object@symmetry) || length(object@symmetry)!=1){stop("'symmetry' must equal TRUE or FALSE")}

    if(!object@symmetry &&
       length(object@margins[!is.na(object@margins)])!= nprods){
      stop("When 'symmetry' is FALSE, all product margins must be supplied")
    }

    return(TRUE)

  }
)


#'@rdname BertrandOther-Classes
#'@export
setClass(
  Class = "LogLin",
  contains="Linear",
  prototype=prototype(
    symmetry=FALSE
  ),
  validity=function(object){




    nprods <- length(object@prices)
    if(any(is.na(object@margins))){
      stop("'margins' cannot contain NA values")
    }
    if(nprods != length(object@priceStart)){
      stop("'priceStart' must have the same length as 'prices'")}
  })


#'@rdname BertrandOther-Classes
#'@export
setClass(
  Class = "AIDS",
  contains="Linear",
  representation=representation(
    priceStart  = "numeric",
    priceDelta       = "numeric",
    mktElast         = "numeric",
    parmStart="numeric",
    insideSize          = "numeric"
  ),
  prototype=prototype(
    insideSize       =  NA_real_,
    priceDelta       =  numeric(),
    mktElast         =  numeric(),
    parmStart        =  numeric(),
    control.slopes = list(
    )
  ),


  validity=function(object){


    if(!length(object@parmStart) %in% c(0,2) || any(object@parmStart > 0,na.rm=TRUE)){stop("'parmStart' must be a length-2 non-positive numeric vector")}

    nprods <- length(object@shares)

    if(!isTRUE(all.equal(rowSums(object@diversion,na.rm=TRUE),rep(0,nprods),check.names=FALSE,tolerance=1e-3))){ stop("'diversions' rows must sum to 0")}

    if(!isTRUE(all.equal(sum(object@shares),1,check.names=FALSE,tolerance=1e-3))){
      stop("The sum of 'shares' values must equal 1")}

    nMargins <- length(object@margins[!is.na(object@margins)])



    if(nMargins<2 && isTRUE(is.na(object@mktElast))){stop("At least 2 elements of 'margins' must not be NA in order to calibrate demand parameters")}
    if(nMargins<1 && !isTRUE(is.na(object@mktElast))){stop("At least 1 element of 'margins' must not be NA in order to calibrate demand parameters")}



    return(NULL)

  }

)


#'@rdname BertrandOther-Classes
#'@export
setClass(
  Class = "PCAIDS",
  contains="AIDS",
  representation=representation(

    knownElast      = "numeric",
    knownElastIndex = "numeric"
  ),

  validity=function(object){




    nprods <- length(object@shares)

    if(length(object@knownElastIndex) != 1 ){stop("'knownElastIndex' must be length 1")}
    if(length(object@knownElast) != 1 ){stop("'knownElast' must be length 1")}
    if(length(object@mktElast) != 1 ){stop("'mktElast' must be length 1")}

    if(!(object@knownElastIndex %in% seq(1,nprods)) ){
      stop("'knownElastIndex' value must be between 1 and the length of 'shares'")}
    if(nprods != length(object@mcDelta)){
      stop("'mcDelta' must have the same length as 'shares'")}
    if(object@knownElast>0 || object@mktElast > 0 ){
      stop("'mktElast', 'knownElast' must be non-positive")}
    if(abs(object@knownElast) < abs(object@mktElast) ){
      stop("'mktElast' must be less  than 'knownElast' in absolute value")}
  }

)


#'@rdname BertrandOther-Classes
#'@export
setClass(
  Class = "PCAIDSNests",
  contains="PCAIDS",
  representation=
    representation(
      nests="factor",
      nestsParms="numeric"),

  validity=function(object){

    nprods <- length(object@shares)

    if(nprods != length(object@nests)){
      stop("'nests' length must equal the number of products")}


    ## Test to see if enough margin info has been supplied to identify all nesting parameters
    nNestParm <- nlevels(object@nests)
    nNestParm <- nNestParm*(nNestParm -1)/2 #calculate the number of nesting parameters
    nMargins  <- length(object@margins[!is.na(object@margins)])

    maxNests <- floor((sqrt(8 * nMargins + 1) + 1)/2) # compute the maximum number of nests that may be used
    # given available margins

    if(!is.vector(object@nestsParms) || nNestParm != length(object@nestsParms)){
      stop(paste("'nestsParmStart' must be a vector of length",nNestParm))}

    if(nNestParm > nMargins){
      stop(paste(
        "Impossible to calibrate nest parameters with the number of margins supplied.\n",
        "The maximum number of nests supported by the supplied margin information is"
        , maxNests,"."))
    }
  }

)
