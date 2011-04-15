## Generate a bunch of generic functions


#setGeneric (
# name= "summary",
# def=function(object){standardGeneric("summary")}
# )


setGeneric (
 name= "ownerToMatrix",
 def=function(object,...){standardGeneric("ownerToMatrix")}
 )

setGeneric (
 name= "calcSlopes",
 def=function(object){standardGeneric("calcSlopes")}
 )

setGeneric (
 name= "calcShares",
 def=function(object,...){standardGeneric("calcShares")}
 )


setGeneric (
 name= "calcMargins",
 def=function(object,...){standardGeneric("calcMargins")}
 )

setGeneric (
 name= "calcPrices",
 def=function(object,...){standardGeneric("calcPrices")}
 )

setGeneric (
 name= "elast",
 def=function(object,...){standardGeneric("elast")}
 )

setGeneric (
 name= "diversion",
 def=function(object,...){standardGeneric("diversion")}
 )

setGeneric (
 name= "CV",
 def=function(object,...){standardGeneric("CV")}
 )

setGeneric (
 name= "cmcr",
 def=function(object,...){standardGeneric("cmcr")}
 )


setGeneric (
 name= "getNestsParms",
 def=function(object,...){standardGeneric("getNestsParms")}
 )


## Create some methods for the Antitrust Class

setMethod(
 f= "ownerToMatrix",
signature= "Antitrust",
definition=function(object,preMerger=TRUE){


    ## transform ownerPre and ownerPost vectors into matrices, when applicable

    if(preMerger) thisOwner <- object@ownerPre
    else{         thisOwner <- object@ownerPost}



    if(is.vector(thisOwner) || is.factor(thisOwner)){

        nprod <- length(object@shares)
        owners <- as.numeric(factor(thisOwner))
        thisOwner <- matrix(0,ncol=nprod,nrow=nprod)


        for( o in unique(owners)){
            thisOwner [owners == o, owners == o] = 1
        }


    }


    return(thisOwner)

}
)

