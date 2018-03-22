setClassUnion("matrixOrVector", c("matrix", "numeric","character","factor"))
setClassUnion("matrixOrList", c("numeric","matrix", "list"))
setClassUnion("characterOrList", c("character", "list"))

setClass(

         Class = "Antitrust",
         representation=representation(
         ownerPre     = "matrixOrVector",
         ownerPost    = "matrixOrVector",
         pricePre     = "numeric",
         pricePost    = "numeric",
         mcPre        = "numeric",
         mcPost       = "numeric",
         labels       = "characterOrList",
         cls        = "character",
         control.slopes = "list",
         control.equ = "list"
         ),
         prototype(
         pricePre  = numeric(),
         pricePost = numeric(),
         mcPre     = numeric(),
         mcPost    = numeric(),
         cls  =    character(),
         ##copied from 'optim' definition:
         control.slopes = list( 
                               reltol = sqrt(.Machine$double.eps) 
                               ),
         ## copied from 'BBsolve' definition:
         ## changed 'maxit' to 2000 from 1500
         ## changed tol to sqrt(.Machine$double.eps) < 1e-07 (default)
         control.equ = list(maxit = 2000, M = c(50, 10), 
                            tol = sqrt(.Machine$double.eps), trace = FALSE, 
                            triter = 10, 
                            noimp = 100, NM = c(TRUE, FALSE))
         ),
         validity=function(object){



           if(is.list(object@labels)){ nprods <- length(object@labels[[1]])}
           else{nprods <- length(object@labels)}


             if(is.matrix(object@ownerPre)){

                 if(nprods != ncol(object@ownerPre)){
                     stop("The number of rows and columns in 'ownerPre' must equal the number of products")}
                 if(nrow(object@ownerPre) != ncol(object@ownerPre)){
                     stop("'ownerPre' must be a square matrix ")}

                 if(
                    any(colSums(unique(object@ownerPre),na.rm=TRUE)>1)
                    ){
                     stop("The columns of the matrix formed from the unique rows of 'ownerPre' must sum to no more than 1")
                     }
             }

             else if (nprods != length(object@ownerPre)) stop("'ownerPre' and 'labels' must be vectors of the same length")
             if(is.matrix(object@ownerPost)){
                 if(nprods != ncol(object@ownerPost)){
                     stop("The number of rows and columns in 'ownerPost' must equal the number of products")}
                 if(nrow(object@ownerPost) != ncol(object@ownerPost)){
                     stop("'ownerPost' must be a square matrix")}
                 if(
                    any(colSums(unique(object@ownerPost))>1,na.rm=TRUE)
                    ){
                     stop("The columns of the matrix formed from the unique rows of 'ownerPost' must sum to no more than 1")
                     }
             }

             else if (nprods != length(object@ownerPost)) stop("'ownerPost' and 'labels' must be vectors of the same length")

             if(identical(object@ownerPre,object@ownerPost)) warning("'ownerPost' and 'ownerPre' are the same")
              return(TRUE)
         }

         )



##
## Antitrust Methods
##

## Generate a bunch of generic functions


setGeneric (
            name= "ownerToMatrix",
            def=function(object,...){standardGeneric("ownerToMatrix")}
 )
setGeneric (
            name= "ownerToVec",
            def=function(object,...){standardGeneric("ownerToVec")}
            )

setGeneric (
 name= "calcPrices",
 def=function(object,...){standardGeneric("calcPrices")}
 )

setGeneric (
            name= "calcPriceDelta",
            def=function(object,...){standardGeneric("calcPriceDelta")}
            )


## print method
setMethod(
 f= "show",
 signature= "Antitrust",
 definition=function(object){

    try(res <- summary(object,market=TRUE), silent = TRUE)
  
  return(NULL)
}
          )



## Method to compute price changes
setMethod(
          f= "calcPriceDelta",
          signature= "Antitrust",
          definition=function(object,levels=FALSE){

              pricePre  <- object@pricePre
              pricePost <- object@pricePost

              if(levels){priceDelta <- pricePost - pricePre}
              else{priceDelta <- pricePost/pricePre - 1}
              #names(priceDelta) <- object@labels

              return(priceDelta)

          }
          )



## create ownership matrix
setMethod(
          f= "ownerToMatrix",
          signature= "Antitrust",
          definition=function(object,preMerger=TRUE){


              ## transform ownerPre/ownerPost vector into matrix, when applicable

              if(preMerger) {thisOwner <- object@ownerPre}
              else{         thisOwner <- object@ownerPost}



              if(is.vector(thisOwner) || is.factor(thisOwner)){

                  nprod <- length(object@labels)
                  owners <- as.numeric(factor(thisOwner))
                  thisOwner <- matrix(0,ncol=nprod,nrow=nprod)


                  for( o in unique(owners)){
                      thisOwner [owners == o, owners == o] = 1
                  }


              }


              return(thisOwner)

          }
          )


## convert ownership matrix to vector
setMethod(
          f= "ownerToVec",
          signature= "Antitrust",
          definition=function(object,preMerger=TRUE){


              ## transform ownerPre/ownerPost matrix into an ownership vector

              if(preMerger) {thisOwner <- object@ownerPre}
              else{         thisOwner <- object@ownerPost}

              if(is.matrix(thisOwner)){

                  thisOwner <- unique(thisOwner)
                  thisOwner <- as.numeric(thisOwner>=0.5) * (1: nrow(thisOwner))
                  thisOwner <- apply(thisOwner,2,max)

              }


    return(as.numeric(thisOwner))
          }

          )
