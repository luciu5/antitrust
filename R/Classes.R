require(Matrix)

setClassUnion("anyMatrix", c("matrix", "Matrix"))
setClassUnion("matrixOrVector", c("matrix", "Matrix", "vector","factor"))
setClassUnion("matrixOrList", c("matrix", "Matrix", "list"))


setClass(

 Class = "Antitrust",
 representation=representation(
 shares       = "vector",
 mcDelta      = "vector",
 slopes       = "matrixOrList",
 ownerPre     = "matrixOrVector",
 ownerPost    = "matrixOrVector",
 labels       = "vector"
 ),
         prototype=prototype(

         slopes          = matrix()

         ),
         validity=function(object){

             ## Sanity Checks

             if(any(object@shares < 0 | object@shares > 1,na.rm=TRUE)){
                 stop("'shares' values must be between 0 and 1")}

             nprods <- length(object@shares)

             ## if(any(object@mcDelta < 0 | object@mcDelta > 1,na.rm=TRUE)){
             ##    stop("'mcDelta' values must be between 0 and 1")}

             if(nprods != length(object@mcDelta)){
                 stop("'mcDelta' must have the same length as 'shares'")}

              if(nprods != length(object@labels)){
                 stop("'labels' must have the same length as 'shares'")}
             if(is.matrix(object@ownerPre) || is(object@ownerPre,"Matrix")){
                 if(nprods != ncol(object@ownerPre)){
                     stop("The number of rows and columns in 'ownerPre' must equal the length of 'shares'")}
                 if(nrow(object@ownerPre) != ncol(object@ownerPre)){
                     stop("'ownerPre' must be a square matrix ")}
             }
             else if (nprods != length(object@ownerPre)) stop("'ownerPre' and shares must be vectors of the same length")
             if(is.matrix(object@ownerPost)|| is(object@ownerPost,"Matrix")){
                 if(nprods != ncol(object@ownerPost)){
                     stop("The number of rows and columns in 'ownerPost' must equal the length of 'shares'")}
                 if(nrow(object@ownerPost) != ncol(object@ownerPost)){
                     stop("'ownerPost' must be a square matrix ")}
             }
             else if (nprods != length(object@ownerPost)) stop("'ownerPost' and shares must be vectors of the same length")
         }

             )

