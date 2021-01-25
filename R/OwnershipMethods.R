#' @title Methods for Manipulating Ownership Matrices
#' @name Ownership-methods
#' @docType methods
#' @aliases ownerToMatrix
#' ownerToVec
#' @author Charles Taragin \email{ctaragin@ftc.gov}
#' @examples showMethods(classes="Antitrust") # show all methods defined for the class
#'
#' @param object An instance of the Antitrust class.
#' @param preMerger The \sQuote{preMerger} takes
#' on a value of TRUE or FALSE, where TRUE invokes the method using the
#' pre-merger values, while FALSE invokes the method using the post-merger ownership structure.
#'
#' @description
#' \code{ownerToMatrix} converts an ownership vector (or factor) to a k x k matrix of
#' 1s and 0s.
#'
#' \code{ownerToVec} converts a k x k  ownership matrix to a length-k
#' vector whose values identify an owner.
#'
#' @include CournotClasses.R
#' @keywords methods
NULL


setGeneric (
  name= "ownerToMatrix",
  def=function(object,...){standardGeneric("ownerToMatrix")}
)
setGeneric (
  name= "ownerToVec",
  def=function(object,...){standardGeneric("ownerToVec")}
)


## create ownership matrix
#'@rdname Ownership-methods
#'@export
setMethod(
  f= "ownerToMatrix",
  signature= "Antitrust",
  definition=function(object,preMerger=TRUE){


    ## transform ownerPre/ownerPost vector into matrix, when applicable

    if(preMerger) {thisOwner <- object@ownerPre}
    else{         thisOwner <- object@ownerPost}



    if(is.vector(thisOwner) || is.factor(thisOwner)){

      if(is.list(object@labels)){nprod <- length(object@labels[[1]])}
      else{nprod <- length(object@labels)}

      owners <- as.numeric(factor(thisOwner, levels= unique(thisOwner)))
      thisOwner <- matrix(0,ncol=nprod,nrow=nprod)


      for( o in unique(owners)){
        thisOwner [owners == o, owners == o] = 1
      }


    }


    return(thisOwner)

  }
)


## create ownership matrix
#'@rdname Ownership-methods
#'@export
setMethod(
  f= "ownerToMatrix",
  signature= "VertBargBertLogit",
  definition=function(object,preMerger=TRUE){
    
    up <- object@up
    down <- object@down
    
   
    
    
    if(preMerger) {thisUpOwner <- up@ownerPre
                   thisDownOwner <- down@ownerPre
                   bargParm <- up@bargpowerPre
                   
    }
    
    else{         thisUpOwner <- up@ownerPost
                  thisDownOwner <- down@ownerPost
                  bargParm <- up@bargpowerPost
                  }
    
    thisUpOwnerMat   <- ownerToMatrix(up,preMerger=preMerger)
    thisDownOwnerMat <- ownerToMatrix(down,preMerger=preMerger)
    
    nprods <- nrow(thisDownOwnerMat)
    
    if(!is.matrix(thisUpOwner) & !is.matrix(thisDownOwner)){
    
      
      
      vertFirms <- intersect(thisUpOwner,thisDownOwner)
        
    ## transform ownerPre/ownerPost vector into matrix, when applicable
    
    thisDownOwnerMatVertical <- matrix(0,nrow=nprods,ncol=nprods)
      
    for( v in vertFirms){
    
    #bargParm[thisUpOwner == v  & thisDownOwner == v] <- 1 
    
    vertrows <- thisUpOwner != v  & thisDownOwner == v
    thisUpOwnerMat[vertrows, thisUpOwner == v] <- -(1-bargParm[vertrows])/bargParm[vertrows]
    }
    
    
    ownerDownLambda <- thisDownOwnerMat * (1-bargParm)/bargParm
    
    for( v in vertFirms){
    
      
      vertrows <-  thisUpOwner == v  & thisDownOwner != v
    
    
    thisDownOwnerMatVertical[thisDownOwner == v, vertrows] <- 1
    #thisDownOwnerMatVertical[thisDownOwner == v, !vertrows] <- 0
    #thisDownOwnerMatVertical[thisDownOwner != v, ] <- 0
    
    
    ownerDownLambda[vertrows, thisDownOwner == v] <- -1
    
    }
    
    if(preMerger){
    object@ownerDownPre <- thisDownOwnerMatVertical
    object@ownerDownLambdaPre <-  ownerDownLambda
    object@ownerUpLambdaPre <- thisUpOwnerMat 
    #object@up@bargpowerPre <- bargParm
    }
    
    else{
      object@ownerDownPost <- thisDownOwnerMatVertical
      object@ownerDownLambdaPost <-  ownerDownLambda
      object@ownerUpLambdaPost <- thisUpOwnerMat
      #object@up@bargpowerPost <- bargParm
      
    }
    #vertical$ownerPostNoSupply.up <-  ids$up.firm == 1  & ids$down.firm != 1
    #vertical$ownerPostNoSupply.down <- ids$up.firm == 1  & ids$down.firm == 1
    
    }
    return(object)
    
    
  }
)

## convert ownership matrix to vector
#'@rdname Ownership-methods
#'@export
setMethod(
  f= "ownerToVec",
  signature= "Antitrust",
  definition=function(object,preMerger=TRUE){


    ## transform ownerPre/ownerPost matrix into an ownership vector

    if(preMerger) {thisOwner <- object@ownerPre}
    else{         thisOwner <- object@ownerPost}

    if(is.matrix(thisOwner)){

      thisOwner <- unique(thisOwner)
      mode(thisOwner) <- "numeric"
      thisOwner <- (thisOwner>=0.5) * (1: nrow(thisOwner))
      thisOwner <- apply(thisOwner,2,max)

    }


    return(as.numeric(thisOwner))
  }

)
