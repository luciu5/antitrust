#' @title Methods For Calculating the Herfindahl-Hirschman Index
#' @name HHI-Methods
#' @docType methods

#' @aliases hhi
#' hhi,ANY-method
#' hhi,Bertrand-method
#' hhi,Cournot-method
#' hhi,VertBargBertLogit-method
#'
#' @description Computes the  Herfindahl-Hirschman Index (HHI) using simulated market
#' shares and either pre- or post-merger ownership information.
#' Outside shares are excluded from the calculation.
#'
#' @param object An instance of one of the classes listed above.
#' @param preMerger If TRUE, returns pre-merger outcome. If
#' FALSE, returns post-merger outcome.  Default is TRUE.
#' @param revenue If TRUE, returns revenues. If FALSE,
#' returns quantities. Default is TRUE.
#' @param insideonly If TRUE, excludes the share of the outside good from the calculation.
#' Default is TRUE.
#'
#' @include HypoMonMethods.R
#' @keywords methods
NULL

setGeneric (
  name= "hhi",
  def=function(object,...){standardGeneric("hhi")}
)


## Method to compute HHI
#'@rdname HHI-Methods
#'@export
setMethod(
  f= "hhi",
  signature= "Bertrand",
  definition=function(object,preMerger=TRUE,revenue=FALSE, insideonly=TRUE){

    if(preMerger){owner <- object@ownerPre}
    else{owner <- object@ownerPost}

    control <- owner>0              #assumes that a firm can set prices
    #on products over which it has partial ownership

    weights <- crossprod(control,owner)
    weights <- t(t(weights)/diag(weights)) # divide each element by its corresponding diagonal

    shares <- calcShares(object,preMerger,revenue)

    if(insideonly) shares <- shares/sum(shares, na.rm=TRUE) # hhi is typically defined over inside goods

    shares <- shares *100

    shares[is.na(shares)] <- 0

    result <- as.vector(shares %*% weights %*% shares)




    return(result)



  }
)


## Method to compute HHI
#'@rdname HHI-Methods
#'@export
setMethod(
  f= "hhi",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE, revenue = FALSE){

    shares <- calcShares(object,preMerger=preMerger,revenue=revenue)
    shares[is.na(shares)] <- 0
    if(preMerger) {owner <- object@ownerPre}
    else{owner <- object@ownerPre}

    hhi <- colSums((owner %*% (shares*100))^2, na.rm =TRUE)

    return(hhi)

  })




## Method to compute HHI
#'@rdname HHI-Methods
#'@export
setMethod(
  f= "hhi",
  signature= "VertBargBertLogit",
  definition=function(object,preMerger=TRUE,revenue=FALSE, insideonly=TRUE){
    
    up <- object@up
    down <- object@down
    
    ownerUp <- ownerToMatrix(up,preMerger=preMerger)
    ownerDown <- ownerToMatrix(down,preMerger=preMerger)
 
    
    if(object@isHorizontal){
      
      if(object@isUpstream) owner <- ownerUp
      else{owner <- ownerDown}
      

     
      
    }
    
    else{ 
      
      if(preMerger) {
        thisDownOwner<- down@ownerPre
        thisUpOwner <- up@ownerPre
      }
      else{ 
        thisDownOwner<- down@ownerPost
        thisUpOwner <- up@ownerPost
        }
      
      vertFirms <- intersect(thisUpOwner,thisDownOwner)

            
      
      for( v in vertFirms){
        
        isParty <- thisUpOwner == v | thisDownOwner== v
        
        ownerDown[isParty,isParty] <- 1
      }
      
      owner <- ownerDown
    }
    

    if(preMerger) down@ownerPre <- owner
    else{down@ownerPost <- owner}
    
    
    result <- hhi(down,preMerger=preMerger,revenue=revenue, insideonly=insideonly)
    
    
    return(result)
    

    
    
    
  }
    
  )
