setClass(
  Class = "Collusion",
  representation=representation(
    stage.function       = "function",
    stage.result         = "Logit"
  )
  
  
  ## compute margins
  setMethod(
    f= "calcMargins",
    signature= "Bertrand",
    definition=function(object,preMerger=TRUE){
      
      
      
      if( preMerger) {
        
        owner  <- object@ownerPre
        revenue<- calcShares(object,preMerger,revenue=TRUE)
        
        elast <-  elast(object,preMerger)
        margins <-  -1 * as.vector(solve(t(elast)*owner) %*% (revenue * diag(owner))) / revenue
        
        
      }
      
      else{
        prices <- object@pricePost
        mc     <- calcMC(object,preMerger)
        
        margins <- 1 - mc/prices
      }
      
      
      names(margins) <- object@labels
      
      return(as.vector(margins))
    }
    
  )
  
  
  
  
  ## Create a method to recover marginal cost using
  ## demand parameters and supplied prices
  setMethod(
    f= "calcMC",
    signature= "Bertrand",
    definition= function(object,preMerger=TRUE){
      
      object@pricePre <- object@prices
      
      
      marginPre <- calcMargins(object,TRUE)
      
      mc <- (1 - marginPre) * object@prices
      
      if(!preMerger){
        mc <- mc*(1+object@mcDelta)
      }
      
      names(mc) <- object@labels
      
      return(as.vector(mc))
    }
  )