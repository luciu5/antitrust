setClass(
  Class = "Collusion",
  representation=representation(
    stage         = "Logit"
  )
  
  
  ## compute all coordinating price games
  setMethod(
    f= "calcCoordinatingEquilibria",
    signature= "Collusion",
    definition=function(object,preMerger=TRUE){
  
      stage  <- object@stage
      nprods <- length(object@labels)
      
      if(preMerger){owner<-stage@ownerPre}
      else{owner<-stage@ownerPost}
      
      firms     <- unique(owner)
      prodOwner <- ownerToVec(object)
       
      ps <- NULL
      for( group in 2:nrow(firms)){
        thisGame <- stage
        theseCombn <- combn(1:nrow(firms),group)
        
        thisPS <- matrix(ncol=ncol(theseCombn),nrow=nprods)
        
        for(c in 1:ncol(theseCombn)){
          
         thisOwner<- owner
         thisOwner[prodOwner %in% theseCombn[c,], prodOwner %in% theseCombn[c,]] <- 1
         thisGame@ownerPre <-  thisOwner
         thisPS[,c] <- calcProducerSurplus(thisGame)
        
        }
        
        ps <- cbind(ps,calcProducerSurplus(thisGame))
        
        
        
      }

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