sttasetClass(
  Class = "Collusion",
  representation=representation(
    stage         = "Logit"
  )
  
  setMethod(
    f= "calcCoalitions",
    signature= "Collusion",
    definition=function(object,preMerger=TRUE){
      
      locations <- unique(object@locations)
      
      
      numLocations <- length(locations)
      
      result <- do.call("cbind",
                        lapply(1:numLocations,
                               function(x){
                                 c=combn(locations,x)
                                 r=matrix(0,ncol=ncol(c),nrow=numLocations)
                                 r[1:nrow(c),]=c
                                 return(r)}
                        )
      )
      
      return(result)
      
    }
  )
  
  
  
  ## compute all coordinating price games
  setMethod(
    f= "calcCoordinationOutcomes",
    signature= "Collusion",
    definition=function(object,preMerger=TRUE){
  
      stage  <- object@stage
      nprods <- length(stage@labels)
      isParty <- rowSums( abs(stage@ownerPost - stage@ownerPre))>0
      
      if(preMerger){owner<-stage@ownerPre}
      else{owner<-stage@ownerPost}
      
      firms     <- unique(owner)
      prodOwner <- ownerToVec(stage)
       
      ps <- list()
      for( group in 2:nrow(firms)){
        thisGame <- stage
        theseCombn <- combn(1:nrow(firms),group)
        
        #thisPS <- matrix(ncol=ncol(theseCombn),nrow=nprods)
         thisPS <- vector("list",ncol(theseCombn))
         names(thisPS)<- apply(theseCombn,2,paste,collapse=",")
        
        for(c in 1:ncol(theseCombn)){
          
         thisOwner<- owner
         thisCoalition <-  prodOwner %in% theseCombn[,c]
         thisOwner[ thisCoalition, thisCoalition]=1
         thisGame@ownerPre <-  thisGame@ownerPost <- thisOwner
         thisGame@ownerPost[isParty,isParty]=1
         thisGame@pricePre <-   calcPrices(thisGame,preMerger=TRUE)
         thisGame@pricePost <-  calcPrices(thisGame,preMerger=FALSE)
         #thisPS[,c] <- calcProducerSurplus(thisGame)
          thisPS[[c]] <- thisGame
        }
        
        #ps <- cbind(ps,calcProducerSurplus(thisGame))
        ps <- c(ps,thisPS)
        
        
      }

      return(ps)
    }
      
  )
  
  
  
  