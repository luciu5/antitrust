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
        theseCoalitions <- combn(1:nrow(firms),group)
        
        #thisPS <- matrix(ncol=ncol(theseCombn),nrow=nprods)
         thisPS <- vector("list",ncol(theseCoalitions))
         names(thisPS)<- apply(theseCoalitions,2,paste,collapse=",")
        
        for(c in 1:ncol(theseCoalitions)){
          
         thisOwner<- owner
         thisCoalition <-  prodOwner %in% theseCoalitions[,c]
         thisOwner[ thisCoalition, thisCoalition]=1
         thisGame@ownerPre <-  thisGame@ownerPost <- thisOwner #All-C scenario for particular coalition
         thisGame@ownerPost[isParty,isParty]=1
         thisGame@pricePre <-   try(calcPrices(thisGame,preMerger=TRUE),silent=TRUE)
         thisGame@pricePost <-  try(calcPrices(thisGame,preMerger=FALSE),silent=TRUE)
         #thisPS[,c] <- calcProducerSurplus(thisGame)
          thisPS[[c]] <- thisGame
        }
        
        #ps <- cbind(ps,calcProducerSurplus(thisGame))
        ps <- c(ps,thisPS)
        
        
      }

      return(ps)
    }
      
  )
  
  
  
  