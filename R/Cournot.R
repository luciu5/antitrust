setClass(
  
  Class = "Cournot",
  contains="Bertrand",
  representation=representation(
    
    intercepts       = "numeric",
    mcfunPre           = "list",
    mcfunPost           = "list",
    vcfunPre           = "list",
    vcfunPost           = "list",
    prices           = "vector",
    quantities       = "matrix",
    margins          = "matrix",
    quantityPre      = "matrix",
    quantityPost     = "matrix",
    capacitiesPre           = "numeric",
    capacitiesPost           = "numeric",
    quantityStart     = "numeric",
    productsPre      = "matrix",
    productsPost     = "matrix",
    demand           = "character",
    cost             = "character",
    mktElast        =  "numeric"
  ), 
  prototype(
    intercepts    =  numeric(),
    mcfunPre =list(),
    mcfunPost=list(),
    vcfunPre =list(),
    vcfunPost=list(),
    mktElast = NA_real_
  ),
  validity=function(object){
    
    nplants <- nrow(object@quantities) # count the number of plants
    nprods  <- length(object@prices)     # count the number of products
  
    
    if(nplants != length(object@capacitiesPre)){
      stop("capacitiesPre' must be a vector whose length equals the number of rows in 'quantities'")}
    
    if( any(object@capacitiesPre<0 ,na.rm=TRUE)){stop("'capacitiesPre' values must be positive, and not NA")}
    
    
    
    if(nplants != length(object@capacitiesPost)){
      stop("capacitiesPre' must be a vector whose length equals the number of rows in 'quantities'")}
    
    if(any(is.na(object@capacitiesPost) |
           object@capacitiesPost<0 ,na.rm=TRUE)){stop("'capacitiesPost' values must be positive, and not NA")}
    
    if(
      !length(object@mcfunPre) %in% c(nplants,0)) {stop("'mcfunPre' must be a list of functions whose length equals the number of plants")}
    if(length(object@mcfunPre) >0 && any(sapply(object@mcfunPre,class) != "function"))
    {stop("'mcfunPre' must be a list of functions")}
    
    if(!identical(length(object@mcfunPre),length(object@vcfunPre) )){stop("'mcfunPre' and 'vcfunPre' should be lists of the same length")}
    
    
    if(
      !length(object@mcfunPost) %in% c(nplants,0)) {stop("'mcfunPost' must be a list of functions whose length equals the number of plants")}
    if(length(object@mcfunPost) >0 && any(sapply(object@mcfunPost,class) != "function"))
    {stop("'mcfunPost' must be a list of functions")}
    
    if( 
      !length(object@vcfunPre) %in% c(nplants,0)) {stop("'vcfunPre' must be a list of functions whose length equals the number of plants")}
    if(length(object@vcfunPre) >0 && any(sapply(object@vcfunPre,class) != "function"))
    {stop("'vcfunPre' must be a list of functions")}
    
    if(
      !length(object@vcfunPost) %in% c(nplants,0)) {stop("'vcfunPost' must be a list of functions whose length equals the number of plants")}
    if(length(object@vcfunPost) >0 && any(sapply(object@vcfunPost,class) != "function"))
    {stop("'vcfunPost' must be a list of functions")}
    
    if(!is.logical(object@productsPre)) stop("'productsPre' must be a logical matrix")
    if(!is.logical(object@productsPost)) stop("'productsPost' must be a logical matrix")
    
    if (!identical(dim(object@quantities), dim(object@margins))) stop("'margins' and 'quantities' must be matrices of the same dimension")
    if (!identical(dim(object@quantities), dim(object@productsPre))) stop("'productsPre' and 'quantities' must be matrices of the same dimension")
    if (!identical(dim(object@quantities), dim(object@productsPost))) stop("'productsPost' and 'quantities' must be matrices of the same dimension")
    
    
    if(!is.list(object@labels)) stop("'labels' must be a list") 
    if (isTRUE(nplants != length(object@labels[[1]]))) stop("'labels' length must be a list whose first element is a vector whose length equals the number of plants")
    
    if (isTRUE(nprods != length(object@labels[[2]]))) stop("'labels' length must be a list whose 2nd element is a vector whose length equals the number of productsS")
    
    
    if(any(object@prices<=0,na.rm=TRUE))             stop("'prices' values must be positive")
    if(any(object@quantities<=0,na.rm=TRUE))          stop("'quantities' values must be positive")
    if(any(object@margins<0 | object@margins>1,na.rm=TRUE)) stop("'margins' values must be between 0 and 1")
    if(any(colSums(object@margins,na.rm=TRUE) == 0  & is.na(object@mktElast))) stop("at least one plant margin or a market elasticity must be supplied for each product")
    
    if(all(!(object@demand %in% c("linear","log")))){stop("'demand' must equal 'linear' or 'log'")}
    if(length(object@demand) != nprods) stop("the length of 'demand' must equal the number of products")
    
    if(all(!(object@cost %in% c("linear","constant")))){stop("'cost' must equal 'linear' or 'constant'")}
    if(length(object@cost) != nplants) stop("the length of 'cost' must equal the number of products")
    
    if(any(rowSums(object@ownerPre) > 0 & !(is.finite(object@capacitiesPre) & is.finite(object@capacitiesPost))  & object@cost == "constant")){ stop("multi-plant firms with 'cost' equal to 'constant' must have finite pre- and post-merger capacity constraints")}
    
    if(length(object@prices) != nprods) stop("the length of 'prices' must equal the number of products")
    
    if(ncol(object@quantities) != nprods) stop("the number of columns in 'quantities' must equal the number of products")
    
    if(any(rowSums(object@quantities,na.rm=TRUE) > object@capacitiesPre)){
      stop("pre-merger plant output must be less than pre-merger capacity constraint")
    }
    
    
    
  }
  
)

setGeneric (
  name= "calcVC",
  def=function(object,...){standardGeneric("calcVC")}
)


##
## Cournot Methods
##


setMethod(
  f= "calcShares",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE,revenue=FALSE){
    
    if(preMerger) quantities <- object@quantityPre
    else{ quantities <- object@quantityPost}
    
    if (revenue){
      if(preMerger){ prices <- object@pricePre}
      else{          prices <- object@pricePost}
      
      totrev <- rowSums(prices*t(quantities), na.rm = TRUE)
      return(t(prices*t(quantities)/totrev))
    }
    
    else{
      totquant <- colSums(quantities,na.rm=TRUE)
      return(t(t(quantities)/totquant))}
  }
)

setMethod(
  f= "elast",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE,market=FALSE){
    
    isLinearD <- object@demand =="linear"
    slopes <- object@slopes
    intercepts <- object@intercepts
    
    
    if(preMerger){
      quantities <- object@quantityPre
      owner <- object@ownerPre
      }
    else{
      quantities <- object@quantityPost
      owner <- object@ownerPost}
    
    quantities[is.na(quantities)] <- 0
    
    quantOwner <- owner %*% quantities
    
    prices <- calcPrices(object,preMerger=preMerger)
    
    mktQuant <-  colSums(quantities,na.rm = TRUE)
    
    ##dPdQ
    partial <- ifelse(isLinearD, 
                      slopes,
                      exp(intercepts)*slopes*mktQuant^(slopes - 1))
    
    ##dQdP
    partial <- 1/partial
    
    
    
      
      elast <- partial*prices/mktQuant
    
    
    if(!market){
      
      sharesOwner <-  t(quantOwner) / mktQuant
      
      
      elast <-  t(elast / sharesOwner)
      
    
      
      dimnames(elast) <- object@labels
    }
    
    return(elast)
    
  }
)

## compute margins
setMethod(
  f= "calcMargins",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){
    
    
    if(preMerger){
      prices <- object@pricePre
    }  
    else{prices <- object@pricePost}
    
    mc <- calcMC(object,preMerger = preMerger)
    prices <- matrix(prices, ncol=length(prices), nrow=length(mc),byrow=TRUE)
    
    
    
    margin <- 1 - mc/prices
    
    
    dimnames(margin) <- object@labels
    return(margin)
  }
  
)


## Create a method to recover marginal cost using
## demand parameters and supplied prices


setMethod(
  f= "calcSlopes",
  signature= "Cournot",
  definition=function(object){
    
    prices <- object@prices
    quantities <- object@quantities
    quantities[is.na(quantities)] <- 0
    margins <- object@margins
    mktElast <- object@mktElast
    
    
    cap <- object@capacitiesPre
    
    mc <- t(t(1 - margins) * prices)
    
    products <- object@productsPre
    demand <- object@demand
    owner <- object@ownerPre
    mcfunPre <- object@mcfunPre
    nprods <- ncol(quantities)
    nplants <- nrow(quantities)
    
    noCosts <- length(mcfunPre) == 0
    isLinearD <- demand=="linear"
    isLinearC <- object@cost=="linear"
    
    quantTot <- colSums(quantities, na.rm = TRUE)
    quantPlants <- rowSums(quantities, na.rm = TRUE)
    quantOwner <- owner %*% quantities
    
    isConstrained <- quantPlants >= cap
    
    if(!noCosts){
      mcPre <- sapply(1:nplants, function(i){object@mcfunPre[[i]](quantities[i,])})
    }
    
    sharesOwner <- t(t(quantOwner)/quantTot)
    
    minDemand <- function(theta){
      
      if(noCosts){
        
        thiscap <- theta[1:nplants]
        theta <- theta[-(1:nplants)]
        mcPre <- ifelse(isLinearC, quantPlants/thiscap, thiscap)
        
      }
      
      thisints <- theta[1:nprods]
      thisslopes <- theta[-(1:nprods)]
      
      thisprices <- ifelse(isLinearD, thisints + thisslopes*quantTot, 
                           exp(thisints)*quantTot^thisslopes)
      
      thisPartial <- ifelse(isLinearD, 
                            thisslopes,
                            exp(thisints)*thisslopes*quantTot^(thisslopes - 1))
      
      
      thisFOC <- (t(quantities) * thisPartial ) %*% owner  + thisprices 
      thisFOC <- t(thisFOC) -   mcPre  
      thisFOC <- thisFOC[!isConstrained,]
      
      dist <- c(thisFOC,thisprices - prices, mktElast * quantTot/prices - 1/thisPartial)
      
      if(noCosts){ dist <- c(dist, mcPre - mc)}
      
      return(sum(dist^2,na.rm=TRUE))
    }
    
    
    margGuess <- margins
    margGuess[is.na(margGuess)] <- -t(t(sharesOwner)/mktElast)[is.na(margGuess)]
    
    bStart      =   ifelse(isLinearD,
                           colMeans(-(prices*margGuess)/(sharesOwner*quantTot),na.rm=TRUE), 
                           colMeans(-margGuess/sharesOwner,na.rm=TRUE))
    intStart    =   ifelse(isLinearD,
                           prices - bStart*quantTot, 
                           log(prices/(quantTot^bStart)))
    intStart    =   abs(intStart)
    
    parmStart   =   c( intStart,bStart)
    
    if(noCosts){
      
      if(isLinearD){margStart <- rowMeans(-(sharesOwner*quantTot)/(prices/bStart),na.rm=TRUE)}
      else{margStart <- rowMeans(-sharesOwner*bStart,na.rm=TRUE)} 
      
      mcStart  <- abs(prices*(margStart - 1)) 
      capStart <- ifelse(isLinearC, quantPlants/mcStart, mcStart)  
      parmStart <- c(capStart,parmStart)
    }
    
    
    bestParms=optim(parmStart,minDemand)$par
    
    
    
    
    ## if no marginal cost functions are supplied 
    ## assume that plant i's marginal cost is
    ## q_i/k_i, where k_i is calculated from FOC
    
    if(noCosts){
      
      mcparm <- bestParms[1:nplants]
      bestParms <- bestParms[-(1:nplants)]
      
      
      mcdef <- ifelse(isLinearC,"function(q,mcparm = %f){ val <- sum(q, na.rm=TRUE) / mcparm; return(val)}",
                      "function(q,mcparm = %f){ val <- mcparm; return(val)}")
      mcdef <- sprintf(mcdef,mcparm)
      mcdef <- lapply(mcdef, function(x){eval(parse(text=x ))})
      
      object@mcfunPre <- mcdef
      names(object@mcfunPre) <- object@labels[[1]]
      
      vcdef <- ifelse(isLinearC,"function(q,mcparm = %f){  val <-  sum(q, na.rm=TRUE)^2 / (mcparm * 2); return(val)}",
                      "function(q,mcparm = %f){  val <-  sum(q, na.rm=TRUE) * mcparm; return(val)}")
      vcdef <- sprintf(vcdef,mcparm)
      vcdef <- lapply(vcdef, function(x){eval(parse(text=x ))})
      
      object@vcfunPre <- vcdef
      names(object@vcfunPre) <- object@labels[[1]]
      
    }
    if(length(object@mcfunPost)==0){
      object@mcfunPost <- object@mcfunPre
      object@vcfunPost <- object@vcfunPre}
    
    intercepts = bestParms[1:nprods]
    slopes = bestParms[-(1:nprods)]
    
    
    object@intercepts <- intercepts
    object@slopes <-     slopes
    
    
    return(object)
    
  })    


setMethod(
  f= "calcMC",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){
    
    if(preMerger){ 
      
      quantity  <- object@quantityPre
      mcfun <- object@mcfunPre
      cap <- object@capacitiesPre
    }
    else{          
      
      quantity  <- object@quantityPost
      mcfun <- object@mcfunPost
      cap <- object@capacitiesPost
    }
    
    plantQuant <- rowSums(quantity,na.rm=TRUE)
    
    
    
    nplants <- nrow(quantity)
    
    mc <- rep(NA, nplants)
    
    for(f in 1:nplants){
      mc[f] <- mcfun[[f]](quantity[f,])
    }
    
    if(!preMerger){mc <- mc*(1 + object@mcDelta)}
    
    mc <- mc + 1/(100*(pmax(cap - plantQuant, 1e-16))) + 1/(100*(pmax(1e-16,plantQuant)))
    #mc <- ifelse(plantQuant <= cap & plantQuant >= 0 , mc, max(mc,na.rm=TRUE) * 1e3)
    
    
    names(mc) <- object@labels[[1]]
    
    return(mc)
  })    


setMethod(
  f= "calcVC",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){
    
    if(preMerger){ 
      
      quantity  <- object@quantityPre
      vcfun <- object@vcfunPre
    }
    else{          
      
      quantity  <- object@quantityPost
      vcfun <- object@vcfunPost
    }
    
    
    
    
    nplants <- nrow(quantity)
    
    vc <- rep(NA, nplants)
    
    for(f in 1:nplants){
      vc[f] <- vcfun[[f]](quantity[f,])
    }
    
    if(!preMerger){vc <- vc*(1 + object@mcDelta)}
    
    names(vc) <- object@labels[[1]]
    
    return(vc)
  })   





## compute producer surplus
setMethod(
  f= "calcProducerSurplus",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){
    
    
    if( preMerger) {
      prices <- object@pricePre
      quantities <- object@quantityPre
      
    }
    else{prices <- object@pricePost
    quantities <- object@quantityPost
    
    }
    
    
    
    vc <- calcVC(object, preMerger= preMerger)
    
    ps <- colSums(prices*t(quantities), na.rm=TRUE) - vc
    names(ps) <- object@labels[[1]]
    
    return(ps)
  }
  
)


setMethod(
  f= "calcQuantities",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE,...){
    
    slopes <- object@slopes
    intercepts <- object@intercepts
    quantityStart <- object@quantityStart
    
    if(preMerger){ owner  <- object@ownerPre
    products  <- object@productsPre
    cap <- object@capacitiesPre
    }
    else{          owner <-  object@ownerPost
    products <-  object@productsPost
    cap <- object@capacitiesPost
    }
    
    nprods <- ncol(products)
    isProducts <- rowSums(products) > 0
    products <- as.vector(products)
    
    isConstrained <- is.finite(cap)
    
    FOC <- function(quantCand){
      
      #quantCand <- quantCand^2 # constrain positive
      
      allquant <- rep(0,length(products))
      allquant[products] <- quantCand
      quantCand <- matrix(allquant,ncol=nprods)
      
      if(preMerger){ object@quantityPre  <- quantCand}
      else{          object@quantityPost <- quantCand}
      
      thisPrice <- calcPrices(object, preMerger= preMerger)
      
      thisMC <- calcMC(object, preMerger= preMerger) 
      
      
      mktQuant <- colSums(quantCand, na.rm=TRUE)
      plantQuant <- rowSums(quantCand, na.rm=TRUE)
      
      thisPartial <- ifelse(object@demand=="linear", 
                            slopes,
                            exp(intercepts)*slopes*mktQuant^(slopes - 1))
      
      
      thisFOC <- (t(quantCand) * thisPartial) %*% owner + thisPrice
      thisFOC <- t(thisFOC) - thisMC
      thisFOC <- t(t(thisFOC)/thisPrice) # rescale
      #thisCons <- (plantQuant - cap)/cap # rescale
      #thisFOC[isConstrained,] <- thisFOC[isConstrained,] + 
      #   thisCons[isConstrained] + 
      #  sqrt(thisFOC[isConstrained,]^2 + 
      #         thisCons[isConstrained]^2)
        thisFOC <- thisFOC[isProducts,]
      return(as.vector(thisFOC))
    }
    
    
    #quantityStart <- sqrt(object@quantityStart[products]) #constrain positive
    quantityStart <- ifelse(quantityStart >= cap, cap-1, quantityStart)
    quantityStart <- quantityStart[products]
    
    
    ## Find price changes that set FOCs equal to 0
    minResult <- BBsolve( quantityStart,FOC, quiet=TRUE,control=object@control.equ,...)
    
    if(minResult$convergence != 0){warning("'calcQuantities' nonlinear solver may not have successfully converged. 'BBsolve' reports: '",minResult$message,"'")}
    
    quantEst        <- rep(0, length(products))
    quantEst[products] <- minResult$par#^2
    quantEst <- matrix(quantEst,ncol = nprods)
    
    dimnames(quantEst) <- object@labels
    
    return(quantEst)
  })


setMethod(
  f= "calcPrices",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){
    
    if(preMerger){
      
      quantities <- object@quantityPre
    }
    else{
      quantities <- object@quantityPost
    }
    
    intercepts <- object@intercepts
    slopes     <- object@slopes
    
    
    mktQuant <- colSums(quantities, na.rm=TRUE)
    
    prices <- ifelse(object@demand == "linear",
                     intercepts + slopes * mktQuant,
                     exp(intercepts) * mktQuant^slopes
    )
    
    
    names(prices) <- object@labels[[2]]
    return(prices)
    
  })  



## Method to compute HHI
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


setMethod(
  f= "cmcr",
  signature= "Cournot",
  definition=function(object){
    
    owner <- object@ownerPre
    isParty <- rowSums( abs(object@ownerPost - object@ownerPre) ) > 0
    
    shares <- calcShares(object,preMerger=TRUE,revenue=FALSE)
    shares[is.na(shares)] <- 0 
    shares <- owner %*% shares
    if(any(isParty)) shares <- unique(shares[isParty,,drop=FALSE])
    if(nrow(shares) == 1 ){shares <- rbind(shares,shares)}
    mktElast <- elast(object, preMerger= TRUE,market=TRUE)
    
    cmcr <- -2 * apply(shares,2,prod) / (mktElast * colSums(shares) )
    
    return(cmcr * 100)
  })


setMethod(
  f= "summary",
  signature= "Cournot",
  definition=function(object,market=FALSE,revenue=FALSE,shares=FALSE,levels=FALSE,parameters=FALSE,digits=2,...){
    
    if(market){nplants <- 1}
    else{ nplants <- nrow(object@quantities) }
    
    curWidth <-  getOption("width")
    curSci  <-  getOption("scipen")
    
    pricePre   <-  object@pricePre
    pricePost  <-  object@pricePost
    priceDelta <- calcPriceDelta(object,levels=levels)
    if(!levels) priceDelta <- priceDelta *100
    
    if(!shares){
      outPre  <-  object@quantityPre
      outPost <-  object@quantityPost
      sumlabels=paste("quantity",c("Pre","Post"),sep="")
      
      if(revenue){
        outPre <- t(pricePre*t(outPre))
        outPost <- t(pricePost*t(outPost))
        sumlabels=paste("revenue",c("Pre","Post"),sep="")
      }
      
    }
    
    else{
      if(!shares){warning("'shares' equals FALSE but 'calcQuantities' not defined. Reporting shares instead of quantities")}
      
      
      outPre  <-  calcShares(object,preMerger=TRUE,revenue=revenue) * 100
      outPost <-  calcShares(object,preMerger=FALSE,revenue=revenue) * 100
      
      
      sumlabels=paste("shares",c("Pre","Post"),sep="")
    }
    
    if(market){
      
      outPre <- colSums(outPre,na.rm=TRUE)
      outPost <- colSums(outPost,na.rm=TRUE)
      ids <- data.frame(plant = 1 ,product= object@labels[[2]])  
    }
    
    
    else{
      
      ids <- expand.grid(plant=object@labels[[1]], product=object@labels[[2]])
    }
    
    
    out <- data.frame(product=ids$product, 
                      plant=ids$plant,outPre=as.vector(t(outPre)), 
                      outPost = as.vector(t(outPost)))
    
    if(market) {out$plant <- NULL}
    else{
      out$isParty <- as.numeric(rowSums( abs(object@ownerPost - object@ownerPre))>0)
      out$isParty <- factor(out$isParty,levels=0:1,labels=c(" ","*"))
    }
    
    mcDelta <- object@mcDelta * 100
    
    if(levels){out$outDelta <- out$outPost - out$outPre}
    else{out$outDelta <- (out$outPost/out$outPre - 1) * 100}
    
    out$pricePre <- rep(pricePre,each=nplants)
    out$pricePost <- rep(pricePost,each=nplants)
    out$priceDelta <- rep(priceDelta, each=nplants)
    
    if(market){
      results <- out[,c("product","pricePre","pricePost","priceDelta","outPre","outPost","outDelta" )]  
    }
    
    else{
      results <- out[, c("isParty","product","plant", "pricePre","pricePost","priceDelta","outPre","outPost","outDelta" )]
    }
    
    colnames(results)[colnames(results) %in% c("outPre","outPost")] <- sumlabels
    
    if(!market && sum(abs(mcDelta))>0) results <- cbind(results,mcDelta=mcDelta)
    
    
    
    sharesPost <- calcShares(object,FALSE,revenue)
    
    cat("\nMerger simulation results under '",class(object),"' demand:\n\n",sep="")
    
    options("width"=100) # this width ensures that everything gets printed on the same line
    options("scipen"=999) # this width ensures that everything gets printed on the same line
    print(format(results,digits=digits),row.names = FALSE)
    options("width"=curWidth) #restore to current width
    options("scipen"=curSci) #restore to current scientific notation
    
    cat("\n\tNotes: '*' indicates merging parties' products.\n ")
    if(levels){cat("\tDeltas are level changes.\n")}
    else{cat("\tDeltas are percent changes.\n")}
    if(revenue){cat("\tOutput is based on revenues.\n")}
    else{cat("\tOutput is based on units sold.\n")}
    
    
    
    ##Only compute cmcr if cmcr method doesn't yield an error
    thisCMCR <- tryCatch(cmcr(object),error=function(e) FALSE)
    if(!is.logical(thisCMCR)){
      cat("\n\nCMCR:\n\n")
      
      cat(format(cmcr(object),digits=digits), fill=TRUE,labels=object@labels[[2]])
    } 
    
    ##Only compute upp if prices are supplied
    #thisUPP <- tryCatch(upp(object),error=function(e) FALSE)
    #if(!is.logical(thisUPP)){
    #  cat("\nShare-Weighted Pricing Pressure:",format(sum(thisUPP*sharesPost[isParty=="*"],na.rm=TRUE)/sum(sharesPost[isParty=="*"],na.rm=TRUE),digits),sep="\t")}
    
    ##Only compute CV if prices  are supplied
    thisCV <- tryCatch(CV(object,...),error=function(e) FALSE)
    if(!is.logical(thisCV)){
      cat("\n\nCompensating Variation (CV):\n\n")
      cat(format(thisCV,digits=digits),fill=TRUE, labels=object@labels[[2]])}
    
    cat("\n\n")
    
    
    if(parameters){
      
      cat("\nDemand Parameter Estimates:\n\n")
      
      print(format(object@slopes,digits=digits))
      
      cat("\n\n")
      
      if(.hasSlot(object,"intercepts")){
        
        cat("\nIntercepts:\n\n")
        print(format(object@intercepts,digits=digits))
        cat("\n\n")
        
      }
      
      
    }
    
    return(invisible(results))
    
  })


setMethod(
  f= "CV",
  signature= "Cournot",
  definition=function(object){
    
    demand <- object@demand
    slopes    <- object@slopes
    intercepts <- object@intercepts
    quantityPre <- colSums(object@quantityPre, na.rm=TRUE)
    quantityPost <- colSums(object@quantityPost, na.rm=TRUE)
    pricePre <- object@pricePre
    pricePost <- object@pricePost
    
    result <- ifelse(demand =="linear",
                     .5*(pricePost - pricePre)*(quantityPre - quantityPost)  ,
                     exp(intercepts)/(slopes + 1) * (quantityPre^(slopes + 1) - quantityPost^(slopes+1)) -  (quantityPre - quantityPost)* pricePre
    )
    
    result <- result  +  (pricePost - pricePre)*quantityPost
    names(result) <-  object@labels[[2]]
    return(result)
  })

setMethod(
  f= "calcPricesHypoMon",
  signature= "Cournot",
  definition=function(object,plantIndex){
    
    nhypoplants <- length(plantIndex)
    intercept <- object@intercepts
    nprods <- length(intercept)
    slopes <- object@slopes
    quantityPre <- as.vector(object@quantityPre)
    demand <- object@demand
    
    
    ## how to deal with multiple products?
    stop("A work in progress!!")
    
    calcMonopolySurplus <- function(quantCand){
      
      
      quantityPre[plantIndex] <- quantCand
      quantCand <- matrix(quantCand,ncol=nprods)
      object@quantityPre <- quantCand
      mktQuant <- colSums(quantCand, na.rm = TRUE)
      
      priceCand <- ifelse(demand == "linear",
                          intercept + slopes * mktQuant,
                          exp(intercept)*mktQuant^slopes)
      
      vcCand <- calcVC(object, preMerger=TRUE)
      vcCand <- vcCand[plantIndex]
      
      revCand <-  colSums(priceCand*t(quantCand[plantIndex,]), na.rm=TRUE)
      
      
      surplus <- sum(revCand - vcCand, na.rm =TRUE)
      
      return(sum(surplus))
    }
    
    if( nhypoplants > 1){
      
      maxResult <- optim(object@quantityPre[plantIndex],
                         calcMonopolySurplus,
                         method="L-BFGS-B",
                         lower = rep(0,nhypoplants)
      )
      
      quantitiesHM <- maxResult$par
    }
    
    
    else{
      
      upperB <- sum(quantityPre,na.rm=TRUE)
      maxResult <- optimize(calcMonopolySurplus,c(0, upperB),maximum = TRUE)
      quantitiesHM <- maxResult$maximum
    }
    
    quantityPre[plantIndex]
    names(pricesHM) <- object@labels[plantIndex]
    
    return(pricesHM)
    
    
  })



cournot <- function(prices,quantities,
                    margins = matrix(NA_real_ , nrow(quantities),ncol(quantities)), 
                    demand = rep("linear",length(prices)),
                    cost   =   rep("linear",nrow(quantities)),
                    mcfunPre=list(),
                    mcfunPost=mcfunPre,
                    vcfunPre=list(),
                    vcfunPost=vcfunPre,
                    capacitiesPre = rep(Inf,nrow(quantities)),
                    capacitiesPost = capacitiesPre,
                    productsPre=!is.na(quantities), 
                    productsPost=productsPre, 
                    ownerPre,ownerPost,
                    mktElast = rep(NA_real_, length(prices)),
                    mcDelta =rep(0,nrow(quantities)),
                    quantityStart=as.vector(quantities),
                    control.slopes,
                    control.equ,
                    labels,
                    ...
){
  
  shares <- as.vector(quantities/sum(quantities))
  
  
  
  if(missing(labels)){
    if(is.null(dimnames(quantities))){ 
      rname <- paste0("O",1:nrow(quantities))
      cname <- paste0("P",1:ncol(quantities))
    }
    else{rname <- rownames(quantities)
    cname <- colnames(quantities)
    }
    labels <- list(rname,cname)
    
  }
  
  result <- new("Cournot",prices=prices, quantities=quantities,margins=margins,
                shares=shares,mcDelta=mcDelta, subset= rep(TRUE,length(shares)), demand = demand, cost=cost,
                mcfunPre=mcfunPre, mcfunPost=mcfunPost,vcfunPre=vcfunPre, vcfunPost=vcfunPost,
                capacitiesPre=capacitiesPre,capacitiesPost=capacitiesPost,
                ownerPre=ownerPre, mktElast = mktElast,productsPre=productsPre,productsPost=productsPost,
                ownerPost=ownerPost, quantityStart=quantityStart,labels=labels)
  
  
  if(!missing(control.slopes)){
    result@control.slopes <- control.slopes
  }
  
  
  ## Convert ownership vectors to ownership matrices
  result@ownerPre  <- ownerToMatrix(result,TRUE)
  result@ownerPost <- ownerToMatrix(result,FALSE)
  
  ## Calculate Demand Slope Coefficients and Intercepts
  result <- calcSlopes(result)
  
  
  result@quantityPre  <- calcQuantities(result, preMerger = TRUE,...)
  result@quantityPost <- calcQuantities(result,preMerger = FALSE,...)
  
  result@pricePre  <- calcPrices(result, preMerger = TRUE)
  result@pricePost <- calcPrices(result,preMerger = FALSE)
  
  return(result)
  
}

