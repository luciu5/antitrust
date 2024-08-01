#' @title Methods for Calculating Margins
#' @name Margins-Methods
#' @docType methods
#'
#' @aliases calcMargins
#' calcMargins,ANY-method
#' calcMargins,AIDS-method
#' calcMargins,Bertrand-method
#' calcMargins,LogitCournot-method
#' calcMargins,VertBargBertLogit-method
#' calcMargins,VertBargBertLogitNests-method
#' calcMargins,LogitCap-method
#' calcMargins,Auction2ndLogit-method
#' calcMargins,Auction2ndLogitNests-method
#' calcMargins,Cournot-method
#' calcMargins,BargainingLogit-method
#' calcMargins,Bargaining2ndLogit-method
#'
#' @description Computes equilibrium product margins assuming that firms are playing a
#' Nash-Bertrand, Cournot, 2nd Score Auction, or Bargaining game. For "LogitCap", assumes firms are
#' playing a Nash-Bertrand or Cournot game with capacity constraints.
#'
#' @param object An instance of one of the classes listed above.
#' @param preMerger If TRUE, returns pre-merger outcome. If
#' FALSE, returns post-merger outcome.  Default is TRUE.
#' @param level IF TRUE, return margins in dollars. If FALSE, returns 
#' margins in proportions. Default for most classes is FALSE.
#' @param exAnte If \sQuote{exAnte} equals TRUE then the
#' \emph{ex ante} expected result for each firm is produced, while FALSE produces the
#' expected result conditional on each firm winning the auction. Default is FALSE.
#'
#' @include CostMethods.R
#' @keywords methods
NULL

setGeneric (
  name= "calcMargins",
  def=function(object,...){standardGeneric("calcMargins")}
)

## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Bertrand",
  definition=function(object,preMerger=TRUE, level=FALSE){

    output <- ifelse(object@output,-1,1)


    if( preMerger) {

      prices <- object@pricePre
      owner  <- object@ownerPre
      
    }

    else{
      prices <- object@pricePost
      owner  <- object@ownerPost

    }

    revenue<- calcShares(object,preMerger,revenue=TRUE)
    
    elast <-  elast(object,preMerger)
    
    margins <-  try(output * as.vector(solve(t(elast)*owner) %*% (revenue * diag(owner))) / revenue,silent=TRUE)
    if(any(class(margins) == "try-error")){margins <- output * as.vector(MASS::ginv(t(elast)*owner) %*% (revenue * diag(owner))) / revenue}
    
    
    
    if(level) {margins <- margins * prices }
    
    names(margins) <- object@labels

    return(as.vector(margins))
  }

)


## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "LogitCournot",
  definition=function(object,preMerger=TRUE, level=FALSE){
    
    
    output <- ifelse(object@output,-1,1)
    alpha <- object@slopes$alpha
    idx   <-  object@normIndex
    
    if( preMerger) {
      
      prices <- object@pricePre
      owner  <- object@ownerPre
      
    }
    
    else{
      prices <- object@pricePost
      owner  <- object@ownerPost
      
    }
    
    shares <- calcShares(object,preMerger=preMerger,revenue=FALSE)
    
    if(is.na(idx)){
      idxShare <- 1 - sum(shares)
    }
    else{
      idxShare <- shares[idx]
     
    }
    sharesFirm <- as.numeric(owner %*% shares) 
    
    
    margins <- output*(1 + sharesFirm/idxShare)/alpha
    
    if(!level) {margins <- margins/prices }
    
    names(margins) <- object@labels
    
    return(as.vector(margins))
  }
  
)


## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Bargaining2ndLogit",
  definition=function(object,preMerger=TRUE,exAnte=FALSE,level=TRUE){
  
    
    
    if( preMerger) {
      
      barg <- object@bargpowerPre
      
      
    }
    
    else{
      
      barg <- object@bargpowerPost
      
    }
      
    (1-barg)*callNextMethod()
    
  })

## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "BargainingLogit",
  definition=function(object,preMerger=TRUE, level=FALSE){
    
    output <- ifelse(object@output,-1,1)
    
    alpha <- object@slopes$alpha
    
    if( preMerger) {
      
      prices <- object@pricePre
      owner  <- object@ownerPre
      barg <- object@bargpowerPre
      
      
    }
    
    else{
      prices <- object@pricePost
      owner  <- object@ownerPost
      barg <- object@bargpowerPost
      
    }
    
    
    barg <- barg/(1-barg) #relative bargaining
    
    nprods <- length(prices)
    
    shares <- calcShares(object,preMerger,revenue=FALSE)
    
    div <- shares/(1-shares)
    
    
    #diag(owner) <- -1*diag(owner)
    
    margins <- -owner * shares
    diag(margins) <- diag(owner) +  diag(margins)
    margins <- solve(t(margins))
    
    margins <-  as.vector(margins %*% ((log(1-shares)*diag(owner))/(-1*output*alpha*(barg*div - 
                                                               log(1-shares)))))
    
    margins <- margins
    
    if(!level) {margins <- margins / prices }
    
    names(margins) <- object@labels
    
    return(as.vector(margins))
  }
  
)


#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "VertBargBertLogit",
  definition=function(object,preMerger=TRUE, level=FALSE){
    
    #check if class is 2nd score Logit
    is2nd <- any(grepl("2nd",class(object)))
    
    up <- object@up
    down <- object@down
    alpha <- down@slopes$alpha
    
    chain_level <- object@chain_level
    
    ## Set downstream margin equal to pre-merger value when chain_level=="wholesaler
    preMerger_chain <- preMerger
    if(chain_level=="wholesaler") {
      preMerger_chain <- TRUE
       if(!preMerger) down@ownerPre <- ownerToMatrix(down,preMerger= TRUE)
    }
  
    #is2nd <- grepl("2nd",class(object))
  
    ownerUp <- ownerToMatrix(up,preMerger = preMerger)
    ownerDown <- ownerToMatrix(down,preMerger= preMerger)
    
    nprods <- nrow(ownerDown)
    
    if( preMerger) {
      bargparm <- up@bargpowerPre
       
      #if(length(up@pricePre) == 0 ){
      #priceUp <- up@prices
      #}
      
      #else{
        priceUp <- up@pricePre
       # }
      #if(length(down@pricePre) == 0 ){
       # priceDown <- down@prices
      #  down@pricePre <- priceDown
      #}
      #else{
        priceDown <- down@pricePre
        
       # }
     
      ownerDownLambda <- object@ownerDownLambdaPre
      ownerUpLambda <- object@ownerUpLambdaPre
      ownerVDown <- object@ownerDownPre
      
      
      #down@mcPre <- down@mcPre + priceUp
      down@ownerPre <- ownerDown
      
    }
    
    else{
      bargparm <- up@bargpowerPost
      
      
      priceUp <- up@pricePost
      priceDown <- down@pricePost
      
      ownerDownLambda <- object@ownerDownLambdaPost
      ownerUpLambda <- object@ownerUpLambdaPost
      ownerVDown <- object@ownerDownPost
      
      down@ownerPost <- ownerDown
      
    }
    
    shareDown <- calcShares(object,preMerger=preMerger, revenue=FALSE)
    elast <-  -alpha*tcrossprod(shareDown)
    diag(elast) <- alpha*shareDown + diag(elast)
    elast.inv <- try(solve(ownerDown * elast),silent=TRUE)
    if(any(class(elast.inv) == "try-error")){elast.inv <- MASS::ginv(ownerDown * elast)}
    
    if(is2nd){
      alpha <- down@slopes$alpha
      down@slopes$meanval <- down@slopes$meanval + alpha *(priceUp - down@priceOutside) 
    }
    
    marginsDown <- calcMargins(down, preMerger=preMerger_chain,level=TRUE)
    
    div <- tcrossprod(1/(1-shareDown),shareDown)*shareDown
    diag(div) <- -shareDown
    div <- as.vector(div)
    
    
    #marginsDown <-  marginsDown - elast.inv %*% ( (ownerVDown * elast) %*% (priceCandUp-mcUp) )
    
    #marginsUp <-  as.vector(solve(ownerUpLambda * div) %*% (((ownerDownLambda * div) %*% (marginsDown)))) 
    
    if(chain_level == "retailer" && !preMerger){ 
      #marginsUp <- calcMargins(object,preMerger=TRUE,level=TRUE)$up
      marginsUp <- up@pricePre - up@mcPost
      }
    else{
    marginsUpPart <-  try(solve(ownerUpLambda * div) %*% (ownerDownLambda * div) ,silent=TRUE)
    if(any(class(marginsUpPart) == "try-error")){
      marginsUpPart <-  MASS::ginv(ownerUpLambda * div) %*% (ownerDownLambda * div) 
    }
    
    marginsUp <- try(solve(diag(nprods) + (marginsUpPart %*% elast.inv %*%  (ownerVDown * elast))),silent=TRUE)
    
    if(any(class(marginsUp) == "try-error")){
    marginsUp <- MASS::ginv(diag(nprods) + (marginsUpPart %*% elast.inv %*%  (ownerVDown * elast)))
    }
    marginsUp <- drop(marginsUp %*% marginsUpPart %*% marginsDown)
    
    }

    if(!is2nd){
      
      if(chain_level == "wholesaler" && !preMerger){ 
        #marginsDown <- calcMargins(object,preMerger=TRUE,level=TRUE)$down
        marginsDown <- down@pricePre - down@mcPost - priceUp
        } 
      else{ marginsDown <-  marginsDown - elast.inv %*% ( (ownerVDown * elast) %*% marginsUp )}
      
      }
  
    if(!level) {
      marginsUp <- marginsUp/priceUp
      marginsDown <- marginsDown/priceDown
    }
     
      
      
    
    names(marginsUp) <- up@labels
    names(marginsDown) <- down@labels
    
    return(list(up=as.vector(marginsUp), down = as.vector(marginsDown))
           )
  }
  
)

#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "VertBargBertLogitNests",
  definition=function(object,preMerger=TRUE, level=FALSE){

    #check if class is 2nd score Logit
    is2nd <- any(grepl("2nd",class(object)))

    up <- object@up
    down <- object@down
    alpha <- down@slopes$alpha
    nestParm <- 1 - down@slopes$sigma #adjust to accomodate below different nesting paramter definition
    
    nests <- down@nests
    nests <- droplevels(nests)
    nestMat <- tcrossprod(model.matrix(~-1+nests))

    
    chain_level <- object@chain_level

    ## Set downstream margin equal to pre-merger value when chain_level=="wholesaler
    preMerger_chain <- preMerger
    if(chain_level=="wholesaler") {
      preMerger_chain <- TRUE
      if(!preMerger) down@ownerPre <- ownerToMatrix(down,preMerger= TRUE)
    }

    #is2nd <- grepl("2nd",class(object))

    ownerUp <- ownerToMatrix(up,preMerger = preMerger)
    ownerDown <- ownerToMatrix(down,preMerger= preMerger)

    nprods <- nrow(ownerDown)

    if( preMerger) {
      bargparm <- up@bargpowerPre

      #if(length(up@pricePre) == 0 ){
      #priceUp <- up@prices
      #}

      #else{
      priceUp <- up@pricePre
      # }
      #if(length(down@pricePre) == 0 ){
      # priceDown <- down@prices
      #  down@pricePre <- priceDown
      #}
      #else{
      priceDown <- down@pricePre

      # }

      ownerDownLambda <- object@ownerDownLambdaPre
      ownerUpLambda <- object@ownerUpLambdaPre
      ownerVDown <- object@ownerDownPre


      #down@mcPre <- down@mcPre + priceUp
      down@ownerPre <- ownerDown

    }

    else{
      bargparm <- up@bargpowerPost


      priceUp <- up@pricePost
      priceDown <- down@pricePost

      ownerDownLambda <- object@ownerDownLambdaPost
      ownerUpLambda <- object@ownerUpLambdaPost
      ownerVDown <- object@ownerDownPost

      down@ownerPost <- ownerDown

    }

    shareDown <- calcShares(object,preMerger=preMerger, revenue=FALSE)
    shareBetweenDown <- as.vector(unname(tapply(shareDown,nests,sum)[nests]))
    shareWithinDown <- shareDown/shareBetweenDown

    #elast <-  -alpha*tcrossprod(shareDown)
    #diag(elast) <- alpha*shareDown + diag(elast)


    elastDiff <- -alpha*tcrossprod(shareDown)
    elastSame <- -alpha*tcrossprod(shareDown*(1 + nestParm[nests] /( 1-nestParm[nests] )*(1/shareBetweenDown)),shareDown)
    elast <- elastSame*nestMat + elastDiff*(1-nestMat)
    diag(elast) <- alpha/(1-nestParm[nests])*shareDown*(1-nestParm[nests]*shareWithinDown - (1-nestParm[nests])*shareDown)
    elast.inv <- try(solve(ownerDown * elast),silent=TRUE)
    if(any(class(elast.inv) == "try-error")){elast.inv <- MASS::ginv(ownerDown * elast)}

    marginsDown = - as.vector(elast.inv %*% shareDown)




    # if(is2nd){
    #   alpha <- down@slopes$alpha
    #   down@slopes$meanval <- down@slopes$meanval + alpha *(priceUp - down@priceOutside)
    # }

    #marginsDown <- calcMargins(down, preMerger=preMerger_chain,level=TRUE)

    # div <- tcrossprod(1/(1-shareDown),shareDown)*shareDown
    # diag(div) <- -shareDown
    # div <- as.vector(div)



    divSameNest <- (1-shareWithinDown)^nestParm[nests]*(1 - shareBetweenDown + shareBetweenDown * (1- shareWithinDown)^(1-nestParm[nests]) )
    divSameNest <-  tcrossprod((1/divSameNest - 1), shareDown)
    ##    matrix(shares * (1/divSameNest - 1), nrow = nprods,ncol=nprods,byrow=TRUE)
    divDiffNest <- shareBetweenDown * (1 - (1- shareWithinDown)^(1-nestParm[nests]) )
    divDiffNest <-  1/divDiffNest - 1
    divDiffNest <-  tcrossprod(1/divDiffNest, shareDown)
    ##   matrix(shares/divDiffNest, nrow=nprods, ncol=nprods, byrow=TRUE)

    div <- divSameNest * nestMat + divDiffNest*(1-nestMat)

    diag(div) <- -shareDown
    div[is.na(div)] <- 0



    #marginsDown <-  marginsDown - elast.inv %*% ( (ownerVDown * elast) %*% (priceCandUp-mcUp) )

    #marginsUp <-  as.vector(solve(ownerUpLambda * div) %*% (((ownerDownLambda * div) %*% (marginsDown))))

    if(chain_level == "retailer" && !preMerger){
      #marginsUp <- calcMargins(object,preMerger=TRUE,level=TRUE)$up
      marginsUp <- up@pricePre - up@mcPost
    }
    else{
      marginsUpPart <-  try(solve(ownerUpLambda * div) %*% (ownerDownLambda * div) ,silent=TRUE)
      if(any(class(marginsUpPart) == "try-error")){
        marginsUpPart <-  MASS::ginv(ownerUpLambda * div) %*% (ownerDownLambda * div)
      }

      marginsUp <- try(solve(diag(nprods) + (marginsUpPart %*% elast.inv %*%  (ownerVDown * elast))),silent=TRUE)

      if(any(class(marginsUp) == "try-error")){
        marginsUp <- MASS::ginv(diag(nprods) + (marginsUpPart %*% elast.inv %*%  (ownerVDown * elast)))
      }
      marginsUp <- drop(marginsUp %*% marginsUpPart %*% marginsDown)

    }

    if(!is2nd){

      if(chain_level == "wholesaler" && !preMerger){
        #marginsDown <- calcMargins(object,preMerger=TRUE,level=TRUE)$down
        marginsDown <- down@pricePre - down@mcPost - priceUp
      }
      else{ marginsDown <-  marginsDown - elast.inv %*% ( (ownerVDown * elast) %*% marginsUp )}

    }

    if(!level) {
      marginsUp <- marginsUp/priceUp
      marginsDown <- marginsDown/priceDown
    }




    names(marginsUp) <- up@labels
    names(marginsDown) <- down@labels

    return(list(up=as.vector(marginsUp), down = as.vector(marginsDown))
    )
  }

)



#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Auction2ndCap",
  definition=function(object,preMerger=TRUE,exAnte=TRUE,level=FALSE){

    result <- calcProducerSurplus(object,preMerger=preMerger,exAnte=exAnte)
    
    if(!level) result <- result / calcPrices(object,preMerger=preMerger,exAnte=exAnte)
    
    return(result)
  })

## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE, level=FALSE){


     output <- object@output
     
    if(preMerger){
      prices <- object@pricePre
    }
    else{prices <- object@pricePost}

    mc <- calcMC(object,preMerger = preMerger)
    prices <- matrix(prices, ncol=length(prices), nrow=length(mc),byrow=TRUE)



    if(output){margin <- prices - mc}
    else{margin <- mc - prices}
    if(!level){margin <- margin/prices}

    dimnames(margin) <- object@labels
    return(margin)
  }

)


#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "AIDS",
  definition=function(object,preMerger=TRUE,level=FALSE){

    priceDelta <- object@priceDelta
    ownerPre   <- object@ownerPre
    shares     <- calcShares(object,TRUE)

    elastPre <-  t(elast(object,TRUE))
    
    elastInv <- try(solve(elastPre * ownerPre),silent=TRUE)
    if(any(class(elastInv)=="try-error")){elastInv <- MASS::ginv(elastPre * ownerPre)}
    
    marginPre <-  -1 * as.vector(elastInv %*% (shares * diag(ownerPre))) / shares

    if(preMerger){
      names(marginPre) <- object@labels
      
      if(level) {marginPre <- marginPre * object@pricePre}
      return(marginPre)}

    else{

      marginPost <- 1 - ((1 + object@mcDelta) * (1 - marginPre) / (priceDelta + 1) )
      names(marginPost) <- object@labels
      if(level) {marginPost <- marginPost * object@pricePost}
      return(marginPost)
    }

  }
)

## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "LogitCap",
  definition=function(object,preMerger=TRUE,level=FALSE){

    margins <- object@margins #capacity-constrained margins not identified -- use supplied margins

    if( preMerger) {
      capacities <- object@capacitiesPre
      prices <- object@pricePre

    }
    else{

      capacities <- object@capacitiesPost
      prices <- object@pricePost
    }



    quantities <- calcQuantities(object, preMerger=TRUE)
    constrained <-  abs(capacities - quantities) < 1e-5

    owner  <- object@ownerPre
    revenue<- calcShares(object,preMerger,revenue=TRUE)[!constrained]
    elast <-  elast(object,preMerger)
    margins[!constrained] <-  -1 * as.vector(MASS::ginv(t(elast*owner)[!constrained,!constrained]) %*% revenue) / revenue



    names(margins) <- object@labels

    if(level) {margins <- margins * prices} 
    return(as.vector(margins))
  }

)



## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Auction2ndLogit",
  definition=function(object,preMerger=TRUE,exAnte=FALSE,level=TRUE){


    output <- ifelse(object@output,1,-1)
    
    nprods <- length(object@shares)
   
   
    margins <- rep(NA,nprods)

    if( preMerger) {
      owner  <- object@ownerPre
      subset <- rep(TRUE,nprods)
      prices <- object@pricePre
    }
    else{
    subset <- object@subset  
    owner  <- object@ownerPost
    prices <- object@pricePost
    }

    owner <- owner[subset,subset]
    alpha <- object@slopes$alpha
    shares <- calcShares(object,preMerger=preMerger,revenue=FALSE)
    shares <- shares[subset]
    firmShares <- drop(owner %*% shares)
    margins[subset] <-  output*log(1-firmShares)/(alpha * firmShares)

    if(exAnte){ margins[subset] <-  margins[subset] * shares}
    if(!level){margins <- margins/prices}
    names(margins) <- object@labels

    return(as.vector(margins))
  }

)


#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Auction2ndLogitNests",
  definition=function(object,preMerger=TRUE,exAnte=FALSE,level=FALSE){
    
    output <- ifelse(object@output,1,-1)
    nprods <- length(object@shares)
    
    
    if( preMerger) {
      owner  <- object@ownerPre
      subset <- rep(TRUE, nprods)
      prices <- object@pricePre}
    else{
      owner  <- object@ownerPost
      subset <- object@subset
      prices <- object@pricePost}
    
    
    nests <- object@nests
    nests <- droplevels(nests[subset])
    nestMat <- tcrossprod(model.matrix(~-1+nests))
    
    margins <- rep(NA,nprods)
    
    
    owner <- owner[subset,subset]
    
    
    alpha <- object@slopes$alpha
    sigma <- object@slopes$sigma
    
    shares <- calcShares(object,preMerger=preMerger,revenue=FALSE)
    
    shares <- shares[subset]
    
    sharesAcross <- as.vector(tapply(shares,nests,sum))[nests]
    sharesIn <- shares/sharesAcross
    
    firmShares <- drop(owner %*% shares)
    
    dupCnt <- rowSums(owner*nestMat) #only include the values in a given nest once
    
    ownerValue <-   1 - (1 - ((owner*nestMat)%*%sharesIn))^sigma[nests]
    ownerValue <-  drop(1 - owner %*%(ownerValue * sharesAcross/dupCnt))
    
    margins[subset] <-  output*log(ownerValue)/(alpha*firmShares) 
    
    
    if(exAnte){ margins[subset] <-  margins[subset] * shares}
    if(!level){margins <- margins/prices}
    names(margins) <- object@labels
    
    return(as.vector(margins))
  }
  
)
