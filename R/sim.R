sim <- function(prices,demand=c("Linear","LogLog","Logit","CES","LogitNests","CESNests"),demand.param,
                     ownerPre,ownerPost,nests=NULL,
                     mcDelta=rep(0,length(prices)),
                     priceStart=prices,
                     labels=paste("Prod",1:length(prices),sep=""),...){

    demand <- match.arg(demand)
    nprods <- length(prices)


    ## Create placeholders values to fill required Class slots
    normIndex   <- 1
    mc <- rep(0,nprods)
    shares <- margins <- mc + .05



    if(!is.null(nests)){nests <- factor(nests)}

    ## Create a function to recover marginal cost using
    ## demand parameters and supplied prices

    calcMC <- function(object){

        object@pricePre <- prices

        shares <- calcShares(object,TRUE)
        elastPre <-  t(elast(object,TRUE))
        marginPre <-  -1 * as.vector(solve(elastPre*object@ownerPre) %*% shares) / shares

        mc <- (1 - marginPre) * prices
        names(mc) <- object@labels

        return(mc)
    }


    ## Create constructors for each demand system specified in the 'demand' paramater

    if(demand == "CESNests"){

        if(!is.list(demand.param) ||
           any(!c("gamma","meanval","sigma") %in% names(demand.param)) ){
            stop("'demand.param' must be a list containing 'gamma' (price coefficient), 'sigma' (nesting parameters) and 'meanval' (mean valuations)")}
        if(is.null(nests) ||
           length(nests)!= nprods ){stop("When 'demand' equals 'CESNests', 'nests' must equal a vector whose length equals the number of products")}

        if(nlevels(nests) != length(demand.param$sigma)){stop("The number of nests in 'nests' does not equal to the number of nesting paramters in 'demand.param$sigma'")}

        ## uncover Numeraire Coefficients
        if(is.null(demand.param$shareInside)){
            warning("'demand.param' does not contain 'shareInside'. Setting demand.param$shareInside=1")
            demand.param$shareInside=1}
        if(demand.param$shareInside < 1 ) {alpha <- 1/demand.param$shareInside - 1}
        else{alpha <- NULL}

        demand.param$alpha <- alpha

         result <- new(demand,prices=prices, shares=shares,margins=margins,
                       mc=mc,mcDelta=mcDelta,
                       ownerPre=ownerPre,
                       ownerPost=ownerPost,
                       nests=nests,
                       normIndex=normIndex,
                       parmsStart=c(demand.param$gamma,demand.param$sigma),
                       priceStart=priceStart,
                       shareInside=demand.param$shareInside,labels=labels)

    }

     else if(demand == "LogitNests"){

        if(!is.list(demand.param) ||
           any(!c("alpha","meanval","sigma") %in% names(demand.param)) ){
            stop("'demand.param' must be a list containing 'alpha' (price coefficient), 'sigma' (nesting parameters) and 'meanval' (mean valuations)")}
        if(is.null(demand.param$shareInside)){
            warning("'demand.param' does not contain 'shareInside'. Setting demand.param$shareInside=1")
            demand.param$shareInside=1}
        if(is.null(nests) ||
           length(nests)!= nprods ){stop("When 'demand' equals 'LogitNests', 'nests' must equal a vector whose length equals the number of products")}
        if(nlevels(nests) != length(demand.param$sigma)){stop("The number of nests in 'nests' does not equal to the number of nesting paramters in 'demand.param$sigma'")}

         result <- new(demand,prices=prices, shares=shares,margins=margins,
                       mc=mc,mcDelta=mcDelta,
                       ownerPre=ownerPre,
                       ownerPost=ownerPost,
                       nests=nests,
                       normIndex=normIndex,
                       parmsStart=c(demand.param$alpha,demand.param$sigma),
                       priceStart=priceStart,
                       shareInside=demand.param$shareInside,labels=labels)

    }

    else if(demand == "CES"){

        if(!is.list(demand.param) ||
           any(!c("gamma","meanval") %in% names(demand.param)) ){
            stop("'demand.param' must be a list containing 'gamma' and 'meanval'")}

        ## uncover Numeraire Coefficients
        if(is.null(demand.param$alpha)){
            if(is.null(demand.param$shareInside)){
                warning("'demand.param' does not contain either 'alpha' or 'shareInside'. Setting demand.param$shareInside=1")
                demand.param$shareInside=1}

            if(demand.param$shareInside < 1 ) {alpha <- 1/demand.param$shareInside - 1}
            else{alpha <- NULL}
        }


        demand.param$alpha <- alpha

        result <- new(demand,prices=prices, shares=shares,
                  margins=margins,
                  normIndex=normIndex,
                  mc=mc,mcDelta=mcDelta,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  priceStart=priceStart,shareInside=demand.param$shareInside,
                  labels=labels)

    }

  else if(demand == "Logit"){

        if(!is.list(demand.param) ||
           any(!c("alpha","meanval") %in% names(demand.param)) ||
           demand.param$alpha>0 ||
           length(demand.param$meanval) != nprods){
            stop("'demand.param' must be a list containing 'alpha' and 'meanval'. 'alpha' must be negative and 'meanval' must be a length-k vector,")}

        if(is.null(demand.param$shareInside)){
            warning("'demand.param' does not contain 'shareInside'. Setting demand.param$shareInside=1")
            demand.param$shareInside=1}

        result <- new(demand,prices=prices, shares=shares,
                  margins=margins,
                  normIndex=normIndex,
                  mc=mc,mcDelta=mcDelta,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  priceStart=priceStart,shareInside=demand.param$shareInside,
                  labels=labels)

    }


 else if(demand == "Linear"){

     if(!(is.matrix(demand.param)    ||
          is(demand.param,"Matrix")) ||
        ncol(demand.param)!=nprods+1 ||
        nrow(demand.param)!=nprods   ||
        any(demand.param[,1]<0)      ||
        any(diag(demand.param[,-1])>0)){
         stop("'demand.param' must be a k x (k+1) matrix whose first column contains the intercepts from a linear demand system and whose remaining columns contain the slope coefficients. The intercepts should all be non-negative and the diagonal of the slopes should be non-positive.")}

      result <- new(demand,prices=prices, quantities=shares,margins=margins,
                   shares=shares,mc=mc,mcDelta=mcDelta,
                   ownerPre=ownerPre,diversion=diag(nprods), symmetry=FALSE,
                   ownerPost=ownerPost, labels=labels)

 }


 else if(demand == "LogLog"){

     if(!(is.matrix(demand.param)      ||
          is(demand.param,"Matrix"))   ||
        ncol(demand.param) != nprods+1 ||
        nrow(demand.param) != nprods   ||
        any(diag(demand.param[,-1])>0)){
         stop("'demand.param' must be a k x (k+1) matrix whose first column contains the intercepts from a linear demand system and whose remaining columns contain the slope coefficients. The diagonal of the slopes should be non-positive.")}

      result <- new(demand,prices=prices, quantities=shares,margins=margins,
                   shares=shares,mc=mc,mcDelta=mcDelta,  symmetry=FALSE, priceStart=priceStart,
                   ownerPre=ownerPre,diversion=diag(nprods),
                   ownerPost=ownerPost, labels=labels)

 }



    result@slopes=demand.param

    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Uncover marginal costs
    result@mc        <- calcMC(result)

    ## Solve Non-Linear System for Price Changes
    result@pricePre  <- calcPrices(result,TRUE,...)
    result@pricePost <- calcPrices(result,FALSE,...)


    return(result)
}
