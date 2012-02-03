sim <- function(prices,demand=c("Linear","AIDS","LogLog","Logit","CES","LogitNests","CESNests"),demand.param,
                     ownerPre,ownerPost,nests=NULL,
                     mcDelta=rep(0,length(prices)),
                     priceStart=prices,
                     labels=paste("Prod",1:length(prices),sep=""),...){

    demand <- match.arg(demand)
    nprods <- length(prices)


    ## Create placeholders values to fill required Class slots

     shares <- margins <- rep(1/nprods,nprods)


    if(!is.null(nests)){nests <- factor(nests)}


    ## general sanity checks
    if(!is.list(demand.param)){stop("'demand.param' must be a list.")}

    ## Sanity checks for discrete choice models
    if(demand %in% c("CESNests","LogitNests","CES","Logit")){



        if(!("normIndex" %in% names(demand.param))){
                warning("'demand.param' does not contain 'normIndex'. Setting normIndex=1.")
                normIndex=1
            }
        else if(!(demand.param$normIndex %in% c(NA,1:nprods))){
            stop("'normIndex' must take on a value between 1 and ",nprods,
                 " OR equal NA.")
            }
        else{
            normIndex <-  demand.param$normIndex
            demand.param$normIndex <- NULL
        }

         if(demand %in% c("CESNests","LogitNests")){

            if(is.null(nests) ||
               length(nests)!= nprods ){stop("When 'demand' equals 'CESNests' or 'LogitNests', 'nests' must equal a vector whose length equals the number of products.")}

            if(nlevels(nests) != length(demand.param$sigma)){stop("The number of nests in 'nests' does not equal to the number of nesting paramters in 'demand.param$sigma'.")}

        }
     }


    ## Sanity checks for Linear-demand style models
    if(demand %in% c("Linear","LogLog","AIDS")){

        if(!("slopes" %in% names(demand.param))){stop("'demand.param' does not contain 'slopes'")}
        if(!("intercepts" %in% names(demand.param))){stop("'demand.param' does not contain 'intercepts'")}

            if(!(is.matrix(demand.param$slopes)    ||
                 is(demand.param$slopes,"Matrix")) ||
               ncol(demand.param$slopes)!=nprods   ||
               nrow(demand.param$slopes)!=nprods   ||
               any(demand.param$intercepts<0)      ||
               any(diag(demand.param$slopes)>0)){
                stop("'demand.param' must be a k x k matrix whose first column contains the intercepts from a linear style demand system and whose remaining columns contain the slope coefficients. The intercepts must all be non-negative and the diagonal of the slopes must be non-positive")}

        if (demand == "AIDS" &&
            !("mktElast" %in% names(demand.param))){
                warning("'demand.param' does not contain 'mktElast'. Setting 'mktElast' equal to -1")
                demand.param$mktElast=-1

            }

        }






    ## Create constructors for each demand system specified in the 'demand' parameter

    if(demand == "CESNests"){

        if(any(!c("gamma","meanval","sigma") %in% names(demand.param)) ){
            stop("'demand.param' must be a list containing 'gamma' (price coefficient), 'sigma' (nesting parameters) and 'meanval' (mean valuations).")}


        ## uncover Numeraire Coefficients
        if(!("alpha" %in% names(demand.param)) &&
           !("shareInside" %in% names(demand.param))){
                warning("'demand.param' does not contain either 'alpha' or 'shareInside'. Setting shareInside=1 and alpha=NULL.")
                shareInside=1
                demand.param$alpha <- NULL}

        else if("shareInside" %in% names(demand.param)){
            shareInside=demand.param$shareInside
            demand.param$shareInside <- NULL

            if(shareInside<1) {demand.param$alpha <- 1/shareInside -1}
            else{ demand.param$alpha <- NULL}
        }




         result <- new(demand,prices=prices, shares=shares,margins=margins,
                       mcDelta=mcDelta,
                       ownerPre=ownerPre,
                       ownerPost=ownerPost,
                       nests=nests,
                       normIndex=normIndex,
                       parmsStart=c(demand.param$gamma,demand.param$sigma),
                       priceStart=priceStart,
                       shareInside=shareInside,labels=labels)

    }

    else if(demand == "LogitNests"){

        if( any(!c("alpha","meanval","sigma") %in% names(demand.param)) ){
            stop("'demand.param' must be a list containing 'alpha' (price coefficient), 'sigma' (nesting parameters) and 'meanval' (mean valuations).")}

        if(!("shareInside" %in% names(demand.param))){
            warning("'demand.param' does not contain 'shareInside'. Setting shareInside=1.")
            shareInside <- 1}
        else{
            shareInside <- demand.param$shareInside
            demand.param$shareInside <- NULL
        }

        result <- new(demand,prices=prices, shares=shares,margins=margins,
                       mcDelta=mcDelta,
                       ownerPre=ownerPre,
                       ownerPost=ownerPost,
                       nests=nests,
                       normIndex=normIndex,
                       parmsStart=c(demand.param$alpha,demand.param$sigma),
                       priceStart=priceStart,
                       shareInside=shareInside,labels=labels)

    }

    else if(demand == "CES"){

        if(any(!c("gamma","meanval") %in% names(demand.param)) ){
            stop("'demand.param' must be a list containing 'gamma' and 'meanval'.")}

        ## uncover Numeraire Coefficients
        if(!("alpha" %in% names(demand.param)) &&
           !("shareInside" %in% names(demand.param))){
                warning("'demand.param' does not contain either 'alpha' or 'shareInside'. Setting shareInside=1 and alpha=NULL.")
                shareInside=1
                demand.param$alpha=NULL}

        else if("shareInside" %in% names(demand.param)){
            shareInside=demand.param$shareInside
            demand.param$shareInside <- NULL

            if(shareInside<1) {demand.param$alpha <- 1/shareInside -1}
            else{ demand.param$alpha <- NULL}
        }


        result <- new(demand,prices=prices, shares=shares,
                  margins=margins,
                  normIndex=normIndex,
                  mcDelta=mcDelta,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  priceStart=priceStart,shareInside=shareInside,
                  labels=labels)

    }

  else if(demand == "Logit"){

        if(any(!c("alpha","meanval") %in% names(demand.param)) ||
           demand.param$alpha>0 ||
           length(demand.param$meanval) != nprods){
            stop("'demand.param' must be a list containing 'alpha' and 'meanval'. 'alpha' must be negative and 'meanval' must be a length-k vector.")}

        if(!("shareInside" %in% names(demand.param))){
            warning("'demand.param' does not contain 'shareInside'. Setting shareInside=1.")
            shareInside=1}
        else{shareInside=demand.param$shareInside}

        result <- new(demand,prices=prices, shares=shares,
                  margins=margins,
                  normIndex=normIndex,
                  mcDelta=mcDelta,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  priceStart=priceStart,shareInside=shareInside,
                  labels=labels)

    }


 else if(demand == "Linear"){



      result <- new(demand,prices=prices, quantities=shares,margins=margins,
                   shares=shares,mcDelta=mcDelta,
                   ownerPre=ownerPre,diversion=-diag(nprods), symmetry=FALSE,
                   ownerPost=ownerPost, labels=labels)

 }

 else if(demand == "AIDS"){

     ## find the market elasticity that best explains user-supplied intercepts and prices


      result <- new(demand,prices=prices, quantities=shares,margins=margins,
                    shares=as.vector(demand.param$intercepts + demand.param$slopes %*% prices), # AIDS needs actual shares for prediction
                    mcDelta=mcDelta, mktElast=demand.param$mktElast,
                    ownerPre=ownerPre,diversion=-diag(nprods),
                    priceDeltaStart=priceStart,
                    ownerPost=ownerPost, labels=labels)

 }



 else if(demand == "LogLog"){


      result <- new(demand,prices=prices, quantities=shares,margins=margins,
                   shares=shares,mcDelta=mcDelta, priceStart=priceStart,
                   ownerPre=ownerPre,diversion=-diag(nprods),
                   ownerPost=ownerPost, labels=labels)

 }


    if(demand %in% c("Linear","LogLog","AIDS")){
        result@slopes <- demand.param$slopes
        result@intercepts <- demand.param$intercepts
    }
    else{result@slopes=demand.param}


    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Uncover marginal costs
    result@mc        <- calcMC(result)

    if(demand == "AIDS"){
        ## Solve Non-Linear System for Price Changes
        result@priceDelta <- calcPriceDelta(result,...)
    }


    ## Solve Non-Linear System for Price Changes
    result@pricePre  <- calcPrices(result,TRUE,...)
    result@pricePost <- calcPrices(result,FALSE,...)


    return(result)
}
