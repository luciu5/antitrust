cmcr.cournot <- function(share,industryElasticity){

    if(!is.vector(share) || length(share)!=2) { stop("'share' argument must be a vector of length 2")}
    if(!is.vector(industryElasticity) || length(industryElasticity)!=1) { stop("'industryElasticity' argument must be a vector of length 1")}
    if(any(share < 0,na.rm=TRUE) ||
       any(share >1,na.rm=TRUE) ) {
        stop("'share' argument must be between 0 and 1")}
    if(any(industryElasticity < 0,na.rm=TRUE)) { stop("'industryElasticity' argument must be positive")}

    mcDelta <- 2*prod(share)/(industryElasticity*sum(share) - sum(share^2))

    return(mcDelta*100)
}


cmcr.bertrand <- function(price, margin, diversion, ownerPre, ownerPost=matrix(1,ncol=length(price), nrow=length(price)), labels=paste("Prod",1:length(price),sep=""))
{


    ## Check to make sure inputs are sane ##

    if(!(is.vector(price) & is.vector(margin))){
        stop("'price' and 'margin' arguments must be vectors")}

    nprod = length(price)


    if(nprod != length(margin)){
        stop("'price'and 'margin' vectors must be the same length")}

    if(any(price < 0,na.rm=TRUE)){ stop("'price' argument must be non-negative")}
    if(any(margin < 0 || margin > 1,na.rm=TRUE) ){ stop("'margin' vector elements  must be between 0 and 1")}

    if(!is.matrix(diversion)){ stop("'diversion' argument must be a matrix")}
    if(!all(diag(diversion) == 1,na.rm=TRUE)){ stop("'diversion' diagonal elements must all equal 1")}
    if(any( abs(diversion) > 1,na.rm=TRUE)){ stop("'diversion' elements must be between -1 and 1")}

    if(!is.matrix(ownerPost)){ stop("'ownerPost' argument must be a matrix")}
    if(any(ownerPost < 0 || ownerPost > 1,na.rm=TRUE)){ stop("'ownerPost' elements must be between 0 and 1")}

     ## transform pre-merger ownership vector into matrix##
    if(is.vector(ownerPre)){

        if(nprod != length(ownerPre)){
            stop("'price'and 'ownerPre' vectors must be the same length")}

        owners <- as.numeric(factor(ownerPre))
        ownerPre <- matrix(0,ncol=nprod,nrow=nprod)

        for( o in unique(owners)){
            ownerPre[owners == o, owners == o] = 1
        }

        rm(owners)
    }




    else if(is.matrix(ownerPre))  {
        if( ncol(ownerPre) != nrow(ownerPre) || any(ownerPre < 0 || ownerPre > 1,na.rm=TRUE)){
            stop("'ownerPre' must be a square matrix whose elements are between 0 and 1")}

        if(length(diag(ownerPre)) != nprod){
            stop("dimensions of 'ownerPre' must be the same as 'price' vector" )}

    }



    ## weight diversion ratios by price ratios and ownership matrices ##
    priceRatio = tcrossprod(1/price, price)
    Bpre =  -1 * diversion * priceRatio * ownerPre;  diag(Bpre) = 1
    Bpost = -1 * diversion * priceRatio * ownerPost; diag(Bpost) = 1



    ##calculate post-merger margin ##
    marginPost = as.vector(solve(Bpost) %*% Bpre %*% margin)

    ## calculate proportional change in marginal cost ##
    mcDelta= (marginPost - margin)/(1 - margin)

    names(mcDelta) <- labels

    return(mcDelta*100)
}
