

upp <- function(price, margin, diversion, ownerPre, efficiency=rep(0,length(price)), labels=paste("Prod",1:length(price),sep="")){


    if(!(is.vector(price) & is.vector(margin)  & is.vector(efficiency))){
        stop("'price' and 'margin' arguments must be vectors")}

    nprod = length(price)


    if(nprod != length(margin) || nprod != length(efficiency) ){
        stop("'price', 'margin', and 'efficiency' vectors must be the same length")}

    if(any(price < 0)){ stop("'price' argument must be non-negative")}
    if(any(margin < 0 || margin > 1) ){ stop("'margin' vector elements  must be between 0 and 1")}
    if(any(efficiency < 0 || efficiency > 1) ){ stop("'efficiency' vector elements  must be between 0 and 1")}

    if(!is.matrix(diversion)){ stop("'diversion' argument must be a matrix")}
    if(!all(diag(diversion) == 1)){ stop("'diversion' diagonal elements must all equal 1")}
    if(any(abs(diversion) > 1)){ stop("'diversion' elements must be between -1 and 1")}


     ## create pre-merger ownership matrix if needed ##
    if(is.vector(ownerPre)){

        if(nprod != length(ownerPre)){
            stop("'price'and 'ownerPre' vectors must be the same length")}

        owners <- as.numeric(factor(ownerPre))
        ownerPre <- matrix(1,ncol=nprod,nrow=nprod)

        for( o in unique(owners)){
            ownerPre[owners == o, owners == o] = NA
        }


        rm(owners)
    }

    else if(is.matrix(ownerPre))  {
        if( ncol(ownerPre) != nrow(ownerPre) || !(ownerPre %in% c(NA,1))){
            stop("'ownerPre' must be a square matrix whose elements are either NA or 1")}

        if(length(diag(ownerPre)) != nprod){
            stop("dimensions of 'ownerPre' must be the same as 'price' vector" )}

    }

    else if(missing(ownerPre)){
        ownerPre <- diag(NA,ncol=length(price), nrow=length(price))+1
    }

    ## compute the implicit "tax" that product i levies on product j ##
    Tax <- t(t(diversion * tcrossprod(1/price,price)) * margin)

    ## subtract any efficiencies ##
    result <- Tax - efficiency * (1-margin)
    result <- result*ownerPre



    dimnames(result) <- list(labels,labels)

    return(result)
}

