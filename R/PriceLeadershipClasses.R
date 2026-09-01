# PriceLeadership subclasses are kept in a separate file so that the
# existing Bertrand RUM class definitions remain unchanged.

#'@rdname BertrandRUM-Classes
#'@export
setClass(
  Class   = "PriceLeadership",
  contains = "Logit",
  slots = list(
    supermarkupPre = "numeric",   # Pre-merger supermarkup m above Bertrand prices
    supermarkupPost = "numeric",  # Post-merger supermarkup m above Bertrand prices
    timingParam = "numeric",      # Firm-specific timing/discount parameters delta in (0,1) - named vector
    coalitionPre = "ANY",         # Pre-merger coalition: vector of indices OR control matrix
    coalitionPost = "ANY",        # Post-merger coalition: vector of indices OR control matrix
    bindingFirm = "numeric",      # Index of firm with binding IC constraint (or NA if none)
    slackValues = "numeric"       # Slack function values for each coalition firm
  ),
  prototype = prototype(
    supermarkupPre = NA_real_,
    supermarkupPost = NA_real_,
    timingParam = numeric(),      # Empty named vector initially
    coalitionPre = numeric(),
    coalitionPost = numeric(),
    bindingFirm = NA_real_,
    slackValues = numeric()
  ),
  validity = function(object){

    # Call parent validity first
    # Parent (Logit) already checks: prices, margins, shares, normIndex, etc.

    nprods <- length(object@prices)

    # PriceLeadership-specific validations only

    # Pre-merger coalition (can be vector or matrix)
    if(length(object@coalitionPre) > 0){
      if(is.matrix(object@coalitionPre)){
        # Matrix form - check dimensions
        if(nrow(object@coalitionPre) != nprods || ncol(object@coalitionPre) != nprods){
          stop("'coalitionPre' matrix must be ", nprods, " x ", nprods)
        }
        # Extract coalition indices to validate
        coalitionIndices <- which(rowSums(object@coalitionPre) > 1)
      } else {
        # Vector form - validate indices
        coalitionIndices <- object@coalitionPre
        if(any(!(coalitionIndices %in% 1:nprods))){
          stop("'coalitionPre' must contain product indices between 1 and ", nprods)
        }
      }

      # Coalition must have at least 2 products
      if(length(coalitionIndices) < 2){
        stop("'coalitionPre' must contain at least 2 products for price leadership")
      }

      # Must have at least one fringe product
      if(length(coalitionIndices) >= nprods){
        stop("At least one fringe product is required for price leadership model")
      }
    }

    # Post-merger coalition (can be vector or matrix)
    if(length(object@coalitionPost) > 0){
      if(is.matrix(object@coalitionPost)){
        # Matrix form - check dimensions
        if(nrow(object@coalitionPost) != nprods || ncol(object@coalitionPost) != nprods){
          stop("'coalitionPost' matrix must be ", nprods, " x ", nprods)
        }
        # Extract coalition indices to validate
        coalitionIndices <- which(rowSums(object@coalitionPost) > 1)
      } else {
        # Vector form - validate indices
        coalitionIndices <- object@coalitionPost
        if(any(!(coalitionIndices %in% 1:nprods))){
          stop("'coalitionPost' must contain product indices between 1 and ", nprods)
        }
      }

      # Coalition must have at least 2 products
      if(length(coalitionIndices) < 2){
        stop("'coalitionPost' must contain at least 2 products for price leadership")
      }

      # Must have at least one fringe product
      if(length(coalitionIndices) >= nprods){
        stop("At least one fringe product is required for price leadership model")
      }
    }

    # Timing parameter must be in (0,1) if specified
    if(length(object@timingParam) > 0){
      timingValues <- object@timingParam[!is.na(object@timingParam)]
      if(length(timingValues) > 0 && isTRUE(any(timingValues <= 0 | timingValues >= 1))){
        stop("'timingParam' values must be between 0 and 1 (exclusive)")
      }
    }

    # Supermarkup should be non-negative if specified
    if(!is.na(object@supermarkupPre)){
      if(object@supermarkupPre < 0){
        stop("'supermarkupPre' must be non-negative")
      }
    }
    if(!is.na(object@supermarkupPost)){
      if(object@supermarkupPost < 0){
        stop("'supermarkupPost' must be non-negative")
      }
    }

    # bindingFirm must be valid product index if specified
    if(!is.na(object@bindingFirm)){
      if(length(object@bindingFirm) != 1 || !(object@bindingFirm %in% 1:nprods)){
        stop("'bindingFirm' must be a single product index between 1 and ", nprods, " or NA")
      }
    }

    return(TRUE)
  }
)

#'@rdname BertrandRUM-Classes
#'@export
setClass(
  Class   = "PriceLeadershipBLP",
  contains = "LogitBLP",
  slots = list(
    supermarkupPre = "numeric",   # Pre-merger supermarkup m above Bertrand prices
    supermarkupPost = "numeric",  # Post-merger supermarkup m above Bertrand prices
    timingParam = "numeric",      # Firm-specific timing/discount parameters delta in (0,1) - named vector
    coalitionPre = "ANY",         # Pre-merger coalition: vector of indices OR control matrix
    coalitionPost = "ANY",        # Post-merger coalition: vector of indices OR control matrix
    bindingFirm = "numeric",      # Index of firm with binding IC constraint (or NA if none)
    slackValues = "numeric"       # Slack function values for each coalition firm
  ),
  prototype = prototype(
    supermarkupPre = NA_real_,
    supermarkupPost = NA_real_,
    timingParam = numeric(),      # Empty named vector initially
    coalitionPre = numeric(),
    coalitionPost = numeric(),
    bindingFirm = NA_real_,
    slackValues = numeric()
  ),
  validity = function(object){

    # Call parent validity first (LogitBLP checks)
    # Parent (LogitBLP -> Logit) already checks: prices, shares, nDraws, etc.

    nprods <- length(object@prices)

    # PriceLeadershipBLP-specific validations

    # Validate coalitionPre is provided
    if(length(object@coalitionPre) == 0){
      stop("'coalitionPre' must be specified (vector of product indices in coordinating coalition)")
    }

    ## Validate that slopes is a list (BLP parameters will be checked more thoroughly)
    if(!is.null(object@slopes) && !is.list(object@slopes)){
      stop("'slopes' must be a list containing pre-calibrated BLP demand parameters")
    }

    # Detailed BLP parameter validation
    if(is.null(object@slopes) || !is.list(object@slopes)){
      stop("BLP demand parameters must be provided in 'slopes' list")
    }

    if(!("sigma" %in% names(object@slopes))){
      stop("'slopes' must contain 'sigma' (std. dev. of random price coefficient)")
    }

    if(!("alpha" %in% names(object@slopes)) && !("alphaMean" %in% names(object@slopes))){
      stop("'slopes' must contain either 'alpha' or 'alphaMean' (mean price coefficient)")
    }

    # Note: Additional validations (like meanval consistency) are deferred to calcSlopes

    # Set default values for optional BLP parameters if not provided
    if(!("sigmaNest" %in% names(object@slopes))){
      object@slopes$sigmaNest <- 1  # Default: no nesting (flat logit)
    }
    if(!("nDemog" %in% names(object@slopes))){
      object@slopes$nDemog <- 0
    }
    if(!("piDemog" %in% names(object@slopes))){
      object@slopes$piDemog <- numeric(0)
    }
    if(!("nDraws" %in% names(object@slopes))){
      if(length(object@nDraws) > 0 && !is.na(object@nDraws)){
        object@slopes$nDraws <- object@nDraws
      } else {
        object@slopes$nDraws <- 500  # Default
      }
    }

    # Pre-merger coalition (can be vector or matrix)
    if(length(object@coalitionPre) > 0){
      if(is.matrix(object@coalitionPre)){
        # Matrix form - check dimensions
        if(nrow(object@coalitionPre) != nprods || ncol(object@coalitionPre) != nprods){
          stop("'coalitionPre' matrix must be ", nprods, " x ", nprods)
        }
        # Extract coalition indices to validate
        coalitionIndices <- which(rowSums(object@coalitionPre) > 1)
      } else {
        # Vector form - validate indices
        coalitionIndices <- object@coalitionPre
        if(any(!(coalitionIndices %in% 1:nprods))){
          stop("'coalitionPre' must contain product indices between 1 and ", nprods)
        }
      }

      # Coalition must have at least 2 products
      if(length(coalitionIndices) < 2){
        stop("'coalitionPre' must contain at least 2 products for price leadership")
      }

      # Must have at least one fringe product
      if(length(coalitionIndices) >= nprods){
        stop("At least one fringe product is required for price leadership model")
      }
    }

    # Post-merger coalition (can be vector or matrix)
    if(length(object@coalitionPost) > 0){
      if(is.matrix(object@coalitionPost)){
        # Matrix form - check dimensions
        if(nrow(object@coalitionPost) != nprods || ncol(object@coalitionPost) != nprods){
          stop("'coalitionPost' matrix must be ", nprods, " x ", nprods)
        }
        # Extract coalition indices to validate
        coalitionIndices <- which(rowSums(object@coalitionPost) > 1)
      } else {
        # Vector form - validate indices
        coalitionIndices <- object@coalitionPost
        if(any(!(coalitionIndices %in% 1:nprods))){
          stop("'coalitionPost' must contain product indices between 1 and ", nprods)
        }
      }

      # Coalition must have at least 2 products
      if(length(coalitionIndices) < 2){
        stop("'coalitionPost' must contain at least 2 products for price leadership")
      }

      # Must have at least one fringe product
      if(length(coalitionIndices) >= nprods){
        stop("At least one fringe product is required for price leadership model")
      }
    }

    # Timing parameter must be in (0,1) if specified
    if(length(object@timingParam) > 0){
      timingValues <- object@timingParam[!is.na(object@timingParam)]
      if(length(timingValues) > 0 && isTRUE(any(timingValues <= 0 | timingValues >= 1))){
        stop("'timingParam' values must be between 0 and 1 (exclusive)")
      }
    }

    # Supermarkup should be non-negative if specified
    if(!is.na(object@supermarkupPre)){
      if(object@supermarkupPre < 0){
        stop("'supermarkupPre' must be non-negative")
      }
    }
    if(!is.na(object@supermarkupPost)){
      if(object@supermarkupPost < 0){
        stop("'supermarkupPost' must be non-negative")
      }
    }

    # bindingFirm must be valid product index if specified
    if(!is.na(object@bindingFirm)){
      if(length(object@bindingFirm) != 1 || !(object@bindingFirm %in% 1:nprods)){
        stop("'bindingFirm' must be a single product index between 1 and ", nprods, " or NA")
      }
    }

    return(TRUE)
  }
)

#'@rdname BertrandRUM-Classes
#'@export
setClass(
  Class   = "LogitCournot",
  contains="Logit"
)
