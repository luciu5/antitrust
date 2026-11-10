#' Advance a solved market without estimating or solving it
#'
#' Copies the solved post-state into the pre-state, preserving physical
#' costs, demand parameters, retention, and the active-product map. Optional
#' entry and quality primitives are then applied using the same transitions
#' as antitrust counterfactuals. Ownership, cost shocks, exit, and policy
#' changes belong to the subsequent simulation and are rejected here.
#'
#' This is an extension point for packages that own additional policy state.
#' It never recalibrates or solves an equilibrium. Constant costs are rebased
#' at the solved effective cost and retention, with zero pending cost shock.
#' The standard Logit second-score auction also advances its cost-relative
#' score index, preserving win probabilities. Quantity-dependent cost-function
#' models and other second-score descendants require their own method.
#'
#' @param model A solved legacy antitrust model.
#' @param step An optional `CounterfactualStep` containing only `entry` and
#'   `quality` changes. Entry support follows the concrete model registry.
#' @return The advanced model, with any supplied environment changes applied.
#' @export
setGeneric("advance_market", function(model, step = NULL) {
    standardGeneric("advance_market")
})

#' @rdname advance_market
#' @export
setMethod("advance_market", "Antitrust", function(model, step = NULL) {
    if (!is.null(step)) {
        if (!methods::is(step, "CounterfactualStep")) {
            stop("'step' must be a CounterfactualStep or NULL")
        }
        unsupported <- setdiff(names(step@changes), c("entry", "quality"))
        if (length(unsupported)) {
            stop("advance_market() accepts only entry and quality changes; simulate ",
                 paste(unsupported, collapse = ", "), " separately")
        }
    }
    if (methods::is(model, "Cournot")) {
        stop("advance_market() requires a model-specific method for quantity-dependent costs")
    }
    if (methods::is(model, "Auction2ndLogit")) {
        if (!identical(class(model)[[1L]], "Auction2ndLogit")) {
            stop("advance_market() requires a model-specific method for this second-score descendant")
        }
        ## .promote_post_to_pre() carries the cost-relative winning index.
    }
    model <- .promote_post_to_pre(model, step)
    model@mcDelta <- rep(0, length(model@mcPre))
    model <- initialize_cost_state(model)
    if (!is.null(step)) model <- .apply_step_environment(model, step)
    model
})
