#' The realized sequence of states from a multi-step Counterfactual
#'
#' `CounterfactualPath` records the legacy S4 result produced by every step
#' of a multi-step `Counterfactual`, together with the ordered
#' `CounterfactualStep` objects that produced them and the originally
#' calibrated legacy S4 state the path began from. The state history is the
#' actual legacy S4 result objects; no parallel economic-state
#' representation is kept.
#'
#' @slot initial The originally calibrated legacy S4 model state.
#' @slot steps The ordered `CounterfactualStep` objects that were solved.
#' @slot results The legacy S4 result after each step, in order.
#' @slot diagnostics Bookkeeping needed to resume simulation from this path.
#' @export
#' @exportClass CounterfactualPath
setClass(
    "CounterfactualPath",
    representation(
        initial = "ANY",
        steps = "list",
        results = "list",
        diagnostics = "list"
    ),
    prototype = prototype(steps = list(), results = list(), diagnostics = list())
)

#' The final solved legacy S4 result of a CounterfactualPath
#'
#' @param path A `CounterfactualPath` object.
#' @return The legacy S4 result after the last step.
#' @export
setGeneric("final_result", function(path) standardGeneric("final_result"))

#' @rdname final_result
#' @export
setMethod("final_result", "CounterfactualPath", function(path) {
    n <- length(path@results)
    if (n == 0L) stop("'path' contains no results")
    path@results[[n]]
})

#' The legacy S4 result at a given step of a CounterfactualPath
#'
#' @param path A `CounterfactualPath` object.
#' @param i A single step index.
#' @return The legacy S4 result after step `i`.
#' @export
setGeneric("result_at", function(path, i) standardGeneric("result_at"))

#' @rdname result_at
#' @export
setMethod("result_at", "CounterfactualPath", function(path, i) {
    n <- length(path@results)
    if (length(i) != 1L || !is.numeric(i) || is.na(i) || i < 1 || i > n || i != as.integer(i)) {
        stop("'i' must be a single valid step index between 1 and ", n)
    }
    path@results[[as.integer(i)]]
})

#' @rdname CounterfactualPath-class
#' @param object A `CounterfactualPath` object.
#' @export
setMethod("show", "CounterfactualPath", function(object) {
    cat("CounterfactualPath\n")
    cat("  steps: ", length(object@steps), "\n", sep = "")
    for (i in seq_along(object@steps)) {
        fields <- names(object@steps[[i]]@changes)
        cat("    [", i, "] ", if (length(fields)) paste(fields, collapse = ", ") else "(none)", "\n", sep = "")
    }
    invisible(NULL)
})
