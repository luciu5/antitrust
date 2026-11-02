#' Print a calibrated structural model
#'
#' A lightweight, inexpensive summary of a fitted model's identity: the
#' demand system, conduct, variant, underlying legacy model class, and
#' construction route. It reads only stored slots and never triggers
#' simulation, contraction, or optimization.
#'
#' @param object An `AntitrustFit`.
#' @return `object`, invisibly.
#' @export
setMethod("show", "AntitrustFit", function(object) {
    spec <- object@spec
    cat("<AntitrustFit>\n")
    cat("  demand:  ", spec$demand, "\n", sep = "")
    cat("  conduct: ", spec$conduct, "\n", sep = "")
    cat("  variant: ", spec$variant, "\n", sep = "")
    cat("  model:   ", class(object@model)[[1]], "\n", sep = "")
    route <- object@diagnostics$route
    cat("  route:   ", if (is.null(route)) "calibrate" else route, "\n", sep = "")
    status <- object@diagnostics$status
    if (!is.null(status)) cat("  status:  ", status, "\n", sep = "")
    invisible(object)
})


#' Summarize a calibrated structural model
#'
#' Returns the stored specification, model class, and diagnostics of an
#' `AntitrustFit` without recomputing anything.
#'
#' @param object An `AntitrustFit`.
#' @param ... Currently unused.
#' @return A list with `spec`, `model_class`, and `diagnostics` elements.
#' @export
setMethod("summary", "AntitrustFit", function(object, ...) {
    list(
        spec = object@spec,
        model_class = class(object@model)[[1]],
        diagnostics = object@diagnostics
    )
})


#' @rdname respecify
#' @export
setMethod("respecify", "StructuralFit", function(object, ...) {
    stop("no respecify() method is defined for objects of class '",
         class(object)[[1]], "'.")
})


#' @rdname simulate
#' @export
setMethod("simulate", "StructuralFit", function(object, ...) {
    stop("no simulate() method is defined for objects of class '",
         class(object)[[1]], "'.")
})


#' @rdname Diagnostics-Methods
#' @export
setMethod("calcDiagnostics", "AntitrustFit", function(object, ...) {
    calcDiagnostics(object@model, ...)
})
