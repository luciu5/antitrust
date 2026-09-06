#' Realize a synthetic market under a model-specific specification
#'
#' The generic is owned by antitrust because synthetic-market realization
#' requires an antitrust model specification. Model-specific methods own all
#' equations, parameter inversion, cost recovery, and equilibrium diagnostics.
#'
#' @param market A `SyntheticMarket` object.
#' @param spec A model-specific specification or lightweight model descriptor.
#' @param ... Additional model-specific arguments.
#' @return A model-specific realized market or fit object.
#' @export
realize_market <- function(market, spec, ...) UseMethod("realize_market")

#' @export
realize_market.default <- function(market, spec, ...) {
    stop("no model-specific realization adapter is registered for this object")
}
