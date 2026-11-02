#' @title Deprecated functions in package \pkg{antitrust}.
#' @description The functions listed below are deprecated and will be defunct in
#'   the near future. When possible, alternative functions with similar
#'   functionality are also mentioned. Help pages for deprecated functions are
#'   available at \code{help("-deprecated")}.
#' @name antitrust-deprecated
#' @keywords internal
NULL

#' Vertical bargaining moved to the vertical package
#'
#' @param ... Ignored.
#' @return This function does not return; it raises a migration error.
#' @export
vertical.barg <- function(...) {
  .Defunct("vertical.barg", package = "vertical",
           msg = "vertical bargaining economics moved to vertical::vertical.barg().")
}
