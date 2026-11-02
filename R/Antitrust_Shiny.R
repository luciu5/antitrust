
#' @title A Shiny Interface to the Antitrust Package
#' @description A Shiny Interface to the Antitrust Package
#' @details \code{antitrust_shiny} launches a shiny interface for
#' the antitrust package. The shiny interface provides users with
#' the ability to calibrate and simulate horizontal mergers using
#' many of the supply and demand models included in the antitrust package.
#' @name antitrust_shiny-deprecated
#' @seealso \code{\link{antitrust-deprecated}}
#' @keywords internal
NULL

#' @rdname antitrust-deprecated
#' @section \code{antitrust_shiny}:
#' For \code{antitrust_shiny}, use \code{\link[competitiontoolbox]{ct_shiny}}.
#'
#'
#' @export
antitrust_shiny <- function() {
  .Deprecated("ct_shiny' in the 'competitiontoolbox' package")
  # requireNamespace("rhandsontable")
  # shiny::runApp(system.file('antitrust_shiny', package='antitrust'))
}

