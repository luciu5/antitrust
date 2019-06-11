#' @import shiny
#' @title A Shiny Interface to the Antitrust Package
#' @name antitrust_shiny
#' @description A Shiny Interface to the Antitrust Package
#' @details \code{antitrust_shiny} launches a shiny interface for
#' the antitrust package. The shiny interface provides users with
#' the ability to calibrate and simulate horizontal mergers using
#' many of the supply and demand models included in the antitrust package.
#' @author Charles Taragin \email{ctaragin@ftc.gov}
#' @examples ## launch shiny interface to antitrust package
#' \dontrun{
#'   antitrust_shiny()
#' }
#' @include SimFunctions.R
NULL

## Function to launch shiny app.
#'@rdname antitrust_shiny
#'@export
antitrust_shiny <- function() {
  requireNamespace("rhandsontable")
  shiny::runApp(system.file('antitrust_shiny', package='antitrust'))
}
