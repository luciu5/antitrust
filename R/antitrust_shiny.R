## Function to launch shiny app.
antitrust_shiny <- function() {
  requireNamespace("rhandsontable")
  shiny::runApp(system.file('antitrust_shiny', package='antitrust'))
  }