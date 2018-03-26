## Function to launch shiny app.
antitrust_shiny <- function() {
  shiny::runApp(system.file('antitrust_shiny', package='antitrust'))
  }