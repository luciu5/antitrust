runShiny <- function() {
  appDir <- system.file("shiny-examples", "atrshiny", package = "antitrust")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `antitrust`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}
