#' @name antitrust-package
#' @aliases antitrust
#' @docType package
#' @title \packageTitle{antitrust}
#' @description \packageDescription{antitrust}
#'
#' @details Choose a complete model with \code{\link{supportedModels}}.
#' Use \code{\link{calibrate}} for observed-market calibration or
#' \code{\link{specify}} when structural parameters are supplied, then
#' \code{\link{counterfactual}} and \code{\link{simulate}} for a scenario.
#' The package vignettes cover workflow, economic models, and extension design.
#' Direct model constructors remain available as compatibility and convenience
#' interfaces. \code{\link{cmcr.bertrand}} and \code{\link{cmcr.cournot}}
#' provide screening measures when a complete market fit is unavailable.
#'
#' \packageDESCRIPTION{antitrust}
#' \packageIndices{antitrust}
#' @author \packageAuthor{antitrust}
#' Maintainer: \packageMaintainer{antitrust}
#' @examples
#' fit <- calibrate("logit", "bertrand", prices = c(2, 2.2, 2.5),
#'   shares = c(.35, .25, .20), margins = c(.40, .35, .30),
#'   ownerPre = c("A", "B", "C"), insideSize = 100)
#' result <- simulate(fit, counterfactual(ownership = c("A", "A", "C")))
#' result@pricePost
#' @include Antitrust_Shiny.R
NULL
