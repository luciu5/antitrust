#' @name antitrust-package
#' @aliases antitrust
#' @docType package
#' @title \packageTitle{antitrust}
#' @description \packageDescription{antitrust}
#'
#' @section Disclaimer
#' The views expressed herein are entirely those of the authors and should
#' not be purported to reflect those of the U.S. Department of Justice.
#' The \code{antitrust} package has been released into the
#' public domain without warranty of any kind, expressed or implied.
#' Address: Economic Analysis Group, Antitrust Division, U.S. Department of
#' Justice, 450 5th St. NW, Washington DC 20530.
#' E-mail: ctaragin@ftc.gov and michael.sandfort@usdoj.gov.
#'
#' @details The DESCRIPTION file:
#' \packageDESCRIPTION{antitrust}
#' \packageIndices{antitrust}
#'
#' To get Started:
#'   \enumerate{
#'     \item Collect data on product prices, shares, margins and diversions
#'     (optional).
#'     \item Specify how firms interact strategically (e.g. Bertrand, Cournot, Auction)
#'     \item If you have data on many/all products in the market consider
#'     calibrating a demand system and simulating a merger with either
#'     a \code{\link{bertrand.alm}},\code{\link{cournot}}, or \code{\link{auction2nd.logit}}.
#'     \item If you only have data on the merging parties' products, consider
#'     using \code{\link{cmcr.bertrand}} or \code{\link{cmcr.cournot}} to
#'     uncover the marginal cost reductions needed to offset a post-merger increase.
#'     }
#' @author \packageAuthor{antitrust}
#' Maintainer: \packageMaintainer{antitrust}
#' @include Antitrust_Shiny.R
