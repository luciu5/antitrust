utils::globalVariables(c(".data"))
#' @title Methods For Calculating Upwards Pricing Pressure Index (Bertrand)
#' @name Plot-Methods
#' @docType methods
#'
#' @aliases plot,Bertrand-method
#'
#' @description Use \code{\link[ggplot2]{ggplot}} to plot pre- and post-merger demand, marginal cost and equilibria.
#' \sQuote{scale} controls the amount above marginal cost and below equilbrium price that is plotted.
#'
#' @param x Used only in plot method. Should always be set equal to object.
#' @param scale The proportion below marginal cost and above equilbrium price that should be plotted. Default is .1.
#'
#' @include OutputMethods.R
#' @keywords methods
NULL

##plot method
#'@rdname Plot-Methods
#'@export
setMethod(
  f = "plot",
  signature = "Bertrand",
  definition = function(x, scale = .1) {
    object = x
    nprods = length(object@labels)
    pricePre = object@pricePre
    pricePost = object@pricePost
    
    if (all(is.na(pricePre)) || all(is.na(pricePost))) {
      stop("'pricePre' or 'pricePost' are all NA")
    }
    
    mcPre = object@mcPre
    mcPost = object@mcPost
    labels = object@labels
    isParty <- rowSums(abs(object@ownerPost - object@ownerPre)) > 0
    isParty <- ifelse(isParty, "*", "")
    labels  <- paste(isParty, labels, sep = "")
    
    missPrices <- all(is.na(object@pricePre))
    
    if (missPrices) {
      outPre  = calcShares(object, preMerger = TRUE)
      outPost = calcShares(object, preMerger = FALSE)
    } else {
      outPre  = calcQuantities(object, preMerger = TRUE)
      outPost = calcQuantities(object, preMerger = FALSE)
    }
    
    prices <- quantPre <- quantPost <- prod <- NULL
    
    plotThis <- function(price, idx, preMerger) {
      thisObj = object
      if (preMerger) {
        thisObj@pricePre[idx] = price
      } else {
        thisObj@pricePost[idx] = price
      }
      
      if (!missPrices) {
        return(calcQuantities(thisObj, preMerger = preMerger)[idx])
      } else {
        return(calcShares(thisObj, preMerger = preMerger)[idx])
      }
    }
    
    for (i in 1:nprods) {
      thesePrices = seq(
        (1 - scale) * min(mcPre[i], mcPost[i], na.rm = TRUE),
        (1 + scale) * max(pricePre[i], pricePost[i], na.rm = TRUE),
        length.out = 100
      )
      quantPre  = c(quantPre, sapply(thesePrices, plotThis, idx = i, preMerger = TRUE))
      quantPost = c(quantPost, sapply(thesePrices, plotThis, idx = i, preMerger = FALSE))
      prices    = c(prices, thesePrices)
      prod      = c(prod, rep(labels[i], length(thesePrices)))
    }
    
    resultPre = data.frame(output = quantPre, price = prices, prod = prod, Demand = "pre-merger")
    resultPost = data.frame(output = quantPost, price = prices, prod = prod, Demand = "post-merger")
    result = rbind(resultPre, resultPost)
    result = result[result$output > 0, ]
    
    equilibria = data.frame(
      output = c(outPre, outPost),
      price  = c(pricePre, pricePost),
      mc     = c(mcPre, mcPost),
      prod   = labels,
      Demand = rep(c("pre-merger", "post-merger"), each = nprods)
    )
    equilibria$Cost = equilibria$Demand
    
    thisPlot = ggplot2::ggplot(
      result,
      ggplot2::aes(
        x = .data$output,
        y = .data$price,
        color = .data$Demand,
        group = .data$Demand
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = ggplot2::element_blank()
      ) +
      ggplot2::facet_wrap(~prod, scales = "free_x") +
      ggplot2::geom_vline(
        data = equilibria[, c("output", "Demand", "prod")],
        ggplot2::aes(
          xintercept = .data$output,
          group = .data$Demand,
          colour = .data$Demand
        ),
        linetype = 3
      ) +
      ggplot2::geom_hline(
        data = equilibria[, c("price", "Demand", "prod")],
        ggplot2::aes(
          yintercept = .data$price,
          group = .data$Demand,
          colour = .data$Demand
        ),
        linetype = 3
      ) +
      ggplot2::geom_point(
        data = equilibria,
        ggplot2::aes(
          x = .data$output,
          y = .data$price,
          color = .data$Demand,
          group = .data$Demand
        )
      ) +
      ggplot2::geom_hline(
        data = equilibria[, c("mc", "Cost", "prod")],
        ggplot2::aes(
          yintercept = .data$mc,
          group = .data$Cost,
          color = .data$Cost
        ),
        show.legend = TRUE
      )
    
    return(thisPlot)
  }
)
