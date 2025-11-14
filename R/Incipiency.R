#' Incipiency Analysis for Bertrand Markets
#'
#' Computes post-estimation incipiency metrics for a Bertrand market object. 
#' Measures how close the market is to monopoly (closeness) and how resilient 
#' it is to additional perturbations (resilience) using hypothetical monopolist 
#' prices and/or consumer welfare.
#'
#' @param object A \code{Bertrand} object (post-estimation).
#' @param outcome Character; either "welfare" or "prices".
#' @param subset Optional vector of product indices to override \code{object@subset}.
#' @param coalition Optional vector of product indices to form the hypothetical monopolist. Default is all products in the subset.
#' @param scenario Character; scenario to implement. Must be one of 
#'   \code{"merger_only"}, \code{"largest_exit"}, \code{"all_other_merge"}.
#'
#' @return A list with components:
#'   \item{ratioOutcome}{Ratio of outcome under scenario to outcome under hypothetical monopoly.}
#'   \item{deltaOutcomeRatio}{Ratio of the change in outcome due to scenario to change under monopoly.}
#'   \item{outcomeMonopoly}{Prices or welfare under hypothetical monopolist.}
#'
#' @export
setMethod(
  f = "incipiency",
  signature = "Bertrand",
  definition = function(object, 
                        outcome = c("welfare", "prices"),
                        subset = NULL,
                        coalition = NULL,
                        scenario = c("merger_only", "largest_exit", "all_other_merge")) {
    
    # Validate arguments
    outcome <- match.arg(outcome)
    scenario <- match.arg(scenario)
    
    # Determine products to include
    useSubset <- if (!is.null(subset)) subset else object@subset
    
    # Default coalition: all products in subset
    if (is.null(coalition)) coalition <- useSubset
    monPrices <- calcPricesHypoMon(object, coalition = coalition)
    
    # Make a copy of object for scenario manipulations
    objScenario <- object
    
    # Determine scenario adjustments
    switch(scenario,
           merger_only = {
             # no change
           },
           largest_exit = {
             largestFirm <- which.max(object@sharesPre[useSubset])
             objScenario@subset <- setdiff(useSubset, largestFirm)
           },
           all_other_merge = {
             objScenario@ownerPost[useSubset] <- objScenario@ownerPost[useSubset[1]]
           })
    
    # Compute the outcome metric
    if (outcome == "prices") {
      outcomePre <- object@pricesPre[useSubset]
      objScenario@pricesPre <- objScenario@pricesPost  # scenario prices
      outcomeScenario <- objScenario@pricesPre[useSubset]
      outcomeMonopoly <- monPrices
    } else if (outcome == "welfare") {
      # For welfare, override pricePre in copy
      objScenario@pricesPre <- objScenario@pricesPost
      outcomePre <- CV(object)
      outcomeScenario <- CV(objScenario)
      # Monopoly welfare: override pricePre with monPrices
      objMon <- object
      objMon@pricesPre <- monPrices
      outcomeMonopoly <- CV(objMon)
    }
    
    ratioOutcome <- outcomeScenario / outcomeMonopoly
    deltaOutcomeRatio <- (outcomeScenario - outcomePre) / (outcomeMonopoly - outcomePre)
    
    list(
      ratioOutcome = ratioOutcome,
      deltaOutcomeRatio = deltaOutcomeRatio,
      outcomeMonopoly = outcomeMonopoly
    )
  }
)
