## Class-aware state-transition primitives for sequential counterfactuals.
##
## The legacy S4 economic model is the only state a Counterfactual path
## carries between steps.  `.promote_post_to_pre()` turns a solved
## post-counterfactual legacy result into the "pre" state for the next
## step (the Markov transition); `.expand_entrant()` grows a model's
## product dimension for a new single-product entrant; `.apply_quality()`
## applies a proportional attractiveness shock under each demand family's
## normalization convention.
## These never call calibrate() and never touch structural parameters
## (alpha/gamma/nests/...), except for an explicitly requested quality shock.

setGeneric(".promote_post_to_pre", function(model, step) {
    standardGeneric(".promote_post_to_pre")
})

## Default: every "*Pre"/"*Post" slot pair on the object is promoted by
## copying the Post value into Pre.  Structural cost state is stored as an
## ordinary attribute by the architecture layer and therefore travels with
## the copied S4 value; `mcDelta` remains the current cumulative shock supplied
## by the simulation loop.  `subset` (the exit mask) is intentionally left as
## is by this default method; sequential exit persistence is handled explicitly
## by the simulate() loop so later steps can extend the mask when products are
## exited or entered.
setMethod(".promote_post_to_pre", "ANY", function(model, step) {
    slots <- methods::slotNames(model)
    post_slots <- grep("Post$", slots, value = TRUE)
    for (post_slot in post_slots) {
        pre_slot <- sub("Post$", "Pre", post_slot)
        if (pre_slot %in% slots) {
            methods::slot(model, pre_slot) <- methods::slot(model, post_slot)
        }
    }
    model
})

## Vertical bargaining stores its economic state in `up`/`down` sub-objects
## rather than directly on the container; promote each side independently.
setMethod(".promote_post_to_pre", "VertBargBertLogit", function(model, step) {
    model@up <- .promote_post_to_pre(model@up, step)
    model@down <- .promote_post_to_pre(model@down, step)
    model
})

## The flat Logit/CES classes verified safe for entry/quality include the
## monopolistic-competition and Cournot descendants (which add no product
## state beyond Logit/CES; conduct is dispatched separately). Every other
## Logit/CES descendant (LogitCap, LogitNests, LogitBLP, Auction2ndLogit*,
## Bargaining*, VertBarg*) is deliberately excluded until individually
## audited.
.entry_supported_classes <- c(
    "Logit", "LogitCournot", "CES", "CESCournot",
    "MonComLogit", "MonComCES"
)

.require_entry_supported <- function(model, action = "entry") {
    if (!(class(model)[[1L]] %in% .entry_supported_classes)) {
        stop("'", action, "' is only supported for models of exact class ",
             paste(.entry_supported_classes, collapse = ", "),
             "; this fit is class '", class(model)[[1L]], "'")
    }
    invisible(model)
}

## Expand a Logit/CES-family model's product dimension for one new
## single-product entrant.  Reuses the existing exported `ownerToMatrix()`
## generic for ownership expansion rather than inventing new ownership
## math; all other product-dimensional slots are appended positionally at
## the end, so pre-existing product indices/labels are never disturbed.
setGeneric(".expand_entrant", function(model, entrant) {
    standardGeneric(".expand_entrant")
})

## A calibrated Bertrand result keeps its inferred constant marginal-cost
## level in the internal cost state.  This avoids re-identifying incumbent
## costs after ownership or demand-state changes.  A direct legacy object
## without that state still follows its historical calcMC() implementation.
## Entry appends the supplied entrant cost to the same state; the entry cost is
## therefore independent of incumbent margins and ownership.
setMethod(".expand_entrant", "Logit", function(model, entrant) {
    .require_entry_supported(model, "entry")
    if (entrant@label %in% model@labels) {
        stop("entrant label '", entrant@label, "' duplicates an existing product label")
    }

    owner_pre <- ownerToMatrix(model, preMerger = TRUE)
    owner_post <- ownerToMatrix(model, preMerger = FALSE)
    n <- nrow(owner_pre)
    expand_owner <- function(owner) {
        expanded <- matrix(0, nrow = n + 1L, ncol = n + 1L)
        expanded[seq_len(n), seq_len(n)] <- owner
        expanded[n + 1L, n + 1L] <- 1
        expanded
    }
    model@ownerPre <- expand_owner(owner_pre)
    model@ownerPost <- expand_owner(owner_post)

    model@labels <- c(model@labels, entrant@label)
    model@pricePre <- c(model@pricePre, entrant@priceStart)
    model@pricePost <- c(model@pricePost, entrant@priceStart)
    ## Entry supplies a new product's marginal-cost primitive directly.  Keep
    ## it in the persistent cost state so incumbents remain tied to their
    ## calibrated costs even though the entrant has no observed margin.
    cost_state <- .cost_state(model)
    if (!is.null(cost_state) && !is.null(cost_state$base)) {
        cost_state$base <- c(cost_state$base, entrant@cost)
        names(cost_state$base) <- model@labels
        model <- .set_cost_state(model, cost_state)
        model@mcPre <- c(model@mcPre, entrant@cost)
        model@mcPost <- c(model@mcPost, entrant@cost)
    } else {
        model@mcPre <- c(model@mcPre, NA_real_)
        model@mcPost <- c(model@mcPost, NA_real_)
    }
    model@mcDelta <- c(model@mcDelta, 0)
    model@subset <- c(model@subset, TRUE)
    model@priceStart <- c(model@priceStart, entrant@priceStart)
    model@prices <- c(model@prices, entrant@priceStart)
    model@margins <- c(model@margins, NA_real_)
    if (length(model@weights)) model@weights <- c(model@weights, 1)
    model@slopes$meanval <- c(model@slopes$meanval, entrant@meanval)
    names(model@slopes$meanval) <- model@labels

    entrant_shares <- calcShares(model, preMerger = TRUE)
    model@shares <- entrant_shares
    model@shareInside <- ifelse(
        isTRUE(all.equal(sum(entrant_shares), 1, check.names = FALSE, tolerance = 1e-3)),
        1, sum(entrant_shares)
    )

    old_diversion <- model@diversion
    expanded_diversion <- matrix(NA_real_, nrow = n + 1L, ncol = n + 1L)
    if (nrow(old_diversion) == n && ncol(old_diversion) == n) {
        expanded_diversion[seq_len(n), seq_len(n)] <- old_diversion
    }
    diag(expanded_diversion) <- -1
    model@diversion <- expanded_diversion

    model
})

## Compatibility helper for direct legacy objects that do not carry the
## internal cost state.  The architecture path uses the supplied entrant cost
## directly whenever that state is available.
setGeneric(".entrant_cost_delta", function(model, entrant) {
    standardGeneric(".entrant_cost_delta")
})

setMethod(".entrant_cost_delta", "Logit", function(model, entrant) {
    idx <- match(entrant@label, model@labels)
    implied_mc_pre <- calcMC(model, preMerger = TRUE)[idx]
    if (!is.finite(implied_mc_pre) || implied_mc_pre <= 0) {
        stop("entry could not derive a consistent pre-entry marginal cost for entrant '",
             entrant@label, "'; check 'meanval' and 'priceStart'")
    }
    entrant@cost / implied_mc_pre - 1
})

## Apply a proportional attractiveness shock for the named products.  The
## demand-specific methods below map that shock into each model's structural
## mean-value normalization.  `meanval` is not a Pre/Post-paired slot, so the
## change persists through promotion and compounds across sequential steps.
setGeneric(".apply_quality", function(model, quality) {
    standardGeneric(".apply_quality")
})

setMethod(".apply_quality", "Logit", function(model, quality) {
    .require_entry_supported(model, "quality")
    if (is.null(names(quality)) || any(!nzchar(names(quality)))) {
        stop("'quality' must be a named numeric vector (product label = proportional change)")
    }
    unknown <- setdiff(names(quality), model@labels)
    if (length(unknown)) {
        stop("'quality' references unknown product label(s): ", paste(unknown, collapse = ", "))
    }
    active_labels <- model@labels[model@subset]
    exited <- setdiff(names(quality), active_labels)
    if (length(exited)) {
        stop("'quality' references product(s) that are not active (exited or not yet entered): ",
             paste(exited, collapse = ", "))
    }
    if (any(quality <= -1)) {
        stop("'quality' values must be greater than -1 for Logit demand")
    }
    meanval <- model@slopes$meanval
    idx <- match(names(quality), model@labels)
    ## Logit mean values are utilities and are only identified up to a
    ## normalization.  A proportional change in positive choice weight is
    ## therefore represented by an additive log(1 + quality) utility shift;
    ## multiplying a possibly negative or zero utility is normalization
    ## dependent and can reverse the intended attractiveness change.
    meanval[idx] <- meanval[idx] + log1p(quality)
    model@slopes$meanval <- meanval
    model
})

## CES mean values enter as positive multiplicative demand weights, so retain
## the proportional convention directly for that demand family.
setMethod(".apply_quality", "CES", function(model, quality) {
    .require_entry_supported(model, "quality")
    if (is.null(names(quality)) || any(!nzchar(names(quality)))) {
        stop("'quality' must be a named numeric vector (product label = proportional change)")
    }
    unknown <- setdiff(names(quality), model@labels)
    if (length(unknown)) {
        stop("'quality' references unknown product label(s): ", paste(unknown, collapse = ", "))
    }
    active_labels <- model@labels[model@subset]
    exited <- setdiff(names(quality), active_labels)
    if (length(exited)) {
        stop("'quality' references product(s) that are not active (exited or not yet entered): ",
             paste(exited, collapse = ", "))
    }
    if (any(!is.finite(quality) | quality <= -1)) {
        stop("'quality' values must be finite and greater than -1")
    }
    meanval <- model@slopes$meanval
    idx <- match(names(quality), model@labels)
    meanval[idx] <- meanval[idx] * (1 + quality)
    model@slopes$meanval <- meanval
    model
})
