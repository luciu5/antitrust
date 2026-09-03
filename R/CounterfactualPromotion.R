## Class-aware state-transition primitives for sequential counterfactuals.
##
## The legacy S4 economic model is the only state a Counterfactual path
## carries between steps.  `.promote_post_to_pre()` turns a solved
## post-counterfactual legacy result into the "pre" state for the next
## step (the Markov transition); `.expand_entrant()` grows a model's
## product dimension for a new single-product entrant; `.apply_quality()`
## multiplies calibrated `meanval` by a product-level percentage shock.
## These never call calibrate() and never touch structural parameters
## (alpha/gamma/nests/...), only Pre/Post-paired economic state.

setGeneric(".promote_post_to_pre", function(model, step) {
    standardGeneric(".promote_post_to_pre")
})

## Default: every "*Pre"/"*Post" slot pair on the object is promoted by
## copying the Post value into Pre. `mcDelta` is deliberately NOT reset here:
## calcMC() for Bertrand- and Cournot-family models always recomputes the
## pre-shock marginal cost fresh from calibration constants (observed
## prices/margins, or cost functions) and treats `mcDelta` as a one-shot
## multiplier relative to that recomputed baseline -- it is never read back
## from a promoted `mcPre`. So the *cumulative* proportional cost change is
## the real persistent cost-environment state, and it must be compounded
## across steps (see `.compound_costs()` in ModelArchitecture.R), not reset
## to zero. `subset` (the exit mask) is intentionally left as-is by this
## default method; sequential exit persistence is handled explicitly by the
## simulate() loop so that later steps can extend the mask when new
## products are exited or entered.
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

## The four leaf classes verified safe for entry/quality: bare Logit, CES,
## and their Cournot-conduct counterparts (which add no new slots beyond
## Logit/CES -- conduct is dispatched separately by simulate(), not by
## class structure). Every other Logit/CES descendant (LogitCap, LogitNests,
## LogitBLP, Auction2ndLogit*, Bargaining*, VertBarg*) is deliberately
## excluded until each is individually audited.
.entry_supported_classes <- c("Logit", "LogitCournot", "CES", "CESCournot")

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

## calcMC() for the Bertrand family never reads mcPre/mcPost directly -- it
## always DERIVES the pre-merger marginal cost from FOC-consistency with
## the calibrated (margin, price, ownership) triple, ignoring whatever is
## assigned to those slots directly. The entrant has no calibrated margin,
## so its cost primitive cannot be set that way; the caller
## (.apply_step_environment(), ModelArchitecture.R) instead computes the
## proportional mcDelta needed so that calcMC()'s implied post-entry cost
## equals entrant@cost, and threads it through the same mcDelta channel
## used for ordinary cost counterfactuals so it is never silently
## overwritten by the step's cost-resolution/compounding logic.
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
    model@mcPre <- c(model@mcPre, NA_real_)
    model@mcPost <- c(model@mcPost, NA_real_)
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

## The proportional mcDelta needed so that calcMC()'s FOC-implied marginal
## cost at the entrant's position equals entrant@cost, evaluated using the
## model's current ownership/elasticity structure (so it stays correct
## whether entry happens against the original baseline or a promoted
## post-merger state).
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

## Multiply calibrated `meanval` by (1 + quality) for the named products.
## `meanval` is a structural (not Pre/Post-paired) slot, so this shock
## persists automatically across `.promote_post_to_pre()` and compounds
## across sequential quality steps for free.
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
    meanval <- model@slopes$meanval
    idx <- match(names(quality), model@labels)
    meanval[idx] <- meanval[idx] * (1 + quality)
    model@slopes$meanval <- meanval
    model
})
