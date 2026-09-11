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

## Entry and quality are implemented for the following concrete descendants.
## The generic Logit/CES methods handle the demand state, while the entrant
## expansion below carries each descendant's additional product primitive.
## BLP remains excluded because its state is not a flat product vector and
## its simulation methods have separate contracts.
.entry_supported_classes <- c(
    "Logit", "LogitCournot", "CES", "CESCournot",
    "MonComLogit", "MonComCES", "LogitNests", "CESNests",
    "LogitCap", "Auction2ndLogit", "BargainingLogit",
    "Bargaining2ndLogit", "Auction2ndCES", "BargainingCES",
    "LogitALM", "LogitCournotALM", "CESALM", "CESCournotALM",
    "LogitNestsALM", "LogitCapALM", "Auction2ndLogitALM",
    "Auction2ndCESALM", "BargainingLogitALM", "BargainingCESALM",
    "Bargaining2ndCES"
)

.quality_supported_classes <- .entry_supported_classes

## The legacy specialized Logit constructors store a positive demand index
## in `meanval`, even though the flat Logit constructor stores utility.  Keep
## their historical multiplicative quality convention while preserving the
## utility shift for flat Logit and its Cournot/MonCom descendants.
.quality_multiplicative_logit_classes <- c(
    "LogitNests", "LogitNestsALM", "LogitCap", "LogitCapALM",
    "Auction2ndLogit", "Auction2ndLogitALM",
    "BargainingLogit", "Bargaining2ndLogit", "BargainingLogitALM",
    "LogitALM", "LogitCournotALM"
)

.require_supported_transition <- function(model, supported, action) {
    if (!(class(model)[[1L]] %in% supported)) {
        stop("'", action, "' is only supported for models of exact class ",
             paste(supported, collapse = ", "),
             "; this fit is class '", class(model)[[1L]], "'")
    }
    invisible(model)
}

.require_entry_supported <- function(model, action = "entry") {
    .require_supported_transition(model, .entry_supported_classes, action)
}

.require_quality_supported <- function(model, action = "quality") {
    .require_supported_transition(model, .quality_supported_classes, action)
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

    ## Nested demand carries product-to-nest membership and a nest-level
    ## curvature vector.  A nested entrant must state its nest explicitly;
    ## a new singleton has the family-specific limiting curvature.
    if (methods::is(model, "LogitNests") || methods::is(model, "CESNests")) {
        nest <- entrant@extras$nest
        if (is.null(nest) || length(nest) != 1L || is.na(nest) ||
            !nzchar(as.character(nest))) {
            stop("entry into nested models requires entrant extra 'nest'")
        }
        nest <- as.character(nest)
        old_nests <- as.character(model@nests)
        model@nests <- factor(c(old_nests, nest),
                              levels = unique(c(levels(model@nests), nest)))
        sigma <- model@slopes$sigma
        if (is.null(sigma)) sigma <- numeric()
        if (!(nest %in% names(sigma))) {
            singleton_sigma <- if (methods::is(model, "CESNests")) 0 else 1
            sigma <- c(sigma, stats::setNames(singleton_sigma, nest))
        }
        model@slopes$sigma <- sigma
    }

    ## Capacity-constrained Logit needs a capacity primitive for the new
    ## product.  Do not infer it from a synthetic margin or price.
    if (methods::is(model, "LogitCap")) {
        capacity <- entrant@extras$capacity
        if (is.null(capacity) || length(capacity) != 1L ||
            !is.numeric(capacity) || is.na(capacity) || capacity < 0 ||
            (is.infinite(capacity) && capacity < 0)) {
            stop("entry into LogitCap requires entrant extra 'capacity'")
        }
        model@capacitiesPre <- c(model@capacitiesPre, as.numeric(capacity))
        model@capacitiesPost <- c(model@capacitiesPost, as.numeric(capacity))
    }

    ## Bargaining classes require a product-level bargaining-power primitive;
    ## the historical entry default is 0.5 and an entrant may override it.
    ## The concrete Logit/CES bargaining classes contain their own bargaining
    ## slots rather than inheriting the virtual `Bargaining` container.
    if (all(c("bargpowerPre", "bargpowerPost") %in% methods::slotNames(model))) {
        bargpower <- entrant@extras$bargpower
        if (is.null(bargpower)) bargpower <- .5
        if (length(bargpower) != 1L || !is.numeric(bargpower) ||
            !is.finite(bargpower) || bargpower < 0 || bargpower > 1) {
            stop("entrant extra 'bargpower' must be a finite number in [0, 1]")
        }
        model@bargpowerPre <- c(model@bargpowerPre, as.numeric(bargpower))
        model@bargpowerPost <- c(model@bargpowerPost, as.numeric(bargpower))
    }

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
    .require_quality_supported(model, "quality")
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
    if (class(model)[[1L]] %in% .quality_multiplicative_logit_classes) {
        ## These specialized legacy constructors expose a positive demand
        ## index under the historical meanval name.
        meanval[idx] <- meanval[idx] * (1 + quality)
    } else {
        ## Flat Logit mean values are utilities and are only identified up to
        ## a normalization.  A proportional change in positive choice weight
        ## is represented by an additive log(1 + quality) utility shift;
        ## multiplying a possibly negative or zero utility can reverse the
        ## intended attractiveness change.
        meanval[idx] <- meanval[idx] + log1p(quality)
    }
    model@slopes$meanval <- meanval
    model
})

## CES mean values enter as positive multiplicative demand weights, so retain
## the proportional convention directly for that demand family.
setMethod(".apply_quality", "CES", function(model, quality) {
    .require_quality_supported(model, "quality")
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
