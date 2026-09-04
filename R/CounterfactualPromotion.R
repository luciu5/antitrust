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

## Quality and entry are supported for the Logit/CES demand family broadly:
## any class that is-a `Logit` and carries a list `@slopes` with a `meanval`
## element (the per-product mean-valuation vector a quality shock scales,
## and the demand primitive an entrant supplies).  This covers bare
## Logit/CES, their Cournot and ALM variants, nested Logit/CES, LogitCap,
## and the second-score-auction and bargaining conduct wrappers -- all
## verified to solve cleanly.  Two families are deliberately excluded:
##   * BLP (`LogitBLP`/`CournotBLP`): its `@slopes` also carries random-
##     coefficient structure (`sigma`, `sigmaNest`), so scaling the mean
##     valuation alone is a different, unverified experiment.
##   * Vertical bargaining (`VertBargBertLogit*`): product-dimensional state
##     lives in `@up`/`@down` sub-objects, so a container-level meanval shock
##     or single-product entrant is ill-defined.
.quality_entry_supported <- function(model) {
    methods::is(model, "Logit") &&
        !methods::is(model, "LogitBLP") &&
        !methods::is(model, "VertBargBertLogit") &&
        is.list(model@slopes) &&
        "meanval" %in% names(model@slopes)
}

.require_quality_entry_supported <- function(model, action = "quality") {
    if (!.quality_entry_supported(model)) {
        stop("'", action, "' is only supported for Logit/CES-family models ",
             "(excluding BLP and vertical bargaining); this fit is class '",
             class(model)[[1L]], "'")
    }
    invisible(model)
}

## Entry additionally requires that the target's product-dimensional slots
## beyond the base Logit/CES set be extensible from primitives an Entrant
## can plausibly supply.  Verified needs beyond label/meanval/cost/
## priceStart:
##   * LogitNests / CESNests / *NestsALM: a nest assignment (`extras$nest`)
##     -- required, no safe default (there is no economically neutral nest
##     to assign; erroring is correct).
##   * LogitCap / LogitCapALM: a capacity (`extras$capacity`) -- required;
##     omitting it leaves `capacitiesPre/Post` shorter than the other
##     product vectors, which breaks the equilibrium solver's dimension
##     checks outright (verified: `Length of fn result <> length of x!`).
##   * Bargaining* / Bargaining2nd*: a bargaining power (`extras$bargpower`)
##     in [0,1] -- optional, defaulting to 0.5, the same default the
##     package's own `bargaining.logit()`/`bargaining.ces()` constructors
##     use for every product, so this is not an invented convention.
##   * Auction2ndLogit / Auction2ndCES / Bargaining2nd*: no extra primitive;
##     appending at the end never disturbs `normIndex`, which is a position
##     into the *existing* product vector and is left untouched by growing
##     the vector at the end.
.require_entrant_extra <- function(entrant, key, action) {
    value <- entrant@extras[[key]]
    if (is.null(value)) {
        stop("entry into this model requires an entrant '", key,
             "' primitive (supply it via entrant(..., ", key, " = ...)); ",
             "none was given for entrant '", entrant@label, "'")
    }
    value
}

## Some legacy calcMC() methods apply the post-counterfactual cost change as
## an ADDITIVE level wedge (`mc + mcDelta`) rather than the multiplicative
## wedge (`mc * (1 + mcDelta)`) used by the Bertrand/Cournot/bargaining
## families.  This distinction determines how sequential cost shocks
## compound and how an entrant's cost primitive is converted into an
## mcDelta.  Only the second-score-auction family (Auction2ndLogit and its
## descendants Auction2ndCES / Bargaining2ndLogit / Bargaining2ndCES) is
## additive; everything else is multiplicative.
setGeneric(".mc_delta_is_additive", function(model) {
    standardGeneric(".mc_delta_is_additive")
})
setMethod(".mc_delta_is_additive", "ANY", function(model) FALSE)
setMethod(".mc_delta_is_additive", "Auction2ndLogit", function(model) TRUE)
setMethod(".mc_delta_is_additive", "Auction2ndCES", function(model) TRUE)

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
## mcDelta needed so that calcMC()'s implied post-entry cost equals
## entrant@cost (via .entrant_cost_delta(), additive or multiplicative per
## .mc_delta_is_additive()), and threads it through the same mcDelta
## channel used for ordinary cost counterfactuals so it is never silently
## overwritten by the step's cost-resolution/compounding logic.
setMethod(".expand_entrant", "Logit", function(model, entrant) {
    .require_quality_entry_supported(model, "entry")
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

    ## LogitCap / LogitCapALM: capacitiesPre/Post are additional product-
    ## dimensional slots with no safe default (verified: omitting them
    ## breaks the equilibrium solver's dimension checks outright), so a
    ## capacity primitive is required from the entrant.
    if (methods::.hasSlot(model, "capacitiesPre")) {
        capacity <- .require_entrant_extra(entrant, "capacity", "entry")
        if (length(capacity) != 1L || !is.finite(capacity) || capacity < 0) {
            stop("entrant 'capacity' must be a single non-negative finite number")
        }
        model@capacitiesPre <- c(model@capacitiesPre, capacity)
        model@capacitiesPost <- c(model@capacitiesPost, capacity)
    }

    ## Bargaining / Bargaining2nd: bargpowerPre/Post are additional product-
    ## dimensional slots. Default to 0.5 -- the same default the package's
    ## own bargaining.logit()/bargaining.ces() constructors use for every
    ## product -- rather than requiring it, since 0.5 (equal bargaining
    ## power) is not an invented economic assumption.
    if (methods::.hasSlot(model, "bargpowerPre")) {
        bargpower <- entrant@extras$bargpower
        if (is.null(bargpower)) bargpower <- 0.5
        if (length(bargpower) != 1L || !is.finite(bargpower) ||
            bargpower < 0 || bargpower > 1) {
            stop("entrant 'bargpower' must be a single number in [0, 1]")
        }
        model@bargpowerPre <- c(model@bargpowerPre, bargpower)
        model@bargpowerPost <- c(model@bargpowerPost, bargpower)
    }

    ## LogitNests / CESNests (and their ALM variants): nests is a factor,
    ## not a plain vector, so appending requires rebuilding it with the
    ## entrant's nest as a valid level. There is no safe default nest, so
    ## this primitive is required. Joining an existing nest reuses that
    ## nest's calibrated sigma; forming a new nest adds a new level, and
    ## the existing single-product-nest normalization (sigma = 1) already
    ## used by calibration for singleton nests applies to it unchanged --
    ## no sigma re-derivation is performed here.
    if (methods::.hasSlot(model, "nests")) {
        nest <- .require_entrant_extra(entrant, "nest", "entry")
        if (length(nest) != 1L || !is.character(nest) || is.na(nest) || !nzchar(nest)) {
            stop("entrant 'nest' must be a single non-empty character value")
        }
        new_levels <- union(levels(model@nests), nest)
        model@nests <- factor(c(as.character(model@nests), nest), levels = new_levels)
        if (!(nest %in% names(model@slopes$sigma))) {
            model@slopes$sigma <- c(model@slopes$sigma, stats::setNames(1, nest))
        }
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

## The mcDelta needed so that calcMC()'s FOC-implied marginal cost at the
## entrant's position equals entrant@cost, evaluated using the model's
## current ownership/elasticity structure (so it stays correct whether
## entry happens against the original baseline or a promoted post-merger
## state). Additive families (second-score auction and its descendants)
## need a level difference; every other family needs a ratio, since
## calcMC() applies mcDelta multiplicatively for them.
setGeneric(".entrant_cost_delta", function(model, entrant) {
    standardGeneric(".entrant_cost_delta")
})

setMethod(".entrant_cost_delta", "Logit", function(model, entrant) {
    idx <- match(entrant@label, model@labels)
    implied_mc_pre <- calcMC(model, preMerger = TRUE)[idx]
    if (!is.finite(implied_mc_pre)) {
        stop("entry could not derive a consistent pre-entry marginal cost for entrant '",
             entrant@label, "'; check 'meanval' and 'priceStart'")
    }
    if (.mc_delta_is_additive(model)) {
        return(entrant@cost - implied_mc_pre)
    }
    if (implied_mc_pre <= 0) {
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
    .require_quality_entry_supported(model, "quality")
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
