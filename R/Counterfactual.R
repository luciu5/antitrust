#' A single simultaneous set of counterfactual changes
#'
#' `CounterfactualStep` stores one bundle of economic-environment changes to
#' be applied and solved as a single equilibrium. Demand, conduct, and model
#' variants are deliberately not counterfactual fields; use `update()` or
#' `respecify()` for those changes.
#'
#' @slot changes A named list of supplied change fields.
#' @export
#' @exportClass CounterfactualStep
setClass(
    "CounterfactualStep",
    representation(changes = "list"),
    prototype = prototype(changes = list())
)

#' An ordered sequence of counterfactual steps
#'
#' `Counterfactual` stores an ordered sequence of `CounterfactualStep`
#' objects. A `Counterfactual` returned by `counterfactual()` has exactly one
#' step; `add_step()` appends additional steps to be solved sequentially,
#' each starting from the equilibrium produced by the previous step.
#'
#' @slot steps A list of `CounterfactualStep` objects.
#' @export
#' @exportClass Counterfactual
setClass(
    "Counterfactual",
    representation(steps = "list"),
    prototype = prototype(steps = list()),
    validity = function(object) {
        if (!length(object@steps)) {
            return("'steps' must contain at least one CounterfactualStep")
        }
        if (!all(vapply(object@steps, methods::is, logical(1), "CounterfactualStep"))) {
            return("all elements of 'steps' must be CounterfactualStep objects")
        }
        TRUE
    }
)

#' A new single-product entrant firm
#'
#' `Entrant` describes a new single-product firm that may be added to a
#' calibrated Logit- or CES-family market via the `entry` counterfactual
#' field. Market-level structural parameters (e.g. the demand price
#' coefficient) are inherited from the current model; only product-level
#' primitives are supplied here.
#'
#' @slot label A unique product label for the entrant.
#' @slot meanval The entrant's mean valuation, on the same scale as the
#'   fitted model's calibrated `meanval` values.
#' @slot cost The entrant's constant marginal cost.
#' @slot priceStart A starting price used to initialize the equilibrium
#'   solver and to compute the entrant's implied share.
#' @slot extras Reserved for model-specific primitives not yet supported.
#' @export
#' @exportClass Entrant
setClass(
    "Entrant",
    representation(
        label = "character",
        meanval = "numeric",
        cost = "numeric",
        priceStart = "numeric",
        extras = "list"
    ),
    prototype = prototype(extras = list()),
    validity = function(object) {
        if (length(object@label) != 1L || is.na(object@label) || !nzchar(object@label)) {
            return("'label' must be a single non-empty character value")
        }
        if (length(object@meanval) != 1L || !is.finite(object@meanval)) {
            return("'meanval' must be a single finite number")
        }
        if (length(object@cost) != 1L || !is.finite(object@cost) || object@cost < 0) {
            return("'cost' must be a single non-negative finite number")
        }
        if (length(object@priceStart) != 1L || !is.finite(object@priceStart) || object@priceStart < 0) {
            return("'priceStart' must be a single non-negative finite number")
        }
        TRUE
    }
)

#' Define a new single-product entrant firm
#'
#' @param label A unique product label for the entrant.
#' @param meanval The entrant's mean valuation.
#' @param cost The entrant's constant marginal cost.
#' @param priceStart A starting price for the entrant.
#' @param ... Reserved for model-specific primitives; currently unused.
#' @return An `Entrant` object.
#' @export
entrant <- function(label, meanval, cost, priceStart, ...) {
    new("Entrant",
        label = as.character(label), meanval = as.numeric(meanval),
        cost = as.numeric(cost), priceStart = as.numeric(priceStart),
        extras = list(...)
    )
}

.normalize_entry_field <- function(entry) {
    if (is.null(entry)) return(NULL)
    if (methods::is(entry, "Entrant")) return(list(entry))
    if (is.list(entry) && length(entry) &&
        all(vapply(entry, methods::is, logical(1), "Entrant"))) {
        return(entry)
    }
    stop("'entry' must be an Entrant object or a list of Entrant objects")
}

.validate_entry_labels <- function(entrants) {
    if (is.null(entrants)) return(invisible(NULL))
    labels <- vapply(entrants, function(e) e@label, character(1))
    if (anyDuplicated(labels)) {
        stop("'entry' contains duplicate entrant labels: ",
             paste(unique(labels[duplicated(labels)]), collapse = ", "))
    }
    invisible(entrants)
}

.validate_quality_field <- function(quality) {
    if (is.null(quality)) return(invisible(NULL))
    if (!is.numeric(quality) || is.null(names(quality)) || any(!nzchar(names(quality)))) {
        stop("'quality' must be a named numeric vector (product label = proportional change)")
    }
    if (anyDuplicated(names(quality))) {
        stop("'quality' contains duplicate product labels")
    }
    invisible(quality)
}

.build_step <- function(ownership = NULL, costs = NULL, exit = NULL,
                        capacity = NULL, tariff = NULL, quota = NULL,
                        bargaining = NULL, leader = NULL, products = NULL,
                        quality = NULL, entry = NULL) {
    entry <- .normalize_entry_field(entry)
    .validate_entry_labels(entry)
    .validate_quality_field(quality)
    changes <- Filter(Negate(is.null), list(
        ownership = ownership, costs = costs, exit = exit,
        capacity = capacity, tariff = tariff, quota = quota,
        bargaining = bargaining, leader = leader, products = products,
        quality = quality, entry = entry
    ))
    new("CounterfactualStep", changes = changes)
}

#' Define a post-calibration economic counterfactual
#'
#' `counterfactual()` stores only supplied changes to the economic
#' environment, bundled as the single `CounterfactualStep` of a
#' `Counterfactual`. Demand, conduct, and model variants are deliberately
#' not counterfactual fields; use `update()` or `respecify()` for those
#' changes. Use `add_step()` to append further, sequential steps.
#'
#' @param ownership Post-counterfactual ownership.
#' @param costs Post-counterfactual marginal-cost changes.
#' @param exit Products to remove, or a logical active-product vector.
#' @param capacity Post-counterfactual capacities.
#' @param tariff A tariff change, unsupported by antitrust models.
#' @param quota A quota change, unsupported by antitrust models.
#' @param bargaining Post-counterfactual bargaining parameters.
#' @param leader Post-counterfactual leader indicators for Stackelberg models.
#' @param products Post-counterfactual product structure for Stackelberg models.
#' @param quality A named numeric vector of proportional changes to
#'   calibrated `meanval`, keyed by product label. `meanval_new <-
#'   meanval_current * (1 + quality)`.
#' @param entry An `Entrant` object, or a list of `Entrant` objects, each
#'   describing a new single-product firm to add to the market.
#' @param ... Reserved; model specification fields are rejected.
#' @return A `Counterfactual` object with exactly one `CounterfactualStep`.
#' @export
counterfactual <- function(ownership = NULL, costs = NULL, exit = NULL,
                           capacity = NULL, tariff = NULL, quota = NULL,
                           bargaining = NULL, leader = NULL, products = NULL,
                           quality = NULL, entry = NULL, ...) {
    extras <- list(...)
    if (length(extras)) {
        stop("counterfactual() accepts economic-environment fields only; use update() or respecify() for model specification changes")
    }
    step <- .build_step(
        ownership = ownership, costs = costs, exit = exit,
        capacity = capacity, tariff = tariff, quota = quota,
        bargaining = bargaining, leader = leader, products = products,
        quality = quality, entry = entry
    )
    new("Counterfactual", steps = list(step))
}

#' Append a sequential counterfactual step
#'
#' `add_step()` appends one new, simultaneous `CounterfactualStep` to a
#' `Counterfactual`. The appended step is solved starting from the
#' equilibrium produced by the immediately preceding step (or, for the
#' first appended step, the fitted baseline); it never replaces or
#' recalibrates prior steps.
#'
#' @param object A `Counterfactual` object.
#' @param ownership Post-counterfactual ownership.
#' @param costs Post-counterfactual marginal-cost changes.
#' @param exit Products to remove, or a logical active-product vector.
#' @param capacity Post-counterfactual capacities.
#' @param tariff A tariff change, unsupported by antitrust models.
#' @param quota A quota change, unsupported by antitrust models.
#' @param bargaining Post-counterfactual bargaining parameters.
#' @param leader Post-counterfactual leader indicators for Stackelberg models.
#' @param products Post-counterfactual product structure for Stackelberg models.
#' @param quality A named numeric vector of proportional changes to
#'   calibrated `meanval`, keyed by product label.
#' @param entry An `Entrant` object, or a list of `Entrant` objects.
#' @param ... Reserved; model specification fields are rejected.
#' @return A `Counterfactual` object with the new step appended.
#' @export
setGeneric("add_step", function(object, ...) standardGeneric("add_step"))

#' @rdname add_step
#' @export
setMethod("add_step", "Counterfactual", function(object, ownership = NULL,
                                                  costs = NULL, exit = NULL,
                                                  capacity = NULL, tariff = NULL,
                                                  quota = NULL, bargaining = NULL,
                                                  leader = NULL, products = NULL,
                                                  quality = NULL, entry = NULL, ...) {
    extras <- list(...)
    if (length(extras)) {
        stop("add_step() accepts economic-environment fields only; use update() or respecify() for model specification changes")
    }
    step <- .build_step(
        ownership = ownership, costs = costs, exit = exit,
        capacity = capacity, tariff = tariff, quota = quota,
        bargaining = bargaining, leader = leader, products = products,
        quality = quality, entry = entry
    )
    object@steps <- c(object@steps, list(step))
    object
})

#' @rdname CounterfactualStep-class
#' @param object A `CounterfactualStep` object.
#' @export
setMethod("show", "CounterfactualStep", function(object) {
    cat("CounterfactualStep\n")
    if (length(object@changes)) {
        cat("  fields: ", paste(names(object@changes), collapse = ", "), "\n", sep = "")
    } else {
        cat("  fields: none\n")
    }
    invisible(NULL)
})

#' @rdname Counterfactual-class
#' @param object A `Counterfactual` object.
#' @export
setMethod("show", "Counterfactual", function(object) {
    cat("Counterfactual\n")
    cat("  steps: ", length(object@steps), "\n", sep = "")
    for (i in seq_along(object@steps)) {
        fields <- names(object@steps[[i]]@changes)
        cat("    [", i, "] ", if (length(fields)) paste(fields, collapse = ", ") else "(none)", "\n", sep = "")
    }
    invisible(NULL)
})

#' @rdname Entrant-class
#' @param object An `Entrant` object.
#' @export
setMethod("show", "Entrant", function(object) {
    cat("Entrant '", object@label, "'\n", sep = "")
    cat("  meanval: ", object@meanval, "\n", sep = "")
    cat("  cost: ", object@cost, "\n", sep = "")
    cat("  priceStart: ", object@priceStart, "\n", sep = "")
    invisible(NULL)
})

#' Combine non-conflicting counterfactual changes into one simultaneous step
#'
#' `combine_counterfactuals()` merges the single steps of several one-step
#' `Counterfactual` objects into one simultaneous `CounterfactualStep`,
#' erroring if any field is supplied with conflicting values. Combining
#' multi-step `Counterfactual` objects is ambiguous and is rejected
#' explicitly; use `add_step()` to sequence changes instead.
#'
#' @param ... `Counterfactual` objects, each with exactly one step.
#' @return A one-step `Counterfactual`.
#' @export
combine_counterfactuals <- function(...) {
    objects <- list(...)
    if (!length(objects) || !all(vapply(objects, methods::is, logical(1), "Counterfactual"))) {
        stop("all arguments must be Counterfactual objects")
    }
    if (any(vapply(objects, function(object) length(object@steps) != 1L, logical(1)))) {
        stop("combine_counterfactuals() only combines single-step Counterfactual objects; use add_step() to sequence multi-step Counterfactuals")
    }
    result <- list()
    for (object in objects) {
        step_changes <- object@steps[[1L]]@changes
        for (name in names(step_changes)) {
            if (!is.null(result[[name]]) &&
                !isTRUE(all.equal(result[[name]], step_changes[[name]]))) {
                stop("conflicting values supplied for counterfactual field '", name, "'")
            }
            result[[name]] <- step_changes[[name]]
        }
    }
    new("Counterfactual", steps = list(new("CounterfactualStep", changes = result)))
}

.counterfactual_subset <- function(exit, n, labels = NULL) {
    if (is.null(exit)) return(rep(TRUE, n))
    if (is.logical(exit)) {
        if (length(exit) != n) stop("'exit' must be a length-k logical vector")
        return(exit)
    }
    if (is.character(exit)) {
        if (is.null(labels) || any(!exit %in% labels)) {
            stop("character 'exit' values must match fitted product labels")
        }
        exit <- match(exit, labels)
    }
    if (!is.numeric(exit) || anyNA(exit) || any(exit < 1) ||
        any(exit > n) || any(exit != as.integer(exit))) {
        stop("'exit' must contain valid product indices or labels")
    }
    active <- rep(TRUE, n)
    active[as.integer(exit)] <- FALSE
    active
}

.validate_counterfactual_step <- function(step, spec) {
    capabilities <- .model_counterfactual_capabilities(spec)
    unsupported <- names(step@changes)[!vapply(names(step@changes), function(name) {
        isTRUE(capabilities[[name]])
    }, logical(1))]
    if (length(unsupported)) {
        stop("model '", spec$id, "' does not support counterfactual field(s): ",
             paste(unsupported, collapse = ", "))
    }
    invisible(step)
}

.validate_counterfactual <- function(cf, spec) {
    for (step in cf@steps) .validate_counterfactual_step(step, spec)
    invisible(cf)
}

## Resolve a possibly-named change vector (quality, costs, ...) against the
## model's current active product labels.  Named vectors are validated
## against the currently-active (non-exited) labels and expanded to a
## full-length vector filled with `default` for untouched products;
## unnamed vectors are returned as-is, preserving today's positional,
## fixed-dimension behavior.
.resolve_named_shock <- function(values, labels, subset, field_name, default = 0) {
    if (is.null(names(values)) || any(!nzchar(names(values)))) {
        stop("'", field_name, "' must be a named vector (product label = value) once the market has changed dimension via entry or exit")
    }
    if (anyDuplicated(names(values))) {
        stop("'", field_name, "' contains duplicate product labels")
    }
    unknown <- setdiff(names(values), labels)
    if (length(unknown)) {
        stop("'", field_name, "' references unknown product label(s): ",
             paste(unknown, collapse = ", "))
    }
    active_labels <- labels[subset]
    exited <- setdiff(names(values), active_labels)
    if (length(exited)) {
        stop("'", field_name, "' references product(s) that are not active (exited or not yet entered): ",
             paste(exited, collapse = ", "))
    }
    full <- stats::setNames(rep(default, length(labels)), labels)
    full[names(values)] <- values
    full
}

.counterfactual_attach <- function(result, fit, cf) {
    attr(result, "counterfactual") <- list(
        model_spec = fit@spec,
        fields = cf
    )
    result
}
