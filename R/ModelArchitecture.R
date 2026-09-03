#' A calibrated structural model
#'
#' `AntitrustFit` is a lightweight wrapper around an existing antitrust S4
#' model.  It records the normalized model identity and calibration metadata
#' while leaving simulation results as the package's existing S4 classes.
#'
#' @slot spec An \code{antitrust_model_spec} object.
#' @slot model The calibrated existing S4 model state.
#' @slot parameters Structural parameters recovered by calibration.
#' @slot observed The observed pre-merger inputs.
#' @slot diagnostics Calibration warnings, messages, and dispatch metadata.
#' @export
#' @exportClass AntitrustFit
setClass(
    Class = "AntitrustFit",
    representation = representation(
        spec = "ANY",
        model = "ANY",
        parameters = "list",
        observed = "list",
        diagnostics = "list"
    ),
    prototype = prototype(
        spec = list(),
        parameters = list(),
        observed = list(),
        diagnostics = list()
    )
)


#' Calibrate a structural antitrust model
#'
#' It delegates calibration to the corresponding model-specific legacy
#' constructor, so existing equations, normalizations, solvers, and warnings
#' remain the behavioral source of truth.
#'
#' @param demand A demand-system name or an object returned by
#'   \code{\link{model_spec}}.
#' @param conduct A conduct name.  Omit it when `demand` is a model
#'   specification object.
#' @param prices A length-k vector of observed product prices.
#' @param shares A length-k vector of observed product shares. Required by
#'   demand systems whose calibration uses shares.
#' @param margins An optional length-k vector of observed product margins;
#'   required by model families whose calibration uses margins.
#' @param ownerPre Pre-merger ownership vector or matrix.
#' @param quantities A length-k vector of observed product quantities for
#'   Linear or LogLin calibration, or an n-by-k plant-quantity matrix for
#'   general Cournot/Stackelberg calibration.
#' @param variant A model-specific calibration variant, such as `"alm"` or
#'   `"auction2nd"` for downstream second-score vertical bargaining.
#' @param knownElast A known own-price elasticity for PCAIDS calibration.
#' @param mktElast A known market own-price elasticity for PCAIDS calibration.
#' @param ... Additional options accepted by the model-specific legacy
#'   calibration constructor. For vertical bargaining, supply the upstream
#'   inputs \code{pricesUp}, \code{marginsUp}, and \code{ownerPreUp} here.
#' @return An \code{AntitrustFit} object.
#' @export
calibrate <- function(demand, conduct = NULL, prices, shares = NULL,
                      margins = NULL,
                      ownerPre, quantities = NULL, variant = "standard",
                      knownElast = NULL, mktElast = NULL, ...) {
    spec <- .architecture_model_spec(demand, conduct, variant)

    entry <- .model_registry_entry(spec$demand, spec$conduct, spec$variant)
    if (!.model_registry_supports(spec, "calibrate")) {
        stop("calibrate() currently supports the registered standard and ALM model variants, including PCAIDS; unsupported structural models retain their legacy constructors.")
    }

    dots <- list(...)
    ## Keep a canonical baseline call so update() can genuinely rerun the
    ## target model's calibration equations.  The legacy constructor later
    ## receives placeholder post-merger state, but that implementation detail
    ## must not become part of the stored observed-data call.
    calibration_args <- c(
        list(demand = spec$demand, conduct = spec$conduct,
             prices = prices, shares = shares, margins = margins,
             ownerPre = ownerPre, quantities = quantities,
             variant = spec$variant, knownElast = knownElast,
             mktElast = mktElast),
        dots
    )
    forbidden <- intersect(names(dots), c("ownerPost", "mcDelta", "subset",
                                           "isLeaderPost"))
    if (length(forbidden)) {
        stop("'", forbidden[[1]], "' is a simulation scenario; supply it to simulate().")
    }
    if (identical(spec$conduct, "vertical_bargaining")) {
        vertical_post <- intersect(names(dots), c("ownerPostUp",
                                                  "ownerPostDown",
                                                  "mcDeltaUp", "mcDeltaDown"))
        if (length(vertical_post)) {
            stop("'", vertical_post[[1]], "' is a vertical simulation scenario; supply it to simulate().")
        }
    }
    if (identical(spec$conduct, "stackelberg")) {
        stack_post <- intersect(names(dots), c("productsPost", "mcfunPost",
                                               "vcfunPost", "dmcfunPost",
                                               "capacitiesPost"))
        if (length(stack_post)) {
            stop("'", stack_post[[1]], "' is a Stackelberg simulation scenario; supply it to simulate().")
        }
    }
    ## These formals are explicit for PCAIDS, but were historically passed
    ## through `...` for other demand systems.  Preserve that forwarding
    ## behavior for existing calibrate() callers.
    if (!is.null(knownElast) &&
        !spec$demand %in% c("pcaids", "pcaids_nests")) {
        dots$knownElast <- knownElast
    }
    if (!is.null(mktElast) &&
        !spec$demand %in% c("pcaids", "pcaids_nests")) {
        dots$mktElast <- mktElast
    }

    observed_extra <- list()

    if (identical(spec$conduct, "vertical_bargaining")) {
        if (is.null(shares) || is.null(margins)) {
            stop("'shares' and 'margins' must be supplied for vertical bargaining calibration.")
        }
        vertical_required <- c("pricesUp", "marginsUp", "ownerPreUp")
        missing_vertical <- vertical_required[
            vapply(vertical_required, function(x) is.null(dots[[x]]), logical(1))
        ]
        if (length(missing_vertical)) {
            stop("'", missing_vertical[[1]], "' must be supplied for vertical bargaining calibration.")
        }
        prices_up <- dots$pricesUp
        margins_up <- dots$marginsUp
        owner_pre_up <- dots$ownerPreUp
        vertical_supply <- if (identical(spec$variant, "auction2nd")) {
            "2nd"
        } else {
            "bertrand"
        }
        if (!is.null(dots$supplyDown)) {
            requested_supply <- match.arg(dots$supplyDown, c("bertrand", "2nd"))
            if (!identical(requested_supply, vertical_supply)) {
                stop("'supplyDown' conflicts with the selected vertical bargaining variant.")
            }
        }
        dots$pricesUp <- NULL
        dots$marginsUp <- NULL
        dots$ownerPreUp <- NULL
        dots$supplyDown <- NULL
        has_nests <- !is.null(dots$nests) && any(!is.na(dots$nests))
        if (identical(spec$demand, "logit_nests") && !has_nests) {
            stop("'nests' must be supplied for nested vertical bargaining calibration.")
        }
        if (!identical(spec$demand, "logit_nests") && has_nests) {
            stop("'nests' is only supported for nested vertical bargaining calibration.")
        }
        constructor_args <- list(
            supplyDown = vertical_supply,
            sharesDown = shares,
            pricesDown = prices,
            marginsDown = margins,
            ownerPreDown = ownerPre,
            ownerPostDown = ownerPre,
            pricesUp = prices_up,
            marginsUp = margins_up,
            ownerPreUp = owner_pre_up,
            ownerPostUp = owner_pre_up,
            mcDeltaDown = rep(0, length(prices)),
            mcDeltaUp = rep(0, length(prices_up))
        )
        observed_extra <- list(
            pricesUp = prices_up,
            marginsUp = margins_up,
            ownerPreUp = owner_pre_up
        )
    } else if (spec$demand %in% c("linear", "loglin") &&
        identical(spec$conduct, "bertrand")) {
        if (is.null(quantities)) {
            stop("'quantities' must be supplied when calibrating Linear or LogLin demand.")
        }
        constructor_args <- list(
            prices = prices,
            quantities = quantities,
            margins = margins,
            ownerPre = ownerPre,
            ## The legacy constructors require both ownership states.
            ## Calibration uses pre-merger ownership for this placeholder;
            ## simulate() replaces it for the counterfactual.
            ownerPost = ownerPre
        )
    } else if (spec$demand == "pcaids") {
        if (is.null(knownElast)) {
            stop("'knownElast' must be supplied when calibrating PCAIDS demand.")
        }
        if (is.null(mktElast)) mktElast <- -1
        constructor_args <- list(
            shares = shares,
            knownElast = knownElast,
            mktElast = mktElast,
            prices = prices,
            ownerPre = ownerPre,
            ## PCAIDS requires both ownership states even though the legacy
            ## constructor performs its own price-change calculation.
            ownerPost = ownerPre
        )
    } else if (spec$demand == "pcaids_nests") {
        if (is.null(knownElast)) {
            stop("'knownElast' must be supplied when calibrating nested PCAIDS demand.")
        }
        if (is.null(margins)) {
            stop("'margins' must be supplied when calibrating nested PCAIDS demand.")
        }
        if (is.null(mktElast)) mktElast <- -1
        if (is.null(dots$nests)) {
            stop("'nests' must be supplied when calibrating nested PCAIDS demand.")
        }
        constructor_args <- list(
            shares = shares,
            margins = margins,
            knownElast = knownElast,
            mktElast = mktElast,
            prices = prices,
            ownerPre = ownerPre,
            ## PCAIDS requires both ownership states even though the legacy
            ## constructor performs its own price-change calculation.
            ownerPost = ownerPre
        )
    } else if (spec$demand == "auction2nd_cap") {
        if (is.null(dots$capacities)) {
            stop("'capacities' must be supplied when calibrating the capacity-constrained auction model.")
        }
        if (is.null(margins)) {
            stop("'margins' must be supplied when calibrating the capacity-constrained auction model.")
        }
        capacities <- dots$capacities
        dots$capacities <- NULL
        constructor_args <- list(
            capacities = capacities,
            margins = margins,
            prices = prices,
            ownerPre = ownerPre,
            ## The auction constructor requires both ownership states;
            ## simulate() replaces this placeholder for the scenario.
            ownerPost = ownerPre
        )
        observed_extra <- list(capacities = capacities)
    } else if (spec$conduct %in% c("cournot", "stackelberg") &&
               spec$demand %in% c("linear", "loglin")) {
        if (is.null(quantities)) {
            stop("'quantities' must be supplied when calibrating a Cournot model.")
        }
        constructor_args <- list(
            prices = prices,
            quantities = quantities,
            margins = margins,
            ownerPre = ownerPre,
            ## The legacy constructor requires both ownership states;
            ## simulate() replaces this placeholder for the scenario.
            ownerPost = ownerPre
        )
        if (spec$conduct == "stackelberg") {
            isLeaderPre <- dots$isLeaderPre
            dots$isLeaderPre <- NULL
            if (is.null(isLeaderPre)) {
                isLeaderPre <- matrix(FALSE, nrow = nrow(quantities),
                                      ncol = ncol(quantities))
            }
            constructor_args$isLeaderPre <- isLeaderPre
            ## Calibration uses the observed leader structure for the
            ## placeholder post-merger state; simulate() supplies the
            ## counterfactual leader structure.
            constructor_args$isLeaderPost <- isLeaderPre
            observed_extra <- list(isLeaderPre = isLeaderPre)
        }
        if (spec$demand == "loglin") {
            constructor_args$demand <- rep("log", length(prices))
        }
    } else {
        constructor_args <- list(
            prices = prices,
            shares = shares,
            margins = margins,
            ownerPre = ownerPre,
            ## The legacy constructors require both ownership states.
            ## Calibration uses pre-merger ownership for this placeholder;
            ## simulate() replaces it for the counterfactual.
            ownerPost = ownerPre
        )
    }
    duplicate_core <- intersect(names(dots), names(constructor_args))
    if (length(duplicate_core)) {
        stop("argument(s) supplied more than once: ", paste(duplicate_core, collapse = ", "))
    }
    constructor_args <- c(constructor_args, dots)

    captured <- .capture_architecture_conditions(
        do.call(.legacy_constructor(entry$calibrator), constructor_args)
    )
    model <- captured$value
    solver <- if (spec$conduct %in% c("bertrand", "bargaining") &&
                  !is.null(dots$solver)) dots$solver else "nleqslv"

    new(
        "AntitrustFit",
        spec = spec,
        model = model,
        parameters = .model_parameters(model),
        observed = c(list(
            prices = prices,
            shares = shares,
            margins = margins,
            ownerPre = ownerPre
        ), observed_extra,
        if (!is.null(quantities)) list(quantities = quantities) else list()),
            diagnostics = list(
                status = "completed",
                model_class = class(model)[[1]],
                solver = solver,
                calibration_args = calibration_args,
                constrain.reserve = if (spec$demand == "auction2nd_cap" &&
                                    !is.null(dots[["constrain.reserve"]])) {
                dots[["constrain.reserve"]]
            } else {
                TRUE
            },
            warnings = captured$warnings,
            messages = captured$messages
        )
    )
}


#' Construct a fitted structural model from supplied parameters
#'
#' @param demand A demand-system name or an object returned by
#'   \code{\link{model_spec}}.
#' @param conduct A conduct name.  Omit it when `demand` is a model
#'   specification object.
#' @param variant A model-specific calibration variant, such as `"alm"` or
#'   `"auction2nd"` for downstream second-score vertical bargaining.
#' @param output Logical indicator for an output (`TRUE`) or input (`FALSE`)
#'   market when the selected model supports both orientations.
#' @param prices A length-k vector of observed product prices.
#' @param parameters A named list of structural demand parameters.
#' @param ownerPre Pre-merger ownership vector or matrix.
#' @param shares Optional observed shares retained as model metadata.
#' @param margins Optional observed margins retained as model metadata.
#' @param quantities Observed quantities for Linear or LogLin demand.
#' @param insideSize Market size passed to the existing `sim()` constructor.
#' @param priceOutside Optional outside-good price.
#' @param priceStart Optional equilibrium price starting values.
#' @param labels Product labels.
#' @param ... Additional options accepted by the legacy parameterized
#'   constructor.
#' @return An \code{AntitrustFit} object.
#' @export
specify <- function(demand, conduct = NULL, prices, parameters, ownerPre,
                    shares = NULL, margins = NULL, quantities = NULL,
                    insideSize = 1,
                    priceOutside, priceStart,
                    labels = paste("Prod", 1:length(prices), sep = ""),
                    variant = "standard", output = NULL, ...) {
    spec <- .architecture_model_spec(demand, conduct, variant)

    if (!.model_registry_supports(spec, "specify")) {
        stop("specify() currently supports the registered standard model variants and BLP parameter loading; ALM variants require their model-specific calibration.")
    }
    if (!is.list(parameters)) {
        stop("'parameters' must be a list.")
    }

    dots <- list(...)
    specification_args <- c(
        list(demand = spec$demand, conduct = spec$conduct,
             prices = prices, parameters = parameters,
             ownerPre = ownerPre, shares = shares, margins = margins,
             quantities = quantities, insideSize = insideSize,
             variant = spec$variant, output = output),
        dots
    )
    forbidden <- intersect(names(dots), c("ownerPost", "mcDelta", "subset"))
    if (length(forbidden)) {
        stop("'", forbidden[[1]], "' is a simulation scenario; supply it to simulate().")
    }
    if (spec$demand %in% c("linear", "loglin") && is.null(quantities)) {
        stop("'quantities' must be supplied when specifying Linear or LogLin demand.")
    }
    sim_shares <- if (spec$demand %in% c("linear", "loglin")) {
        quantities / sum(quantities)
    } else {
        shares
    }
    constructor_args <- list(
        prices = prices,
        shares = sim_shares,
        margins = margins,
        supply = spec$conduct,
        demand = switch(spec$demand,
                        linear = "Linear",
                        aids = "AIDS",
                        loglin = "LogLin",
                        logit = "Logit",
                        ces = "CES",
                        logit_nests = "LogitNests",
                        ces_nests = "CESNests",
                        logit_cap = "LogitCap",
                        blp = "BLP"),
        demand.param = parameters,
        ownerPre = ownerPre,
        ownerPost = ownerPre,
        insideSize = insideSize,
        labels = labels
    )
    if (!missing(priceOutside)) constructor_args$priceOutside <- priceOutside
    if (!missing(priceStart)) constructor_args$priceStart <- priceStart
    if (!is.null(output)) constructor_args$output <- output
    duplicate_core <- intersect(names(dots), names(constructor_args))
    if (length(duplicate_core)) {
        stop("argument(s) supplied more than once: ", paste(duplicate_core, collapse = ", "))
    }

    captured <- .capture_architecture_conditions(
        do.call(.legacy_sim_constructor, c(constructor_args, dots))
    )
    model <- captured$value
    if (spec$demand %in% c("linear", "loglin")) {
        ## The legacy parameterized constructor exposes its quantity input
        ## through the historical `shares` argument and validates that slot as
        ## a share vector.  Restore the explicit quantity scale for the new
        ## API, then refresh cost state using the supplied structural demand.
        model@quantities <- quantities
        model@shares <- quantities / sum(quantities)
        model@pricePre <- prices
        model@mcPre <- calcMC(model, preMerger = TRUE)
        model@mcPost <- calcMC(model, preMerger = FALSE)
        model@pricePost <- calcPrices(model, preMerger = FALSE,
                                      subset = rep(TRUE, length(prices)))
    }
    new(
        "AntitrustFit",
        spec = spec,
        model = model,
        parameters = .model_parameters(model),
        observed = list(
            prices = prices,
            shares = shares,
            margins = margins,
            quantities = quantities,
            ownerPre = ownerPre
        ),
        diagnostics = list(
            status = "completed",
            source = "specified",
            route = "specify",
            model_class = class(model)[[1]],
            solver = if (identical(spec$conduct, "bargaining") &&
                        !is.null(dots$solver)) dots$solver else "nleqslv",
            specification_args = specification_args,
            warnings = captured$warnings,
            messages = captured$messages
        )
    )
}


## Apply the quality/entry environmental changes of one CounterfactualStep to
## `model`, before that step's ownership/cost/capacity/bargaining/exit
## changes are solved.  Entry expands the product dimension first, so
## same-step ownership/cost vectors are validated against the expanded
## model exactly like any other sequential-dimension-changing field.
.apply_step_environment <- function(model, step) {
    entrants <- step@changes$entry
    if (!is.null(entrants)) {
        for (one_entrant in entrants) {
            model <- .expand_entrant(model, one_entrant)
        }
    }
    quality <- step@changes$quality
    if (!is.null(quality)) {
        model <- .apply_quality(model, quality)
    }
    model
}

## `costs` may be supplied as a plain positional vector (existing,
## fixed-dimension behavior) or, once the market has changed dimension via
## entry/exit, as a named vector resolved against the model's current
## active labels.  Vertical bargaining's `costs` is always a list with
## `up`/`down` numeric vectors (not itself a "named shock" vector) and is
## passed through unchanged; vertical models are not entry/quality-
## supported, so dimension-changing named resolution never applies to them.
.resolve_step_costs <- function(model, costs, is_vertical) {
    if (is.null(costs) || is_vertical || is.null(names(costs))) return(costs)
    .resolve_named_shock(costs, model@labels, model@subset, "costs", default = 0)
}

## calcMC() for every supported model family always recomputes the
## pre-shock marginal cost fresh from immutable calibration constants
## (observed prices/margins, or cost functions) and applies `mcDelta` as a
## one-shot multiplicative wedge relative to that fixed baseline -- it never
## reads a promoted `mcPre`.  So the cumulative proportional cost change
## (not a single step's shock) is the real persistent cost-environment
## state, and successive steps' cost shocks must compound multiplicatively:
## (1 + cumulative) <- (1 + cumulative) * (1 + this step's shock).
.compound_costs <- function(current, step_shock) {
    if (is.null(step_shock)) return(current)
    if (is.null(current)) current <- rep(0, length(step_shock))
    (1 + current) * (1 + step_shock) - 1
}

.compound_vertical_costs <- function(current, step_shock) {
    if (is.null(step_shock)) return(current)
    list(
        up = .compound_costs(current$up, step_shock$up),
        down = .compound_costs(current$down, step_shock$down)
    )
}

## Solve one equilibrium: apply ownership/cost/capacity/bargaining/exit
## changes to the (already environment-adjusted) `model` and return the
## legacy S4 result.  This is the extracted core of the original one-step
## `simulate()`; it never reads or writes `fit@model` directly so it can be
## called repeatedly with a promoted current state for sequential paths.
.model_simulate_step <- function(fit, model, ownerPost = NULL, mcDelta = NULL,
                                 exit = NULL, subset = NULL, priceStart = NULL,
                                 capacitiesPost = NULL, bargpowerPost = NULL,
                                 solver = NULL, isMax = FALSE, dots = list()) {
    is_vertical <- methods::is(model, "VertBargBertLogit")
    if (is.null(ownerPost)) {
        ownerPost <- if (is_vertical) {
            list(up = model@up@ownerPre, down = model@down@ownerPre)
        } else {
            model@ownerPre
        }
    }
    nprods <- if (is_vertical) length(model@down@prices) else length(model@prices)
    has_subset_slot <- is_vertical || "subset" %in% methods::slotNames(model)
    current_subset <- if (is_vertical) {
        model@down@subset
    } else if (has_subset_slot) {
        model@subset
    } else {
        rep(TRUE, nprods)
    }
    if (length(current_subset) != nprods) current_subset <- rep(TRUE, nprods)
    if (!is.null(exit)) {
        labels <- if (is_vertical) model@down@labels else model@labels
        subset <- current_subset & .counterfactual_subset(exit, nprods, labels)
    }
    if (is.null(subset)) subset <- current_subset
    nmc <- if (is_vertical) {
        NA_integer_
    } else if (methods::is(model, "Cournot")) {
        nrow(model@quantities)
    } else {
        nprods
    }
    if (is_vertical) {
        if (is.null(mcDelta)) {
            mcDelta <- list(
                up = rep(0, length(model@up@prices)),
                down = rep(0, length(model@down@prices))
            )
        }
        if (!is.list(mcDelta) ||
            !all(c("up", "down") %in% names(mcDelta)) ||
            !is.numeric(mcDelta$up) || !is.numeric(mcDelta$down) ||
            length(mcDelta$up) != length(model@up@prices) ||
            length(mcDelta$down) != length(model@down@prices) ||
            anyNA(mcDelta$up) || anyNA(mcDelta$down)) {
            stop("For vertical bargaining, 'mcDelta' must be a list with numeric 'up' and 'down' vectors matching the fitted markets and no NAs.")
        }
    } else {
        if (is.null(mcDelta)) mcDelta <- rep(0, nmc)
        if (!is.numeric(mcDelta) || length(mcDelta) != nmc || anyNA(mcDelta)) {
            stop("'mcDelta' must be a numeric vector with the same length as the fitted cost units and no NAs.")
        }
    }
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
        stop("'subset' must be a logical vector the same length as the fitted prices with at least one TRUE value.")
    }
    if ((!is_vertical && any(mcDelta > 0, na.rm = TRUE)) ||
        (is_vertical && (any(mcDelta$up > 0, na.rm = TRUE) ||
                         any(mcDelta$down > 0, na.rm = TRUE)))) {
        warning("positive values of 'mcDelta' imply an INCREASE in marginal costs")
    }

    if (is_vertical) {
        if (!is.list(ownerPost) ||
            !all(c("up", "down") %in% names(ownerPost))) {
            stop("For vertical bargaining, 'ownerPost' must be a list with 'up' and 'down' ownership vectors.")
        }
        if (length(ownerPost$up) != length(model@up@prices) ||
            length(ownerPost$down) != length(model@down@prices)) {
            stop("Vertical 'ownerPost' ownership vectors must match the fitted upstream and downstream markets.")
        }
        model@up@ownerPost <- ownerPost$up
        model@down@ownerPost <- ownerPost$down
        model@up@mcDelta <- mcDelta$up
        model@down@mcDelta <- mcDelta$down
        model@down@subset <- subset
        ## The vertical constructor fixes bargaining power at one for any
        ## newly integrated upstream/downstream product.  Recreate that
        ## post-merger state from the calibrated bargaining powers without
        ## rerunning the pre-merger calibration optimizer.
        model@up@bargpowerPost <- model@up@bargpowerPre
        integrated_post <- ownerPost$up == ownerPost$down
        model@up@bargpowerPost[integrated_post] <- 1
        pre_vertical <- model@up@ownerPre == model@down@ownerPre
        post_vertical <- model@up@ownerPost == model@down@ownerPost
        model@isHorizontal <- !any(!pre_vertical & post_vertical)
        is_upstream_horizontal <- !isTRUE(all.equal(
            model@up@ownerPre, model@up@ownerPost, check.attributes = FALSE
        ))
        model@isUpstream <- model@isHorizontal && is_upstream_horizontal
        model <- ownerToMatrix(model, preMerger = FALSE)
        mc_post <- calcMC(model, preMerger = FALSE)
        model@up@mcPost <- mc_post$up
        model@down@mcPost <- mc_post$down
        prices_post <- do.call(calcPrices, c(list(model, preMerger = FALSE), dots))
        model@up@pricePost <- prices_post$up
        model@down@pricePost <- prices_post$down
        return(model)
    }

    if (methods::is(model, "Auction2ndCap")) {
        if (is.matrix(ownerPost) || length(ownerPost) != nprods) {
            stop("'ownerPost' must be a length-k ownership vector for the capacity-constrained auction model.")
        }
        model@ownerPost <- ownerPost
    } else {
        model@ownerPost <- ownerPost
        model@ownerPost <- ownerToMatrix(model, preMerger = FALSE)
    }
    model@mcDelta <- mcDelta
    if (!methods::is(model, "Auction2ndCap")) {
        model@subset <- subset
    }
    if (!is.null(capacitiesPost)) {
        if (methods::is(model, "Stackelberg")) {
            if (length(capacitiesPost) != nrow(model@quantities)) {
                stop("'capacitiesPost' must have the same length as the fitted plants.")
            }
        } else {
            if (!methods::is(model, "LogitCap")) {
                stop("'capacitiesPost' is only supported for LogitCap or Stackelberg fits.")
            }
            if (length(capacitiesPost) != nprods) {
                stop("'capacitiesPost' must have the same length as the fitted prices.")
            }
        }
        model@capacitiesPost <- capacitiesPost
    }
    if (!is.null(bargpowerPost)) {
        if (!(fit@spec$conduct %in% c("bargaining", "bargaining2nd"))) {
            stop("'bargpowerPost' is only supported for bargaining fits.")
        }
        if (!is.numeric(bargpowerPost) || length(bargpowerPost) != nprods ||
            any(!is.finite(bargpowerPost)) || any(bargpowerPost < 0) ||
            any(bargpowerPost > 1)) {
            stop("'bargpowerPost' must be a numeric vector of values in [0, 1] with the same length as the fitted prices.")
        }
        model@bargpowerPost <- bargpowerPost
    }
    if (!is.null(priceStart)) {
        model@priceStart <- priceStart
    }
    if (methods::is(model, "Stackelberg")) {
        scenario_fields <- c("isLeaderPost", "productsPost", "mcfunPost",
                             "vcfunPost", "dmcfunPost")
        if (is.null(dots$isLeaderPost)) {
            model@isLeaderPost <- model@isLeaderPre
        } else {
            model@isLeaderPost <- dots$isLeaderPost
        }
        if (!is.null(dots$productsPost)) {
            model@productsPost <- dots$productsPost
        }
        if (!is.null(dots$mcfunPost)) {
            model@mcfunPost <- dots$mcfunPost
        }
        if (!is.null(dots$vcfunPost)) {
            model@vcfunPost <- dots$vcfunPost
        }
        if (!is.null(dots$dmcfunPost)) {
            model@dmcfunPost <- dots$dmcfunPost
        }
    }
    if (methods::is(model, "Auction2ndCap")) {
        if (!is.null(solver)) {
            stop("The capacity-constrained auction uses its existing equilibrium methods; 'solver' cannot be overridden.")
        }
        if (any(!subset)) {
            stop("'subset' is not supported for the capacity-constrained auction model.")
        }
        if (isTRUE(fit@diagnostics$constrain.reserve)) {
            model@reservePost <- model@reservePre
        } else {
            model@reservePost <- calcOptimalReserve(model, preMerger = FALSE)
        }
        model@mcPost <- calcMC(model, preMerger = FALSE, exAnte = FALSE)
        model@pricePost <- calcPrices(model, preMerger = FALSE, exAnte = FALSE)
        return(model)
    }
    if (methods::is(model, "Stackelberg")) {
        if (!is.null(solver)) {
            stop("Stackelberg uses its existing nonlinear quantity solver; 'solver' cannot be overridden.")
        }
        scenario_fields <- c("isLeaderPost", "productsPost", "mcfunPost",
                             "vcfunPost", "dmcfunPost")
        solver_args <- dots[setdiff(names(dots), scenario_fields)]
        model@quantityPost <- do.call(
            calcQuantities,
            c(list(object = model, preMerger = FALSE), solver_args)
        )
        ## The legacy Stackelberg constructor leaves the marginal-cost slots
        ## empty; downstream calcMC() computes them from equilibrium state.
        model@pricePost <- calcPrices(model, preMerger = FALSE)
        return(model)
    }
    if (methods::is(model, "Cournot")) {
        if (!is.null(solver)) {
            stop("Cournot uses its existing nonlinear quantity solver; 'solver' cannot be overridden.")
        }
        model@quantityPost <- do.call(calcQuantities, c(list(model, preMerger = FALSE), dots))
        ## The legacy Cournot constructor leaves the marginal-cost slots
        ## empty; downstream calcMC() computes them from the equilibrium
        ## quantities.  Preserve that result-object convention.
        model@pricePost <- calcPrices(model, preMerger = FALSE)
        return(model)
    }
    if (fit@spec$demand %in% c("aids", "pcaids", "pcaids_nests")) {
        if (!is.null(solver) && !identical(solver, "nleqslv")) {
            stop("AIDS uses its existing nonlinear price solver; 'solver' cannot be overridden.")
        }
        model@priceDelta <- do.call(calcPriceDelta, c(
            list(model, isMax = isMax, subset = subset), dots
        ))
    }
    model@mcPost <- calcMC(model, preMerger = FALSE)

    if (identical(fit@spec$conduct, "bertrand")) {
        if (is.null(solver)) solver <- fit@diagnostics$solver
        solver <- match.arg(solver, c("nleqslv", "ag"))
        if (!(fit@spec$demand %in% c("logit", "ces", "logit_nests",
                                      "ces_nests", "logit_cap", "blp",
                                      "pcaids", "pcaids_nests")) &&
            identical(solver, "ag")) {
            stop("This model uses its existing nonlinear price solver; 'solver = \"ag\"' is not supported.")
        }
        if (identical(fit@spec$demand, "blp") && identical(solver, "ag")) {
            stop("BLP uses its existing nonlinear price solver; 'solver = \"ag\"' is not supported.")
        }
        if (identical(solver, "ag")) {
            model@pricePost <- calcPricesAG(model, preMerger = FALSE, isMax = isMax)
        } else if (fit@spec$demand %in% c("linear", "loglin", "aids",
                                           "pcaids", "pcaids_nests")) {
            ## These legacy price methods do not expose `isMax`; passing it
            ## through `...` would alter their optimizer call signature.
            model@pricePost <- do.call(calcPrices, c(
                list(model, preMerger = FALSE, subset = subset), dots
            ))
        } else {
            model@pricePost <- do.call(calcPrices, c(
                list(model, preMerger = FALSE, subset = subset, isMax = isMax), dots
            ))
        }
    } else if (identical(fit@spec$conduct, "cournot")) {
        if (!is.null(solver) && !identical(solver, "nleqslv")) {
            stop("Logit-Cournot uses its existing nonlinear price solver; 'solver' cannot be overridden.")
        }
        model@pricePost <- calcPrices(model, preMerger = FALSE)
    } else if (identical(fit@spec$conduct, "bargaining")) {
        if (is.null(solver)) solver <- fit@diagnostics$solver
        solver <- match.arg(solver, c("nleqslv", "ag"))
        if (identical(solver, "ag")) {
            model@pricePost <- calcPricesAG(model, preMerger = FALSE,
                                             isMax = isMax, subset = subset)
        } else {
            model@pricePost <- do.call(calcPrices, c(
                list(model, preMerger = FALSE, subset = subset, isMax = isMax), dots
            ))
        }
    } else {
        if (!is.null(solver)) {
            stop("This model uses its existing price solver; 'solver' cannot be overridden.")
        }
        ## Auction and second-score bargaining models retain their direct or
        ## model-specific equilibrium methods.  The legacy second-score
        ## methods intentionally do not receive a subset argument.
        model@pricePost <- calcPrices(model, preMerger = FALSE)
    }

    model
}


## Solve every step of a multi-step Counterfactual in order, promoting the
## solved state between steps so step t+1 begins from step t's equilibrium
## (never resetting to the originally calibrated baseline).  Returns a list
## of legacy S4 results, one per step.
.model_simulate_steps <- function(fit, initial_model, steps, initial_costs = NULL) {
    state <- initial_model
    cumulative_costs <- initial_costs
    results <- vector("list", length(steps))
    for (i in seq_along(steps)) {
        step <- steps[[i]]
        state <- .apply_step_environment(state, step)
        changes <- step@changes
        is_vertical <- methods::is(state, "VertBargBertLogit")
        step_shock <- .resolve_step_costs(state, changes$costs, is_vertical)
        cumulative_costs <- if (is_vertical) {
            .compound_vertical_costs(cumulative_costs, step_shock)
        } else {
            .compound_costs(cumulative_costs, step_shock)
        }
        dots <- list()
        if (!is.null(changes$leader)) dots$isLeaderPost <- changes$leader
        if (!is.null(changes$products)) dots$productsPost <- changes$products
        result <- .model_simulate_step(
            fit, state,
            ownerPost = changes$ownership,
            mcDelta = cumulative_costs,
            exit = changes$exit,
            capacitiesPost = changes$capacity,
            bargpowerPost = changes$bargaining,
            dots = dots
        )
        results[[i]] <- result
        state <- .promote_post_to_pre(result, step)
    }
    list(results = results, cumulative_costs = cumulative_costs)
}

#' Simulate a counterfactual from a calibrated structural model
#'
#' For a one-step `Counterfactual` (or the equivalent direct legacy scenario
#' arguments), `simulate()` returns the existing S4 simulation-result object
#' unchanged, exactly as before. For a multi-step `Counterfactual` (built
#' with \code{\link{add_step}}), each step is solved in turn, promoting the
#' previous step's solved equilibrium into the next step's starting state,
#' and a `CounterfactualPath` recording every step's result is returned.
#' `fit` may also be a `CounterfactualPath`, in which case simulation
#' resumes from that path's final solved state and the returned path's
#' history includes the earlier steps.
#'
#' @param fit An \code{AntitrustFit} returned by \code{\link{calibrate}}, or
#'   a `CounterfactualPath` to resume from.
#' @param ownerPost Post-counterfactual ownership vector or matrix, or a
#' `Counterfactual` object. For
#'   vertical bargaining, use a list with \code{up} and \code{down} ownership
#'   vectors.
#' @param mcDelta A model-specific vector of proportional cost changes, with
#'   one element per product for most models and one element per plant for
#'   general Cournot models. For vertical bargaining, use a list with numeric
#'   \code{up} and \code{down} vectors.
#' @param subset A length-k logical vector selecting products in the
#'   counterfactual equilibrium.
#' @param priceStart Optional price starting values.
#' @param capacitiesPost Optional post-counterfactual capacities for LogitCap
#'   or plant capacities for Stackelberg.
#' @param bargpowerPost Optional post-counterfactual bargaining powers for
#'   bargaining models.
#' @param solver Optional solver override.  The calibration solver is reused
#'   by default for Logit-Bertrand; Logit-Cournot retains its legacy solver.
#' @param isMax Whether to run the existing local profit-maximum check.
#' @param ... Additional arguments passed to the existing price solver.
#'   For Stackelberg fits, \code{isLeaderPost}, \code{productsPost},
#'   \code{mcfunPost}, \code{vcfunPost}, and \code{dmcfunPost} may be supplied as counterfactual
#'   structure; the remaining arguments are passed to the existing quantity
#'   solver.
#' @return For a one-step counterfactual, an existing S4 simulation-result
#'   object, such as \code{Logit} or \code{LogitCournot}. For a multi-step
#'   counterfactual, a `CounterfactualPath`.
#' @export
simulate <- function(fit, ownerPost = NULL,
                     mcDelta = NULL,
                     subset = NULL,
                     priceStart, capacitiesPost = NULL,
                     bargpowerPost = NULL,
                     solver = NULL, isMax = FALSE, ...) {
    resume_path <- if (methods::is(fit, "CounterfactualPath")) fit else NULL
    if (!is.null(resume_path)) fit <- NULL

    dots <- list(...)
    cf <- if (methods::is(ownerPost, "Counterfactual")) ownerPost else NULL

    if (!is.null(resume_path)) {
        if (is.null(cf)) {
            stop("simulate() on a CounterfactualPath requires a Counterfactual as its second argument.")
        }
        base_fit <- resume_path@diagnostics$fit
        if (is.null(base_fit)) {
            stop("'fit' CounterfactualPath does not retain enough diagnostics to resume simulation.")
        }
        .validate_counterfactual(cf, base_fit@spec)
        last_step <- resume_path@steps[[length(resume_path@steps)]]
        state <- .promote_post_to_pre(final_result(resume_path), last_step)
        solved <- .model_simulate_steps(base_fit, state, cf@steps,
                                        initial_costs = resume_path@diagnostics$cumulative_costs)
        return(new(
            "CounterfactualPath",
            initial = resume_path@initial,
            steps = c(resume_path@steps, cf@steps),
            results = c(resume_path@results, solved$results),
            diagnostics = list(fit = base_fit, cumulative_costs = solved$cumulative_costs)
        ))
    }

    if (!methods::is(fit, "AntitrustFit")) {
        stop("'fit' must be an AntitrustFit returned by calibrate() or specify(), or a CounterfactualPath.")
    }
    if (!.model_registry_supports(fit@spec, "simulate")) {
        stop("simulate() currently supports Linear, LogLin, AIDS, and PCAIDS Bertrand models; Logit/CES Bertrand, Cournot, Stackelberg, auction, and bargaining models; nested Logit/CES Bertrand models; LogitCap-Bertrand; and BLP simulations.")
    }

    model <- fit@model
    is_vertical <- methods::is(model, "VertBargBertLogit")

    if (!is.null(cf)) {
        .validate_counterfactual(cf, fit@spec)
        conflicts <- c(
            if (!is.null(mcDelta)) "mcDelta",
            if (!is.null(subset)) "subset",
            if (!is.null(capacitiesPost)) "capacitiesPost",
            if (!is.null(bargpowerPost)) "bargpowerPost",
            if (!missing(priceStart)) "priceStart",
            intersect(names(dots), c("isLeaderPost", "productsPost",
                                     "mcfunPost", "vcfunPost", "dmcfunPost"))
        )
        if (length(conflicts)) {
            stop("cannot combine a Counterfactual with legacy scenario argument(s): ",
                 paste(unique(conflicts), collapse = ", "))
        }

        if (length(cf@steps) > 1L) {
            solved <- .model_simulate_steps(fit, model, cf@steps)
            return(new(
                "CounterfactualPath",
                initial = model,
                steps = cf@steps,
                results = solved$results,
                diagnostics = list(fit = fit, cumulative_costs = solved$cumulative_costs)
            ))
        }

        step <- cf@steps[[1L]]
        model <- .apply_step_environment(model, step)
        changes <- step@changes
        if (!is.null(changes$leader)) dots$isLeaderPost <- changes$leader
        if (!is.null(changes$products)) dots$productsPost <- changes$products
        result <- .model_simulate_step(
            fit, model,
            ownerPost = changes$ownership,
            mcDelta = .resolve_step_costs(model, changes$costs, is_vertical),
            exit = changes$exit,
            priceStart = if (!missing(priceStart)) priceStart,
            capacitiesPost = changes$capacity,
            bargpowerPost = changes$bargaining,
            solver = solver, isMax = isMax, dots = dots
        )
        return(.counterfactual_attach(result, fit, cf))
    }

    if (is.null(ownerPost)) {
        if (missing(ownerPost)) stop("'ownerPost' must be supplied for simulate().")
        stop("'ownerPost' must not be NULL for simulate().")
    }
    fields <- list(ownership = ownerPost, costs = mcDelta,
                   capacity = capacitiesPost, bargaining = bargpowerPost)
    fields <- fields[!vapply(fields, is.null, logical(1))]
    cf <- do.call(counterfactual, fields)
    result <- .model_simulate_step(
        fit, model,
        ownerPost = ownerPost, mcDelta = mcDelta, subset = subset,
        priceStart = if (!missing(priceStart)) priceStart,
        capacitiesPost = capacitiesPost, bargpowerPost = bargpowerPost,
        solver = solver, isMax = isMax, dots = dots
    )
    .counterfactual_attach(result, fit, cf)
}


#' Recalibrate a fitted antitrust model
#'
#' `update()` rebuilds the complete calibration call from the baseline inputs
#' stored in an `AntitrustFit`, applies the supplied replacements, and invokes
#' the target model's ordinary calibration routine.  It is intentionally not a
#' class conversion and is available for fits created by `calibrate()`.
#'
#' @param object An `AntitrustFit` returned by `calibrate()`.
#' @param ... Baseline data, model-specification arguments, or model-specific
#'   calibration options to replace.
#' @param evaluate If `FALSE`, return the reconstructed calibration call.
#' @return A newly calibrated `AntitrustFit`, or a call when `evaluate` is
#'   `FALSE`.
#' @export
#' @exportS3Method stats::update AntitrustFit
update.AntitrustFit <- function(object, ..., evaluate = TRUE) {
    if (!methods::is(object, "AntitrustFit")) {
        stop("'object' must be an AntitrustFit returned by calibrate().")
    }
    calibration_args <- object@diagnostics$calibration_args
    if (!is.list(calibration_args) || is.null(names(calibration_args))) {
        if (identical(object@diagnostics$route, "respecify") ||
            !is.null(object@diagnostics$source_calibration_args)) {
            stop("update() requires a fit whose current specification was created by calibrate(); this fit was created by respecify()")
        }
        stop("this fit does not retain a calibration call; update() requires a fit created by calibrate()")
    }

    replacements <- list(...)
    if (length(replacements)) {
        if (is.null(names(replacements)) || any(!nzchar(names(replacements)))) {
            stop("update() arguments must be named calibration or model-specification arguments")
        }
        calibration_args[names(replacements)] <- replacements
    }
    if (!isTRUE(evaluate)) {
        return(as.call(c(list(quote(calibrate)), calibration_args)))
    }
    do.call(calibrate, calibration_args)
}


#' Respecify a fitted antitrust model
#'
#' `respecify()` applies an explicitly registered model transition. Same-demand
#' conduct transitions carry portable structural parameters. Demand transitions
#' use the taxonomy stored in the transition registry: structural restrictions,
#' algebraic or conditional translations, and canonical first-order Linear or
#' LogLin representations. Target mean values, intercepts, and slopes are
#' solved analytically from the fitted baseline, or target structural primitives
#' are required explicitly. No transition calibrates from source margins or
#' minimizes an elasticity-distance objective.
#'
#' @param fit An `AntitrustFit` returned by `calibrate()` or `specify()`.
#' @param demand Optional target demand-system name.
#' @param conduct Optional target conduct name.
#' @param variant Optional target model variant.
#' @param ... Transition-specific target primitives. Depending on the
#'   registered transition, these may include `alpha`, `gamma`, `nests`, and
#'   `sigma`.
#' @return A newly constructed `AntitrustFit` under the target specification.
#' @seealso [`specify()`], [`update.AntitrustFit()`]
#' @export
respecify <- function(fit, demand = NULL, conduct = NULL,
                      variant = NULL, ...) {
    if (!methods::is(fit, "AntitrustFit")) {
        stop("'fit' must be an AntitrustFit returned by calibrate() or specify().")
    }
    extras <- list(...)
    if (length(extras) &&
        (is.null(names(extras)) || any(!nzchar(names(extras))))) {
        stop("respecify() transition arguments must be named")
    }
    allowed_extras <- c("alpha", "gamma", "nests", "sigma")
    unsupported_extras <- setdiff(names(extras), allowed_extras)
    if (length(unsupported_extras)) {
        stop("unsupported respecify() transition argument(s): ",
             paste(unsupported_extras, collapse = ", "))
    }

    source <- fit@spec
    target <- model_spec(
        demand = if (is.null(demand)) source$demand else demand,
        conduct = if (is.null(conduct)) source$conduct else conduct,
        variant = if (is.null(variant)) source$variant else variant
    )
    if (identical(source$id, target$id)) {
        stop("respecify() requires a different registered model specification")
    }
    transition <- .model_transition_entry(source, target)
    if (identical(transition$handler, "portable") && length(extras)) {
        stop("respecify() transition from '", source$id, "' to '",
             target$id, "' does not accept transition-specific demand arguments")
    }

    if (identical(transition$handler, "demand")) {
        translated <- .translate_antitrust_demand(fit, target, transition,
                                                  extras)
        result <- translated$fit
        result@parameters <- .model_parameters(result@model)
        result@observed <- list(
            prices = translated$state$prices,
            shares = translated$shares,
            margins = fit@observed$margins,
            quantities = translated$state$quantities,
            ownerPre = translated$state$owner
        )
        result@diagnostics$source <- "respecify"
        result@diagnostics$route <- "respecify"
        result@diagnostics$transition <- list(
            from = source$id,
            to = target$id,
            kind = transition$kind,
            required_arguments = transition$required_arguments,
            retained = transition$retain,
            derived = transition$derived,
            discarded = transition$discarded,
            recomputed = transition$recompute,
            invalidated = transition$invalidate,
            calibration_required = transition$calibration_required
        )
        result@diagnostics$translation <- translated$diagnostics
        ## Retain the older diagnostic name for callers of the pilot API.
        result@diagnostics$local_translation <- translated$diagnostics
        result@diagnostics$source_calibration_args <-
            fit@diagnostics$calibration_args

        ## A translated fit was not calibrated under its target specification.
        ## Preserve the source call only as provenance; never turn it into a
        ## synthetic target calibration call.
        result@diagnostics$calibration_args <- NULL
        return(result)
    }

    portable <- fit@parameters[intersect(transition$retain,
                                         names(fit@parameters))]
    missing_portable <- setdiff(transition$retain, names(portable))
    if (length(missing_portable)) {
        for (name in missing_portable) {
            if (.hasSlot(fit@model, name)) {
                portable[[name]] <- methods::slot(fit@model, name)
            }
        }
    }
    missing_portable <- setdiff(transition$retain, names(portable))
    if (length(missing_portable)) {
        stop("source fit does not contain portable parameter(s): ",
             paste(missing_portable, collapse = ", "))
    }

    baseline <- fit@diagnostics$calibration_args
    if (!is.list(baseline)) baseline <- fit@diagnostics$specification_args
    if (!is.list(baseline)) {
        stop("this fit does not retain the baseline inputs needed for respecify()")
    }
    prices <- baseline$prices
    owner_pre <- baseline$ownerPre
    if (is.null(prices) || is.null(owner_pre)) {
        stop("source fit does not retain prices and ownerPre needed for respecify()")
    }

    specify_args <- list(
        demand = target$demand,
        conduct = target$conduct,
        variant = target$variant,
        prices = prices,
        parameters = portable,
        ownerPre = owner_pre
    )
    ## These are observed baseline quantities used to construct the target
    ## state, not identifying margins.  In particular, margins are deliberately
    ## omitted from the respecification call.
    for (name in c("shares", "quantities", "insideSize", "priceOutside",
                   "priceStart", "labels")) {
        if (!is.null(baseline[[name]])) specify_args[[name]] <- baseline[[name]]
    }

    result <- do.call(specify, specify_args)
    result@parameters <- portable
    result@observed <- fit@observed
    result@diagnostics$source <- "respecify"
    result@diagnostics$route <- "respecify"
        result@diagnostics$transition <- list(
            from = source$id,
            to = target$id,
            kind = transition$kind,
            required_arguments = transition$required_arguments,
            retained = transition$retain,
            derived = transition$derived,
            discarded = transition$discarded,
            recomputed = transition$recompute,
            invalidated = transition$invalidate,
            calibration_required = transition$calibration_required
    )
    ## A portable respecification is not a calibration under the target
    ## specification.  Keep any source call as provenance only.
    result@diagnostics$source_calibration_args <-
        fit@diagnostics$calibration_args
    result@diagnostics$calibration_args <- NULL
    result
}


.architecture_model_spec <- function(demand, conduct, variant = "standard") {
    if (inherits(demand, "antitrust_model_spec")) {
        if (!is.null(conduct)) {
            stop("'conduct' must be omitted when 'demand' is a model specification.")
        }
        requested_variant <- .normalize_variant_name(variant)
        if (!identical(requested_variant, "standard") &&
            !identical(requested_variant, demand$variant)) {
            stop("'variant' conflicts with the supplied model specification.")
        }
        demand
    } else {
        model_spec(demand = demand, conduct = conduct, variant = variant)
    }
}


.legacy_constructor <- function(name) {
    get(name, envir = asNamespace("antitrust"), inherits = FALSE)
}


.legacy_sim_constructor <- function(...) {
    get(".sim_legacy", envir = asNamespace("antitrust"), inherits = FALSE)(...)
}


.model_parameters <- function(model) {
    if (methods::is(model, "Auction2ndCap")) {
        list(sellerCostParms = model@sellerCostParms,
             buyerValuation = model@buyerValuation,
             reservePre = model@reservePre,
             reservePost = model@reservePost)
    } else if (methods::is(model, "PCAIDSNests")) {
        list(slopes = model@slopes, intercepts = model@intercepts,
             mktElast = model@mktElast, knownElast = model@knownElast,
             knownElastIndex = model@knownElastIndex,
             nests = model@nests, nestsParms = model@nestsParms)
    } else if (methods::is(model, "PCAIDS")) {
        list(slopes = model@slopes, intercepts = model@intercepts,
             mktElast = model@mktElast, knownElast = model@knownElast,
             knownElastIndex = model@knownElastIndex)
    } else if (methods::is(model, "AIDS")) {
        list(slopes = model@slopes, intercepts = model@intercepts,
             mktElast = model@mktElast)
    } else if (methods::is(model, "LogLin")) {
        list(slopes = model@slopes, intercepts = model@intercepts)
    } else if (methods::is(model, "Linear")) {
        list(slopes = model@slopes, intercepts = model@intercepts)
    } else if (methods::is(model, "Stackelberg")) {
        list(slopes = model@slopes, intercepts = model@intercepts,
             mktElast = model@mktElast, isLeaderPre = model@isLeaderPre)
    } else if (methods::is(model, "VertBargBertLogit")) {
        list(down = model@down@slopes,
             bargpowerPre = model@up@bargpowerPre,
             bargpowerPost = model@up@bargpowerPost)
    } else if (methods::is(model, "Cournot")) {
        list(slopes = model@slopes, intercepts = model@intercepts,
             mktElast = model@mktElast)
    } else if (methods::is(model, "Bertrand") && is.list(model@slopes)) {
        model@slopes
    } else if (methods::is(model, "Bertrand")) {
        list(slopes = model@slopes)
    } else {
        list()
    }
}


.capture_architecture_conditions <- function(expr) {
    warnings <- character()
    messages <- character()
    value <- withCallingHandlers(
        force(expr),
        warning = function(condition) {
            message_text <- conditionMessage(condition)
            warnings <<- c(warnings, message_text)
            ## Calibration intentionally uses the pre-merger ownership for the
            ## legacy constructor's required post-merger slot.  Do not print
            ## its no-merger warning, but retain it in diagnostics.
            if (grepl("'ownerPost' and 'ownerPre' are the same", message_text,
                      fixed = TRUE)) {
                invokeRestart("muffleWarning")
            }
        },
        message = function(condition) {
            message_text <- conditionMessage(condition)
            messages <<- c(messages, message_text)
            invokeRestart("muffleMessage")
        }
    )
    list(value = value, warnings = warnings, messages = messages)
}
