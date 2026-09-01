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
#' @param variant A model-specific calibration variant, such as `"alm"`.
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
        dots$pricesUp <- NULL
        dots$marginsUp <- NULL
        dots$ownerPreUp <- NULL
        if (!is.null(dots$nests) && any(!is.na(dots$nests))) {
            stop("Nested vertical bargaining is not yet supported by calibrate().")
        }
        constructor_args <- list(
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
#' @param variant A model-specific calibration variant, such as `"alm"`.
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
                    variant = "standard", ...) {
    spec <- .architecture_model_spec(demand, conduct, variant)

    if (!.model_registry_supports(spec, "specify")) {
        stop("specify() currently supports the registered standard model variants and BLP parameter loading; ALM variants require their model-specific calibration.")
    }
    if (!is.list(parameters)) {
        stop("'parameters' must be a list.")
    }

    dots <- list(...)
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
            model_class = class(model)[[1]],
            solver = if (identical(spec$conduct, "bargaining") &&
                        !is.null(dots$solver)) dots$solver else "nleqslv",
            warnings = captured$warnings,
            messages = captured$messages
        )
    )
}


#' Simulate a counterfactual from a calibrated structural model
#'
#' @param fit An \code{AntitrustFit} returned by \code{\link{calibrate}}.
#' @param ownerPost Post-counterfactual ownership vector or matrix. For
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
#' @return An existing S4 simulation-result object, such as \code{Logit} or
#'   \code{LogitCournot}.
#' @export
simulate <- function(fit, ownerPost,
                     mcDelta = NULL,
                     subset = NULL,
                     priceStart, capacitiesPost = NULL,
                     bargpowerPost = NULL,
                     solver = NULL, isMax = FALSE, ...) {
    if (!methods::is(fit, "AntitrustFit")) {
        stop("'fit' must be an AntitrustFit returned by calibrate() or specify().")
    }
    if (missing(ownerPost)) {
        stop("'ownerPost' must be supplied for simulate().")
    }
    if (!.model_registry_supports(fit@spec, "simulate")) {
        stop("simulate() currently supports Linear, LogLin, AIDS, and PCAIDS Bertrand models; Logit/CES Bertrand, Cournot, Stackelberg, auction, and bargaining models; nested Logit/CES Bertrand models; LogitCap-Bertrand; and BLP simulations.")
    }

    model <- fit@model
    is_vertical <- methods::is(model, "VertBargBertLogit")
    nprods <- if (is_vertical) length(model@down@prices) else length(model@prices)
    if (is.null(subset)) subset <- rep(TRUE, nprods)
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
        prices_post <- calcPrices(model, preMerger = FALSE, ...)
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
    dots <- list(...)
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
    if (!missing(priceStart)) {
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
        model@quantityPost <- calcQuantities(model, preMerger = FALSE, ...)
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
        model@priceDelta <- calcPriceDelta(model, isMax = isMax,
                                            subset = subset, ...)
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
            model@pricePost <- calcPrices(model, preMerger = FALSE,
                                           subset = subset, ...)
        } else {
            model@pricePost <- calcPrices(model, preMerger = FALSE,
                                           subset = subset, isMax = isMax, ...)
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
            model@pricePost <- calcPrices(model, preMerger = FALSE,
                                           subset = subset, isMax = isMax, ...)
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
