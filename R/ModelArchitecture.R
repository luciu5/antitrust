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
        model = "Antitrust",
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
#'   Linear or LogLin calibration.
#' @param variant A model-specific calibration variant, such as `"alm"`.
#' @param knownElast A known own-price elasticity for PCAIDS calibration.
#' @param mktElast A known market own-price elasticity for PCAIDS calibration.
#' @param ... Additional options accepted by the model-specific legacy
#'   calibration constructor.
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
    forbidden <- intersect(names(dots), c("ownerPost", "mcDelta", "subset"))
    if (length(forbidden)) {
        stop("'", forbidden[[1]], "' is a simulation scenario; supply it to simulate().")
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

    if (spec$demand %in% c("linear", "loglin")) {
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
        ), if (!is.null(quantities)) list(quantities = quantities) else list()),
        diagnostics = list(
            status = "completed",
            model_class = class(model)[[1]],
            solver = solver,
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
#' @param ownerPost Post-counterfactual ownership vector or matrix.
#' @param mcDelta A length-k vector of proportional marginal-cost changes.
#' @param subset A length-k logical vector selecting products in the
#'   counterfactual equilibrium.
#' @param priceStart Optional price starting values.
#' @param capacitiesPost Optional post-counterfactual capacities for LogitCap.
#' @param bargpowerPost Optional post-counterfactual bargaining powers for
#'   bargaining models.
#' @param solver Optional solver override.  The calibration solver is reused
#'   by default for Logit-Bertrand; Logit-Cournot retains its legacy solver.
#' @param isMax Whether to run the existing local profit-maximum check.
#' @param ... Additional arguments passed to the existing price solver.
#' @return An existing S4 simulation-result object, such as \code{Logit} or
#'   \code{LogitCournot}.
#' @export
simulate <- function(fit, ownerPost,
                     mcDelta = rep(0, length(fit@model@prices)),
                     subset = rep(TRUE, length(fit@model@prices)),
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
        stop("simulate() currently supports Linear, LogLin, AIDS, and PCAIDS Bertrand models; Logit/CES Bertrand, Cournot, auction, and bargaining models; nested Logit/CES Bertrand models; LogitCap-Bertrand; and BLP simulations.")
    }

    model <- fit@model
    nprods <- length(model@prices)
    if (!is.numeric(mcDelta) || length(mcDelta) != nprods || anyNA(mcDelta)) {
        stop("'mcDelta' must be a numeric vector with the same length as the fitted prices and no NAs.")
    }
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
        stop("'subset' must be a logical vector the same length as the fitted prices with at least one TRUE value.")
    }
    if (any(mcDelta > 0, na.rm = TRUE)) {
        warning("positive values of 'mcDelta' imply an INCREASE in marginal costs")
    }

    model@ownerPost <- ownerPost
    model@ownerPost <- ownerToMatrix(model, preMerger = FALSE)
    model@mcDelta <- mcDelta
    model@subset <- subset
    if (!is.null(capacitiesPost)) {
        if (!methods::is(model, "LogitCap")) {
            stop("'capacitiesPost' is only supported for LogitCap fits.")
        }
        if (length(capacitiesPost) != nprods) {
            stop("'capacitiesPost' must have the same length as the fitted prices.")
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
    if (methods::is(model, "PCAIDSNests")) {
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
