#' A calibrated structural model
#'
#' `AntitrustFit` is a lightweight wrapper around an existing antitrust S4
#' model.  It records the normalized model identity and calibration metadata
#' while leaving simulation results as the package's existing S4 classes.
#'
#' @slot spec An code{antitrust_model_spec} object.
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
#' The initial implementation supports Logit-Bertrand and Logit-Cournot.  It
#' delegates calibration to the corresponding legacy constructor, so their
#' model-specific equations, normalizations, solvers, and warnings remain the
#' behavioral source of truth.
#'
#' @param demand A demand-system name or an object returned by
#'   code{link{model_spec}}.
#' @param conduct A conduct name.  Omit it when `demand` is a model
#'   specification object.
#' @param prices A length-k vector of observed product prices.
#' @param shares A length-k vector of observed product shares.
#' @param margins A length-k vector of observed product margins.
#' @param ownerPre Pre-merger ownership vector or matrix.
#' @param ... Additional options accepted by the model-specific legacy
#'   calibration constructor.
#' @return An code{AntitrustFit} object.
#' @export
calibrate <- function(demand, conduct = NULL, prices, shares, margins,
                      ownerPre, ...) {
    spec <- if (inherits(demand, "antitrust_model_spec")) {
        if (!is.null(conduct)) {
            stop("'conduct' must be omitted when 'demand' is a model specification.")
        }
        demand
    } else {
        model_spec(demand = demand, conduct = conduct)
    }

    if (!identical(spec$demand, "logit") ||
        !spec$conduct %in% c("bertrand", "cournot")) {
        stop("calibrate() currently supports Logit-Bertrand and Logit-Cournot.")
    }

    dots <- list(...)
    forbidden <- intersect(names(dots), c("ownerPost", "mcDelta", "subset"))
    if (length(forbidden)) {
        stop("'", forbidden[[1]], "' is a simulation scenario; supply it to simulate().")
    }

    entry <- .model_registry_entry(spec$demand, spec$conduct)
    constructor_args <- list(
        prices = prices,
        shares = shares,
        margins = margins,
        ownerPre = ownerPre,
        ## The legacy constructors require both ownership states.  Calibration
        ## uses pre-merger ownership for this placeholder; simulate() replaces
        ## it before recovering post-merger costs and solving the counterfactual.
        ownerPost = ownerPre
    )
    duplicate_core <- intersect(names(dots), names(constructor_args))
    if (length(duplicate_core)) {
        stop("argument(s) supplied more than once: ", paste(duplicate_core, collapse = ", "))
    }
    constructor_args <- c(constructor_args, dots)

    captured <- .capture_architecture_conditions(
        do.call(.legacy_constructor(entry$calibrator), constructor_args)
    )
    model <- captured$value
    solver <- if (identical(spec$conduct, "bertrand") &&
                  !is.null(dots$solver)) dots$solver else "nleqslv"

    new(
        "AntitrustFit",
        spec = spec,
        model = model,
        parameters = .model_parameters(model),
        observed = list(
            prices = prices,
            shares = shares,
            margins = margins,
            ownerPre = ownerPre
        ),
        diagnostics = list(
            status = "completed",
            model_class = class(model)[[1]],
            solver = solver,
            warnings = captured$warnings,
            messages = captured$messages
        )
    )
}


#' Simulate a counterfactual from a calibrated structural model
#'
#' @param fit An code{AntitrustFit} returned by code{link{calibrate}}.
#' @param ownerPost Post-counterfactual ownership vector or matrix.
#' @param mcDelta A length-k vector of proportional marginal-cost changes.
#' @param subset A length-k logical vector selecting products in the
#'   counterfactual equilibrium.
#' @param priceStart Optional price starting values.
#' @param solver Optional solver override.  The calibration solver is reused
#'   by default for Logit-Bertrand; Logit-Cournot retains its legacy solver.
#' @param isMax Whether to run the existing local profit-maximum check.
#' @param ... Additional arguments passed to the existing price solver.
#' @return An existing S4 simulation-result object, such as code{Logit} or
#'   code{LogitCournot}.
#' @export
simulate <- function(fit, ownerPost,
                     mcDelta = rep(0, length(fit@model@prices)),
                     subset = rep(TRUE, length(fit@model@prices)),
                     priceStart, solver = NULL, isMax = FALSE, ...) {
    if (!is(fit, "AntitrustFit")) {
        stop("'fit' must be an AntitrustFit returned by calibrate() or specify().")
    }
    if (missing(ownerPost)) {
        stop("'ownerPost' must be supplied for simulate().")
    }
    if (!identical(fit@spec$demand, "logit") ||
        !fit@spec$conduct %in% c("bertrand", "cournot")) {
        stop("simulate() currently supports Logit-Bertrand and Logit-Cournot.")
    }

    model <- fit@model
    nprods <- length(model@prices)
    if (!is.numeric(mcDelta) || length(mcDelta) != nprods || anyNA(mcDelta)) {
        stop("'mcDelta' must be a numeric vector with the same length as the fitted prices and no NAs.")
    }
    if (!is.logical(subset) || length(subset) != nprods || !any(subset)) {
        stop("'subset' must be a logical vector the same length as the fitted prices with at least one TRUE value.")
    }

    model@ownerPost <- ownerPost
    model@ownerPost <- ownerToMatrix(model, preMerger = FALSE)
    model@mcDelta <- mcDelta
    model@subset <- subset
    if (!missing(priceStart)) {
        model@priceStart <- priceStart
    }
    model@mcPost <- calcMC(model, preMerger = FALSE)

    if (identical(fit@spec$conduct, "bertrand")) {
        if (is.null(solver)) solver <- fit@diagnostics$solver
        solver <- match.arg(solver, c("nleqslv", "ag"))
        if (identical(solver, "ag")) {
            model@pricePost <- calcPricesAG(model, preMerger = FALSE, isMax = isMax)
        } else {
            model@pricePost <- calcPrices(model, preMerger = FALSE,
                                           isMax = isMax, ...)
        }
    } else {
        if (!is.null(solver) && !identical(solver, "nleqslv")) {
            stop("Logit-Cournot uses its existing nonlinear price solver; 'solver' cannot be overridden.")
        }
        model@pricePost <- calcPrices(model, preMerger = FALSE)
    }

    model
}


.legacy_constructor <- function(name) {
    get(name, envir = asNamespace("antitrust"), inherits = FALSE)
}


.model_parameters <- function(model) {
    if (is(model, "Bertrand") && is.list(model@slopes)) {
        model@slopes
    } else if (is(model, "Bertrand")) {
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
