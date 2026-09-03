#' Describe a demand/conduct model combination
#'
#' @param demand A human-readable demand-system name, such as \code{"logit"}
#'   or \code{"ces"}.
#' @param conduct A human-readable conduct name, such as \code{"bertrand"}
#'   or \code{"cournot"}.
#' @param variant A model-specific calibration variant. The default is
#'   \code{"standard"}; \code{"alm"} selects the existing unknown-market-
#'   elasticity calibration where supported, and \code{"auction2nd"} selects
#'   downstream second-score vertical bargaining.
#' @return A small object of class \code{antitrust_model_spec} containing
#'   normalized model names.
#' @export
model_spec <- function(demand, conduct, variant = "standard") {
    if (missing(demand) || length(demand) != 1L || is.na(demand)) {
        stop("'demand' must be a single model name.")
    }
    if (missing(conduct) || length(conduct) != 1L || is.na(conduct)) {
        stop("'conduct' must be a single conduct name.")
    }

    demand_text <- tolower(trimws(as.character(demand)))
    demand_text <- gsub("[[:space:].-]+", "_", demand_text)
    variant_aliases <- c(
        logit_alm = "logit",
        logitalm = "logit",
        ces_alm = "ces",
        cesalm = "ces"
    )
    if (demand_text %in% names(variant_aliases)) {
        requested_variant <- .normalize_variant_name(variant)
        if (!missing(variant) && !identical(requested_variant, "standard") &&
            !identical(requested_variant, "alm")) {
            stop("variant is inconsistent with the ALM demand alias.")
        }
        demand <- variant_aliases[[demand_text]]
        variant <- "alm"
    } else {
        demand <- .normalize_demand_name(demand)
        variant <- .normalize_variant_name(variant)
    }
    conduct <- .normalize_conduct_name(conduct)
    entry <- .model_registry_entry(demand, conduct, variant)
    if (is.null(entry)) {
        stop(demand, " / ", conduct, " currently not supported.")
    }

    structure(
        list(
            demand = demand,
            conduct = conduct,
            variant = variant,
            id = entry$id
        ),
        class = c("antitrust_model_spec", "list")
    )
}


#' List demand/conduct combinations supported by antitrust
#'
#' @return A data frame with normalized demand names, conduct names, existing
#'   S4 result classes, and legacy constructor mappings.
#' @export
supportedModels <- function() {
    registry <- .model_registry()
    do.call(rbind, lapply(registry, function(entry) {
        data.frame(
            demand = entry$demand,
            conduct = entry$conduct,
            variant = if (is.null(entry$variant)) "standard" else entry$variant,
            class = entry$class,
            calibrator = entry$calibrator,
            calibrate = entry$calibrate,
            specify = entry$specify,
            simulate = entry$simulate,
            stringsAsFactors = FALSE,
            row.names = entry$id
        )
    }))
}


## Explicit fitted-model transitions.  These are deliberately narrower than
## the model registry: a transition is allowed only when the source and target
## specifications have portable structural primitives and a tested target
## construction path.  This is not a general S4 coercion table.
.model_transition_registry <- local({
    entries <- list(
        list(from = "logit", to = "logit_nests",
             kind = "conditional-translation", required_arguments = c("nests", "sigma"),
             retain = c("alpha", "prices", "quantities", "ownership", "conduct"),
             derived = c("meanval"), discarded = character(),
             recompute = c("marginal costs", "target supply state")),
        list(from = "ces", to = "ces_nests",
             kind = "conditional-translation", required_arguments = c("nests", "sigma"),
             retain = c("gamma", "prices", "quantities", "ownership", "conduct"),
             derived = c("meanval"), discarded = character(),
             recompute = c("marginal costs", "target supply state")),
        list(from = "logit_nests", to = "logit",
             kind = "structural-restriction", required_arguments = character(),
             retain = c("alpha", "prices", "quantities", "ownership", "conduct"),
             derived = c("meanval"), discarded = c("nests", "sigma"),
             recompute = c("marginal costs", "target supply state")),
        list(from = "ces_nests", to = "ces",
             kind = "structural-restriction", required_arguments = character(),
             retain = c("gamma", "prices", "quantities", "ownership", "conduct"),
             derived = c("meanval"), discarded = c("nests", "sigma"),
             recompute = c("marginal costs", "target supply state")),
        list(from = "logit", to = "ces",
             kind = "algebraic-translation", required_arguments = c("gamma"),
             retain = c("prices", "quantities", "ownership", "conduct"),
             derived = c("meanval"), discarded = c("alpha"),
             recompute = c("marginal costs", "target supply state")),
        list(from = "ces", to = "logit",
             kind = "algebraic-translation", required_arguments = c("alpha"),
             retain = c("prices", "quantities", "ownership", "conduct"),
             derived = c("meanval"), discarded = c("gamma"),
             recompute = c("marginal costs", "target supply state")),
        list(from = "logit_nests", to = "ces_nests",
             kind = "algebraic-translation", required_arguments = c("gamma"),
             retain = c("nests", "lambda", "prices", "quantities", "ownership", "conduct"),
             derived = c("sigma", "meanval"), discarded = c("alpha"),
             recompute = c("marginal costs", "target supply state")),
        list(from = "ces_nests", to = "logit_nests",
             kind = "algebraic-translation", required_arguments = c("alpha"),
             retain = c("nests", "sigmaCES", "prices", "quantities", "ownership", "conduct"),
             derived = c("sigma", "meanval"), discarded = c("gamma"),
             recompute = c("marginal costs", "target supply state")),
        list(from = "logit_nests", to = "ces",
             kind = "algebraic-translation", required_arguments = c("gamma"),
             retain = c("prices", "quantities", "ownership", "conduct"),
             derived = c("meanval"), discarded = c("alpha", "nests", "sigma"),
             recompute = c("marginal costs", "target supply state")),
        list(from = "ces_nests", to = "logit",
             kind = "algebraic-translation", required_arguments = c("alpha"),
             retain = c("prices", "quantities", "ownership", "conduct"),
             derived = c("meanval"), discarded = c("gamma", "nests", "sigma"),
             recompute = c("marginal costs", "target supply state")),
        list(from = "logit", to = "ces_nests",
             kind = "conditional-translation", required_arguments = c("gamma", "nests", "sigma"),
             retain = c("prices", "quantities", "ownership", "conduct"),
             derived = c("meanval"), discarded = c("alpha"),
             recompute = c("marginal costs", "target supply state")),
        list(from = "ces", to = "logit_nests",
             kind = "conditional-translation", required_arguments = c("alpha", "nests", "sigma"),
             retain = c("prices", "quantities", "ownership", "conduct"),
             derived = c("meanval"), discarded = c("gamma"),
             recompute = c("marginal costs", "target supply state"))
    )
    ## Linear and LogLin are canonical local representations of any supported
    ## differentiable demand for which the target constructor accepts supplied
    ## slopes and intercepts.  Add these declarative records to the same
    ## registry rather than hiding the family expansion in the dispatcher.
    first_order_demands <- c("logit", "logit_nests", "ces", "ces_nests",
                             "aids", "linear", "loglin")
    first_order <- lapply(setdiff(first_order_demands, "linear"),
                          function(from) list(
        from = from, to = "linear",
        kind = "first-order-linearization",
        required_arguments = character(),
        retain = c("prices", "quantities", "ownership", "conduct"),
        derived = c("slopes", "intercepts"), discarded = character(),
        recompute = c("marginal costs", "target supply state")
    ))
    first_order <- c(first_order, lapply(
        setdiff(first_order_demands, "loglin"), function(from) list(
            from = from, to = "loglin",
            kind = "first-order-loglinearization",
            required_arguments = character(),
            retain = c("prices", "quantities", "ownership", "conduct"),
            derived = c("slopes", "intercepts"), discarded = character(),
            recompute = c("marginal costs", "target supply state")
        )))
    entries <- c(entries, first_order)
    function() entries
})

.model_transition_entry <- function(from, to) {
    if (!identical(from$conduct, to$conduct) ||
        !identical(from$variant, to$variant)) {
        ## Existing conduct transitions are structural restrictions that retain
        ## the portable demand primitives.  They remain explicit below.
        entries <- list(
            list(from = "logit", to = "logit", kind = "structural-restriction",
                 required_arguments = character(), retain = c("alpha", "meanval"),
                 derived = character(), discarded = character(),
                 recompute = c("marginal costs", "target supply state")),
            list(from = "ces", to = "ces", kind = "structural-restriction",
                 required_arguments = character(), retain = c("gamma", "alpha", "meanval", "shareInside"),
                 derived = character(), discarded = character(),
                 recompute = c("marginal costs", "target supply state")))
    } else {
        entries <- .model_transition_registry()
    }
    matches <- Filter(function(entry) {
        identical(entry$from, from$demand) && identical(entry$to, to$demand)
    }, entries)
    target_entry <- .model_registry_entry(to$demand, to$conduct, to$variant)
    if (length(matches) && !is.null(target_entry) &&
        isTRUE(target_entry$specify)) {
        entry <- matches[[1L]]
        entry$from <- from$id
        entry$to <- to$id
        if (is.null(entry$calibration_required)) {
            entry$calibration_required <- FALSE
        }
        entry$handler <- if (identical(from$demand, to$demand)) {
            "portable"
        } else {
            "demand"
        }
        return(entry)
    }
    stop("respecify() transition from '", from$id, "' to '",
         to$id, "' is not supported; use update() to recalibrate the target model")
}


#' @export
print.antitrust_model_spec <- function(x, ...) {
    cat("antitrust model specification\n")
    cat("  demand:  ", x$demand, "\n", sep = "")
    cat("  conduct: ", x$conduct, "\n", sep = "")
    cat("  variant: ", x$variant, "\n", sep = "")
    invisible(x)
}


.normalize_demand_name <- function(demand) {
    demand <- tolower(trimws(as.character(demand)))
    demand <- gsub("[[:space:].-]+", "_", demand)
    aliases <- c(
        logitblp = "blp",
        cournotblp = "blp",
        logitnests = "logit_nests",
        logit_nests = "logit_nests",
        cesnests = "ces_nests",
        ces_nests = "ces_nests",
        logitcap = "logit_cap",
        logit_cap = "logit_cap",
        pcaidsnests = "pcaids_nests",
        pcaids_nests = "pcaids_nests",
        auction2ndcap = "auction2nd_cap",
        auction2nd_cap = "auction2nd_cap",
        loglin = "loglin"
    )
    if (demand %in% names(aliases)) aliases[[demand]] else demand
}


.normalize_conduct_name <- function(conduct) {
    conduct <- tolower(trimws(as.character(conduct)))
    conduct <- gsub("[[:space:]._-]+", "", conduct)
    aliases <- c(
        bertrand = "bertrand",
        cournot = "cournot",
        auction2nd = "auction2nd",
        auction2 = "auction2nd",
        bargaining = "bargaining",
        bargaining2nd = "bargaining2nd",
        bargaining2 = "bargaining2nd",
        stackelberg = "stackelberg",
        stack = "stackelberg",
        vertical = "vertical_bargaining",
        verticalbarg = "vertical_bargaining",
        verticalbargaining = "vertical_bargaining"
    )
    if (conduct %in% names(aliases)) aliases[[conduct]] else conduct
}


.normalize_variant_name <- function(variant) {
    if (length(variant) != 1L || is.na(variant)) {
        stop("'variant' must be a single model variant name.")
    }
    variant <- tolower(trimws(as.character(variant)))
    variant <- gsub("[[:space:]._-]+", "", variant)
    aliases <- c(
        default = "standard", standard = "standard", alm = "alm",
        auction2nd = "auction2nd", auction2 = "auction2nd",
        second = "auction2nd", secondscore = "auction2nd"
    )
    if (variant %in% names(aliases)) aliases[[variant]] else variant
}


.model_registry_key <- function(demand, conduct, variant = "standard") {
    variant <- .normalize_variant_name(variant)
    if (identical(variant, "standard")) {
        paste(demand, conduct, sep = "::")
    } else {
        paste(demand, conduct, variant, sep = "::")
    }
}


.model_registry <- local({
    entries <- list(
        list(id = "linear::bertrand", demand = "linear", conduct = "bertrand",
             class = "Linear", calibrator = "linear", calibrate = TRUE,
             specify = TRUE, simulate = TRUE),
        list(id = "aids::bertrand", demand = "aids", conduct = "bertrand",
             class = "AIDS", calibrator = "aids", calibrate = TRUE,
             specify = TRUE, simulate = TRUE),
        list(id = "loglin::bertrand", demand = "loglin", conduct = "bertrand",
             class = "LogLin", calibrator = "loglinear", calibrate = TRUE,
             specify = TRUE, simulate = TRUE),
        list(id = "logit::bertrand", demand = "logit", conduct = "bertrand",
             class = "Logit", calibrator = "logit", calibrate = TRUE,
             specify = TRUE, simulate = TRUE),
        list(id = "logit::cournot", demand = "logit", conduct = "cournot",
             class = "LogitCournot", calibrator = "logit.cournot", calibrate = TRUE,
             specify = TRUE, simulate = TRUE),
        list(id = "logit::auction2nd", demand = "logit", conduct = "auction2nd",
             class = "Auction2ndLogit", calibrator = "auction2nd.logit", calibrate = TRUE,
             specify = TRUE, simulate = TRUE),
        list(id = "logit::bargaining", demand = "logit", conduct = "bargaining",
             class = "BargainingLogit", calibrator = "bargaining.logit", calibrate = TRUE,
             specify = TRUE, simulate = TRUE),
        list(id = "logit::bargaining2nd", demand = "logit", conduct = "bargaining2nd",
             class = "Bargaining2ndLogit", calibrator = "bargaining2nd.logit", calibrate = TRUE,
             specify = TRUE, simulate = TRUE),
        list(id = "ces::bertrand", demand = "ces", conduct = "bertrand",
             class = "CES", calibrator = "ces", calibrate = TRUE,
             specify = TRUE, simulate = TRUE),
        list(id = "ces::cournot", demand = "ces", conduct = "cournot",
             class = "CESCournot", calibrator = "ces.cournot", calibrate = TRUE,
             specify = TRUE, simulate = TRUE),
        list(id = "ces::auction2nd", demand = "ces", conduct = "auction2nd",
             class = "Auction2ndCES", calibrator = "auction2nd.ces", calibrate = TRUE,
             specify = TRUE, simulate = TRUE),
        list(id = "ces::bargaining", demand = "ces", conduct = "bargaining",
             class = "BargainingCES", calibrator = "bargaining.ces", calibrate = TRUE,
             specify = TRUE, simulate = TRUE),
        list(id = "ces::bargaining2nd", demand = "ces", conduct = "bargaining2nd",
             class = "Bargaining2ndCES", calibrator = "bargaining2nd.ces", calibrate = TRUE,
             specify = TRUE, simulate = TRUE),
        list(id = "logit_nests::bertrand", demand = "logit_nests", conduct = "bertrand",
             class = "LogitNests", calibrator = "logit.nests", calibrate = TRUE,
             specify = TRUE, simulate = TRUE),
        list(id = "ces_nests::bertrand", demand = "ces_nests", conduct = "bertrand",
             class = "CESNests", calibrator = "ces.nests", calibrate = TRUE,
             specify = TRUE, simulate = TRUE),
        list(id = "logit_cap::bertrand", demand = "logit_cap", conduct = "bertrand",
             class = "LogitCap", calibrator = "logit.cap", calibrate = TRUE,
             specify = TRUE, simulate = TRUE),
        list(id = "pcaids::bertrand", demand = "pcaids", conduct = "bertrand",
             class = "PCAIDS", calibrator = "pcaids", calibrate = TRUE,
             specify = FALSE, simulate = TRUE),
        list(id = "pcaids_nests::bertrand", demand = "pcaids_nests",
             conduct = "bertrand", class = "PCAIDSNests",
             calibrator = "pcaids.nests", calibrate = TRUE,
             specify = FALSE, simulate = TRUE),
        list(id = "blp::bertrand", demand = "blp", conduct = "bertrand",
             class = "LogitBLP", calibrator = "sim", calibrate = FALSE,
             specify = TRUE, simulate = TRUE),
        list(id = "blp::cournot", demand = "blp", conduct = "cournot",
             class = "CournotBLP", calibrator = "sim", calibrate = FALSE,
             specify = TRUE, simulate = TRUE),
        list(id = "auction2nd_cap::auction2nd", demand = "auction2nd_cap",
             conduct = "auction2nd", class = "Auction2ndCap",
             calibrator = "auction2nd.cap", calibrate = TRUE,
             specify = FALSE, simulate = TRUE),
        list(id = "linear::cournot", demand = "linear", conduct = "cournot",
             class = "Cournot", calibrator = "cournot", calibrate = TRUE,
             specify = FALSE, simulate = TRUE),
        list(id = "loglin::cournot", demand = "loglin", conduct = "cournot",
             class = "Cournot", calibrator = "cournot", calibrate = TRUE,
             specify = FALSE, simulate = TRUE),
        list(id = "linear::stackelberg", demand = "linear", conduct = "stackelberg",
             class = "Stackelberg", calibrator = "stackelberg", calibrate = TRUE,
             specify = FALSE, simulate = TRUE),
        list(id = "loglin::stackelberg", demand = "loglin", conduct = "stackelberg",
             class = "Stackelberg", calibrator = "stackelberg", calibrate = TRUE,
             specify = FALSE, simulate = TRUE),
        list(id = "logit::vertical_bargaining", demand = "logit",
             conduct = "vertical_bargaining", class = "VertBargBertLogit",
             calibrator = "vertical.barg", calibrate = TRUE,
             specify = FALSE, simulate = TRUE),
        list(id = "logit_nests::vertical_bargaining", demand = "logit_nests",
             conduct = "vertical_bargaining", class = "VertBargBertLogitNests",
             calibrator = "vertical.barg", calibrate = TRUE,
             specify = FALSE, simulate = TRUE),
        list(id = "logit::vertical_bargaining::auction2nd", demand = "logit",
             conduct = "vertical_bargaining", variant = "auction2nd",
             class = "VertBarg2ndLogit", calibrator = "vertical.barg",
             calibrate = TRUE, specify = FALSE, simulate = TRUE),
        list(id = "logit::bertrand::alm", demand = "logit", conduct = "bertrand",
             variant = "alm", class = "LogitALM", calibrator = "logit.alm",
             calibrate = TRUE, specify = FALSE, simulate = TRUE),
        list(id = "logit::cournot::alm", demand = "logit", conduct = "cournot",
             variant = "alm", class = "LogitCournotALM",
             calibrator = "logit.cournot.alm", calibrate = TRUE,
             specify = FALSE, simulate = TRUE),
        list(id = "ces::bertrand::alm", demand = "ces", conduct = "bertrand",
             variant = "alm", class = "CESALM", calibrator = "ces.alm",
             calibrate = TRUE, specify = FALSE, simulate = TRUE),
        list(id = "ces::cournot::alm", demand = "ces", conduct = "cournot",
             variant = "alm", class = "CESCournotALM",
             calibrator = "ces.cournot.alm", calibrate = TRUE,
             specify = FALSE, simulate = TRUE),
        list(id = "logit_nests::bertrand::alm", demand = "logit_nests",
             conduct = "bertrand", variant = "alm", class = "LogitNestsALM",
             calibrator = "logit.nests.alm", calibrate = TRUE,
             specify = FALSE, simulate = TRUE),
        list(id = "logit_cap::bertrand::alm", demand = "logit_cap",
             conduct = "bertrand", variant = "alm", class = "LogitCapALM",
             calibrator = "logit.cap.alm", calibrate = TRUE,
             specify = FALSE, simulate = TRUE),
        list(id = "logit::auction2nd::alm", demand = "logit",
             conduct = "auction2nd", variant = "alm",
             class = "Auction2ndLogitALM", calibrator = "auction2nd.logit.alm",
             calibrate = TRUE, specify = FALSE, simulate = TRUE),
        list(id = "ces::auction2nd::alm", demand = "ces",
             conduct = "auction2nd", variant = "alm",
             class = "Auction2ndCESALM", calibrator = "auction2nd.ces.alm",
             calibrate = TRUE, specify = FALSE, simulate = TRUE),
        list(id = "logit::bargaining::alm", demand = "logit",
             conduct = "bargaining", variant = "alm",
             class = "BargainingLogitALM", calibrator = "bargaining.logit.alm",
             calibrate = TRUE, specify = FALSE, simulate = TRUE),
        list(id = "ces::bargaining::alm", demand = "ces",
             conduct = "bargaining", variant = "alm",
             class = "BargainingCESALM", calibrator = "bargaining.ces.alm",
             calibrate = TRUE, specify = FALSE, simulate = TRUE)
    )
    registry <- stats::setNames(entries, vapply(entries, `[[`, character(1), "id"))
    function() registry
})


.model_registry_entry <- function(demand, conduct, variant = "standard") {
    .model_registry()[[.model_registry_key(demand, conduct, variant)]]
}


.model_registry_supports <- function(spec, operation) {
    entry <- .model_registry_entry(spec$demand, spec$conduct, spec$variant)
    !is.null(entry) && isTRUE(entry[[operation]])
}

## Quality and entry are verified only for the four leaf classes whose
## product-dimensional slots are byte-identical Logit/CES vectors: bare
## Logit, CES, and their Cournot-conduct counterparts (which add no new
## slots -- conduct is dispatched separately from class structure). Every
## other Logit/CES descendant (LogitCap, LogitNests, LogitBLP,
## Auction2ndLogit*, Bargaining*, VertBarg*) is excluded until individually
## audited; all non-Logit/CES demands (Linear, LogLin, AIDS, PCAIDS*,
## Cournot, Stackelberg) are excluded outright.
.entry_quality_supported_classes <- c("Logit", "LogitCournot", "CES", "CESCournot")

.model_counterfactual_capabilities <- function(spec) {
    entry <- .model_registry_entry(spec$demand, spec$conduct, spec$variant)
    if (is.null(entry)) return(stats::setNames(logical(), character()))
    cls <- entry$class
    c(
        ownership = TRUE,
        costs = TRUE,
        exit = cls != "Auction2ndCap",
        capacity = cls %in% c("LogitCap", "LogitCapALM", "Stackelberg"),
        bargaining = spec$conduct %in% c("bargaining", "bargaining2nd"),
        leader = cls == "Stackelberg",
        products = cls == "Stackelberg",
        tariff = FALSE,
        quota = FALSE,
        quality = cls %in% .entry_quality_supported_classes,
        entry = cls %in% .entry_quality_supported_classes
    )
}
