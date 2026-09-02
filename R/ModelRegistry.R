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
        list(from = "logit::bertrand", to = "logit::cournot",
             retain = c("alpha", "meanval"),
             recompute = c("marginal costs", "Cournot supply state"),
             invalidate = c("Bertrand margins"),
             calibration_required = FALSE),
        list(from = "logit::cournot", to = "logit::bertrand",
             retain = c("alpha", "meanval"),
             recompute = c("marginal costs", "Bertrand supply state"),
             invalidate = c("Cournot margins"),
             calibration_required = FALSE),
        list(from = "ces::bertrand", to = "ces::cournot",
             retain = c("gamma", "alpha", "meanval", "shareInside"),
             recompute = c("marginal costs", "Cournot supply state"),
             invalidate = c("Bertrand margins"),
             calibration_required = FALSE),
        list(from = "ces::cournot", to = "ces::bertrand",
             retain = c("gamma", "alpha", "meanval", "shareInside"),
             recompute = c("marginal costs", "Bertrand supply state"),
             invalidate = c("Cournot margins"),
             calibration_required = FALSE)
    )
    function() entries
})

.model_transition_entry <- function(from, to) {
    entries <- .model_transition_registry()
    matches <- Filter(function(entry) {
        identical(entry$from, from$id) && identical(entry$to, to$id)
    }, entries)
    if (!length(matches)) {
        stop("respecify() transition from '", from$id, "' to '",
             to$id, "' is not supported; use update() to recalibrate the target model")
    }
    matches[[1L]]
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
    registry <- setNames(entries, vapply(entries, `[[`, character(1), "id"))
    function() registry
})


.model_registry_entry <- function(demand, conduct, variant = "standard") {
    .model_registry()[[.model_registry_key(demand, conduct, variant)]]
}


.model_registry_supports <- function(spec, operation) {
    entry <- .model_registry_entry(spec$demand, spec$conduct, spec$variant)
    !is.null(entry) && isTRUE(entry[[operation]])
}
