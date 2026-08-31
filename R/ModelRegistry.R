#' Describe a demand/conduct model combination
#'
#' @param demand A human-readable demand-system name, such as code{"logit"}
#'   or code{"ces"}.
#' @param conduct A human-readable conduct name, such as code{"bertrand"}
#'   or code{"cournot"}.
#' @return A small object of class code{antitrust_model_spec} containing
#'   normalized model names.
#' @export
model_spec <- function(demand, conduct) {
    if (missing(demand) || length(demand) != 1L || is.na(demand)) {
        stop("'demand' must be a single model name.")
    }
    if (missing(conduct) || length(conduct) != 1L || is.na(conduct)) {
        stop("'conduct' must be a single conduct name.")
    }

    demand <- .normalize_demand_name(demand)
    conduct <- .normalize_conduct_name(conduct)
    entry <- .model_registry_entry(demand, conduct)
    if (is.null(entry)) {
        stop(demand, " / ", conduct, " currently not supported.")
    }

    structure(
        list(
            demand = demand,
            conduct = conduct,
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


#' @export
print.antitrust_model_spec <- function(x, ...) {
    cat("antitrust model specification\n")
    cat("  demand:  ", x$demand, "\n", sep = "")
    cat("  conduct: ", x$conduct, "\n", sep = "")
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
        bargaining2 = "bargaining2nd"
    )
    if (conduct %in% names(aliases)) aliases[[conduct]] else conduct
}


.model_registry_key <- function(demand, conduct) {
    paste(demand, conduct, sep = "::")
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
        list(id = "blp::bertrand", demand = "blp", conduct = "bertrand",
             class = "LogitBLP", calibrator = "sim", calibrate = FALSE,
             specify = TRUE, simulate = TRUE),
        list(id = "blp::cournot", demand = "blp", conduct = "cournot",
             class = "CournotBLP", calibrator = "sim", calibrate = FALSE,
             specify = TRUE, simulate = TRUE)
    )
    registry <- setNames(entries, vapply(entries, `[[`, character(1), "id"))
    function() registry
})


.model_registry_entry <- function(demand, conduct) {
    .model_registry()[[.model_registry_key(demand, conduct)]]
}


.model_registry_supports <- function(spec, operation) {
    entry <- .model_registry_entry(spec$demand, spec$conduct)
    !is.null(entry) && isTRUE(entry[[operation]])
}
