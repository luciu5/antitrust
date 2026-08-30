# Shared QA helpers.  The warning recorder deliberately muffles only after it
# has recorded the condition; no warning emitted by a formal test is silently
# discarded.

source(testthat::test_path("fixtures", "synthetic-fixtures.R"))

.qa_state <- new.env(parent = emptyenv())
.qa_state$warnings <- data.frame(
    context = character(), class = character(), message = character(),
    stringsAsFactors = FALSE
)
.qa_state$messages <- data.frame(
    context = character(), message = character(), stringsAsFactors = FALSE
)

qa_classify_warning <- function(message) {
    patterns <- c(
        "solver-nonconvergence" = "nonlinear solver.*(not|may not|did not).*converg|may not have successfully converged",
        "calibration-fit" = "optimizer may not have found a good solution|optimizer may not have successfully|Calibration routine may not have converged",
        "boundary-normalization" = "outside share is close to|Normalizing relative to largest good",
        "second-order-condition" = "Hessian.*(positive definite|maximize profits)",
        "unidentified-nest" = "nest.*(not identified|singleton)",
        "default-assumption" = "Defaulting to|Setting .* equal to|does not contain either",
        "ignored-input" = "is only used for BLP|are not used for BLP|Ignoring",
        "economic-input" = "positive values of 'mcDelta'|Negative .*marginal costs|intercepts are negative|own-slope coefficients are positive|Matrix of demand slopes coefficients is not symmetric",
        "diagnostic-fallback" = "Reporting shares instead of quantities|Calculating CV as|calcQuantities.*yielded all NAs",
        "no-merger-identity" = "ownerPost.*ownerPre.*same"
    )
    hit <- names(patterns)[vapply(patterns, function(pattern) {
        grepl(pattern, message, ignore.case = TRUE, perl = TRUE)
    }, logical(1))]
    if (length(hit)) hit[[1]] else "unclassified-review-required"
}

qa_record_warning <- function(context, message, disposition = NULL) {
    cls <- if (is.null(disposition)) qa_classify_warning(message) else disposition
    row <- data.frame(
        context = as.character(context), class = as.character(cls),
        message = as.character(message), stringsAsFactors = FALSE
    )
    .qa_state$warnings <- rbind(.qa_state$warnings, row)
    log_path <- Sys.getenv("ANTITRUST_QA_WARNING_LOG", "")
    if (nzchar(log_path)) {
        dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
        write.table(row, log_path, sep = ",", row.names = FALSE,
                    col.names = !file.exists(log_path), append = file.exists(log_path),
                    qmethod = "double")
    }
    invisible(cls)
}

qa_capture <- function(expr, context, disposition = NULL) {
    warnings <- character()
    messages <- character()
    value <- withCallingHandlers(
        force(expr),
        warning = function(condition) {
            message_text <- conditionMessage(condition)
            warnings <<- c(warnings, message_text)
            qa_record_warning(context, message_text, disposition)
            invokeRestart("muffleWarning")
        },
        message = function(condition) {
            message_text <- conditionMessage(condition)
            messages <<- c(messages, message_text)
            .qa_state$messages <- rbind(
                .qa_state$messages,
                data.frame(context = as.character(context), message = message_text,
                           stringsAsFactors = FALSE)
            )
            invokeRestart("muffleMessage")
        }
    )
    list(value = value, warnings = warnings, messages = messages)
}

qa_value <- function(expr, context, disposition = NULL) {
    qa_capture(expr, context = context, disposition = disposition)$value
}

qa_warnings <- function() .qa_state$warnings
qa_messages <- function() .qa_state$messages

qa_expect_error <- function(expr, pattern, context = "validation") {
    captured <- tryCatch(
        qa_capture(expr, context),
        error = function(condition) condition
    )
    testthat::expect_s3_class(captured, "condition")
    testthat::expect_match(conditionMessage(captured), pattern)
    invisible(captured)
}

qa_assert_finite <- function(x, message = "value must be finite") {
    testthat::expect_true(all(is.finite(x)), info = message)
    invisible(x)
}

qa_assert_close <- function(actual, expected, tolerance = 1e-7, message = NULL) {
    testthat::expect_equal(unname(actual), unname(expected), tolerance = tolerance,
                          info = message)
    invisible(actual)
}
