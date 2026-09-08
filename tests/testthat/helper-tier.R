# Test-tier controls.  The fast tier is the default for ordinary development
# and pull-request CI.  The extended tier is an additive nightly/manual gate
# for expensive legacy-parity, solver-matrix, and compound-scenario checks.

qa_test_tier <- function() {
    tier <- tolower(trimws(Sys.getenv("ANTITRUST_TEST_TIER", "fast")))
    if (!nzchar(tier)) tier <- "fast"
    if (!tier %in% c("fast", "extended")) {
        stop("ANTITRUST_TEST_TIER must be 'fast' or 'extended'.")
    }
    tier
}


qa_extended_enabled <- function() identical(qa_test_tier(), "extended")


qa_skip_if_not_extended <- function() {
    if (!qa_extended_enabled()) {
        testthat::skip(
            "extended parity test; rerun with ANTITRUST_TEST_TIER=extended"
        )
    }
    invisible(TRUE)
}
