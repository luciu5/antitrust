# Test-tier controls.  Tiers are additive: an extended or nightly run includes
# every lower-tier check.  Keep the routing here so block-level gates use one
# validated vocabulary across the suite.

qa_test_tier <- function() {
    tier <- tolower(trimws(Sys.getenv("ANTITRUST_TEST_TIER", "fast")))
    if (!nzchar(tier)) tier <- "fast"
    if (!tier %in% c("fast", "extended", "nightly")) {
        stop("ANTITRUST_TEST_TIER must be 'fast', 'extended', or 'nightly'.")
    }
    tier
}


qa_tier_enabled <- function(required) {
    required <- tolower(required)
    if (length(required) != 1L ||
        !required %in% c("fast", "extended", "nightly")) {
        stop("required test tier must be 'fast', 'extended', or 'nightly'.")
    }
    match(qa_test_tier(), c("fast", "extended", "nightly")) >=
        match(required, c("fast", "extended", "nightly"))
}


qa_extended_enabled <- function() qa_tier_enabled("extended")


qa_nightly_enabled <- function() identical(qa_test_tier(), "nightly")


qa_skip_unless_tier <- function(required, reason = NULL) {
    required <- tolower(required)
    if (!qa_tier_enabled(required)) {
        if (is.null(reason)) {
            reason <- paste0(required,
                             "-tier test; rerun with ANTITRUST_TEST_TIER=",
                             required)
        }
        testthat::skip(reason)
    }
    invisible(TRUE)
}


qa_skip_if_not_extended <- function() {
    qa_skip_unless_tier(
        "extended",
        "extended-tier test; rerun with ANTITRUST_TEST_TIER=extended"
    )
}


qa_skip_if_not_nightly <- function() {
    qa_skip_unless_tier(
        "nightly",
        "nightly-tier test; rerun with ANTITRUST_TEST_TIER=nightly"
    )
}
