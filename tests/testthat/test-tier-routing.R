test_that("test tiers are ordered and validated centrally", {
    previous <- Sys.getenv("ANTITRUST_TEST_TIER", unset = NA_character_)
    on.exit({
        if (is.na(previous)) {
            Sys.unsetenv("ANTITRUST_TEST_TIER")
        } else {
            Sys.setenv(ANTITRUST_TEST_TIER = previous)
        }
    }, add = TRUE)

    Sys.setenv(ANTITRUST_TEST_TIER = "fast")
    expect_identical(qa_test_tier(), "fast")
    expect_true(qa_tier_enabled("fast"))
    expect_false(qa_tier_enabled("extended"))
    expect_false(qa_tier_enabled("nightly"))

    Sys.setenv(ANTITRUST_TEST_TIER = "extended")
    expect_true(qa_extended_enabled())
    expect_false(qa_nightly_enabled())
    expect_true(qa_tier_enabled("extended"))
    expect_false(qa_tier_enabled("nightly"))

    Sys.setenv(ANTITRUST_TEST_TIER = "nightly")
    expect_true(qa_extended_enabled())
    expect_true(qa_nightly_enabled())
    expect_true(qa_tier_enabled("nightly"))

    Sys.setenv(ANTITRUST_TEST_TIER = "unknown")
    expect_error(qa_test_tier(), "fast.*extended.*nightly")
})
