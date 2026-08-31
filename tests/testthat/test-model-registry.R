test_that("the registry covers every currently supported sim combination", {
    registry <- supportedModels()
    expected <- data.frame(
        demand = c(
            "linear", "aids", "loglin", "logit", "logit", "logit",
            "logit", "logit", "ces", "ces", "ces", "ces", "ces",
            "logit_nests", "ces_nests", "logit_cap", "blp", "blp"
        ),
        conduct = c(
            "bertrand", "bertrand", "bertrand", "bertrand", "cournot",
            "auction2nd", "bargaining", "bargaining2nd", "bertrand",
            "cournot", "auction2nd", "bargaining", "bargaining2nd",
            "bertrand", "bertrand", "bertrand", "bertrand", "cournot"
        ),
        class = c(
            "Linear", "AIDS", "LogLin", "Logit", "LogitCournot",
            "Auction2ndLogit", "BargainingLogit", "Bargaining2ndLogit",
            "CES", "CESCournot", "Auction2ndCES", "BargainingCES",
            "Bargaining2ndCES", "LogitNests", "CESNests", "LogitCap",
            "LogitBLP", "CournotBLP"
        ),
        calibrator = c(
            "linear", "aids", "loglinear", "logit", "logit.cournot",
            "auction2nd.logit", "bargaining.logit", "bargaining2nd.logit",
            "ces", "ces.cournot", "auction2nd.ces", "bargaining.ces",
            "bargaining2nd.ces", "logit.nests", "ces.nests", "logit.cap",
            "sim", "sim"
        ),
        calibrate = c(rep(TRUE, 16), FALSE, FALSE),
        specify = rep(TRUE, 18),
        simulate = rep(TRUE, 18),
        stringsAsFactors = FALSE,
        row.names = row.names(registry)
    )
    expect_equal(registry, expected)
})

test_that("model specifications normalize names and reject unsupported combinations", {
    spec <- model_spec("Logit", "Cournot")
    expect_s3_class(spec, "antitrust_model_spec")
    expect_equal(spec$demand, "logit")
    expect_equal(spec$conduct, "cournot")
    expect_equal(spec$id, "logit::cournot")

    expect_equal(model_spec("CESNests", "Bertrand")$demand, "ces_nests")
    expect_equal(model_spec("logit-cap", "bertrand")$demand, "logit_cap")
    expect_equal(model_spec("logit", "auction_2nd")$conduct, "auction2nd")

    qa_expect_error(
        model_spec("linear", "cournot"),
        "currently not supported",
        "unsupported model specification"
    )
})
