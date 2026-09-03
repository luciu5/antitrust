test_that("the registry covers generalized and migrated model combinations", {
    registry <- supportedModels()
    expected <- data.frame(
        demand = c(
            "linear", "aids", "loglin", "logit", "logit", "logit",
            "logit", "logit", "ces", "ces", "ces", "ces", "ces",
            "logit_nests", "ces_nests", "logit_cap", "pcaids", "pcaids_nests", "blp", "blp",
            "auction2nd_cap", "linear", "loglin", "linear", "loglin", "logit",
            "logit_nests", "logit", "logit", "logit", "ces", "ces", "logit_nests", "logit_cap",
            "logit", "ces", "logit", "ces"
        ),
        conduct = c(
            "bertrand", "bertrand", "bertrand", "bertrand", "cournot",
            "auction2nd", "bargaining", "bargaining2nd", "bertrand",
            "cournot", "auction2nd", "bargaining", "bargaining2nd",
            "bertrand", "bertrand", "bertrand", "bertrand", "bertrand", "bertrand", "cournot",
            "auction2nd", "cournot", "cournot", "stackelberg", "stackelberg", "vertical_bargaining", "vertical_bargaining", "vertical_bargaining",
            "bertrand", "cournot", "bertrand", "cournot", "bertrand",
            "bertrand", "auction2nd", "auction2nd", "bargaining", "bargaining"
        ),
        variant = c(rep("standard", 27), "auction2nd", rep("alm", 10)),
        class = c(
            "Linear", "AIDS", "LogLin", "Logit", "LogitCournot",
            "Auction2ndLogit", "BargainingLogit", "Bargaining2ndLogit",
            "CES", "CESCournot", "Auction2ndCES", "BargainingCES",
            "Bargaining2ndCES", "LogitNests", "CESNests", "LogitCap",
            "PCAIDS", "PCAIDSNests", "LogitBLP", "CournotBLP", "Auction2ndCap",
            "Cournot", "Cournot", "Stackelberg", "Stackelberg", "VertBargBertLogit", "VertBargBertLogitNests", "VertBarg2ndLogit",
            "LogitALM", "LogitCournotALM",
            "CESALM", "CESCournotALM", "LogitNestsALM", "LogitCapALM",
            "Auction2ndLogitALM", "Auction2ndCESALM",
            "BargainingLogitALM", "BargainingCESALM"
        ),
        calibrator = c(
            "linear", "aids", "loglinear", "logit", "logit.cournot",
            "auction2nd.logit", "bargaining.logit", "bargaining2nd.logit",
            "ces", "ces.cournot", "auction2nd.ces", "bargaining.ces",
            "bargaining2nd.ces", "logit.nests", "ces.nests", "logit.cap",
            "pcaids", "pcaids.nests", "sim", "sim", "auction2nd.cap", "cournot", "cournot",
            "stackelberg", "stackelberg",
            "vertical.barg", "vertical.barg", "vertical.barg",
            "logit.alm", "logit.cournot.alm", "ces.alm",
            "ces.cournot.alm", "logit.nests.alm", "logit.cap.alm",
            "auction2nd.logit.alm", "auction2nd.ces.alm",
            "bargaining.logit.alm", "bargaining.ces.alm"
        ),
        calibrate = c(rep(TRUE, 18), FALSE, FALSE, rep(TRUE, 18)),
        specify = c(rep(TRUE, 16), FALSE, FALSE, TRUE, TRUE, rep(FALSE, 18)),
        simulate = rep(TRUE, 38),
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
    expect_equal(model_spec("PCAIDS.Nests", "Bertrand")$demand, "pcaids_nests")
    expect_equal(model_spec("auction2nd.cap", "auction2nd")$id,
                 "auction2nd_cap::auction2nd")
    expect_equal(model_spec("logit", "auction_2nd")$conduct, "auction2nd")
    expect_equal(model_spec("linear", "stack")$conduct, "stackelberg")
    expect_equal(model_spec("logit", "vertical")$conduct, "vertical_bargaining")
    expect_equal(model_spec("logit_nests", "vertical")$id,
                 "logit_nests::vertical_bargaining")
    expect_equal(model_spec("logit", "vertical_bargaining", variant = "second_score")$id,
                 "logit::vertical_bargaining::auction2nd")
    expect_equal(model_spec("logit", "bertrand", variant = "ALM")$variant, "alm")
    expect_equal(model_spec("LogitALM", "Bertrand")$id,
                 "logit::bertrand::alm")

    qa_expect_error(
        model_spec("linear", "bargaining"),
        "currently not supported",
        "unsupported model specification"
    )
    qa_expect_error(
        model_spec("logit", "bertrand", variant = "unknown"),
        "currently not supported",
        "unsupported model variant"
    )
})

test_that("demand transition registry is explicit and non-duplicated", {
    transitions <- getFromNamespace(
        ".model_transition_registry", "antitrust"
    )()
    keys <- vapply(
        transitions,
        function(entry) paste(entry$from, entry$to, sep = "->"),
        character(1)
    )
    kinds <- vapply(transitions, `[[`, character(1), "kind")

    expect_false(anyDuplicated(keys) > 0L)
    expect_true(all(kinds %in% c(
        "structural-restriction", "algebraic-translation",
        "conditional-translation", "first-order-linearization",
        "first-order-loglinearization"
    )))
    expect_true(all(vapply(
        transitions,
        function(entry) is.character(entry$required_arguments),
        logical(1)
    )))
    expect_true(any(keys == "aids->linear"))
    expect_true(any(keys == "aids->loglin"))
    expect_true(any(keys == "logit_nests->logit"))
})
