test_that("the registry covers generalized and migrated model combinations", {
    registry <- supportedModels()
    expected <- data.frame(
        demand = c(
            "linear", "aids", "loglin", "logit", "logit", "logit", "logit",
            "logit", "logit", "ces", "ces", "ces", "ces", "ces", "ces",
            "logit_nests", "ces_nests", "logit_cap", "pcaids", "pcaids_nests", "blp", "blp", "blp", "blp", "blp",
            "auction2nd_cap", "linear", "loglin", "linear", "loglin",
            "logit", "logit", "ces", "ces", "logit_nests", "logit_cap",
            "logit", "ces", "logit", "ces"
        ),
        conduct = c(
            "bertrand", "bertrand", "bertrand", "bertrand", "moncom", "cournot",
            "auction2nd", "bargaining", "bargaining2nd", "bertrand", "moncom",
            "cournot", "auction2nd", "bargaining", "bargaining2nd",
            "bertrand", "bertrand", "bertrand", "bertrand", "bertrand", "bertrand", "moncom", "cournot",
            "auction2nd", "bargaining", "auction2nd", "cournot", "cournot", "stackelberg", "stackelberg",
            "bertrand", "cournot", "bertrand", "cournot", "bertrand",
            "bertrand", "auction2nd", "auction2nd", "bargaining", "bargaining"
        ),
        variant = c(rep("standard", 30), rep("alm", 10)),
        class = c(
            "Linear", "AIDS", "LogLin", "Logit", "MonComLogit", "LogitCournot",
            "Auction2ndLogit", "BargainingLogit", "Bargaining2ndLogit",
            "CES", "MonComCES", "CESCournot", "Auction2ndCES", "BargainingCES",
            "Bargaining2ndCES", "LogitNests", "CESNests", "LogitCap",
            "PCAIDS", "PCAIDSNests", "LogitBLP", "MonComBLP", "CournotBLP", "Auction2ndBLP",
            "BargainingBLP", "Auction2ndCap",
            "Cournot", "Cournot", "Stackelberg", "Stackelberg",
            "LogitALM", "LogitCournotALM",
            "CESALM", "CESCournotALM", "LogitNestsALM", "LogitCapALM",
            "Auction2ndLogitALM", "Auction2ndCESALM",
            "BargainingLogitALM", "BargainingCESALM"
        ),
        calibrator = c(
            "linear", "aids", "loglinear", "logit", "moncom.logit", "logit.cournot",
            "auction2nd.logit", "bargaining.logit", "bargaining2nd.logit",
            "ces", "moncom.ces", "ces.cournot", "auction2nd.ces", "bargaining.ces",
            "bargaining2nd.ces", "logit.nests", "ces.nests", "logit.cap",
            "pcaids", "pcaids.nests", "blp", "blp", "blp", "blp", "blp", "auction2nd.cap", "cournot", "cournot",
            "stackelberg", "stackelberg",
            "logit.alm", "logit.cournot.alm", "ces.alm",
            "ces.cournot.alm", "logit.nests.alm", "logit.cap.alm",
            "auction2nd.logit.alm", "auction2nd.ces.alm",
            "bargaining.logit.alm", "bargaining.ces.alm"
        ),
        calibrate = rep(TRUE, 40),
        specify = c(rep(TRUE, 18), rep(FALSE, 2), rep(TRUE, 5), rep(FALSE, 15)),
        simulate = rep(TRUE, 40),
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
    expect_error(model_spec("logit", "vertical"), "vertical::model_spec")
    expect_error(model_spec("logit", "vertical_bargaining"),
                 "vertical::model_spec")
    expect_error(calibrate("logit", "vertical"), "vertical::calibrate")
    expect_error(specify("logit", "vertical"), "vertical::model_spec")
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

test_that("model_transition exposes the registered transition record", {
    from <- model_spec("logit", "bertrand")
    to <- model_spec("ces", "bertrand")

    transition <- model_transition(from, to)

    expect_type(transition, "list")
    expect_equal(transition$from, from$id)
    expect_equal(transition$to, to$id)
    expect_equal(transition$kind, "algebraic-translation")
    expect_equal(transition$required_arguments, "gamma")
    expect_equal(
        transition,
        getFromNamespace(".model_transition_entry", "antitrust")(from, to)
    )
})

test_that("model_transition requires valid antitrust model specifications", {
    spec <- model_spec("logit", "bertrand")

    expect_error(model_transition(list(), spec), "antitrust_model_spec")
    expect_error(model_transition(spec, list()), "antitrust_model_spec")
    malformed <- structure(
        list(demand = "logit", conduct = "bertrand",
             variant = "standard", id = "wrong"),
        class = c("antitrust_model_spec", "list")
    )
    expect_error(model_transition(malformed, spec), "valid antitrust_model_spec")
})
