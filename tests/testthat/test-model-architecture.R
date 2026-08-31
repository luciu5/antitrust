refactor_logit_args <- function(conduct = "bertrand", ...) {
    c(
        list(
            prices = c(2, 2.2, 2.5),
            shares = c(.35, .25, .20),
            margins = c(.40, .35, .30),
            ownerPre = c("A", "B", "C"),
            ownerPost = c("A", "A", "C")
        ),
        list(...)
    )
}

refactor_logit_constructor <- function(conduct = "bertrand", ...) {
    args <- refactor_logit_args(conduct, ...)
    if (identical(conduct, "bertrand")) {
        do.call(logit, args)
    } else {
        do.call(logit.cournot, args)
    }
}

expect_logit_result_parity <- function(actual, expected) {
    expect_equal(class(actual), class(expected))
    expect_equal(actual@slopes, expected@slopes, tolerance = 1e-9)
    expect_equal(actual@mcPre, expected@mcPre, tolerance = 1e-9)
    expect_equal(actual@mcPost, expected@mcPost, tolerance = 1e-9)
    expect_equal(actual@pricePre, expected@pricePre, tolerance = 1e-9)
    expect_equal(actual@pricePost, expected@pricePost, tolerance = 1e-9)
    expect_equal(actual@ownerPre, expected@ownerPre, tolerance = 1e-12)
    expect_equal(actual@ownerPost, expected@ownerPost, tolerance = 1e-12)
    expect_equal(calcShares(actual, TRUE), calcShares(expected, TRUE), tolerance = 1e-9)
    expect_equal(calcShares(actual, FALSE), calcShares(expected, FALSE), tolerance = 1e-9)
    expect_equal(calcMargins(actual, TRUE), calcMargins(expected, TRUE), tolerance = 1e-9)
    expect_equal(calcMargins(actual, FALSE), calcMargins(expected, FALSE), tolerance = 1e-9)
    expect_equal(elast(actual, TRUE), elast(expected, TRUE), tolerance = 1e-9)
    expect_equal(elast(actual, FALSE), elast(expected, FALSE), tolerance = 1e-9)
    expect_equal(diversion(actual, TRUE), diversion(expected, TRUE), tolerance = 1e-9)
    expect_equal(diversion(actual, FALSE), diversion(expected, FALSE), tolerance = 1e-9)
    expect_equal(upp(actual), upp(expected), tolerance = 1e-9)
    expect_equal(cmcr(actual), cmcr(expected), tolerance = 1e-7)
    expect_equal(CV(actual), CV(expected), tolerance = 1e-9)
    invisible(actual)
}

test_that("calibrate and simulate preserve Logit-Bertrand behavior", {
    old <- qa_value(refactor_logit_constructor(), "legacy Logit Bertrand")
    fit <- qa_value(calibrate(
        demand = "logit", conduct = "bertrand",
        prices = c(2, 2.2, 2.5), shares = c(.35, .25, .20),
        margins = c(.40, .35, .30), ownerPre = c("A", "B", "C")
    ), "calibrate Logit Bertrand")
    actual <- qa_value(simulate(
        fit, ownerPost = c("A", "A", "C")
    ), "simulate Logit Bertrand")

    expect_s4_class(fit, "AntitrustFit")
    expect_equal(fit@spec$demand, "logit")
    expect_equal(fit@spec$conduct, "bertrand")
    expect_equal(fit@parameters, old@slopes, tolerance = 1e-9)
    expect_equal(fit@diagnostics$status, "completed")
    expect_logit_result_parity(actual, old)
})

test_that("calibrate and simulate preserve Logit-Cournot behavior", {
    old <- qa_value(refactor_logit_constructor("cournot"), "legacy Logit Cournot")
    fit <- qa_value(calibrate(
        model_spec("Logit", "Cournot"),
        prices = c(2, 2.2, 2.5), shares = c(.35, .25, .20),
        margins = c(.40, .35, .30), ownerPre = c("A", "B", "C")
    ), "calibrate Logit Cournot")
    actual <- qa_value(simulate(
        fit, ownerPost = c("A", "A", "C")
    ), "simulate Logit Cournot")

    expect_s4_class(actual, "LogitCournot")
    expect_equal(fit@parameters, old@slopes, tolerance = 1e-9)
    expect_logit_result_parity(actual, old)
})

test_that("a calibrated Logit can be reused for distinct counterfactuals", {
    common <- list(
        prices = c(2, 2.2, 2.5), shares = c(.35, .25, .20),
        margins = c(.40, .35, .30), ownerPre = c("A", "B", "C")
    )
    fit <- qa_value(do.call(calibrate, c(
        list(demand = "logit", conduct = "bertrand"), common
    )), "reusable Logit calibration")

    for (owners in list(c("A", "A", "C"), c("A", "B", "B"))) {
        actual <- qa_value(simulate(fit, ownerPost = owners), "reused Logit simulation")
        old <- qa_value(do.call(logit, c(common, list(ownerPost = owners))),
                        "legacy repeated Logit simulation")
        expect_logit_result_parity(actual, old)
    }
})

test_that("Logit-Bertrand scenario controls and AG solver retain parity", {
    common <- list(
        prices = c(2, 2.2, 2.5), shares = c(.35, .25, .20),
        margins = c(.40, .35, .30), ownerPre = c("A", "B", "C"),
        ownerPost = c("A", "A", "C"), mcDelta = c(.1, 0, 0),
        subset = c(TRUE, TRUE, FALSE)
    )
    old <- qa_value(do.call(logit, c(common, list(solver = "ag"))),
                    "legacy Logit AG scenario")
    fit <- qa_value(do.call(calibrate, c(
        list(demand = "logit", conduct = "bertrand", solver = "ag"),
        common[setdiff(names(common), c("ownerPost", "mcDelta", "subset"))]
    )), "calibrate Logit AG scenario")
    actual <- qa_value(simulate(
        fit, ownerPost = common$ownerPost, mcDelta = common$mcDelta,
        subset = common$subset
    ), "simulate Logit AG scenario")

    expect_logit_result_parity(actual, old)
})

test_that("specify loads supplied Logit parameters into the shared pipeline", {
    parameters <- list(alpha = -1.2, meanval = c(.5, .3, .1))
    common <- list(
        prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .20),
        margins = c(.40, .35, .30),
        ownerPre = c("A", "B", "C"),
        ownerPost = c("A", "A", "C"),
        insideSize = 100
    )
    legacy_sim <- getFromNamespace(".sim_legacy", "antitrust")
    old <- qa_value(do.call(legacy_sim, c(common, list(
        supply = "bertrand", demand = "Logit", demand.param = parameters
    ))), "legacy supplied-parameter Logit")
    fit <- qa_value(do.call(specify, c(
        list(demand = "logit", conduct = "bertrand", parameters = parameters),
        common[names(common) %in% c("prices", "shares", "margins", "ownerPre", "insideSize")]
    )), "specify Logit Bertrand")
    actual <- qa_value(simulate(fit, ownerPost = common$ownerPost),
                       "simulate specified Logit Bertrand")

    expect_equal(fit@diagnostics$source, "specified")
    expect_logit_result_parity(actual, old)
})

test_that("legacy sim remains a wrapper for migrated Logit combinations", {
    common <- list(
        prices = c(2, 2.2, 2.5), shares = c(.35, .25, .20),
        margins = c(.40, .35, .30), ownerPre = c("A", "B", "C"),
        ownerPost = c("A", "A", "C"), insideSize = 100,
        demand.param = list(alpha = -1.2, meanval = c(.5, .3, .1))
    )
    legacy_sim <- getFromNamespace(".sim_legacy", "antitrust")
    old <- qa_value(do.call(legacy_sim, c(common, list(
        supply = "bertrand", demand = "Logit"
    ))), "legacy sim Bertrand")
    actual <- qa_value(do.call(sim, c(common, list(
        supply = "bertrand", demand = "Logit"
    ))), "compatibility sim Bertrand")
    expect_logit_result_parity(actual, old)
})
