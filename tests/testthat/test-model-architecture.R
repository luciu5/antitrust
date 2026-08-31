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

refactor_ces_result_parity <- function(actual, expected) {
    expect_equal(class(actual), class(expected))
    expect_equal(actual@slopes, expected@slopes, tolerance = 1e-8)
    expect_equal(actual@mcPre, expected@mcPre, tolerance = 1e-8)
    expect_equal(actual@mcPost, expected@mcPost, tolerance = 1e-8)
    expect_equal(actual@pricePre, expected@pricePre, tolerance = 1e-8)
    expect_equal(actual@pricePost, expected@pricePost, tolerance = 1e-8)
    expect_equal(calcShares(actual, TRUE), calcShares(expected, TRUE), tolerance = 1e-8)
    expect_equal(calcShares(actual, FALSE), calcShares(expected, FALSE), tolerance = 1e-8)
    expect_equal(calcMargins(actual, TRUE), calcMargins(expected, TRUE), tolerance = 1e-8)
    expect_equal(calcMargins(actual, FALSE), calcMargins(expected, FALSE), tolerance = 1e-8)
    expect_equal(elast(actual, TRUE), elast(expected, TRUE), tolerance = 1e-8)
    expect_equal(elast(actual, FALSE), elast(expected, FALSE), tolerance = 1e-8)
    expect_equal(diversion(actual, TRUE), diversion(expected, TRUE), tolerance = 1e-8)
    expect_equal(diversion(actual, FALSE), diversion(expected, FALSE), tolerance = 1e-8)
    expect_equal(CV(actual), CV(expected), tolerance = 1e-8)
    invisible(actual)
}

test_that("CES Bertrand and Cournot retain distinct model-specific parity", {
    common <- list(
        prices = c(1.5, 1.8, 2), shares = c(.30, .20, .10),
        margins = c(.30, .25, .20), ownerPre = c("A", "B", "C"),
        ownerPost = c("A", "A", "C"), insideSize = 100
    )
    for (conduct in c("bertrand", "cournot")) {
        old <- qa_value(do.call(
            if (conduct == "bertrand") ces else ces.cournot,
            common
        ), paste("legacy CES", conduct))
        fit <- qa_value(do.call(calibrate, c(
            list(demand = "ces", conduct = conduct),
            common[names(common) != "ownerPost"]
        )), paste("calibrate CES", conduct))
        actual <- qa_value(simulate(fit, ownerPost = common$ownerPost),
                           paste("simulate CES", conduct))
        refactor_ces_result_parity(actual, old)
    }
})

test_that("legacy sim routes CES Bertrand and Cournot through the fit pipeline", {
    common <- list(
        prices = c(1.5, 1.8, 2), shares = c(.30, .20, .10),
        margins = c(.30, .25, .20), ownerPre = c("A", "B", "C"),
        ownerPost = c("A", "A", "C"), insideSize = 100,
        demand.param = list(gamma = 2, meanval = c(1, .8, .6), shareInside = .7)
    )
    legacy_sim <- getFromNamespace(".sim_legacy", "antitrust")
    for (conduct in c("bertrand", "cournot")) {
        old <- qa_value(do.call(legacy_sim, c(
            common, list(supply = conduct, demand = "CES")
        )), paste("legacy supplied CES", conduct))
        actual <- qa_value(do.call(sim, c(
            common, list(supply = conduct, demand = "CES")
        )), paste("compatibility supplied CES", conduct))
        refactor_ces_result_parity(actual, old)
    }
})

test_that("nested Logit and CES models retain model-specific calibration parity", {
    common <- list(
        prices = c(1.8, 2, 2.2, 2.5),
        shares = c(.30, .20, .10, .05),
        margins = c(.30, .25, .20, .15),
        ownerPre = c("A", "B", "C", "D"),
        ownerPost = c("A", "A", "C", "D"),
        nests = c("N1", "N1", "N2", "N2"),
        insideSize = 100
    )
    cases <- list(
        list(demand = "logit_nests", constructor = logit.nests,
             parameters = list(alpha = -1, meanval = c(.4, .2, .1, .05), sigma = .7),
             parmsStart = c(-1, .7)),
        list(demand = "ces_nests", constructor = ces.nests,
             parameters = list(gamma = 2, meanval = c(1, .8, .6, .5), sigma = 2.5,
                               shareInside = .7),
             parmsStart = c(1.5, 2.5))
    )
    for (case in cases) {
        old <- qa_value(do.call(case$constructor, c(
            common, list(parmsStart = case$parmsStart)
        )), paste("legacy", case$demand))
        fit <- qa_value(do.call(calibrate, c(
            list(demand = case$demand, conduct = "bertrand"),
            common[names(common) != "ownerPost"],
            list(parmsStart = case$parmsStart)
        )), paste("calibrate", case$demand))
        actual <- qa_value(simulate(fit, ownerPost = common$ownerPost),
                           paste("simulate", case$demand))
        expect_equal(class(actual), class(old))
        expect_equal(actual@slopes, old@slopes, tolerance = 1e-8)
        expect_equal(actual@mcPre, old@mcPre, tolerance = 1e-8)
        expect_equal(actual@mcPost, old@mcPost, tolerance = 1e-8)
        expect_equal(actual@pricePre, old@pricePre, tolerance = 1e-8)
        expect_equal(actual@pricePost, old@pricePost, tolerance = 1e-8)
        expect_equal(calcShares(actual, FALSE), calcShares(old, FALSE), tolerance = 1e-8)
        expect_equal(calcMargins(actual, FALSE), calcMargins(old, FALSE), tolerance = 1e-8)
    }
})

test_that("legacy sim routes nested Bertrand models through specify and simulate", {
    common <- list(
        prices = c(1.8, 2, 2.2, 2.5),
        shares = c(.30, .20, .10, .05),
        margins = c(.30, .25, .20, .15),
        ownerPre = c("A", "B", "C", "D"),
        ownerPost = c("A", "A", "C", "D"),
        nests = c("N1", "N1", "N2", "N2"),
        insideSize = 100
    )
    cases <- list(
        list(demand = "LogitNests", parameters = list(
            alpha = -1, meanval = c(.4, .2, .1, .05), sigma = .7
        ), parmsStart = c(-1, .7)),
        list(demand = "CESNests", parameters = list(
            gamma = 2, meanval = c(1, .8, .6, .5), sigma = c(.7, .8),
            shareInside = .7
        ))
    )
    legacy_sim <- getFromNamespace(".sim_legacy", "antitrust")
    for (case in cases) {
        args <- c(common, list(
            supply = "bertrand", demand = case$demand,
            demand.param = case$parameters
        ))
        old <- qa_value(do.call(legacy_sim, args), paste("legacy", case$demand))
        actual <- qa_value(do.call(sim, args), paste("compatibility", case$demand))
        expect_equal(class(actual), class(old))
        expect_equal(actual@pricePost, old@pricePost, tolerance = 1e-8)
        expect_equal(actual@mcPost, old@mcPost, tolerance = 1e-8)
    }
})

test_that("LogitCap calibration and supplied parameters use the shared simulation boundary", {
    prices <- c(2, 2.2, 2.5)
    shares <- c(.25, .20, .15)
    margins <- c(.30, .28, .25)
    capacities <- c(25, 20, 15)
    ownerPre <- c("A", "B", "C")
    ownerPost <- c("A", "A", "C")

    old <- qa_value(logit.cap(
        prices = prices, shares = shares, margins = margins,
        capacitiesPre = capacities, capacitiesPost = capacities,
        insideSize = 100, ownerPre = ownerPre, ownerPost = ownerPost
    ), "legacy LogitCap")
    fit <- qa_value(calibrate(
        demand = "logit_cap", conduct = "bertrand",
        prices = prices, shares = shares, margins = margins,
        ownerPre = ownerPre, capacitiesPre = capacities,
        capacitiesPost = capacities, insideSize = 100
    ), "calibrate LogitCap")
    actual <- qa_value(simulate(
        fit, ownerPost = ownerPost, capacitiesPost = capacities
    ), "simulate LogitCap")

    expect_s4_class(actual, "LogitCap")
    expect_equal(actual@slopes, old@slopes, tolerance = 1e-8)
    expect_equal(actual@mcPre, old@mcPre, tolerance = 1e-8)
    expect_equal(actual@mcPost, old@mcPost, tolerance = 1e-8)
    expect_equal(actual@pricePre, old@pricePre, tolerance = 1e-8)
    expect_equal(actual@pricePost, old@pricePost, tolerance = 1e-8)
    expect_equal(calcQuantities(actual, FALSE), calcQuantities(old, FALSE), tolerance = 1e-8)

    parameters <- list(
        alpha = old@slopes$alpha,
        meanval = old@slopes$meanval,
        mktSize = 100
    )
    legacy_sim <- getFromNamespace(".sim_legacy", "antitrust")
    supplied_old <- qa_value(legacy_sim(
        prices = prices, demand = "LogitCap", supply = "bertrand",
        demand.param = parameters, capacities = capacities, margins = margins,
        ownerPre = ownerPre, ownerPost = ownerPost
    ), "legacy supplied LogitCap")
    supplied_actual <- qa_value(sim(
        prices = prices, demand = "LogitCap", supply = "bertrand",
        demand.param = parameters, capacities = capacities, margins = margins,
        ownerPre = ownerPre, ownerPost = ownerPost
    ), "compatibility supplied LogitCap")
    expect_equal(supplied_actual@pricePost, supplied_old@pricePost, tolerance = 1e-8)
    expect_equal(supplied_actual@mcPost, supplied_old@mcPost, tolerance = 1e-8)
})

test_that("BLP supplied parameters reuse the existing Bertrand and Cournot solvers", {
    parameters <- qa_fixture_blp_parameters()
    common <- list(
        prices = c(2, 2.2, 2.5), shares = c(.35, .25, .20),
        ownerPre = c("A", "B", "C"), ownerPost = c("A", "A", "C"),
        insideSize = 100
    )
    legacy_sim <- getFromNamespace(".sim_legacy", "antitrust")
    for (conduct in c("bertrand", "cournot")) {
        old <- qa_value(do.call(legacy_sim, c(
            common, list(supply = conduct, demand = "BLP",
                         demand.param = parameters)
        )), paste("legacy BLP", conduct))
        fit <- qa_value(do.call(specify, c(
            list(demand = "blp", conduct = conduct,
                 parameters = parameters),
            common[names(common) != "ownerPost"]
        )), paste("specify BLP", conduct))
        actual <- qa_value(simulate(fit, ownerPost = common$ownerPost),
                           paste("simulate BLP", conduct))
        expect_equal(class(actual), class(old))
        expect_equal(actual@slopes, old@slopes, tolerance = 1e-8)
        expect_equal(actual@mcPre, old@mcPre, tolerance = 1e-8)
        expect_equal(actual@mcPost, old@mcPost, tolerance = 1e-8)
        expect_equal(actual@pricePre, old@pricePre, tolerance = 1e-8)
        expect_equal(actual@pricePost, old@pricePost, tolerance = 1e-8)
    }
})
