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

test_that("legacy BLP aliases preserve their conduct mapping", {
    parameters <- qa_fixture_blp_parameters()
    common <- list(
        prices = c(2, 2.2, 2.5), shares = c(.35, .25, .20),
        ownerPre = c("A", "B", "C"), ownerPost = c("A", "A", "C"),
        insideSize = 100
    )
    explicit <- suppressWarnings(sim(
        prices = common$prices, shares = common$shares,
        demand = "BLP", demand.param = parameters, supply = "cournot",
        ownerPre = common$ownerPre, ownerPost = common$ownerPost,
        insideSize = common$insideSize
    ))
    alias <- suppressWarnings(sim(
        prices = common$prices, shares = common$shares,
        demand = "CournotBLP", demand.param = parameters,
        ownerPre = common$ownerPre, ownerPost = common$ownerPost,
        insideSize = common$insideSize
    ))
    expect_s4_class(alias, "CournotBLP")
    expect_equal(alias@pricePost, explicit@pricePost, tolerance = 1e-8)
})

test_that("legacy BLP simulation messages remain visible through sim", {
    parameters <- qa_fixture_blp_parameters()
    common <- list(
        prices = c(2, 2.2, 2.5), shares = c(.35, .25, .20),
        ownerPre = c("A", "B", "C"), ownerPost = c("A", "A", "C"),
        insideSize = 100
    )
    nonempty <- function(x) x[nzchar(trimws(x))]
    legacy_messages <- capture.output(suppressWarnings(getFromNamespace(
        ".sim_legacy", "antitrust")(
            prices = common$prices, shares = common$shares,
            demand = "BLP", demand.param = parameters,
            supply = "bertrand", ownerPre = common$ownerPre,
            ownerPost = common$ownerPost, insideSize = common$insideSize
        )
    ), type = "message")
    wrapper_messages <- capture.output(suppressWarnings(sim(
        prices = common$prices, shares = common$shares,
        demand = "BLP", demand.param = parameters, supply = "bertrand",
        ownerPre = common$ownerPre, ownerPost = common$ownerPost,
        insideSize = common$insideSize
    )), type = "message")
    expect_equal(nonempty(wrapper_messages), nonempty(legacy_messages))
})

test_that("ALM model variants calibrate and simulate through their legacy methods", {
    p <- c(2, 2.5, 3)
    s <- c(.40, .35, .25)
    m <- c(.45, .40, .35)
    owner_pre <- c("A", "B", "C")
    owner_post <- c("A", "A", "C")
    common <- list(
        prices = p, shares = s, margins = m,
        ownerPre = owner_pre, ownerPost = owner_post
    )
    cases <- list(
        logit_bertrand = list(
            demand = "logit", conduct = "bertrand",
            constructor = logit.alm,
            options = list(parmsStart = c(-.5, .1))
        ),
        logit_cournot = list(
            demand = "logit", conduct = "cournot",
            constructor = logit.cournot.alm,
            options = list(parmsStart = c(-.5, .1))
        ),
        ces_bertrand = list(
            demand = "ces", conduct = "bertrand",
            constructor = ces.alm,
            options = list(parmsStart = c(1, .1))
        ),
        ces_cournot = list(
            demand = "ces", conduct = "cournot",
            constructor = ces.cournot.alm,
            options = list(parmsStart = c(1, .1))
        ),
        logit_auction2nd = list(
            demand = "logit", conduct = "auction2nd",
            constructor = auction2nd.logit.alm,
            options = list(parmsStart = c(-.5, .1))
        ),
        ces_auction2nd = list(
            demand = "ces", conduct = "auction2nd",
            constructor = auction2nd.ces.alm,
            options = list()
        ),
        logit_bargaining = list(
            demand = "logit", conduct = "bargaining",
            constructor = bargaining.logit.alm,
            options = list(parmsStart = c(-.5, .1),
                           bargpowerPre = rep(.5, 3),
                           bargpowerPost = rep(.5, 3))
        ),
        ces_bargaining = list(
            demand = "ces", conduct = "bargaining",
            constructor = bargaining.ces.alm,
            options = list(parmsStart = c(1, .1),
                           bargpowerPre = rep(.5, 3),
                           bargpowerPost = rep(.5, 3))
        )
    )
    cases$logit_nests <- list(
        demand = "logit_nests", conduct = "bertrand",
        constructor = logit.nests.alm,
        options = list(nests = factor(c("N1", "N1", "N2", "N2")),
                       parmsStart = c(-.5, .5, .5)),
        prices = c(1.8, 2, 2.2, 2.5),
        shares = c(.35, .30, .20, .15),
        margins = c(.45, .40, .35, .30),
        ownerPre = c("A", "B", "C", "D"),
        ownerPost = c("A", "A", "C", "D")
    )
    cases$logit_cap <- list(
        demand = "logit_cap", conduct = "bertrand",
        constructor = logit.cap.alm,
        options = list(capacitiesPre = c(60, 60, 60),
                       capacitiesPost = c(60, 60, 60),
                       insideSize = 180,
                       parmsStart = c(-.5, .1))
    )

    for (case_name in names(cases)) {
        case <- cases[[case_name]]
        prices <- if (is.null(case$prices)) common$prices else case$prices
        shares <- if (is.null(case$shares)) common$shares else case$shares
        margins <- if (is.null(case$margins)) common$margins else case$margins
        owner_pre <- if (is.null(case$ownerPre)) common$ownerPre else case$ownerPre
        owner_post <- if (is.null(case$ownerPost)) common$ownerPost else case$ownerPost
        legacy_args <- c(
            list(prices = prices, shares = shares, margins = margins,
                 ownerPre = owner_pre, ownerPost = owner_post),
            case$options
        )
        old <- suppressWarnings(do.call(case$constructor, legacy_args))
        fit_args <- c(
            list(demand = model_spec(case$demand, case$conduct, "alm"),
                 prices = prices, shares = shares, margins = margins,
                 ownerPre = owner_pre),
            case$options
        )
        fit <- suppressWarnings(do.call(calibrate, fit_args))
        actual <- suppressWarnings(simulate(fit, ownerPost = owner_post))
        expect_true(is(actual, class(old)[[1]]), info = case_name)
        expect_equal(actual@slopes, old@slopes, tolerance = 1e-7,
                     info = case_name)
        expect_equal(actual@mcPre, old@mcPre, tolerance = 1e-7,
                     info = case_name)
        expect_equal(actual@mcPost, old@mcPost, tolerance = 1e-7,
                     info = case_name)
        expect_equal(actual@pricePre, old@pricePre, tolerance = 1e-7,
                     info = case_name)
        expect_equal(actual@pricePost, old@pricePost, tolerance = 1e-7,
                     info = case_name)
    }
})

test_that("second-score auction calibration and simulation retain model-specific parity", {
    cases <- list(
        list(
            demand = "logit",
            constructor = auction2nd.logit,
            prices = c(2, 2.2, 2.5),
            shares = c(.35, .25, .20),
            margins = c(.40, .35, .30)
        ),
        list(
            demand = "ces",
            constructor = auction2nd.ces,
            prices = c(1.5, 1.5, 1.5),
            shares = c(.50, .30, .20),
            margins = c(.60, .40, .30)
        )
    )
    owner_pre <- c("A", "B", "C")
    owner_post <- c("A", "A", "C")

    for (case in cases) {
        common <- list(
            prices = case$prices,
            shares = case$shares,
            margins = case$margins,
            ownerPre = owner_pre,
            ownerPost = owner_post
        )
        old <- qa_value(do.call(case$constructor, common),
                        paste("legacy auction", case$demand))
        fit <- qa_value(do.call(calibrate, c(
            list(demand = case$demand, conduct = "auction2nd"),
            common[names(common) != "ownerPost"]
        )), paste("calibrate auction", case$demand))
        actual <- qa_value(simulate(fit, ownerPost = owner_post),
                           paste("simulate auction", case$demand))

        expect_equal(fit@parameters, old@slopes, tolerance = 1e-8)
        expect_equal(class(actual), class(old))
        expect_equal(actual@mcPre, old@mcPre, tolerance = 1e-8)
        expect_equal(actual@mcPost, old@mcPost, tolerance = 1e-8)
        expect_equal(actual@pricePre, old@pricePre, tolerance = 1e-8)
        expect_equal(actual@pricePost, old@pricePost, tolerance = 1e-8)
        expect_equal(calcShares(actual, TRUE), calcShares(old, TRUE), tolerance = 1e-8)
        expect_equal(calcShares(actual, FALSE), calcShares(old, FALSE), tolerance = 1e-8)
        expect_equal(calcMargins(actual, TRUE), calcMargins(old, TRUE), tolerance = 1e-8)
        expect_equal(calcMargins(actual, FALSE), calcMargins(old, FALSE), tolerance = 1e-8)
        expect_equal(elast(actual, TRUE), elast(old, TRUE), tolerance = 1e-8)
        expect_equal(diversion(actual, TRUE), diversion(old, TRUE), tolerance = 1e-8)
    }
})

test_that("legacy sim routes second-score auction models through specify and simulate", {
    legacy_sim <- getFromNamespace(".sim_legacy", "antitrust")
    cases <- list(
        list(
            demand = "Logit",
            parameters = list(alpha = -1.2, meanval = c(.5, .3, .1)),
            prices = c(2, 2.2, 2.5),
            shares = c(.35, .25, .20),
            margins = c(.40, .35, .30)
        ),
        list(
            demand = "CES",
            parameters = list(gamma = 2, meanval = c(1, .8, .6), shareInside = .7),
            prices = c(1.5, 1.5, 1.5),
            shares = c(.50, .30, .20),
            margins = c(.60, .40, .30)
        )
    )
    for (case in cases) {
        common <- list(
            prices = case$prices,
            shares = case$shares,
            margins = case$margins,
            demand = case$demand,
            demand.param = case$parameters,
            supply = "auction2nd",
            ownerPre = c("A", "B", "C"),
            ownerPost = c("A", "A", "C")
        )
        old <- qa_value(do.call(legacy_sim, common),
                        paste("legacy supplied auction", case$demand))
        actual <- qa_value(do.call(sim, common),
                           paste("compatibility supplied auction", case$demand))
        expect_equal(class(actual), class(old))
        expect_equal(actual@slopes, old@slopes, tolerance = 1e-8)
        expect_equal(actual@mcPost, old@mcPost, tolerance = 1e-8)
        expect_equal(actual@pricePost, old@pricePost, tolerance = 1e-8)
    }
})

test_that("bargaining calibration and simulation retain model-specific parity", {
    cases <- list(
        list(
            demand = "logit",
            constructor = bargaining.logit,
            prices = c(2, 2.2, 2.5),
            shares = c(.35, .25, .20),
            margins = c(.40, .35, .30)
        ),
        list(
            demand = "ces",
            constructor = bargaining.ces,
            prices = c(1.5, 1.5, 1.5),
            shares = c(.50, .30, .20),
            margins = c(.30, .20, .15)
        ),
        list(
            demand = "logit",
            constructor = bargaining2nd.logit,
            prices = c(2, 2.2, 2.5),
            shares = c(.35, .25, .20),
            margins = c(.40, .35, .30)
        ),
        list(
            demand = "ces",
            constructor = bargaining2nd.ces,
            prices = c(1.5, 1.5, 1.5),
            shares = c(.50, .30, .20),
            margins = c(.30, .20, .15)
        )
    )
    barg_pre <- c(.5, .4, .6)
    barg_post <- c(.4, .3, .5)
    owner_pre <- c("A", "B", "C")
    owner_post <- c("A", "A", "C")

    for (case in cases) {
        common <- list(
            prices = case$prices,
            shares = case$shares,
            margins = case$margins,
            ownerPre = owner_pre,
            ownerPost = owner_post,
            bargpowerPre = barg_pre,
            bargpowerPost = barg_post,
            mcDelta = c(-.05, 0, 0)
        )
        old <- qa_value(do.call(case$constructor, common),
                        paste("legacy", case$demand, "bargaining"))
        fit_args <- list(
            demand = case$demand,
            conduct = "bargaining",
            prices = case$prices,
            shares = case$shares,
            margins = case$margins,
            ownerPre = owner_pre,
            bargpowerPre = barg_pre,
            bargpowerPost = barg_pre
        )
        ## The constructor is stored as a function object, so identify the
        ## conduct from its class name rather than relying on call-site names.
        fit_args$conduct <- if (identical(case$constructor, bargaining2nd.logit) ||
                                identical(case$constructor, bargaining2nd.ces)) {
            "bargaining2nd"
        } else {
            "bargaining"
        }
        fit <- qa_value(do.call(calibrate, fit_args),
                        paste("calibrate", case$demand, fit_args$conduct))
        actual <- qa_value(simulate(
            fit, ownerPost = owner_post, mcDelta = c(-.05, 0, 0),
            bargpowerPost = barg_post
        ), paste("simulate", case$demand, fit_args$conduct))

        expect_equal(fit@parameters, old@slopes, tolerance = 1e-8)
        expect_equal(class(actual), class(old))
        expect_equal(actual@mcPre, old@mcPre, tolerance = 1e-8)
        expect_equal(actual@mcPost, old@mcPost, tolerance = 1e-8)
        expect_equal(actual@pricePre, old@pricePre, tolerance = 1e-8)
        expect_equal(actual@pricePost, old@pricePost, tolerance = 1e-8)
        expect_equal(calcShares(actual, TRUE), calcShares(old, TRUE), tolerance = 1e-8)
        expect_equal(calcShares(actual, FALSE), calcShares(old, FALSE), tolerance = 1e-8)
        expect_equal(calcMargins(actual, TRUE), calcMargins(old, TRUE), tolerance = 1e-8)
        expect_equal(calcMargins(actual, FALSE), calcMargins(old, FALSE), tolerance = 1e-8)
        expect_equal(elast(actual, TRUE), elast(old, TRUE), tolerance = 1e-8)
        expect_equal(diversion(actual, TRUE), diversion(old, TRUE), tolerance = 1e-8)
    }
})

test_that("bargaining AG solver selection is retained by calibrate and simulate", {
    common <- list(
        prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .20),
        margins = c(.40, .35, .30),
        ownerPre = c("A", "B", "C"),
        ownerPost = c("A", "A", "C"),
        bargpowerPre = rep(.5, 3),
        bargpowerPost = rep(.4, 3),
        solver = "ag"
    )
    old <- qa_value(do.call(bargaining.logit, common),
                    "legacy bargaining AG")
    fit <- qa_value(do.call(calibrate, c(
        list(demand = "logit", conduct = "bargaining"),
        common[names(common) != "ownerPost"]
    )), "calibrate bargaining AG")
    actual <- qa_value(simulate(
        fit, ownerPost = common$ownerPost,
        bargpowerPost = common$bargpowerPost
    ), "simulate bargaining AG")
    expect_equal(fit@diagnostics$solver, "ag")
    expect_equal(actual@pricePost, old@pricePost, tolerance = 1e-8)
})

test_that("legacy sim routes bargaining models through specify and simulate", {
    legacy_sim <- getFromNamespace(".sim_legacy", "antitrust")
    cases <- list(
        list(
            demand = "Logit",
            supply = "bargaining",
            parameters = list(alpha = -1.2, meanval = c(.5, .3, .1)),
            prices = c(2, 2.2, 2.5),
            shares = c(.35, .25, .20),
            margins = c(.40, .35, .30)
        ),
        list(
            demand = "CES",
            supply = "bargaining2nd",
            parameters = list(gamma = 2, meanval = c(1, .8, .6), shareInside = .7),
            prices = c(1.5, 1.5, 1.5),
            shares = c(.50, .30, .20),
            margins = c(.30, .20, .15)
        )
    )
    for (case in cases) {
        common <- list(
            prices = case$prices,
            shares = case$shares,
            margins = case$margins,
            demand = case$demand,
            supply = case$supply,
            demand.param = case$parameters,
            ownerPre = c("A", "B", "C"),
            ownerPost = c("A", "A", "C"),
            bargpowerPre = rep(.5, 3),
            bargpowerPost = rep(.4, 3)
        )
        old <- qa_value(do.call(legacy_sim, common),
                        paste("legacy supplied", case$supply))
        actual <- qa_value(do.call(sim, common),
                           paste("compatibility supplied", case$supply))
        expect_equal(class(actual), class(old))
        expect_equal(actual@slopes, old@slopes, tolerance = 1e-8)
        expect_equal(actual@mcPost, old@mcPost, tolerance = 1e-8)
        expect_equal(actual@pricePost, old@pricePost, tolerance = 1e-8)
    }
})

test_that("Linear, LogLin, and AIDS calibration preserve their Bertrand equations", {
    prices <- c(2, 2.2, 2.5)
    quantities <- c(40, 35, 25)
    shares <- quantities / sum(quantities)
    margins <- c(.40, .35, .30)
    owner_pre <- c("A", "B", "C")
    owner_post <- c("A", "A", "C")
    diversions <- matrix(c(
        -1, .5, .5,
        .5, -1, .5,
        .5, .5, -1
    ), 3, byrow = TRUE)

    cases <- list(
        list(
            demand = "linear",
            constructor = linear,
            args = list(prices = prices, quantities = quantities,
                        margins = margins, diversions = diversions),
            calibration = list(shares = NULL, quantities = quantities,
                               diversions = diversions)
        ),
        list(
            demand = "loglin",
            constructor = loglinear,
            args = list(prices = prices, quantities = quantities,
                        margins = margins, diversions = diversions),
            calibration = list(shares = NULL, quantities = quantities,
                               diversions = diversions)
        ),
        list(
            demand = "aids",
            constructor = aids,
            args = list(shares = shares, prices = prices,
                        margins = margins, diversions = diversions,
                        mktElast = -1.2, priceStart = rep(.2, 3)),
            calibration = list(shares = shares, diversions = diversions,
                               mktElast = -1.2, priceStart = rep(.2, 3))
        )
    )

    for (case in cases) {
        old <- qa_value(do.call(case$constructor, c(
            case$args, list(ownerPre = owner_pre, ownerPost = owner_post)
        )), paste("legacy", case$demand, "calibration"))
        fit <- qa_value(do.call(calibrate, c(
            list(demand = case$demand, conduct = "bertrand",
                 prices = prices, margins = margins,
                 ownerPre = owner_pre),
            case$calibration
        )), paste("calibrate", case$demand))
        actual <- qa_value(simulate(fit, ownerPost = owner_post),
                           paste("simulate", case$demand))

        expect_equal(class(actual), class(old))
        expect_equal(fit@parameters,
                     getFromNamespace(".model_parameters", "antitrust")(old),
                     tolerance = 1e-8)
        expect_equal(actual@mcPre, old@mcPre, tolerance = 1e-8)
        expect_equal(actual@mcPost, old@mcPost, tolerance = 1e-8)
        expect_equal(actual@pricePre, old@pricePre, tolerance = 1e-8)
        expect_equal(actual@pricePost, old@pricePost, tolerance = 1e-8)
        expect_equal(calcShares(actual, TRUE), calcShares(old, TRUE), tolerance = 1e-8)
        expect_equal(calcShares(actual, FALSE), calcShares(old, FALSE), tolerance = 1e-8)
        expect_equal(calcMargins(actual, TRUE), calcMargins(old, TRUE), tolerance = 1e-8)
        expect_equal(calcMargins(actual, FALSE), calcMargins(old, FALSE), tolerance = 1e-8)
    }
})

test_that("supplied Linear, LogLin, and AIDS parameters use the shared simulation boundary", {
    prices <- c(2, 2.2, 2.5)
    quantities <- c(.40, .35, .25)
    margins <- c(.40, .35, .30)
    owner_pre <- c("A", "B", "C")
    owner_post <- c("A", "A", "C")
    legacy_sim <- getFromNamespace(".sim_legacy", "antitrust")

    linear_parameters <- list(
        slopes = {
            value <- matrix(.1, 3, 3)
            diag(value) <- -2
            value
        },
        intercepts = c(5, 5, 5)
    )
    loglin_parameters <- linear_parameters
    aids_fit <- qa_value(aids(
        shares = quantities, prices = prices, margins = margins,
        diversions = qa_fixture_diversions(), ownerPre = owner_pre,
        ownerPost = owner_post, mktElast = -1.2, priceStart = rep(.2, 3)
    ), "AIDS supplied-parameter source")
    cases <- list(
        list(demand = "Linear", parameters = linear_parameters,
             quantities = quantities),
        list(demand = "LogLin", parameters = loglin_parameters,
             quantities = quantities),
        list(demand = "AIDS", parameters = list(
            slopes = aids_fit@slopes, intercepts = aids_fit@intercepts,
            mktElast = aids_fit@mktElast
        ), quantities = NULL)
    )

    for (case in cases) {
        common <- list(
            prices = prices,
            shares = if (is.null(case$quantities)) quantities else case$quantities,
            margins = margins,
            demand = case$demand,
            supply = "bertrand",
            demand.param = case$parameters,
            ownerPre = owner_pre,
            ownerPost = owner_post,
            priceStart = if (case$demand == "AIDS") rep(.2, 3) else prices
        )
        old <- qa_value(do.call(legacy_sim, common),
                        paste("legacy supplied", case$demand))
        specify_args <- list(
            demand = tolower(case$demand), conduct = "bertrand",
            prices = prices, parameters = case$parameters,
            ownerPre = owner_pre, margins = margins
        )
        if (is.null(case$quantities)) {
            specify_args$shares <- quantities
            specify_args$priceStart <- rep(.2, 3)
        } else {
            specify_args$quantities <- case$quantities
        }
        fit <- qa_value(do.call(specify, specify_args),
                        paste("specify", case$demand))
        actual <- qa_value(simulate(fit, ownerPost = owner_post),
                           paste("simulate specified", case$demand))
        expect_equal(class(actual), class(old))
        expect_equal(actual@mcPost, old@mcPost, tolerance = 1e-8)
        expect_equal(actual@pricePost, old@pricePost, tolerance = 1e-8)
    }
})
