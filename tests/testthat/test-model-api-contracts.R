api_logit_observed_market <- function() {
    list(
        prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .20),
        margins = c(.40, .35, .30),
        ownerPre = c("A", "B", "C"),
        ownerPost = c("A", "A", "C"),
        insideSize = 100
    )
}


test_that("calibrate and simulate retain the canonical Logit public contract", {
    market <- api_logit_observed_market()
    legacy <- qa_value(do.call(logit, market), "legacy Logit contract")
    fit <- qa_value(do.call(calibrate, c(
        list(demand = "logit", conduct = "bertrand"),
        market[names(market) != "ownerPost"]
    )), "calibrate Logit contract")
    result <- qa_value(simulate(fit, ownerPost = market$ownerPost),
                       "simulate Logit contract")

    expect_s4_class(fit, "AntitrustFit")
    expect_equal(fit@spec$id, "logit::bertrand")
    expect_equal(fit@diagnostics$status, "completed")
    expect_s4_class(result, "Logit")
    expect_equal(result@pricePost, legacy@pricePost, tolerance = 1e-9)
    expect_equal(result@mcPost, legacy@mcPost, tolerance = 1e-9)
})


test_that("specify and legacy sim share the supplied-parameter boundary", {
    market <- api_logit_observed_market()
    parameters <- list(alpha = -1.2, meanval = c(.5, .3, .1))
    fit <- qa_value(do.call(specify, c(
        list(demand = "logit", conduct = "bertrand", parameters = parameters),
        market[names(market) %in% c("prices", "shares", "margins", "ownerPre", "insideSize")]
    )), "specify Logit contract")
    result <- qa_value(simulate(fit, ownerPost = market$ownerPost),
                       "simulate specified Logit contract")
    legacy <- qa_value(sim(
        prices = market$prices, shares = market$shares,
        margins = market$margins, ownerPre = market$ownerPre,
        ownerPost = market$ownerPost, insideSize = market$insideSize,
        supply = "bertrand", demand = "Logit", demand.param = parameters
    ), "legacy supplied Logit contract")

    expect_equal(fit@diagnostics$source, "specified")
    expect_equal(result@pricePost, legacy@pricePost, tolerance = 1e-9)
    expect_equal(result@mcPost, legacy@mcPost, tolerance = 1e-9)
})


test_that("calibration and specification reject post-state fields", {
    market <- api_logit_observed_market()
    expect_error(
        do.call(calibrate, c(
            list(demand = "logit", conduct = "bertrand"),
            market[names(market) != "ownerPost"],
            list(mcfunPost = list(function(q) q))
        )),
        "simulation scenario"
    )
    expect_error(
        do.call(specify, c(
            list(demand = "logit", conduct = "bertrand",
                 parameters = list(alpha = -1.2, meanval = c(.5, .3, .1))),
            market[names(market) %in% c("prices", "shares", "ownerPre", "insideSize")],
            list(mcfunPost = list(function(q) q))
        )),
        "simulation scenario"
    )
})


test_that("update recalibrates the same model while respecify changes conduct", {
    market <- api_logit_observed_market()
    fit <- qa_value(do.call(calibrate, c(
        list(demand = "logit", conduct = "bertrand"),
        market[names(market) != "ownerPost"]
    )), "source fit for update/respecify contract")
    revised_margins <- c(.30, .28, .24)
    updated <- qa_value(update(fit, margins = revised_margins),
                        "updated same-model calibration")
    direct <- qa_value(calibrate(
        demand = "logit", conduct = "bertrand", prices = market$prices,
        shares = market$shares, margins = revised_margins,
        ownerPre = market$ownerPre, insideSize = market$insideSize
    ), "direct same-model calibration")
    respecified <- qa_value(respecify(fit, conduct = "cournot"),
                             "respecified Cournot fit")

    expect_equal(updated@parameters, direct@parameters, tolerance = 1e-9)
    expect_equal(respecified@parameters, fit@parameters, tolerance = 0)
    expect_equal(updated@spec$id, fit@spec$id)
    expect_error(update(fit, conduct = "cournot"), "same demand, conduct, and variant")
    expect_error(update(respecified), "created by calibrate.*respecify")
})

test_that("respecify retains source baseline metadata", {
    fit <- calibrate(
        demand = "logit", conduct = "bertrand",
        prices = c(2, 2.2, 2.5), shares = c(.35, .25, .20),
        margins = c(.40, .35, .30), ownerPre = c("A", "B", "C"),
        labels = c("X", "Y", "Z"), priceOutside = .7,
        priceStart = c(2.1, 2.3, 2.6), insideSize = 100
    )
    target <- respecify(fit, demand = "ces", gamma = 2)

    for (name in c("labels", "priceOutside", "priceStart", "output",
                   "insideSize")) {
        expect_equal(target@observed[[name]], fit@observed[[name]],
                     info = paste("metadata field", name))
    }
})
