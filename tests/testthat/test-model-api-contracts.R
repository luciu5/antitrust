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


test_that("update recalibrates while respecify preserves fitted demand", {
    market <- api_logit_observed_market()
    fit <- qa_value(do.call(calibrate, c(
        list(demand = "logit", conduct = "bertrand"),
        market[names(market) != "ownerPost"]
    )), "source fit for update/respecify contract")
    updated <- qa_value(update(fit, conduct = "cournot"),
                        "updated Cournot calibration")
    direct <- qa_value(do.call(calibrate, c(
        list(demand = "logit", conduct = "cournot"),
        market[names(market) != "ownerPost"]
    )), "direct Cournot calibration")
    respecified <- qa_value(respecify(fit, conduct = "cournot"),
                             "respecified Cournot fit")

    expect_equal(updated@parameters, direct@parameters, tolerance = 1e-9)
    expect_equal(respecified@parameters, fit@parameters, tolerance = 0)
    expect_false(isTRUE(all.equal(updated@parameters, respecified@parameters)))
    expect_equal(fit@spec$id, "logit::bertrand")
    expect_error(update(respecified), "created by calibrate.*respecify")
})
