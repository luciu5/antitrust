test_that("registered core counterfactual capabilities have behavior-backed anchors", {
    capabilities <- getFromNamespace(
        ".model_counterfactual_capabilities", "antitrust"
    )

    logit_spec <- model_spec("logit", "bertrand")
    logit_capabilities <- capabilities(logit_spec)
    expect_true(all(logit_capabilities[c(
        "ownership", "costs", "exit", "quality", "entry"
    )]))
    expect_false(any(logit_capabilities[c(
        "capacity", "bargaining", "leader", "products", "tariff", "quota"
    )]))

    fit <- calibrate(
        "logit", "bertrand", prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .2), margins = c(.4, .35, .3),
        ownerPre = c("A", "B", "C"), insideSize = 100
    )
    quality <- simulate(fit, counterfactual(quality = c(Prod1 = .1)))
    expect_true(all(is.finite(quality@pricePost)))
    expect_equal(quality@mcPost, fit@model@mcPre, tolerance = 1e-8)

    ## A demand transition is a structural construction. Changing the
    ## source's observed margins after calibration cannot alter its target
    ## primitives or invoke a second target calibration.
    translated <- respecify(fit, demand = "ces", gamma = 2)
    altered <- fit
    altered@observed$margins <- c(.9, .8, .7)
    translated_altered <- respecify(altered, demand = "ces", gamma = 2)
    expect_equal(translated@parameters, translated_altered@parameters,
                 tolerance = 0)
    expect_equal(translated@model@pricePre,
                 translated_altered@model@pricePre, tolerance = 0)
    expect_null(translated@diagnostics$calibration_args)

    entered <- simulate(fit, counterfactual(entry = entrant(
        "E1", meanval = .1, cost = 1.5, priceStart = 2.2
    )))
    expect_true(all(is.finite(entered@pricePost)))
    expect_equal(length(entered@labels), length(fit@model@labels) + 1L)

    cap_spec <- model_spec("logit_cap", "bertrand")
    cap_capabilities <- capabilities(cap_spec)
    expect_true(isTRUE(cap_capabilities[["capacity"]]))
    expect_true(isTRUE(cap_capabilities[["quality"]]))
    expect_true(isTRUE(cap_capabilities[["entry"]]))

    cap_fit <- suppressWarnings(specify(
        "logit_cap", "bertrand", prices = c(2, 2.2, 2.5),
        parameters = list(alpha = -1, meanval = c(.4, .2, .1),
                          mktSize = 200),
        ownerPre = c("A", "B", "C"), insideSize = 200,
        capacities = c(50, 50, 50), margins = c(.4, .35, .3)
    ))
    cap_quality <- simulate(cap_fit,
                            counterfactual(quality = c(Prod1 = .1)))
    expect_equal(cap_quality@slopes$meanval[[1]], .4 * 1.1,
                 tolerance = 1e-12)
    cap_entry <- simulate(cap_fit, counterfactual(entry = entrant(
        "E1", meanval = .2, cost = 1, priceStart = 2, capacity = 50
    )))
    expect_true("E1" %in% cap_entry@labels)
    expect_equal(cap_entry@capacitiesPre[[4]], 50, tolerance = 0)

    bargaining_spec <- model_spec("logit", "bargaining")
    bargaining_capabilities <- capabilities(bargaining_spec)
    expect_true(isTRUE(bargaining_capabilities[["bargaining"]]))
    expect_true(isTRUE(bargaining_capabilities[["quality"]]))

    bargaining <- respecify(fit, conduct = "bargaining",
                            bargpowerPre = rep(.5, length(fit@model@labels)))
    bargaining_quality <- simulate(
        bargaining, counterfactual(quality = c(Prod1 = .1))
    )
    expect_equal(bargaining_quality@slopes$meanval[[1]],
                 fit@model@slopes$meanval[[1]] * 1.1, tolerance = 1e-12)

    stackelberg_spec <- model_spec("linear", "stackelberg")
    stackelberg_capabilities <- capabilities(stackelberg_spec)
    expect_true(all(stackelberg_capabilities[c(
        "capacity", "leader", "products"
    )]))
})
