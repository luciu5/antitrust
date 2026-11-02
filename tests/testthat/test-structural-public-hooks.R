test_that("persistent structural cost state is initialized without changing the fit", {
    fit <- calibrate(
        "logit", "bertrand", prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .2), margins = c(.4, .35, .3),
        ownerPre = c("A", "B", "C"), insideSize = 100
    )
    initialized <- initialize_cost_state(fit@model)

    expect_s4_class(initialized, "Logit")
    expect_equal(initialized@pricePre, fit@model@pricePre, tolerance = 0)
    state <- attr(initialized, "antitrust_cost_state", exact = TRUE)
    expect_type(state, "list")
    expect_equal(state$base, unname(fit@model@mcPre), tolerance = 1e-12)
    expect_identical(state$mode, "multiplicative")
})

test_that("supplied flat demand can bind an observed baseline without solving", {
    prices <- c(2, 2.2, 2.5)
    owner <- c("A", "A", "B")
    parameters <- list(alpha = -1.2, meanval = c(.3, .1, -.2))

    fit <- suppressWarnings(specify(
        "logit", "bertrand", prices = prices, parameters = parameters,
        ownerPre = owner, insideSize = 100, baseline = "observed"
    ))

    expect_s4_class(fit, "AntitrustFit")
    expect_identical(fit@diagnostics$baseline, "observed")
    expect_identical(fit@diagnostics$baseline_equilibrium, "supplied_observed")
    expect_identical(fit@diagnostics$specification_args$baseline, "observed")
    expect_equal(unname(fit@model@pricePre), prices, tolerance = 0)
    expect_equal(unname(fit@model@pricePost), prices, tolerance = 0)
    expect_true(all(is.finite(fit@model@mcPre)))
    expect_true(all(is.finite(fit@model@mcPost)))

    expect_error(
        specify("linear", "bertrand", prices = prices,
                parameters = list(slopes = diag(-1, 3), intercepts = prices),
                ownerPre = owner, baseline = "observed"),
        "supported only for standard supplied-parameter Logit/CES"
    )

    bad_owner <- fit
    bad_owner@model@ownerPost[1, 1] <- 0
    expect_error(
        initialize_baseline_state(bad_owner, baseline = "observed"),
        "neutral pre/post ownership"
    )

    bad_cost <- fit
    bad_cost@model@mcDelta[1] <- .1
    expect_error(
        initialize_baseline_state(bad_cost, baseline = "observed"),
        "neutral cost-shock"
    )

    bad_subset <- fit
    bad_subset@model@subset[1] <- FALSE
    expect_error(
        initialize_baseline_state(bad_subset, baseline = "observed"),
        "all products active"
    )
})

test_that("cost shocks preserve multiplicative and second-score additive conventions", {
    fit <- calibrate(
        "logit", "bertrand", prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .2), margins = c(.4, .35, .3),
        ownerPre = c("A", "B", "C"), insideSize = 100
    )
    expect_equal(
        compound_cost_shocks(fit@model, c(.10, 0, -.05), c(.20, .10, 0)),
        c(.32, .10, -.05), tolerance = 1e-12
    )

    auction <- auction2nd.logit(
        prices = c(2, 2.2, 2.5), shares = c(.35, .25, .2),
        margins = c(.4, .35, .3), ownerPre = c("A", "B", "C"),
        ownerPost = c("A", "B", "C")
    )
    expect_equal(
        compound_cost_shocks(auction, c(.10, 0, -.05), c(.20, .10, 0)),
        c(.30, .10, -.05), tolerance = 1e-12
    )
})

test_that("counterfactual validation dispatches through StructuralFit", {
    fit <- calibrate(
        "logit", "bertrand", prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .2), margins = c(.4, .35, .3),
        ownerPre = c("A", "B", "C"), insideSize = 100
    )
    cf <- counterfactual(ownership = c("A", "A", "C"))
    expect_identical(validate_counterfactual(fit, cf), cf)

    if (!methods::isClass("StructuralHookTestFit")) {
        methods::setClass("StructuralHookTestFit", contains = "StructuralFit")
    }
    unknown <- methods::new("StructuralHookTestFit")
    expect_error(
        validate_counterfactual(unknown, cf),
        "no validate_counterfactual\\(\\) method"
    )
})

test_that("CounterfactualPath resume uses the public validation generic", {
    fit <- calibrate(
        "logit", "bertrand", prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .2), margins = c(.4, .35, .3),
        ownerPre = c("A", "B", "C"), insideSize = 100
    )
    cf <- counterfactual(ownership = c("A", "A", "C"))
    cf <- add_step(cf, costs = c(-.01, 0, 0))
    path <- simulate(fit, cf)
    resumed <- simulate(path, counterfactual(ownership = c("A", "A", "A")))

    expect_s4_class(path, "CounterfactualPath")
    expect_s4_class(resumed, "CounterfactualPath")
    expect_length(resumed@results, 3L)
})

test_that("the old antitrust vertical constructor is an explicit migration error", {
    expect_error(vertical.barg(), "vertical::vertical.barg")
})
