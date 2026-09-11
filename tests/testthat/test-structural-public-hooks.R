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
