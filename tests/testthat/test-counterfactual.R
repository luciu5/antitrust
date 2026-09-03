test_that("Counterfactual translates legacy antitrust scenario arguments", {
    fit <- calibrate(
        "logit", "bertrand", prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .2), margins = c(.4, .35, .3),
        ownerPre = c("A", "B", "C"), insideSize = 100
    )
    owner <- c("A", "A", "C")
    costs <- c(0, -.1, 0)
    old <- simulate(fit, ownerPost = owner, mcDelta = costs)
    cf <- counterfactual(ownership = owner, costs = costs)
    new <- simulate(fit, cf)

    expect_equal(new@pricePost, old@pricePost, tolerance = 1e-8)
    expect_equal(new@mcPost, old@mcPost, tolerance = 1e-8)
    expect_equal(attr(new, "counterfactual")$fields, cf)
    expect_equal(fit@model@pricePost, fit@model@pricePre, tolerance = 1e-12)
})

test_that("an empty Counterfactual preserves the baseline equilibrium", {
    fit <- calibrate(
        "logit", "bertrand", prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .2), margins = c(.4, .35, .3),
        ownerPre = c("A", "B", "C"), insideSize = 100
    )
    result <- simulate(fit, counterfactual())
    expect_equal(result@pricePost, fit@model@pricePre, tolerance = 1e-8)
})

test_that("Counterfactual supports exit and rejects model changes", {
    fit <- specify(
        "logit", "bertrand", prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .2), ownerPre = c("A", "B", "C"),
        parameters = list(alpha = -1, meanval = c(.4, .2, .1)),
        insideSize = 100
    )
    source <- fit@model@pricePre
    result <- simulate(fit, counterfactual(exit = 3))
    logical_result <- simulate(fit, counterfactual(
        exit = c(TRUE, TRUE, FALSE)
    ))
    expect_equal(fit@model@pricePre, source, tolerance = 0)
    expect_error(counterfactual(demand = "ces"), "update|respecify")
    expect_error(simulate(fit, counterfactual(quota = 1)), "does not support")
    expect_error(simulate(fit, counterfactual(exit = 3), subset = rep(TRUE, 3)),
                 "cannot combine")
    expect_true(inherits(result, "Logit"))
    expect_true(inherits(logical_result, "Logit"))
})

test_that("Counterfactuals compose without mutation", {
    first <- counterfactual(ownership = c("A", "A"))
    second <- counterfactual(costs = c(0, -.1))
    combined <- combine_counterfactuals(first, second)
    expect_equal(names(combined), c("ownership", "costs"))
    expect_error(combine_counterfactuals(
        first, counterfactual(ownership = c("A", "B"))
    ), "conflicting")
})
