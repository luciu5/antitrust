# Regression tests for the S4 lifecycle dispatch contract: respecify() and
# simulate() are S4 generics owned by antitrust, with explicit AntitrustFit
# methods that reproduce the pre-conversion economic behavior exactly, and
# explicit failure (never silent economics) for unsupported StructuralFit
# subclasses. See R/ModelArchitecture.R and R/StructuralFitMethods.R.

s4_lifecycle_fit <- function(conduct = "bertrand") {
    calibrate(
        demand = "logit", conduct = conduct,
        prices = c(2, 2.2, 2.5), shares = c(.35, .25, .20),
        margins = c(.40, .35, .30), ownerPre = c("A", "B", "C"),
        insideSize = 100
    )
}


test_that("AntitrustFit satisfies the StructuralFit hierarchy", {
    fit <- s4_lifecycle_fit()
    expect_true(is(fit, "AntitrustFit"))
    expect_true(is(fit, "StructuralFit"))
})


## ---- respecify() dispatch and invariants -----------------------------------

test_that("respecify() is an S4 generic with an explicit AntitrustFit method", {
    expect_true(methods::isGeneric("respecify"))
    expect_true(methods::existsMethod("respecify", "AntitrustFit"))
})


test_that("respecify() Logit Bertrand to Cournot preserves observed baseline", {
    fit <- s4_lifecycle_fit("bertrand")
    target <- suppressWarnings(respecify(fit, conduct = "cournot"))

    expect_s4_class(target, "AntitrustFit")
    expect_equal(target@spec$conduct, "cournot")
    expect_equal(target@spec$demand, "logit")
    expect_equal(unname(target@model@pricePre), unname(fit@model@pricePre),
                 tolerance = 1e-10)
    expect_equal(target@model@labels, fit@model@labels)
    expect_equal(target@model@priceOutside, fit@model@priceOutside)
    expect_equal(target@model@insideSize, fit@model@insideSize)
    expect_equal(target@model@ownerPre, fit@model@ownerPre)
    expect_equal(target@diagnostics$route, "respecify")
    expect_equal(target@diagnostics$transition$to, "logit::cournot")
})


test_that("respecify() Logit Bertrand to MonCom preserves observed baseline", {
    fit <- s4_lifecycle_fit("bertrand")
    target <- suppressWarnings(respecify(fit, conduct = "moncom"))

    expect_s4_class(target, "AntitrustFit")
    expect_equal(target@spec$conduct, "moncom")
    expect_equal(unname(target@model@pricePre), unname(fit@model@pricePre),
                 tolerance = 1e-10)
    expect_equal(target@model@labels, fit@model@labels)
    expect_equal(target@model@ownerPre, fit@model@ownerPre)
    expect_equal(target@model@insideSize, fit@model@insideSize)
})


test_that("respecify() Logit Bertrand to second-score auction preserves invariants", {
    fit <- s4_lifecycle_fit("bertrand")
    target <- suppressWarnings(respecify(fit, conduct = "auction2nd"))

    expect_s4_class(target, "AntitrustFit")
    expect_equal(target@spec$conduct, "auction2nd")
    expect_equal(unname(target@model@pricePre), unname(fit@model@pricePre),
                 tolerance = 1e-10)
    expect_equal(target@model@labels, fit@model@labels)
    expect_equal(target@model@ownerPre, fit@model@ownerPre)
})


test_that("respecify() a CES demand transition preserves observed baseline", {
    fit <- s4_lifecycle_fit("bertrand")
    target <- suppressWarnings(respecify(fit, demand = "ces", gamma = 2))

    expect_s4_class(target, "AntitrustFit")
    expect_equal(target@spec$demand, "ces")
    expect_s4_class(target@model, "CES")
    expect_equal(unname(target@model@pricePre), unname(fit@model@pricePre),
                 tolerance = 1e-10)
    expect_equal(target@model@labels, fit@model@labels)
    expect_equal(target@model@ownerPre, fit@model@ownerPre)
    expect_equal(target@observed$insideSize, fit@observed$insideSize)
    expect_equal(target@diagnostics$transition$to, "ces::bertrand")
})


test_that("respecify() fails clearly on unsupported/unregistered transitions", {
    fit <- s4_lifecycle_fit("bertrand")
    expect_error(respecify(fit), "different registered")
    expect_error(respecify(fit, demand = "ces"), "requires explicit target primitive")
})


## ---- simulate() dispatch and invariants ------------------------------------

test_that("simulate() is an S4 generic with AntitrustFit and CounterfactualPath methods", {
    expect_true(methods::isGeneric("simulate"))
    expect_true(methods::existsMethod("simulate", "AntitrustFit"))
    expect_true(methods::existsMethod("simulate", "CounterfactualPath"))
})


test_that("simulate(AntitrustFit, ownerPost=) matches the legacy result", {
    fit <- s4_lifecycle_fit("bertrand")
    result <- simulate(fit, ownerPost = c("A", "A", "C"))

    expect_s4_class(result, "Logit")
    expect_length(result@pricePost, 3)
    expect_true(all(is.finite(result@pricePost)))
})


test_that("simulate() with an empty counterfactual is a baseline identity", {
    fit <- s4_lifecycle_fit("bertrand")
    baseline <- simulate(fit, counterfactual())

    expect_equal(unname(baseline@pricePost), unname(baseline@pricePre),
                 tolerance = 1e-6)
})


test_that("an ownership-only counterfactual does not reidentify structural costs", {
    fit <- s4_lifecycle_fit("bertrand")
    result <- simulate(fit, counterfactual(ownership = c("A", "A", "C")))

    expect_equal(unname(result@mcPost), unname(fit@model@mcPre), tolerance = 1e-10)
})


test_that("simulate() sequential counterfactual state invariants remain intact", {
    fit <- s4_lifecycle_fit("bertrand")
    cf <- add_step(
        counterfactual(ownership = c("A", "A", "C")),
        costs = c(-.05, 0, 0)
    )
    path <- simulate(fit, cf)

    expect_s4_class(path, "CounterfactualPath")
    expect_length(path@steps, 2L)
    expect_length(path@results, 2L)

    resumed <- simulate(path, counterfactual(costs = c(0, -.05, 0)))
    expect_s4_class(resumed, "CounterfactualPath")
    expect_length(resumed@steps, 3L)
    expect_equal(resumed@results[[1]]@pricePost, path@results[[1]]@pricePost,
                 tolerance = 1e-10)
    expect_equal(resumed@results[[2]]@pricePost, path@results[[2]]@pricePost,
                 tolerance = 1e-10)
})


## ---- Unsupported StructuralFit subclass safety -----------------------------

test_that("an unsupported StructuralFit subclass fails explicitly, never running economics", {
    setClass("DummyFit", contains = "StructuralFit")
    on.exit(removeClass("DummyFit"), add = TRUE)
    dummy <- methods::new("DummyFit")

    expect_true(is(dummy, "StructuralFit"))
    expect_false(is(dummy, "AntitrustFit"))

    expect_error(
        simulate(dummy, ownerPost = c("A", "B")),
        "no simulate\\(\\) method is defined"
    )
    expect_error(
        respecify(dummy, conduct = "cournot"),
        "no respecify\\(\\) method is defined"
    )
})


## ---- stats::simulate() is not hijacked -------------------------------------

test_that("attaching antitrust does not break ordinary stats::simulate() dispatch", {
    fit_lm <- lm(mpg ~ wt, data = mtcars)

    set.seed(1)
    named <- simulate(fit_lm, nsim = 2)
    set.seed(1)
    positional <- simulate(fit_lm, 2)
    set.seed(1)
    reference <- stats::simulate(fit_lm, nsim = 2)

    expect_equal(named, reference)
    expect_equal(positional, reference)
    expect_equal(dim(simulate(fit_lm, seed = 99)), c(32L, 1L))
})


## ---- Namespace / dispatch table --------------------------------------------

test_that("registered simulate() and respecify() signatures are exactly as intended", {
    simulate_signatures <- sort(methods::showMethods("simulate", printTo = FALSE))
    respecify_signatures <- sort(methods::showMethods("respecify", printTo = FALSE))

    expect_true(methods::existsMethod("simulate", "AntitrustFit"))
    expect_true(methods::existsMethod("simulate", "CounterfactualPath"))
    expect_true(methods::existsMethod("simulate", "StructuralFit"))
    expect_true(methods::existsMethod("simulate", "ANY"))

    expect_true(methods::existsMethod("respecify", "AntitrustFit"))
    expect_true(methods::existsMethod("respecify", "StructuralFit"))
})
