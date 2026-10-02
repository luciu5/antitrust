test_that("StructuralFit exposes the shared extension contract", {
    slots <- c("spec", "model", "parameters", "observed", "diagnostics")
    expect_true(methods::getClassDef("StructuralFit")@virtual)
    expect_equal(methods::slotNames("StructuralFit"), slots)
    expect_true(methods::extends("AntitrustFit", "StructuralFit"))

    fit <- calibrate(
        "logit", "bertrand", prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .2), margins = c(.4, .35, .3),
        ownerPre = c("A", "B", "C")
    )
    result <- simulate(fit, ownerPost = c("A", "A", "C"))
    expect_s4_class(fit, "StructuralFit")
    expect_s4_class(result, "Logit")
    expect_true(all(is.finite(result@pricePost)))
})

test_that("unsupported StructuralFit subclasses fail before running economics", {
    methods::setClass("StructuralFitNoFallbackDummy", contains = "StructuralFit")
    on.exit(methods::removeClass("StructuralFitNoFallbackDummy"), add = TRUE)
    dummy <- methods::new("StructuralFitNoFallbackDummy")
    expect_error(simulate(dummy, ownerPost = c("A", "B")),
                 "no simulate\\(\\) method is defined")
    expect_error(respecify(dummy, conduct = "cournot"),
                 "no respecify\\(\\) method is defined")
})

test_that("antitrust keeps ordinary stats::simulate dispatch", {
    lm_fit <- lm(mpg ~ wt, data = mtcars)
    set.seed(1)
    actual <- simulate(lm_fit, nsim = 2)
    set.seed(1)
    reference <- stats::simulate(lm_fit, nsim = 2)
    expect_equal(actual, reference)
})
