test_that("StructuralFit is the shared virtual fit contract", {
    expected_slots <- c("spec", "model", "parameters", "observed",
                        "diagnostics")
    structural <- methods::getClassDef("StructuralFit")

    expect_true(methods::isClass("StructuralFit"))
    expect_true(structural@virtual)
    expect_equal(names(structural@slots), expected_slots)
    expect_true(methods::extends("AntitrustFit", "StructuralFit"))
    expect_equal(methods::slotNames("AntitrustFit"), expected_slots)
})


test_that("AntitrustFit construction retains the shared slot contract", {
    fit <- methods::new(
        "AntitrustFit",
        spec = list(id = "test::fit"),
        model = list(state = TRUE),
        parameters = list(alpha = -1.2),
        observed = list(prices = c(1, 2)),
        diagnostics = list(status = "constructed")
    )

    expect_s4_class(fit, "AntitrustFit")
    expect_s4_class(fit, "StructuralFit")
    expect_equal(fit@spec$id, "test::fit")
    expect_equal(fit@parameters$alpha, -1.2)
})


test_that("named simulate calls preserve AntitrustFit dispatch", {
    fit <- calibrate(
        demand = "logit", conduct = "bertrand",
        prices = c(2, 2.2, 2.5), shares = c(.35, .25, .20),
        margins = c(.40, .35, .30), ownerPre = c("A", "B", "C")
    )
    result <- simulate(fit, ownerPost = c("A", "A", "C"))

    expect_s4_class(fit, "StructuralFit")
    expect_s4_class(result, "Logit")
    expect_length(result@pricePost, 3)
})


test_that("AntitrustFit survives an RDS round trip", {
    fit <- methods::new(
        "AntitrustFit",
        spec = list(id = "test::rds"),
        model = list(state = TRUE),
        parameters = list(alpha = -1.2),
        observed = list(prices = c(1, 2)),
        diagnostics = list(status = "constructed")
    )
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    saveRDS(fit, path)
    restored <- readRDS(path)

    expect_s4_class(restored, "AntitrustFit")
    expect_s4_class(restored, "StructuralFit")
    expect_equal(restored, fit)
})


test_that("StructuralFit has no broad simulation fallback", {
    expect_true(methods::isGeneric("simulate"))

    setClass("StructuralFitNoFallbackDummy", contains = "StructuralFit")
    on.exit(removeClass("StructuralFitNoFallbackDummy"), add = TRUE)
    dummy <- methods::new("StructuralFitNoFallbackDummy")

    expect_error(
        simulate(dummy, ownerPost = c("A", "A")),
        "no simulate\\(\\) method is defined"
    )
    expect_error(
        respecify(dummy),
        "no respecify\\(\\) method is defined"
    )
    expect_error(
        methods::new("StructuralFit"),
        "virtual class"
    )
})
