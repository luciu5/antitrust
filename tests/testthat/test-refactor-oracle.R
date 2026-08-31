## Phase 0: stable behavioral anchors for the calibration/simulation refactor.
## These values were recorded from master with R 4.6.1.  They intentionally
## exercise separate Bertrand and Cournot supply implementations that happen
## to use the same Logit demand family.

refactor_oracle_market <- function(conduct) {
    args <- list(
        prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .20),
        margins = c(.40, .35, .30),
        ownerPre = c("A", "B", "C"),
        ownerPost = c("A", "A", "C")
    )
    if (identical(conduct, "bertrand")) {
        do.call(logit, args)
    } else {
        do.call(logit.cournot, args)
    }
}

test_that("master Logit-Bertrand behavior is anchored before refactoring", {
    fit <- qa_value(refactor_oracle_market("bertrand"), "Phase 0 Logit Bertrand")

    expect_s4_class(fit, "Logit")
    expect_equal(unname(fit@slopes$alpha), -1.80501870371948, tolerance = 1e-10)
    expect_equal(unname(fit@slopes$meanval),
                 c(4.16965319537438, 4.19418469949706, 4.51254675929870),
                 tolerance = 1e-10)
    expect_equal(unname(fit@mcPre),
                 c(1.14767556962632, 1.46131882700948, 1.80748640032139),
                 tolerance = 1e-10)
    expect_equal(unname(fit@pricePre), c(2, 2.2, 2.5), tolerance = 1e-12)
    expect_equal(unname(fit@pricePost),
                 c(2.22868077358827, 2.54232403097142, 2.54263243829094),
                 tolerance = 1e-10)
    expect_equal(unname(calcShares(fit, FALSE)),
                 c(.308192520661549, .179311445518731, .246393436039781),
                 tolerance = 1e-10)
    expect_equal(unname(calcMargins(fit, FALSE)),
                 c(.485042639175656, .425203550425570, .289127923870463),
                 tolerance = 1e-10)
    expect_equal(unname(elast(fit, TRUE)), matrix(c(
        -2.34652431483532, 1.26351309260363, 1.26351309260363,
        .992760287045714, -2.97828086113714, .992760287045714,
        .902509351859740, .902509351859740, -3.61003740743896
    ), 3), tolerance = 1e-10)
    expect_equal(unname(diversion(fit, TRUE)), matrix(c(
        -1, .466666666666667, .4375,
        .384615384615385, -1, .3125,
        .307692307692308, .266666666666667, -1
    ), 3), tolerance = 1e-10)
    expect_equal(unname(upp(fit)), c(.106060606060606, .148076923076923),
                 tolerance = 1e-10)
    expect_equal(unname(cmcr(fit)), c(46.415797554791, 44.2303222555082),
                 tolerance = 1e-9)
    expect_equal(CV(fit), .158205880758702, tolerance = 1e-10)
})

test_that("master Logit-Cournot behavior is anchored before refactoring", {
    fit <- qa_value(refactor_oracle_market("cournot"), "Phase 0 Logit Cournot")

    expect_s4_class(fit, "LogitCournot")
    expect_equal(unname(fit@slopes$alpha), -2.89957783198386, tolerance = 1e-10)
    expect_equal(unname(fit@slopes$meanval),
                 c(6.23584219105484, 6.53324858827096, 7.19777427805901),
                 tolerance = 1e-10)
    expect_equal(unname(fit@mcPre),
                 c(1.12139938050303, 1.45275405005246, 1.82744800690037),
                 tolerance = 1e-10)
    expect_equal(unname(fit@pricePre), c(2, 2.2, 2.5), tolerance = 1e-12)
    expect_equal(unname(fit@pricePost),
                 c(2.08971000307171, 2.42106467262115, 2.5),
                 tolerance = 1e-10)
    expect_equal(unname(calcShares(fit, TRUE)),
                 c(.331784513198520, .250128936866548, .203696037938389),
                 tolerance = 1e-10)
    expect_equal(unname(calcShares(fit, FALSE)),
                 c(.317502162956362, .163547321016086, .252837992557633),
                 tolerance = 1e-10)
    expect_equal(unname(calcMargins(fit, FALSE)),
                 c(.463370812781458, .399952398443181, .269020797239854),
                 tolerance = 1e-10)
    expect_equal(unname(elast(fit, TRUE)), matrix(c(
        -3.87508562503575, 1.92407003893197, 1.92407003893197,
        1.59559030504706, -4.78348092531743, 1.59559030504706,
        1.47658129017274, 1.47658129017274, -5.77236328978691
    ), 3), tolerance = 1e-10)
    expect_equal(unname(upp(fit)), c(.163465347926664, .232562136626279),
                 tolerance = 1e-10)
    expect_equal(unname(cmcr(fit)), c(35.8808971264294, 36.7386914271462),
                 tolerance = 1e-9)
    expect_equal(CV(fit), .0745350109396821, tolerance = 1e-10)
})
