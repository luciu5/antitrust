translation_assertions <- function(source, target, target_demand,
                                    tolerance = 1e-7) {
    source_quantities <- calcQuantities(source@model, preMerger = TRUE)
    target_is_ces <- target_demand %in% c("ces", "ces_nests")
    source_shares <- if (target_demand %in% c("linear", "loglin")) {
        source_quantities / sum(source_quantities)
    } else {
        calcShares(source@model, preMerger = TRUE,
                   revenue = target_is_ces)
    }
    target_shares <- calcShares(target@model, preMerger = TRUE,
                                revenue = target_is_ces)
    expect_equal(unname(target@model@pricePre),
                 unname(source@model@pricePre), tolerance = tolerance)
    expect_equal(unname(target_shares), unname(source_shares),
                 tolerance = tolerance)
    expect_equal(unname(calcQuantities(target@model, preMerger = TRUE)),
                 unname(source_quantities), tolerance = 10 * tolerance)

    diagnostics <- target@diagnostics$translation
    expect_true(is.list(diagnostics))
    expect_true(is.finite(diagnostics$baseline_price_discrepancy))
    expect_lte(diagnostics$baseline_price_discrepancy, tolerance)
    expect_lte(diagnostics$baseline_share_discrepancy, tolerance)
    expect_lte(diagnostics$baseline_quantity_discrepancy, 10 * tolerance)
    expect_true(is.finite(diagnostics$elasticity_rmse))
    expect_true(is.finite(diagnostics$maximum_absolute_elasticity_difference))
    invisible(target)
}

flat_logit_fit <- function(conduct = "bertrand") {
    specify(
        "logit", conduct, prices = c(1.8, 2, 2.2),
        shares = c(.30, .20, .10), ownerPre = c("A", "B", "C"),
        parameters = list(alpha = -2,
                          meanval = c(.6, .4, .2)),
        insideSize = 100
    )
}

flat_ces_fit <- function(conduct = "bertrand") {
    specify(
        "ces", conduct, prices = c(1.8, 2, 2.2),
        shares = c(.30, .20, .10), ownerPre = c("A", "B", "C"),
        parameters = list(gamma = 1.5,
                          meanval = c(1, .8, .6)),
        insideSize = 100
    )
}

nested_logit_fit <- function() {
    specify(
        "logit_nests", "bertrand", prices = c(1.8, 2, 2.2, 2.5),
        shares = rep(.25, 4), ownerPre = c("A", "B", "C", "D"),
        parameters = list(alpha = -2,
                          meanval = c(0, -.1, .2, .05),
                          sigma = c(.7, .8)),
        nests = c("N1", "N1", "N2", "N2"), insideSize = 100
    )
}

nested_ces_fit <- function() {
    specify(
        "ces_nests", "bertrand", prices = c(1.8, 2, 2.2, 2.5),
        shares = rep(.25, 4), ownerPre = c("A", "B", "C", "D"),
        parameters = list(gamma = .8,
                          meanval = c(1, .6, .8, .5),
                          sigma = c(.7, .8)),
        nests = c("N1", "N1", "N2", "N2"), insideSize = 200
    )
}

test_that("flat Logit and CES translations require and use target primitives", {
    source <- flat_logit_fit()
    source_parameters <- source@parameters

    expect_error(respecify(source, demand = "ces"),
                 "requires explicit target primitive.*gamma")
    target <- respecify(source, demand = "ces", gamma = .8)
    translation_assertions(source, target, "ces")
    expect_s4_class(target@model, "CES")
    expect_equal(target@parameters$gamma, .8, tolerance = 0)
    expect_false(isTRUE(all.equal(source@model@mcPre, target@model@mcPre)))

    reverse <- respecify(target, demand = "logit", alpha = -1.5)
    translation_assertions(target, reverse, "logit")
    expect_s4_class(reverse@model, "Logit")
    expect_equal(reverse@parameters$alpha, -1.5, tolerance = 0)
    expect_equal(source@parameters, source_parameters, tolerance = 0)
})

test_that("flat Logit/CES translations retain Cournot conduct", {
    source <- flat_logit_fit("cournot")
    target <- respecify(source, demand = "ces", gamma = .8)
    translation_assertions(source, target, "ces")
    expect_s4_class(target@model, "CESCournot")
    expect_equal(target@spec$conduct, "cournot")
})

test_that("nested-to-flat transitions are structural restrictions", {
    logit_source <- nested_logit_fit()
    logit_target <- respecify(logit_source, demand = "logit")
    translation_assertions(logit_source, logit_target, "logit")
    expect_s4_class(logit_target@model, "Logit")
    expect_equal(logit_target@parameters$alpha,
                 logit_source@parameters$alpha, tolerance = 0)
    expect_equal(logit_target@diagnostics$translation$transition_kind,
                 "structural-restriction")

    ces_source <- nested_ces_fit()
    ces_target <- respecify(ces_source, demand = "ces")
    translation_assertions(ces_source, ces_target, "ces")
    expect_s4_class(ces_target@model, "CES")
    expect_equal(ces_target@parameters$gamma,
                 ces_source@parameters$gamma, tolerance = 0)
})

test_that("nested Logit/CES translations preserve nests and derive curvature", {
    logit_source <- nested_logit_fit()
    ces_target <- respecify(logit_source, demand = "ces_nests", gamma = .8)
    translation_assertions(logit_source, ces_target, "ces_nests")
    expect_s4_class(ces_target@model, "CESNests")
    expect_equal(as.character(ces_target@model@nests),
                 as.character(logit_source@model@nests))
    expect_equal(unname(ces_target@parameters$sigma),
                 unname(1 + (.8 - 1) / logit_source@parameters$sigma),
                 tolerance = 1e-12)

    logit_target <- respecify(ces_target, demand = "logit_nests", alpha = -1.5)
    translation_assertions(ces_target, logit_target, "logit_nests")
    expect_s4_class(logit_target@model, "LogitNests")
    expect_equal(as.character(logit_target@model@nests),
                 as.character(ces_target@model@nests))
    expect_equal(unname(logit_target@parameters$sigma),
                 unname((.2) / (1 - ces_target@parameters$sigma)),
                 tolerance = 1e-12)
})

test_that("flat-to-nested transitions require all new structural primitives", {
    prices <- c(1.8, 2, 2.2, 2.5)
    owner <- c("A", "B", "C", "D")
    nests <- c("N1", "N1", "N2", "N2")
    logit_source <- specify(
        "logit", "bertrand", prices = prices, shares = rep(.25, 4),
        ownerPre = owner,
        parameters = list(alpha = -2, meanval = c(.6, .4, .2, .1)),
        insideSize = 100
    )
    expect_error(respecify(logit_source, demand = "logit_nests"),
                 "requires explicit target primitive.*nests.*sigma")
    target <- respecify(logit_source, demand = "logit_nests",
                        nests = nests, sigma = .7)
    translation_assertions(logit_source, target, "logit_nests")
    expect_equal(as.character(target@model@nests), nests)

    ces_source <- specify(
        "ces", "bertrand", prices = prices, shares = rep(.25, 4),
        ownerPre = owner,
        parameters = list(gamma = 1.5, meanval = c(1, .8, .6, .4)),
        insideSize = 100
    )
    ces_target <- respecify(ces_source, demand = "ces_nests", gamma = .8,
                            nests = nests, sigma = .7)
    translation_assertions(ces_source, ces_target, "ces_nests")
    expect_equal(as.character(ces_target@model@nests), nests)
})

test_that("first-order Linear and LogLin translations are canonical", {
    source <- flat_logit_fit()
    linear <- respecify(source, demand = "linear")
    loglin <- respecify(source, demand = "loglin")
    translation_assertions(source, linear, "linear")
    translation_assertions(source, loglin, "loglin")
    expect_lte(linear@diagnostics$translation$jacobian_rmse, 1e-8)
    expect_lte(loglin@diagnostics$translation$elasticity_rmse, 1e-8)

    ces_source <- flat_ces_fit()
    ces_linear <- respecify(ces_source, demand = "linear")
    translation_assertions(ces_source, ces_linear, "linear")
    expect_lte(ces_linear@diagnostics$translation$jacobian_rmse, 1e-8)

    round_trip <- respecify(loglin, demand = "linear")
    translation_assertions(loglin, round_trip, "linear")
})

test_that("AIDS can be translated to tangent Linear and LogLin demand", {
    prices <- c(1.8, 2, 2.2)
    shares <- c(.4, .3, .3)
    fit <- calibrate(
        "aids", "bertrand", prices = prices, shares = shares,
        margins = c(.2, .25, .3), ownerPre = c("A", "B", "C"),
        mktElast = -1, insideSize = 100
    )
    linear <- respecify(fit, demand = "linear")
    loglin <- respecify(fit, demand = "loglin")
    translation_assertions(fit, linear, "linear")
    translation_assertions(fit, loglin, "loglin")
    expect_lte(linear@diagnostics$translation$jacobian_rmse, 1e-8)
    expect_lte(loglin@diagnostics$translation$elasticity_rmse, 1e-8)
})

test_that("respecification is deterministic, explicit, and non-mutating", {
    source <- flat_logit_fit()
    source_parameters <- source@parameters
    source_model <- source@model

    expect_error(respecify(source, demand = "aids"),
                 "not supported.*update")
    expect_error(respecify(source, conduct = "cournot", alpha = -1),
                 "does not accept transition-specific")
    expect_error(respecify(source, demand = "logit_nests",
                           nests = c("N1", "N1", "N2")),
                 "requires explicit target primitive.*sigma")

    target <- respecify(source, demand = "ces", gamma = .8)
    expect_equal(source@parameters, source_parameters, tolerance = 0)
    expect_equal(source@model, source_model, tolerance = 0)
    expect_equal(target@diagnostics$translation$required_arguments, "gamma")
    expect_false("objective" %in% names(target@diagnostics$translation))
})
