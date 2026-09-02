translation_assertions <- function(source, target, target_demand) {
    source_quantities <- calcQuantities(source@model, preMerger = TRUE)
    source_shares <- calcShares(
        source@model, preMerger = TRUE,
        revenue = identical(target_demand, "ces") ||
            identical(target_demand, "ces_nests")
    )
    target_shares <- calcShares(
        target@model, preMerger = TRUE,
        revenue = identical(target_demand, "ces") ||
            identical(target_demand, "ces_nests")
    )

    expect_equal(unname(target@model@pricePre),
                 unname(source@model@pricePre), tolerance = 1e-8)
    expect_equal(unname(target_shares), unname(source_shares), tolerance = 1e-8)
    expect_equal(unname(calcQuantities(target@model, preMerger = TRUE)),
                 unname(source_quantities), tolerance = 1e-7)
    expect_true(is.finite(target@diagnostics$local_translation$elasticity_rmse))
    expect_true(is.finite(target@diagnostics$local_translation$
                          maximum_absolute_elasticity_difference))
    expect_lte(target@diagnostics$local_translation$baseline_share_discrepancy,
               1e-8)
    expect_lte(target@diagnostics$local_translation$baseline_quantity_discrepancy,
               1e-6)
    expect_equal(target@diagnostics$local_translation$transition_type,
                 "local-demand-translation")
    invisible(target)
}

test_that("flat Logit and CES respecification preserves the fitted baseline", {
    prices <- c(1.8, 2, 2.2, 2.5)
    owner <- c("A", "B", "C", "D")
    source <- specify(
        "logit", "bertrand", prices = prices,
        shares = c(.30, .20, .10, .05), ownerPre = owner,
        parameters = list(alpha = -1,
                          meanval = c(.4, .2, .1, .05)),
        insideSize = 100
    )
    source_parameters <- source@parameters
    target <- respecify(source, demand = "ces")
    translation_assertions(source, target, "ces")
    expect_s4_class(target@model, "CES")
    expect_false(isTRUE(all.equal(source@model@mcPre, target@model@mcPre)))
    expect_equal(source@parameters, source_parameters, tolerance = 0)
    expect_equal(target@diagnostics$local_translation$source_demand, "logit")
    expect_equal(target@diagnostics$local_translation$target_demand, "ces")

    reverse <- respecify(target, demand = "logit")
    translation_assertions(target, reverse, "logit")
    expect_s4_class(reverse@model, "Logit")
})

test_that("flat Logit/CES translation retains Cournot conduct", {
    prices <- c(1.8, 2, 2.2)
    owner <- c("A", "B", "C")
    source <- specify(
        "logit", "cournot", prices = prices,
        shares = c(.30, .20, .10), ownerPre = owner,
        parameters = list(alpha = -1, meanval = c(.4, .2, .1)),
        insideSize = 100
    )
    target <- respecify(source, demand = "ces")
    translation_assertions(source, target, "ces")
    expect_s4_class(target@model, "CESCournot")
    expect_equal(target@spec$conduct, "cournot")
})

test_that("nested Logit and nested CES respecification preserves nests and shares", {
    prices <- c(1.8, 2, 2.2, 2.5)
    owner <- c("A", "B", "C", "D")
    nests <- c("N1", "N1", "N2", "N2")
    source <- specify(
        "logit_nests", "bertrand", prices = prices,
        shares = c(.30, .20, .10, .05), ownerPre = owner,
        parameters = list(alpha = -1,
                          meanval = c(.4, .2, .1, .05), sigma = .7),
        nests = nests, insideSize = 100
    )
    source_parameters <- source@parameters
    target <- respecify(source, demand = "ces_nests")
    translation_assertions(source, target, "ces_nests")
    expect_s4_class(target@model, "CESNests")
    expect_equal(as.character(target@model@nests),
                 as.character(source@model@nests))
    expect_false(isTRUE(all.equal(source@model@mcPre, target@model@mcPre)))
    expect_equal(source@parameters, source_parameters, tolerance = 0)

    reverse <- respecify(target, demand = "logit_nests")
    translation_assertions(target, reverse, "logit_nests")
    expect_s4_class(reverse@model, "LogitNests")
    expect_equal(as.character(reverse@model@nests), nests)
})

test_that("nested CES to nested Logit translates supplied target parameters", {
    prices <- c(1.8, 2, 2.2, 2.5)
    owner <- c("A", "B", "C", "D")
    nests <- c("N1", "N1", "N2", "N2")
    source <- specify(
        "ces_nests", "bertrand", prices = prices,
        shares = rep(.25, 4), ownerPre = owner,
        parameters = list(gamma = .8,
                          meanval = c(1, .6, .8, .5),
                          sigma = c(.7, .8)),
        nests = nests, insideSize = 200
    )
    target <- respecify(source, demand = "logit_nests")
    translation_assertions(source, target, "logit_nests")
    expect_s4_class(target@model, "LogitNests")
    expect_true(all(target@model@slopes$sigma > 0 &
                    target@model@slopes$sigma <= 1))
})
