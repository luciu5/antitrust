test_that("all documented sim demand/supply combinations have formal coverage", {
    p <- c(2, 2.2, 2.5)
    op <- c("A", "B", "C")
    oo <- c("A", "A", "C")
    common <- list(prices = p, ownerPre = op, ownerPost = oo)
    results <- list()

    results$bertrand_linear <- qa_value(sim(
        c(2, 2.2, 2.5), supply = "bertrand", demand = "Linear",
        demand.param = qa_fixture_linear_parameters(),
        ownerPre = op, ownerPost = oo
    ), "sim bertrand Linear")
    results$bertrand_loglin <- qa_value(sim(
        p, supply = "bertrand", demand = "LogLin",
        demand.param = qa_fixture_linear_parameters(),
        ownerPre = op, ownerPost = oo
    ), "sim bertrand LogLin")
    results$bertrand_aids <- qa_value(sim(
        p, supply = "bertrand", demand = "AIDS",
        demand.param = qa_fixture_aids_parameters(),
        ownerPre = op, ownerPost = oo
    ), "sim bertrand AIDS")

    for (supply in c("bertrand", "cournot", "auction2nd", "bargaining", "bargaining2nd")) {
        results[[paste("logit", supply, sep = "_")]] <- qa_value(sim(
            p, supply = supply, demand = "Logit",
            demand.param = list(alpha = -1, meanval = c(.4, .2, .1)),
            ownerPre = op, ownerPost = oo, insideSize = 100
        ), paste("sim Logit", supply))
        results[[paste("ces", supply, sep = "_")]] <- qa_value(sim(
            p, supply = supply, demand = "CES",
            demand.param = list(gamma = 2, meanval = c(1, .8, .6),
                                shareInside = .7),
            ownerPre = op, ownerPost = oo, insideSize = 100
        ), paste("sim CES", supply))
    }

    p4 <- c(1.8, 2, 2.2, 2.5)
    op4 <- c("A", "B", "C", "D")
    oo4 <- c("A", "A", "C", "D")
    nests4 <- c("N1", "N1", "N2", "N2")
    results$bertrand_logit_nests <- qa_value(sim(
        p4, supply = "bertrand", demand = "LogitNests",
        demand.param = list(alpha = -1, meanval = c(.4, .2, .1, .05), sigma = .7),
        nests = nests4, ownerPre = op4, ownerPost = oo4
    ), "sim bertrand LogitNests")
    results$bertrand_ces_nests <- qa_value(sim(
        p4, supply = "bertrand", demand = "CESNests",
        demand.param = list(gamma = 2, meanval = c(1, .8, .6, .5),
                            sigma = c(3, 4), shareInside = .7),
        nests = nests4, ownerPre = op4, ownerPost = oo4
    ), "sim bertrand CESNests")
    results$bertrand_logit_cap <- qa_value(sim(
        p, supply = "bertrand", demand = "LogitCap",
        demand.param = list(alpha = -1, meanval = c(.4, .2, .1), mktSize = 100),
        capacities = c(25, 20, 15), ownerPre = op, ownerPost = oo
    ), "sim bertrand LogitCap")

    blp <- qa_fixture_blp_parameters()
    results$bertrand_blp <- qa_value(sim(
        p, shares = c(.35, .25, .20), supply = "bertrand", demand = "BLP",
        demand.param = blp, ownerPre = op, ownerPost = oo, insideSize = 100
    ), "sim bertrand BLP")
    results$cournot_blp <- qa_value(sim(
        p, shares = c(.35, .25, .20), supply = "cournot", demand = "BLP",
        demand.param = blp, ownerPre = op, ownerPost = oo, insideSize = 100
    ), "sim cournot BLP")

    expected_classes <- c(
        bertrand_linear = "Linear", bertrand_loglin = "LogLin", bertrand_aids = "AIDS",
        logit_bertrand = "Logit", logit_cournot = "LogitCournot",
        logit_auction2nd = "Auction2ndLogit", logit_bargaining = "BargainingLogit",
        logit_bargaining2nd = "Bargaining2ndLogit", ces_bertrand = "CES",
        ces_cournot = "CESCournot", ces_auction2nd = "Auction2ndCES",
        ces_bargaining = "BargainingCES", ces_bargaining2nd = "Bargaining2ndCES",
        bertrand_logit_nests = "LogitNests", bertrand_ces_nests = "CESNests",
        bertrand_logit_cap = "LogitCap", bertrand_blp = "LogitBLP",
        cournot_blp = "CournotBLP"
    )
    for (nm in names(expected_classes)) {
        testthat::expect_true(is(results[[nm]], expected_classes[[nm]]), info = nm)
        qa_assert_finite(results[[nm]]@pricePost, paste(nm, "post prices"))
    }
    ## Structural LogitCap parameters must survive simulation.  Previously
    ## calcSlopes() recalibrated alpha from capacity-bound observations and
    ## replaced the supplied demand system with an unidentified near-zero
    ## coefficient.
    testthat::expect_equal(unname(results$bertrand_logit_cap@slopes$alpha), -1,
                           tolerance = 1e-12)
    testthat::expect_equal(unname(results$bertrand_logit_cap@slopes$meanval),
                           c(.4, .2, .1), tolerance = 1e-12)
    testthat::expect_true(all(calcQuantities(results$bertrand_logit_cap, TRUE) <=
                              results$bertrand_logit_cap@capacitiesPre + 1e-8))

    ## Nesting-parameter accessors are public methods and should preserve the
    ## one-column sigma representation for both nested demand families.
    logit_nest_parms <- qa_value(getNestsParms(results$bertrand_logit_nests),
                                 "LogitNests nesting parameters")
    ces_nest_parms <- qa_value(getNestsParms(results$bertrand_ces_nests),
                               "CESNests nesting parameters")
    testthat::expect_equal(ncol(logit_nest_parms), 1L)
    testthat::expect_equal(ncol(ces_nest_parms), 1L)
    qa_assert_finite(logit_nest_parms, "LogitNests nesting parameters")
    qa_assert_finite(ces_nest_parms, "CESNests nesting parameters")
})

test_that("LogitCap solves binding-capacity KKT conditions", {
    prices <- c(2, 2.2, 2.5)
    capacities <- c(25, 20, 15)
    shares <- capacities / 100
    meanval <- prices + log(shares / (1 - sum(shares)))

    captured <- qa_capture(sim(
        prices, supply = "bertrand", demand = "LogitCap",
        demand.param = list(alpha = -1, meanval = meanval, mktSize = 100),
        capacities = capacities, margins = rep(.30, 3),
        ownerPre = c("A", "B", "C"), ownerPost = c("A", "A", "C")
    ), "LogitCap binding KKT", disposition = "capacity-boundary")
    fit <- captured$value

    testthat::expect_equal(unname(fit@slopes$alpha), -1, tolerance = 1e-12)
    testthat::expect_true(all(fit@pricePost > fit@pricePre))
    testthat::expect_true(all(calcQuantities(fit, FALSE) <= capacities + 1e-7))
    testthat::expect_true(max(abs(calcQuantities(fit, FALSE) - capacities)) > 1e-3)
})

test_that("unsupported sim combinations fail explicitly", {
    qa_expect_error(sim(
        c(2, 2.2), supply = "cournot", demand = "Linear",
        demand.param = qa_fixture_linear_parameters()[1:2],
        ownerPre = c("A", "B"), ownerPost = c("A", "A")
    ), "currently not supported", "unsupported Linear/Cournot combination")
    qa_expect_error(sim(
        c(2, 2.2), supply = "auction2nd", demand = "BLP",
        demand.param = qa_fixture_blp_parameters(),
        ownerPre = c("A", "B"), ownerPost = c("A", "A"),
        shares = c(.3, .2)
    ), "currently not supported|only compatible", "unsupported BLP/auction combination")
})
