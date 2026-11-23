zeta_fixture <- function(demand, conduct, output = TRUE) {
    prices <- c(2, 2.2, 2.4, 2.6)
    owner <- c("A", "A", "B", "C")
    params <- switch(demand,
        logit = list(alpha = if (output) -1.6 else 1.6,
                     meanval = if (output) c(1, .8, .6, .4) else c(-4, -4.2, -4.5, -5)),
        ces = list(alpha = 1.5, gamma = if (output) 4 else -2, meanval = c(.95, .8, .6, .4)),
        blp = list(alphaMean = -1.6, sigma = .2, meanval = c(1, .8, .6, .4),
                   integration = "provided", draws = c(-1, 0, 1),
                   drawWeights = c(.15, .55, .30)))
    args <- list(demand, conduct, prices = prices, ownerPre = owner,
                 parameters = params, insideSize = 100, output = output)
    if (demand == "blp") {
        args$s0 <- .35
        args$shares <- c(.2, .18, .15, .12)
    } else {
        args$baseline <- "observed"
        args$priceOutside <- if (demand == "ces") 1 else .2
    }
    do.call(specify, args)
}

test_that("post-merger price updates preserve conduct-specific roots with exit and unequal retention", {
    for (demand in c("logit", "ces", "blp")) {
        conducts <- if (demand == "blp") c("bertrand", "moncom") else
            c("bertrand", "cournot", "moncom")
        for (conduct in conducts) {
            fit <- zeta_fixture(demand, conduct)
            scenario <- list(fit, c("AB", "AB", "AB", "C"),
                             mcDelta = c(-.04, 0, -.02, 0),
                             subset = c(TRUE, FALSE, TRUE, TRUE),
                             revenueRetentionPost = c(.94, 1, .97, 1))
            ref <- do.call(simulate, scenario)
            obj <- do.call(simulate, c(scenario, list(solver = "zeta")))
            price <- .zeta_prices(obj)
            label <- paste(demand, conduct)
            expect_false(is.null(price), info = label)
            expect_equal(obj@pricePost, price, tolerance = 1e-7, info = label)
            expect_equal(obj@pricePost, ref@pricePost, tolerance = 1e-6, info = label)
            expect_true(is.na(obj@pricePost[2]), info = label)
            if (conduct == "bertrand") {
                i <- 1L
                profit <- function(p) {
                    obj@pricePost[i] <- p
                    q <- calcQuantities(obj, FALSE)
                    retention <- getRetention(obj, FALSE)
                    sum((obj@pricePost - obj@mcPost) * retention * q *
                        obj@ownerPost[i, ], na.rm = TRUE)
                }
                expect_lt(abs(numDeriv::grad(profit, obj@pricePost[i])), 1e-5)
            }
        }
    }
})

test_that("input-market price updates use markdowns under Logit and CES demand", {
    for (demand in c("logit", "ces")) {
        for (conduct in c("bertrand", "moncom")) {
            fit <- zeta_fixture(demand, conduct, output = FALSE)
            args <- list(fit, c("A", "A", "A", "C"), mcDelta = rep(-.01, 4))
            if (demand == "logit") args$price_domain <- "real"
            result <- do.call(simulate, c(args, list(solver = "zeta")))
            ref <- do.call(simulate, args)
            expect_equal(.zeta_prices(result), ref@pricePost, tolerance = 1e-6)
            expect_equal(result@pricePost, ref@pricePost, tolerance = 1e-6)
        }
    }
})

test_that("a failed fixed-point iteration restarts at the canonical price", {
    fit <- zeta_fixture("logit", "bertrand")
    post <- c("AB", "AB", "AB", "C")
    # The legacy solver retains its diagnostic if both attempts exhaust their budget.
    fit@model@control.equ$maxit <- 1L
    expect_warning(reference <- simulate(fit, post), "iterations exceeded")
    expect_null(.zeta_prices(reference))
    expect_warning(retry <- simulate(fit, post, solver = "zeta"), "iterations exceeded")
    expect_equal(retry@pricePost, reference@pricePost)
})

test_that("zeta preserves the chosen solver and compounded costs on resumed paths", {
    fit <- zeta_fixture("logit", "bertrand")
    cf <- add_step(counterfactual(ownership = c("A", "A", "A", "B")),
                   costs = c(-.05, 0, 0, 0))
    path <- simulate(fit, cf, solver = "zeta")
    step <- counterfactual(costs = rep(-.02, 4))
    resumed <- simulate(path, step)
    ref <- simulate(simulate(fit, cf), step)
    expect_equal(resumed@diagnostics$solver, "zeta")
    expect_equal(final_result(resumed)@pricePost, final_result(ref)@pricePost,
                 tolerance = 1e-7)
})

test_that("zeta preserves default signed rates and the explicit positive-domain error", {
    fit <- specify("logit", "bertrand", prices = c(.2,.1,.3), output = FALSE,
        parameters = list(alpha = .8, meanval = log(c(.405,.315,.18)/.1) -
            .8 * (c(.2,.1,.3) - 1)), ownerPre = letters[1:3],
        priceOutside = 1, insideSize = 100, baseline = "observed")
    post <- c("a","a","c")
    expect_error(simulate(fit, post, mcDelta = rep(-.2,3), solver = "zeta", price_domain = "positive"),
                 class = "antitrust_price_domain_error")
    result <- simulate(fit, post, mcDelta = rep(-.2,3), solver = "zeta",
                       price_domain = "real")
    expect_lt(min(result@pricePost), 0)
    expect_equal(.zeta_prices(result), result@pricePost, tolerance = 1e-7)
})

test_that("a CES monopoly without outside demand retains its no-equilibrium diagnosis", {
    fit <- specify("ces", "bertrand", prices = c(1.5, 1.8),
        parameters = list(gamma = 2, meanval = c(1, .7), shareInside = 1),
        ownerPre = c("A", "B"), insideSize = 100, baseline = "observed")
    expect_error(simulate(fit, c("A", "A"), solver = "zeta"),
                 "no finite Bertrand equilibrium")
})

test_that("bargaining updates preserve product-specific disagreement and retention", {
    for (demand in c("logit", "ces")) for (output in c(TRUE, FALSE)) {
        base <- zeta_fixture(demand, "bertrand", output)
        fit <- specify(demand, "bargaining", prices = base@model@prices,
            parameters = base@parameters, bargpowerPre = c(.15,.35,.5,.65),
            ownerPre = c("A","A","B","C"), insideSize = 100,
            priceOutside = base@model@priceOutside, baseline = "observed", output = output)
        args <- list(fit, c("AB","AB","AB","C"), mcDelta = c(-.04,0,-.02,0),
            bargpowerPost = c(.05,.21,.48,.57), revenueRetentionPost = c(.94,1,.97,1),
            subset = c(TRUE,FALSE,TRUE,TRUE))
        ref <- do.call(simulate, args)
        fast <- do.call(simulate, c(args, list(solver = "zeta")))
        price <- .zeta_prices(fast)
        expect_false(is.null(price), info = paste(demand, output))
        expect_equal(price, ref@pricePost, tolerance = 1e-6)
        expect_equal(fast@pricePost, price, tolerance = 1e-7)
    }
})

test_that("mixed-sign BLP nodes preserve signed Bertrand and MonCom markdowns", {
    for (conduct in c("bertrand", "moncom")) {
        fit <- suppressWarnings(specify("blp", conduct, prices = c(.2,.1,.3),
            parameters = list(alphaMean = .8, sigma = .1, meanval = c(0,.1,.2)),
            ownerPre = LETTERS[1:3], insideSize = 100, output = FALSE))
        expect_true(any(fit@model@slopes$alphas < 0))
        fast <- simulate(fit, c("A","A","C"), mcDelta = rep(-.2,3), solver = "zeta")
        ref <- simulate(fit, c("A","A","C"), mcDelta = rep(-.2,3))
        expect_lt(min(fast@pricePost), 0)
        expect_equal(.zeta_prices(fast), ref@pricePost, tolerance = 1e-7)
        expect_equal(fast@pricePost, ref@pricePost, tolerance = 1e-7)
    }
})

test_that("BLP bargaining integrates disagreement before solving product prices", {
    for (output in c(TRUE, FALSE)) {
        fit <- suppressWarnings(specify("blp", "bargaining", prices = c(1.7,2,2.2,2.5),
            parameters = list(alphaMean = if (output) -1.3 else .8, sigma = .1,
                meanval = if (output) c(.3,.1,.2,.4) else c(-2,-2.3,-2.4,-2.2),
                integration = "provided", draws = c(-1,0,1), drawWeights = c(.2,.5,.3)),
            bargpowerPre = c(.1,.2,.45,.6), ownerPre = c("A","A","B","C"),
            insideSize = 100, output = output))
        args <- list(fit, c("AB","AB","AB","C"), bargpowerPost = c(.1,.4,.3,.55),
            revenueRetentionPost = c(.94,1,.97,1), subset = c(TRUE,FALSE,TRUE,TRUE))
        ref <- do.call(simulate,args)
        fast <- do.call(simulate,c(args,list(solver="zeta")))
        expect_equal(.zeta_prices(fast),ref@pricePost,tolerance=1e-6)
        expect_equal(fast@pricePost,ref@pricePost,tolerance=1e-6)
    }
})
