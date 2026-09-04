blp_quadrature_market <- function(demand.param) {
    qa_value(sim(
        prices = c(2.00, 2.20, 2.50),
        shares = c(.35, .25, .20),
        supply = "bertrand", demand = "BLP",
        demand.param = demand.param,
        ownerPre = c("A", "B", "C"),
        ownerPost = c("A", "A", "C"),
        insideSize = 100
    ), "BLP quadrature market")
}


test_that("Gauss-Hermite nodes integrate standard-normal moments", {
    rule <- antitrust:::.blp_normal_nodes(10)
    expect_equal(sum(rule$weights), 1, tolerance = 1e-14)
    expect_equal(sum(rule$weights * rule$nodes), 0, tolerance = 1e-14)
    expect_equal(sum(rule$weights * rule$nodes^2), 1, tolerance = 1e-13)
    expect_equal(sum(rule$weights * rule$nodes^3), 0, tolerance = 1e-13)
    expect_equal(sum(rule$weights * rule$nodes^4), 3, tolerance = 1e-12)
    expect_equal(sum(rule$weights * rule$nodes^6), 15, tolerance = 1e-11)

    reference <- stats::integrate(function(z) exp(.3 * z) * stats::dnorm(z),
                                  -Inf, Inf)$value
    expect_equal(sum(rule$weights * exp(.3 * rule$nodes)), reference,
                 tolerance = 1e-10)
})


test_that("price-only BLP defaults to deterministic Gauss-Hermite", {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(20260903)
    before <- .Random.seed
    first <- blp_quadrature_market(list(alpha = -1, sigma = .1))
    after <- .Random.seed
    second <- blp_quadrature_market(list(alpha = -1, sigma = .1))

    expect_identical(before, after)
    expect_identical(first@slopes$integration, "gauss-hermite")
    expect_length(first@slopes$consDraws, 31L)
    expect_equal(sum(first@slopes$drawWeights), 1, tolerance = 1e-14)
    expect_identical(first@slopes$consDraws, second@slopes$consDraws)
    expect_identical(first@pricePost, second@pricePost)
})


test_that("one demographic with zero price sigma defaults to Gauss-Hermite", {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(20260904)
    before <- .Random.seed
    first <- blp_quadrature_market(list(
        alpha = -1, sigma = 0, piDemog = .1,
        demogMean = .3, demogCov = matrix(.25, nrow = 1, ncol = 1),
        meanval = c(.5, .3, .1)
    ))
    after <- .Random.seed
    second <- blp_quadrature_market(list(
        alpha = -1, sigma = 0, piDemog = .1,
        demogMean = .3, demogCov = matrix(.25, nrow = 1, ncol = 1),
        meanval = c(.5, .3, .1)
    ))

    expect_identical(before, after)
    expect_identical(first@slopes$integration, "gauss-hermite")
    expect_length(first@slopes$consDraws, 31L)
    expect_equal(
        first@slopes$demogDraws[, 1],
        .3 + .5 * first@slopes$consDraws,
        tolerance = 1e-14
    )
    expect_equal(
        first@slopes$alphas,
        -1 + .1 * (first@slopes$demogDraws[, 1] - .3),
        tolerance = 1e-14
    )
    expect_identical(first@slopes$consDraws, second@slopes$consDraws)
    expect_identical(first@pricePost, second@pricePost)
})


test_that("BLP integration mode distinguishes nodes, draws, and supplied points", {
    legacy_mc <- blp_quadrature_market(list(alpha = -1, sigma = .1, nDraws = 12))
    expect_identical(legacy_mc@slopes$integration, "monte-carlo")
    expect_length(legacy_mc@slopes$consDraws, 12L)
    expect_error(
        blp_quadrature_market(list(alpha = -1, sigma = .1,
                                   integration = "gauss-hermite", nDraws = 12)),
        "only valid with integration = 'monte-carlo'"
    )

    quadrature <- blp_quadrature_market(list(
        alpha = -1, sigma = .1, integration = "gauss-hermite", nNodes = 15
    ))
    expect_identical(quadrature@slopes$integration, "gauss-hermite")
    expect_length(quadrature@slopes$consDraws, 15L)

    set.seed(4711)
    mc_first <- blp_quadrature_market(list(
        alpha = -1, sigma = .1, integration = "monte-carlo", nDraws = 12
    ))
    set.seed(4711)
    mc_second <- blp_quadrature_market(list(
        alpha = -1, sigma = .1, integration = "monte-carlo", nDraws = 12
    ))
    expect_identical(mc_first@slopes$integration, "monte-carlo")
    expect_identical(mc_first@slopes$consDraws, mc_second@slopes$consDraws)
    expect_equal(mc_first@slopes$drawWeights, rep(1 / 12, 12))

    supplied <- blp_quadrature_market(list(
        alpha = -1, sigma = .1, consDraws = c(-1, 0, 2),
        integrationWeights = c(.2, .3, .5)
    ))
    expect_identical(supplied@slopes$integration, "provided")
    expect_equal(supplied@slopes$drawWeights, c(.2, .3, .5))
    draw_shares <- calcShares(supplied, aggregate = FALSE)
    expect_equal(
        unname(calcShares(supplied)),
        unname(as.vector(draw_shares %*% c(.2, .3, .5))),
        tolerance = 1e-12
    )
})


test_that("weighted BLP derivatives match a finite difference", {
    fit <- blp_quadrature_market(list(
        alpha = -1, sigma = .1, consDraws = c(-1, 0, 2),
        integrationWeights = c(.2, .3, .5)
    ))
    partials <- elast(fit, partial = TRUE)
    bumped <- fit
    step <- 1e-6
    bumped@pricePre[1] <- bumped@pricePre[1] + step
    finite_difference <- (calcShares(bumped, preMerger = TRUE) -
        calcShares(fit, preMerger = TRUE)) / step
    expect_equal(unname(partials[, 1]), unname(finite_difference),
                 tolerance = 1e-5)
})


test_that("weighted quantiles retain equal-weight compatibility", {
    x <- c(-1, 0, 2)
    expect_equal(
        antitrust:::.blp_weighted_quantile(x, c(.25, .75), c(.2, .3, .5)),
        c(0, 2)
    )
    expect_equal(
        antitrust:::.blp_weighted_quantile(x, c(.25, .75), rep(1 / 3, 3)),
        as.numeric(stats::quantile(x, c(.25, .75), names = FALSE, type = 7))
    )
})


test_that("multidimensional BLP retains Monte Carlo fallback", {
    multidimensional <- blp_quadrature_market(list(
        alpha = -1, sigma = .1, piDemog = .05,
        integration = "monte-carlo", nDraws = 12
    ))
    expect_identical(multidimensional@slopes$integration, "monte-carlo")
    expect_error(
        blp_quadrature_market(list(
            alpha = -1, sigma = .1, piDemog = .05,
            integration = "gauss-hermite", nNodes = 12
        )),
        "only for one-dimensional"
    )
})


test_that("PriceLeadershipBLP reuses the shared quadrature rule", {
    shares <- c(.35, .25, .25, .15)
    result <- qa_value(ple.blp(
        prices = c(.93, .88, 1.10, 1.02),
        shares = shares,
        ownerPre = c("Bank1", "Bank2", "Bank3", "Fringe"),
        ownerPost = c("Bank1", "Bank2", "Bank3", "Fringe"),
        coalitionPre = 1:3,
        coalitionPost = 1:3,
        insideSize = 1000,
        slopes = list(
            alphaMean = -5.767013,
            alpha = -5.767013,
            sigma = .5,
            meanval = c(0, log(shares[-1] / shares[1]) -
                (-5.767013) * (c(.88, 1.10, 1.02) - .93)),
            sigmaNest = 1
        )
    ), "PriceLeadershipBLP quadrature")
    expect_identical(result@slopes$integration, "gauss-hermite")
    expect_length(result@slopes$consDraws, 31L)
    expect_equal(sum(result@slopes$drawWeights), 1, tolerance = 1e-14)
})
