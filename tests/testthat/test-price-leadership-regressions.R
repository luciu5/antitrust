test_that("price leadership examples construct finite equilibria", {
    result <- qa_value(ple(
        prices = c(10, 12, 11, 9),
        shares = c(.25, .20, .18, .17),
        margins = c(.40, NA, .35, .25),
        ownerPre = c("A", "A", "B", "C"),
        ownerPost = c("A", "A", "B", "C"),
        coalitionPre = 1:3,
        insideSize = 1000
    ), "standard price leadership regression")

    testthat::expect_s4_class(result, "PriceLeadership")
    qa_assert_finite(result@pricePre, "price leadership pre-merger prices")
    qa_assert_finite(result@pricePost, "price leadership post-merger prices")
    qa_assert_finite(calcMargins(result, TRUE), "price leadership margins")
    qa_assert_finite(calcShares(result, TRUE), "price leadership shares")
    qa_assert_finite(qa_value(calcSlack(result, TRUE),
                              "price leadership slack values"),
                     "price leadership slack values")
    qa_assert_finite(calcSupermarkup(result, FALSE),
                     "price leadership supermarkup")
    testthat::expect_length(result@timingParam, 2L)
})

test_that("BLP price leadership handles omitted margins", {
    prices <- c(.93, .88, 1.10, 1.02)
    shares <- c(.35, .25, .25, .15)
    owner <- c("Bank1", "Bank2", "Bank3", "Fringe")
    slopes <- list(
        alphaMean = -5.767013,
        alpha = -5.767013,
        sigma = .5,
        meanval = c(0, log(shares[-1] / shares[1]) -
            (-5.767013) * (prices[-1] - prices[1])),
        sigmaNest = 1,
        nDraws = 64,
        consDraws = seq(-1, 1, length.out = 64),
        nDemog = 0,
        piDemog = numeric(0)
    )

    result <- qa_value(ple.blp(
        prices = prices,
        shares = shares,
        ownerPre = owner,
        ownerPost = owner,
        coalitionPre = 1:3,
        coalitionPost = 1:3,
        insideSize = 1000,
        slopes = slopes
    ), "BLP price leadership without margins")

    testthat::expect_s4_class(result, "PriceLeadershipBLP")
    testthat::expect_true(all(is.na(result@margins)))
    qa_assert_finite(result@pricePre, "BLP price leadership pre-merger prices")
    qa_assert_finite(result@pricePost, "BLP price leadership post-merger prices")
    qa_assert_finite(result@supermarkupPost,
                     "BLP price leadership post-merger supermarkup")
})

test_that("LogitCap direct simulation accepts structurally consistent meanval", {
    capacities <- c(30, 25, 20)
    market_size <- 150
    shares <- capacities / market_size
    prices <- c(10, 11, 9)
    meanval <- log(shares / (1 - sum(shares))) - (-1) * prices

    result <- qa_value(sim(
        prices = prices,
        margins = c(.30, .25, .20),
        demand = "LogitCap",
        demand.param = list(alpha = -1, meanval = meanval,
                            mktSize = market_size),
        ownerPre = c("A", "B", "C"),
        ownerPost = c("A", "A", "C"),
        capacities = capacities
    ), "LogitCap direct simulation")

    testthat::expect_s4_class(result, "LogitCap")
    qa_assert_close(result@shares, shares, tolerance = 1e-12,
                    message = "LogitCap shares")
    testthat::expect_equal(result@insideSize, sum(capacities),
                           tolerance = 1e-12)
    qa_assert_finite(result@pricePost, "LogitCap post-merger prices")
})
