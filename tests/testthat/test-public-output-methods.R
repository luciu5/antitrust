test_that("public output, diagnostic, and ownership methods work on Bertrand families", {
    f <- qa_fixture_market()
    fit <- qa_value(logit(
        f$prices, f$shares, f$margins, qa_fixture_diversions(),
        f$ownerPre, f$ownerPost, labels = f$labels
    ), "public Logit methods")

    outputs <- list(
        prices = calcPrices(fit, TRUE), shares = calcShares(fit, TRUE),
        margins = calcMargins(fit, TRUE), mc = calcMC(fit, TRUE),
        quantities = calcQuantities(fit, TRUE), revenues = calcRevenues(fit, TRUE),
        slopes = calcSlopes(fit), priceDelta = calcPriceDelta(fit),
        elast = elast(fit, TRUE), diversion = diversion(fit, TRUE),
        hhi = hhi(fit, TRUE), upp = upp(fit), cmcr = cmcr(fit),
        cv = CV(fit)
    )
    testthat::expect_true(all(vapply(outputs, function(x) length(x) > 0, logical(1))))
    qa_assert_finite(outputs$prices, "Logit prices")
    qa_assert_finite(outputs$shares, "Logit shares")
    qa_assert_finite(outputs$elast, "Logit elasticities")
    testthat::expect_true(is.list(getParms(fit)))
    testthat::expect_true(is.list(calcDiagnostics(fit)))
    owner_pre_matrix <- qa_value(ownerToMatrix(fit, TRUE),
                                 "Bertrand ownership matrix")
    owner_pre_vector <- qa_value(ownerToVec(fit, TRUE),
                                 "Bertrand ownership vector")
    testthat::expect_equal(dim(owner_pre_matrix), c(3L, 3L))
    testthat::expect_equal(length(owner_pre_vector), 3L)
    testthat::expect_true(all(diag(owner_pre_matrix) == 1))
    testthat::expect_true(length(qa_value(capture.output(summary(fit)),
                                           "summary method")) > 0)
    testthat::expect_true(length(qa_value(capture.output(show(fit)),
                                           "show method")) > 0)
})

test_that("capacity, auction, and bargaining-specific methods have formal coverage", {
    f <- qa_fixture_market()
    cap <- qa_value(logit.cap(
        prices = f$prices, shares = f$shares, margins = f$margins,
        ownerPre = f$ownerPre, ownerPost = f$ownerPost,
        capacitiesPre = c(50, 40, 30), capacitiesPost = c(55, 40, 30),
        insideSize = 100
    ), "Logit capacity methods")
    qa_assert_finite(calcPrices(cap, TRUE), "capacity prices")
    qa_assert_finite(calcMargins(cap, TRUE), "capacity margins")
    qa_assert_finite(calcMC(cap, TRUE), "capacity marginal costs")
    qa_assert_finite(calcQuantities(cap, TRUE), "capacity quantities")

    auction_cap <- qa_value(auction2nd.cap(
        capacities = c(.65, .30, .05),
        margins = c(.228, .209, .197), prices = c(3.89, 3.79, 3.74),
        reserve = NA, shareInside = .67, sellerCostCDF = "punif",
        ownerPre = c("A", "B", "C"), ownerPost = c("A", "A", "C")
    ), "capacity auction methods")
    cap_values <- c(
        calcBuyerExpectedCost(auction_cap), calcBuyerValuation(auction_cap),
        calcExpectedLowestCost(auction_cap), calcExpectedPrice(auction_cap),
        calcOptimalReserve(auction_cap), cdfG(auction_cap),
        calcProducerSurplus(auction_cap)
    )
    qa_assert_finite(cap_values, "capacity-auction outputs")

    bargaining <- qa_value(bargaining.logit(
        f$prices, f$shares, f$margins, f$ownerPre, f$ownerPost,
        bargpowerPre = rep(.5, 3)
    ), "bargaining Logit methods")
    qa_assert_finite(calcPrices(bargaining, TRUE), "bargaining prices")
    qa_assert_finite(calcMargins(bargaining, TRUE), "bargaining margins")
    qa_assert_finite(qa_value(calcProducerSurplus(bargaining, TRUE),
                               "bargaining producer surplus"),
                     "bargaining producer surplus")
})

test_that("Cournot, Stackelberg, and vertical bargaining paths expose stable outputs", {
    n <- 3
    cap <- c(.5, .6, .7)
    intercept <- 10
    slope <- -.25
    bmat <- matrix(slope, nrow = n, ncol = n)
    diag(bmat) <- 2 * slope - 1 / cap
    quantities <- rowSums(solve(bmat) * -intercept)
    prices <- intercept + slope * sum(quantities)
    margins <- 1 - (quantities / cap) / prices
    owner_pre <- diag(n)
    owner_post <- owner_pre
    owner_post[1, 2] <- owner_post[2, 1] <- 1

    cournot_fit <- qa_value(cournot(
        prices = prices, quantities = matrix(quantities, ncol = 1),
        margins = matrix(margins, ncol = 1), ownerPre = owner_pre,
        ownerPost = owner_post
    ), "Cournot methods")
    stack_fit <- qa_value(stackelberg(
        prices = prices, quantities = matrix(quantities, ncol = 1),
        margins = matrix(margins, ncol = 1), ownerPre = owner_pre,
        ownerPost = owner_post, isLeaderPre = matrix(c(FALSE, FALSE, TRUE), ncol = 1),
        isLeaderPost = matrix(TRUE, nrow = n, ncol = 1)
    ), "Stackelberg methods")
    for (fit in list(cournot_fit, stack_fit)) {
        qa_assert_finite(fit@pricePost, "quantity-game prices")
        qa_assert_finite(calcQuantities(fit, TRUE), "quantity-game quantities")
        qa_assert_finite(calcMargins(fit, TRUE), "quantity-game margins")
        qa_assert_finite(calcRevenues(fit, TRUE), "quantity-game revenues")
        qa_assert_finite(hhi(fit, TRUE), "quantity-game HHI")
    }

    pD <- c(63.08158, 50.70465, 95.82960, 83.45267)
    sD <- c(.1293482, .1422541, .4631014, .2152962)
    pU <- c(58.109, 53.31135, 58.109, 53.31135)
    mD <- c(13.04232, 13.04233, 29.53958, 29.53958) / pD
    mU <- c(23.31, 14.78715, 23.31, 14.78715) / pU
    down_owner <- paste0("D", rep(c(1, 2), each = 2))
    up_owner <- paste0("U", rep(c(1, 2), 2))
    vertical_fit <- qa_value(vertical.barg(
        sharesDown = sD, pricesDown = pD, marginsDown = mD,
        ownerPreDown = down_owner, ownerPostDown = down_owner,
        pricesUp = pU, marginsUp = mU, ownerPreUp = up_owner,
        ownerPostUp = rep("U1", 4), priceOutside = 10
    ), "vertical bargaining methods")
    qa_assert_finite(vertical_fit@down@pricePost, "vertical downstream prices")
    qa_assert_finite(vertical_fit@up@pricePost, "vertical upstream prices")
    qa_assert_finite(hhi(vertical_fit, TRUE), "vertical HHI")
    vertical_owner <- qa_value(ownerToMatrix(vertical_fit, TRUE),
                               "vertical ownership transformation")
    testthat::expect_true(is(vertical_owner, "VertBargBertLogit"))
    testthat::expect_equal(dim(vertical_owner@ownerDownPre), c(4L, 4L))
    testthat::expect_true(is.list(calcDiagnostics(vertical_fit)))
})
