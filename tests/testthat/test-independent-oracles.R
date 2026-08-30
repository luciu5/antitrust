test_that("Logit shares, margins, and accounting satisfy hand-coded oracles", {
    f <- qa_fixture_market()
    fit <- qa_value(logit(
        prices = f$prices, shares = f$shares, margins = f$margins,
        ownerPre = f$ownerPre, ownerPost = f$ownerPost,
        labels = f$labels
    ), "oracle Logit calibration")

    alpha <- unname(fit@slopes$alpha)
    outside <- 1 - sum(f$shares)
    meanval <- log(f$shares / outside) - alpha * f$prices
    expected <- exp(meanval + alpha * fit@pricePre)
    expected <- expected / (1 + sum(expected))
    qa_assert_close(calcShares(fit, preMerger = TRUE), expected, 1e-8,
                    "Logit share oracle")

    expected_margins <- -1 / (alpha * f$prices * (1 - f$shares))
    qa_assert_close(calcMargins(fit, preMerger = TRUE), expected_margins, 1e-7,
                    "Logit Bertrand margin oracle")
    qa_assert_close(calcRevenues(fit, preMerger = TRUE) / f$prices,
                    calcQuantities(fit, preMerger = TRUE), 1e-8,
                    "revenue = price * quantity")
    qa_assert_close(calcMC(fit, preMerger = TRUE),
                    f$prices * (1 - expected_margins), 1e-7,
                    "margin and marginal-cost accounting")

    ## For a multi-product Bertrand firm, the post-merger Logit FOC has a
    ## closed form in terms of the firm's aggregate share.  This checks the
    ## equilibrium residual independently of the package's price solver.
    post_shares <- calcShares(fit, preMerger = FALSE)
    firm_share <- as.numeric(fit@ownerPost %*% post_shares)
    expected_post_markup <- -1 / (alpha * (1 - firm_share))
    qa_assert_close(fit@pricePost - fit@mcPost, expected_post_markup,
                    1e-6, "post-merger Logit FOC residual")
})

test_that("CES revenue and quantity shares satisfy the independent CES oracle", {
    fit <- qa_value(ces(
        prices = c(1.5, 1.8, 2.1), shares = c(.35, .25, .20),
        margins = c(.30, .25, .20), ownerPre = c("A", "B", "C"),
        ownerPost = c("A", "A", "C"), insideSize = 100
    ), "oracle CES calibration")

    gamma <- fit@slopes$gamma
    meanval <- fit@slopes$meanval
    prices <- fit@pricePre
    raw <- meanval * prices^(1 - gamma)
    outside_raw <- fit@priceOutside^(1 - gamma)
    expected_revenue <- raw / (sum(raw) + outside_raw)
    expected_quantity <- (expected_revenue / prices) /
        ((1 - sum(expected_revenue)) / fit@priceOutside +
             sum(expected_revenue / prices))

    qa_assert_close(calcShares(fit, TRUE, revenue = TRUE), expected_revenue,
                    1e-8, "CES revenue-share oracle")
    qa_assert_close(calcShares(fit, TRUE, revenue = FALSE), expected_quantity,
                    1e-8, "CES quantity-share oracle")
    qa_assert_close(calcQuantities(fit, TRUE),
                    calcRevenues(fit, TRUE) / prices, 1e-8,
                    "CES revenue/quantity accounting")
})

test_that("Linear calibration recovers the hand-coded diversion-to-slope map", {
    f <- qa_fixture_market()
    q <- c(30, 22, 16)
    d <- qa_fixture_diversions()
    fit <- qa_value(linear(
        prices = f$prices, quantities = q, margins = f$margins,
        diversions = d, symmetry = FALSE,
        ownerPre = f$ownerPre, ownerPost = f$ownerPost
    ), "oracle Linear calibration")

    owner <- diag(length(q))
    expected <- matrix(f$margins * f$prices, nrow = length(q),
                       ncol = length(q), byrow = TRUE)
    expected <- diag(owner) / rowSums(expected * d * owner) * q
    expected <- -t(expected * d)
    qa_assert_close(fit@slopes, expected, 1e-8,
                    "Linear diversion-to-slope oracle")
    qa_assert_close(calcMargins(fit, TRUE), f$margins, 1e-8,
                    "Linear pre-merger FOC oracle")
})

test_that("AIDS and PCAIDS use documented price-share transformations", {
    p <- c(2, 2.2, 2.5)
    s <- c(.40, .35, .25)
    d <- qa_fixture_diversions()
    op <- c("A", "B", "C")
    oo <- c("A", "A", "C")

    aids_fit <- qa_value(aids(
        shares = s, margins = c(.40, .35, .30), prices = p,
        diversions = d, ownerPre = op, ownerPost = oo,
        mktElast = -1.2, priceStart = rep(.2, 3)
    ), "oracle AIDS calibration")
    expected_pre <- aids_fit@shares
    qa_assert_close(calcShares(aids_fit, TRUE, revenue = TRUE), expected_pre,
                    1e-7, "AIDS pre-merger share oracle")
    qa_assert_close(calcQuantities(aids_fit, TRUE),
                    calcRevenues(aids_fit, TRUE) / calcPrices(aids_fit, TRUE),
                    1e-8, "AIDS revenue/quantity accounting")

    pcaids_fit <- qa_value(pcaids(
        shares = s, knownElast = -2, mktElast = -1.2, prices = p,
        diversions = d, ownerPre = op, ownerPost = oo,
        priceStart = rep(.2, 3)
    ), "oracle PCAIDS calibration")
    testthat::expect_true(is(pcaids_fit, "PCAIDS"))
    testthat::expect_true(length(getParms(pcaids_fit)) > 0)
})
