test_that("no-merger and zero-efficiency identities hold", {
    f <- qa_fixture_market()
    no_merger <- qa_value(logit(
        prices = f$prices, shares = f$shares, margins = f$margins,
        ownerPre = f$ownerPre, ownerPost = f$ownerPre,
        mcDelta = rep(0, 3)
    ), "no-merger identity")
    qa_assert_close(no_merger@pricePost, no_merger@pricePre, 1e-8,
                    "no-merger prices")
    qa_assert_close(calcPriceDelta(no_merger), rep(0, 3), 1e-8,
                    "no-merger price changes")

    zero_efficiency <- qa_value(logit(
        prices = f$prices, shares = f$shares, margins = f$margins,
        ownerPre = f$ownerPre, ownerPost = f$ownerPost,
        mcDelta = rep(0, 3)
    ), "zero-efficiency identity")
    with_efficiency <- qa_value(logit(
        prices = f$prices, shares = f$shares, margins = f$margins,
        ownerPre = f$ownerPre, ownerPost = f$ownerPost,
        mcDelta = c(-.05, 0, 0)
    ), "efficiency comparison")
    qa_assert_close(zero_efficiency@mcPost, zero_efficiency@mcPre, 1e-8,
                    "zero-efficiency marginal costs")
    testthat::expect_lt(with_efficiency@pricePost[1], zero_efficiency@pricePost[1])
})

test_that("permutation and ownership-label invariance hold", {
    f <- qa_fixture_market()
    fit <- qa_value(logit(
        f$prices, f$shares, f$margins, qa_fixture_diversions(),
        f$ownerPre, f$ownerPost
    ), "baseline permutation invariant")
    perm <- c(3, 1, 2)
    fit_perm <- qa_value(logit(
        f$prices[perm], f$shares[perm], f$margins[perm],
        qa_fixture_diversions()[perm, perm], f$ownerPre[perm], f$ownerPost[perm]
    ), "product permutation invariant")
    qa_assert_close(fit_perm@pricePost, fit@pricePost[perm], 1e-7,
                    "product permutation")
    qa_assert_close(calcPriceDelta(fit_perm), calcPriceDelta(fit)[perm], 1e-7,
                    "price-delta permutation")

    relabeled <- qa_value(logit(
        f$prices, f$shares, f$margins, qa_fixture_diversions(),
        c("X", "Y", "Z"), c("X", "X", "Z")
    ), "ownership relabeling invariant")
    qa_assert_close(relabeled@pricePost, fit@pricePost, 1e-7,
                    "ownership relabeling")
})

test_that("unit scaling preserves Logit percentage effects", {
    f <- qa_fixture_market()
    fit <- qa_value(logit(
        f$prices, f$shares, f$margins,
        ownerPre = f$ownerPre, ownerPost = f$ownerPost
    ), "unit-scale baseline")
    scaled <- qa_value(logit(
        100 * f$prices, f$shares, f$margins,
        ownerPre = f$ownerPre, ownerPost = f$ownerPost,
        priceOutside = 0
    ), "unit-scale prices")
    qa_assert_close(scaled@pricePost / 100, fit@pricePost, 1e-3,
                    "price unit scaling")
    qa_assert_close(calcPriceDelta(scaled), calcPriceDelta(fit), 2e-3,
                    "price-change unit scaling")
})

test_that("subset simulation excludes products from post-merger outputs", {
    f <- qa_fixture_market()
    f$subset <- c(TRUE, TRUE, FALSE)
    fit <- qa_value(logit(
        f$prices, f$shares, f$margins, ownerPre = f$ownerPre,
        ownerPost = f$ownerPost, subset = f$subset
    ), "subset identity")
    testthat::expect_true(is.na(fit@pricePost[3]))
    testthat::expect_equal(unname(calcShares(fit, FALSE)[3]), 0)
    testthat::expect_true(is.na(calcQuantities(fit, FALSE)[3]))
    testthat::expect_true(is.na(calcRevenues(fit, FALSE)[3]))
})

test_that("standalone HHI has an independently calculated ownership oracle", {
    shares <- c(.20, .30, .10, .40)
    owner <- c("A", "A", "B", "C")
    # A owns the first two products, so its combined share is .50.
    expected <- (50^2 + 10^2 + 40^2)
    testthat::expect_equal(HHI(shares, owner), expected, tolerance = 1e-10)
    testthat::expect_equal(HHI(shares, owner, control = owner), expected,
                           tolerance = 1e-10)
})
