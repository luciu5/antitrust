test_that("fixed BLP draws are repeatable and preserve caller RNG state", {
    f <- qa_fixture_market()
    params <- qa_fixture_blp_parameters()
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(20260830)
    before <- .Random.seed
    first <- qa_value(sim(
        f$prices, shares = f$shares, demand = "BLP", demand.param = params,
        ownerPre = f$ownerPre, ownerPost = f$ownerPost, insideSize = 100
    ), "deterministic BLP first run")
    after <- .Random.seed
    second <- qa_value(sim(
        f$prices, shares = f$shares, demand = "BLP", demand.param = params,
        ownerPre = f$ownerPre, ownerPost = f$ownerPost, insideSize = 100
    ), "deterministic BLP second run")

    testthat::expect_identical(before, after)
    testthat::expect_identical(first@pricePost, second@pricePost)
    testthat::expect_identical(first@slopes$consDraws, params$consDraws)
    qa_assert_finite(first@pricePost, "BLP deterministic post prices")
})

test_that("seeded generated BLP draws are repeatable", {
    f <- qa_fixture_market()
    run <- function() {
        set.seed(918273)
        qa_value(sim(
            f$prices, shares = f$shares, demand = "BLP",
            demand.param = list(alpha = -1, sigma = .1, sigmaNest = 1, nDraws = 12),
            ownerPre = f$ownerPre, ownerPost = f$ownerPost, insideSize = 100
        ), "seeded generated BLP")
    }
    first <- run()
    second <- run()
    testthat::expect_identical(first@slopes$consDraws, second@slopes$consDraws)
    testthat::expect_identical(first@pricePost, second@pricePost)
})
