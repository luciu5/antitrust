test_that("standalone helpers reject malformed dimensions and values", {
    qa_expect_error(HHI(c(.2, 1.1), c("A", "B")), "between 0 and 1",
                     "HHI share bounds")
    qa_expect_error(HHI(c(.2, .3), c("A")), "same length",
                     "HHI owner dimensions")
    qa_expect_error(HHI(c(.2, .3), matrix(1, 3, 3)), "square matrix|dimensions",
                     "HHI ownership matrix dimensions")
    qa_expect_error(upp.bertrand(c(1, 2), c(.2), matrix(1, 2, 2), c("A", "B")),
                    "same length|length", "UPP dimensions")
})

test_that("constructors reject malformed core inputs", {
    f <- qa_fixture_market()
    qa_expect_error(logit(
        f$prices, f$shares, f$margins[1:2], ownerPre = f$ownerPre,
        ownerPost = f$ownerPost
    ), "same length", "Logit margin dimensions")
    qa_expect_error(logit(
        f$prices, c(1.2, .2, .2), f$margins, ownerPre = f$ownerPre,
        ownerPost = f$ownerPost
    ), "between 0 and 1", "Logit share bounds")
    qa_expect_error(logit(
        f$prices, f$shares, c(.4, -.1, .3), ownerPre = f$ownerPre,
        ownerPost = f$ownerPost
    ), "positive", "Logit margin bounds")
    qa_expect_error(sim(
        f$prices, demand = "BLP", demand.param = list(alpha = -1, sigma = .1),
        ownerPre = f$ownerPre, ownerPost = f$ownerPost
    ), "meanval|shares", "BLP identification inputs")
    qa_expect_error(sim(
        f$prices, demand = "Logit", demand.param = list(alpha = -1, meanval = c(1, 2)),
        ownerPre = f$ownerPre, ownerPost = f$ownerPost
    ), "length-k|length", "Logit meanval dimensions")
})

test_that("BLP dimensions, nesting, capacities, and stochastic controls are validated", {
    f <- qa_fixture_market()
    blp <- qa_fixture_blp_parameters()
    blp_bad <- blp
    blp_bad$nDraws <- 0
    qa_expect_error(sim(
        f$prices, shares = f$shares, demand = "BLP",
        demand.param = blp_bad, ownerPre = f$ownerPre,
        ownerPost = f$ownerPost
    ), "positive scalar", "BLP draw count")
    qa_expect_error(sim(
        f$prices, shares = f$shares, demand = "BLP",
        demand.param = c(blp, list(prodChar = matrix(1, 2, 2), beta = c(1, 1))),
        ownerPre = f$ownerPre, ownerPost = f$ownerPost
    ), "k x L|number of characteristics", "BLP product characteristics")
    qa_expect_error(sim(
        f$prices, demand = "LogitNests",
        demand.param = list(alpha = -1, meanval = f$shares, sigma = .7),
        nests = c("N1", "N2"), ownerPre = f$ownerPre,
        ownerPost = f$ownerPost
    ), "length equals|length-k", "nested Logit dimensions")
    qa_expect_error(sim(
        f$prices, demand = "LogitCap",
        demand.param = list(alpha = -1, meanval = f$shares, mktSize = 10),
        capacities = c(8, 8, 1), ownerPre = f$ownerPre,
        ownerPost = f$ownerPost
    ), "cannot exceed", "capacity market size")
})
