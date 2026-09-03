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
    qa_expect_error(sim(
        c(1.5, 1.8), demand = "CES",
        demand.param = list(gamma = 2, meanval = c(1, .7)),
        ownerPre = c("A", "B"), ownerPost = c("A", "A")
    ), "no finite Bertrand equilibrium", "CES no-outside monopoly")
    qa_expect_error(ces(
        prices = c(1.5, 1.8), shares = c(.60, .40),
        margins = c(.20, .20), ownerPre = c("A", "B"),
        ownerPost = c("A", "A"), solver = "ag"
    ), "no finite Bertrand equilibrium", "CES AG no-outside monopoly")
    qa_expect_error(sim(
        c(2, 2.2), demand = "LogLin",
        demand.param = list(
            slopes = matrix(c(-1.2, .1, .1, -1), 2, byrow = TRUE),
            intercepts = c(3, 3)
        ), ownerPre = c("A", "B"), ownerPost = c("A", "A")
    ), "strictly below -1", "LogLin unit-elastic boundary")
    qa_expect_error(sim(
        c(2, 2.2), demand = "LogitCap",
        demand.param = list(alpha = -1, mktSize = 100),
        capacities = c(25, 20), margins = c(.3, .3),
        ownerPre = c("A", "B"), ownerPost = c("A", "A")
    ), "meanval", "LogitCap demand identification")
    qa_expect_error(sim(
        c(2, 2.2), demand = "LogitCap",
        demand.param = list(alpha = -1, meanval = c(.2, .1), mktSize = 100),
        capacities = c(Inf, 20), margins = c(.3, .3),
        ownerPre = c("A", "B"), ownerPost = c("A", "A")
    ), "non-negative|finite", "LogitCap capacity finiteness")
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

test_that("LogitCap accepts +Inf as the unbounded-capacity sentinel", {
    object <- new("LogitCap",
        prices = c(2, 2.2),
        shares = c(.2, .3),
        margins = c(.3, .3),
        capacitiesPre = c(Inf, 20),
        capacitiesPost = c(Inf, 20),
        mcDelta = c(0, 0),
        insideSize = 30,
        ownerPre = c("A", "B"),
        ownerPost = c("A", "A"),
        subset = c(TRUE, TRUE),
        priceStart = c(2, 2.2),
        shareInside = .5,
        labels = c("P1", "P2")
    )
    testthat::expect_true(methods::validObject(object, test = TRUE))

    testthat::expect_error(new("LogitCap",
        prices = c(2, 2.2),
        shares = c(.2, .3),
        margins = c(.3, .3),
        capacitiesPre = c(-Inf, 20),
        capacitiesPost = c(Inf, 20),
        mcDelta = c(0, 0),
        insideSize = 30,
        ownerPre = c("A", "B"),
        ownerPost = c("A", "A"),
        subset = c(TRUE, TRUE),
        priceStart = c(2, 2.2),
        shareInside = .5,
        labels = c("P1", "P2")
    ), "non-negative")
})
