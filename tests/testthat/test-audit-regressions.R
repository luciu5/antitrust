test_that("sim initializes BLP, CES, and bargaining2nd CES correctly", {
  owners_pre <- c("A", "B")
  owners_post <- c("A", "A")

  fit_blp <- qa_value(sim(
    prices = c(10, 12),
    shares = c(0.30, 0.20),
    demand = "BLP",
    demand.param = list(alpha = -2, sigma = 0.1, sigmaNest = 1, nDraws = 20),
    ownerPre = owners_pre,
    ownerPost = owners_post
  ), "sim BLP initialization")
  expect_true(is(fit_blp, "LogitBLP"))
  expect_true(is.na(fit_blp@normIndex))
  expect_true(all(is.finite(fit_blp@pricePost)))

  fit_ces <- qa_value(sim(
    prices = c(1.5, 1.8),
    margins = c(0.30, 0.25),
    demand = "CES",
    demand.param = list(gamma = 2, meanval = c(1, 0.7)),
    ownerPre = owners_pre,
    ownerPost = owners_post
  ), "sim CES initialization")
  expect_true(is(fit_ces, "CES"))
  expect_equal(fit_ces@normIndex, 1)
  expect_true(isTRUE(fit_ces@output))
  expect_true(is.finite(fit_ces@mktSize))

  fit_bargaining <- qa_value(sim(
    prices = c(1.5, 1.8),
    margins = c(0.30, 0.25),
    demand = "CES",
    supply = "bargaining",
    demand.param = list(gamma = 2, meanval = c(1, 0.7)),
    ownerPre = owners_pre,
    ownerPost = owners_post
  ), "sim bargaining CES initialization")
  expect_true(is(fit_bargaining, "BargainingCES"))

  fit_bargaining2nd <- qa_value(sim(
    prices = c(1.5, 1.8),
    margins = c(0.30, 0.25),
    demand = "CES",
    supply = "bargaining2nd",
    demand.param = list(gamma = 2, meanval = c(1, 0.7)),
    ownerPre = owners_pre,
    ownerPost = owners_post
  ), "sim second-score CES initialization")
  expect_true(is(fit_bargaining2nd, "Bargaining2ndCES"))
  expect_true(all(is.finite(fit_bargaining2nd@pricePost)))
})

test_that("bargaining CES ALM and second-score calibration respect bargaining power", {
  prices <- c(2, 2.5, 3)
  shares <- c(0.4, 0.35, 0.25)
  ownerPre <- c("A", "B", "C")
  ownerPost <- c("A", "A", "C")

  fit_ces <- qa_value(ces.alm(
    prices = prices, shares = shares, margins = c(0.45, 0.40, 0.35),
    ownerPre = ownerPre, ownerPost = ownerPost
  ), "CES ALM calibration")
  bargaining_margins <- calcMargins(fit_ces, preMerger = TRUE, level = FALSE) * 0.5
  fit_bargaining_alm <- qa_value(bargaining.ces.alm(
    prices = prices, shares = shares, margins = bargaining_margins,
    bargpowerPre = rep(0.5, 3), ownerPre = ownerPre, ownerPost = ownerPost
  ), "bargaining CES ALM calibration")
  expect_true(is.finite(fit_bargaining_alm@slopes$gamma))
  expect_equal(fit_bargaining_alm@slopes$gamma, fit_ces@slopes$gamma, tolerance = 0.01)

  common <- list(
    prices = c(1.5, 1.5, 1.5), shares = c(0.5, 0.3, 0.2),
    ownerPre = ownerPre, ownerPost = ownerPost
  )
  fit_auction <- qa_value(
    do.call(auction2nd.ces, c(common, list(margins = c(0.60, 0.40, 0.30)))),
    "second-score CES calibration"
  )
  fit_bargaining2nd <- qa_value(
    do.call(bargaining2nd.ces, c(common, list(
      margins = c(0.30, 0.20, 0.15),
      bargpowerPre = rep(0.5, 3),
      mcDeltaOutside = 0.7
    ))),
    "second-score bargaining CES calibration"
  )
  expect_equal(fit_bargaining2nd@slopes$gamma, fit_auction@slopes$gamma, tolerance = 1e-8)
  expect_equal(fit_bargaining2nd@priceOutside, 0.7)
})

test_that("AG solver honors post-merger subsets", {
  args <- list(
    prices = c(2, 2.2, 2.5),
    shares = c(0.35, 0.25, 0.20),
    margins = c(0.40, 0.35, 0.30),
    ownerPre = c("A", "B", "C"),
    ownerPost = c("A", "A", "C"),
    subset = c(TRUE, TRUE, FALSE)
  )
  fit_standard <- qa_value(do.call(logit, c(args, list(solver = "nleqslv"))), "subset nleqslv")
  fit_ag <- qa_value(do.call(logit, c(args, list(solver = "ag"))), "subset AG")
  expect_equal(fit_ag@pricePost, fit_standard@pricePost, tolerance = 1e-4)
  expect_true(is.na(fit_ag@pricePost[3]))
})

test_that("AG CES and second-score wrappers preserve model-specific equations", {
  common <- list(
    prices = c(1.5, 1.8, 2), shares = c(0.30, 0.20, 0.10),
    margins = c(0.30, 0.25, 0.20), ownerPre = c("A", "B", "C"),
    ownerPost = c("A", "A", "C"), output = FALSE
  )
  ces_std <- qa_value(do.call(ces, c(common, list(solver = "nleqslv"))), "CES nleqslv")
  ces_ag <- qa_value(do.call(ces, c(common, list(solver = "ag"))), "CES AG")
  expect_equal(ces_ag@pricePost, ces_std@pricePost, tolerance = 1e-4)

  auction_args <- common
  auction_args$margins <- c(0.60, 0.40, 0.30)
  auction <- qa_value(do.call(auction2nd.ces, auction_args), "second-score CES")
  bargaining_args <- common
  bargaining_args$margins <- c(0.30, 0.20, 0.15)
  bargaining_args$bargpowerPre <- rep(0.5, 3)
  bargaining_args$mcDeltaOutside <- 0.7
  bargaining <- qa_value(do.call(bargaining2nd.ces, bargaining_args), "second-score bargaining CES")
  expect_equal(calcPricesAG(bargaining), calcPrices(bargaining), tolerance = 1e-4)
  expect_equal(bargaining@priceOutside, 0.7)
  expect_true(is.finite(auction@slopes$gamma))
})

test_that("labels remain positional-compatible after solver addition", {
  fit <- qa_value(do.call(logit, list(
    c(1, 1), c(0.4, 0.3), c(0.3, 0.2), matrix(NA_real_, 2, 2),
    c("A", "B"), c("A", "A"), TRUE, rep(1, 2), NA_integer_, rep(0, 2),
    rep(TRUE, 2), NA_real_, 0, c(1, 1), FALSE,
    list(reltol = 1e-6), list(tol = 1e-8), c("X", "Y")
  )), "positional logit arguments")
  expect_equal(fit@labels, c("X", "Y"))
})
