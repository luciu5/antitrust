test_that("sim initializes BLP, CES, and bargaining2nd CES correctly", {
  owners_pre <- c("A", "B")
  owners_post <- c("A", "A")

  fit_blp <- suppressWarnings(suppressMessages(sim(
    prices = c(10, 12),
    shares = c(0.30, 0.20),
    demand = "BLP",
    demand.param = list(alpha = -2, sigma = 0.1, sigmaNest = 1, nDraws = 20),
    ownerPre = owners_pre,
    ownerPost = owners_post
  )))
  expect_true(is(fit_blp, "LogitBLP"))
  expect_true(is.na(fit_blp@normIndex))
  expect_true(all(is.finite(fit_blp@pricePost)))

  fit_ces <- suppressWarnings(sim(
    prices = c(1.5, 1.8),
    margins = c(0.30, 0.25),
    demand = "CES",
    demand.param = list(gamma = 2, meanval = c(1, 0.7)),
    ownerPre = owners_pre,
    ownerPost = owners_post
  ))
  expect_true(is(fit_ces, "CES"))
  expect_equal(fit_ces@normIndex, 1)
  expect_true(isTRUE(fit_ces@output))
  expect_true(is.finite(fit_ces@mktSize))

  fit_bargaining <- suppressWarnings(sim(
    prices = c(1.5, 1.8),
    margins = c(0.30, 0.25),
    demand = "CES",
    supply = "bargaining",
    demand.param = list(gamma = 2, meanval = c(1, 0.7)),
    ownerPre = owners_pre,
    ownerPost = owners_post
  ))
  expect_true(is(fit_bargaining, "BargainingCES"))

  fit_bargaining2nd <- suppressWarnings(sim(
    prices = c(1.5, 1.8),
    margins = c(0.30, 0.25),
    demand = "CES",
    supply = "bargaining2nd",
    demand.param = list(gamma = 2, meanval = c(1, 0.7)),
    ownerPre = owners_pre,
    ownerPost = owners_post
  ))
  expect_true(is(fit_bargaining2nd, "Bargaining2ndCES"))
  expect_true(all(is.finite(fit_bargaining2nd@pricePost)))
})

test_that("bargaining CES ALM and second-score calibration respect bargaining power", {
  prices <- c(2, 2.5, 3)
  shares <- c(0.4, 0.35, 0.25)
  ownerPre <- c("A", "B", "C")
  ownerPost <- c("A", "A", "C")

  fit_ces <- suppressWarnings(ces.alm(
    prices = prices, shares = shares, margins = c(0.45, 0.40, 0.35),
    ownerPre = ownerPre, ownerPost = ownerPost
  ))
  bargaining_margins <- calcMargins(fit_ces, preMerger = TRUE, level = FALSE) * 0.5
  fit_bargaining_alm <- suppressWarnings(bargaining.ces.alm(
    prices = prices, shares = shares, margins = bargaining_margins,
    bargpowerPre = rep(0.5, 3), ownerPre = ownerPre, ownerPost = ownerPost
  ))
  expect_true(is.finite(fit_bargaining_alm@slopes$gamma))
  expect_equal(fit_bargaining_alm@slopes$gamma, fit_ces@slopes$gamma, tolerance = 0.01)

  common <- list(
    prices = c(1.5, 1.5, 1.5), shares = c(0.5, 0.3, 0.2),
    ownerPre = ownerPre, ownerPost = ownerPost
  )
  fit_auction <- suppressWarnings(do.call(auction2nd.ces, c(
    common, list(margins = c(0.60, 0.40, 0.30))
  )))
  fit_bargaining2nd <- suppressWarnings(do.call(bargaining2nd.ces, c(
    common, list(
      margins = c(0.30, 0.20, 0.15),
      bargpowerPre = rep(0.5, 3),
      mcDeltaOutside = 0.7
    )
  )))
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
  fit_standard <- suppressWarnings(do.call(logit, c(args, list(solver = "nleqslv"))))
  fit_ag <- suppressWarnings(do.call(logit, c(args, list(solver = "ag"))))
  expect_equal(fit_ag@pricePost, fit_standard@pricePost, tolerance = 1e-4)
  expect_true(is.na(fit_ag@pricePost[3]))
})

test_that("AG CES and second-score wrappers preserve model-specific equations", {
  common <- list(
    prices = c(1.5, 1.8, 2), shares = c(0.30, 0.20, 0.10),
    margins = c(0.30, 0.25, 0.20), ownerPre = c("A", "B", "C"),
    ownerPost = c("A", "A", "C"), output = FALSE
  )
  ces_std <- suppressWarnings(do.call(ces, c(common, list(solver = "nleqslv"))))
  ces_ag <- suppressWarnings(do.call(ces, c(common, list(solver = "ag"))))
  expect_equal(ces_ag@pricePost, ces_std@pricePost, tolerance = 1e-4)

  auction_args <- common
  auction_args$margins <- c(0.60, 0.40, 0.30)
  auction <- suppressWarnings(do.call(auction2nd.ces, auction_args))
  bargaining_args <- common
  bargaining_args$margins <- c(0.30, 0.20, 0.15)
  bargaining_args$bargpowerPre <- rep(0.5, 3)
  bargaining_args$mcDeltaOutside <- 0.7
  bargaining <- suppressWarnings(do.call(bargaining2nd.ces, bargaining_args))
  expect_equal(calcPricesAG(bargaining), calcPrices(bargaining), tolerance = 1e-4)
  expect_equal(bargaining@priceOutside, 0.7)
  expect_true(is.finite(auction@slopes$gamma))
})

test_that("labels remain positional-compatible after solver addition", {
  fit <- suppressWarnings(do.call(logit, list(
    c(1, 1), c(0.4, 0.3), c(0.3, 0.2), matrix(NA_real_, 2, 2),
    c("A", "B"), c("A", "A"), TRUE, rep(1, 2), NA_integer_, rep(0, 2),
    rep(TRUE, 2), NA_real_, 0, c(1, 1), FALSE,
    list(reltol = 1e-6), list(tol = 1e-8), c("X", "Y")
  )))
  expect_equal(fit@labels, c("X", "Y"))
})
