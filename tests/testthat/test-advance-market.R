# A permanent bridge tests: a downstream policy package must advance a solved
# physical-cost market without calling calibration or an equilibrium solver.
test_that("A: advancing a tariff equilibrium retains physical costs and entry primitives without a solve", {
  source <- specify("logit", "bertrand", prices = c(2, 2.2, 2.4),
    ownerPre = c("A", "A", "B"), insideSize = 100,
    parameters = list(alpha = -1.5, meanval = c(2, 1.5, 1)),
    labels = c("a", "b", "c"))
  solved <- suppressWarnings(simulate(source, c("A", "A", "B"),
    mcDelta = c(.1, 0, 0), revenueRetentionPost = c(.9, .7, .8)))
  before <- serialize(solved, NULL)
  physical <- calcMC(solved, FALSE) * getRetention(solved, FALSE)
  step <- counterfactual(entry = entrant("e", 1, .8, 2))@steps[[1L]]
  testthat::local_mocked_bindings(
    calibrate = function(...) stop("unexpected recalibration"),
    specify = function(...) stop("unexpected specification"),
    calcPrices = function(...) stop("unexpected equilibrium solve"),
    calcSlopes = function(...) stop("unexpected demand estimation"),
    .package = "antitrust")
  advanced <- advance_market(solved, step)
  expect_equal(as.numeric(calcMC(advanced, TRUE) * getRetention(advanced, TRUE)),
               c(as.numeric(physical), .8), tolerance = 1e-10)
  expect_equal(advanced@pricePre[1:3], solved@pricePost, tolerance = 1e-12)
  expect_equal(advanced@slopes$meanval[1:3], solved@slopes$meanval,
               ignore_attr = TRUE)
  expect_equal(advanced@mcDelta, rep(0, 4))
  expect_identical(serialize(solved, NULL), before)
  expect_error(advance_market(solved, counterfactual(costs = rep(0, 3))@steps[[1L]]),
               "only entry and quality")
})

# A permanent exit->entry regression: NA prices for an exited product must
# not make the new inside-market mass undefined during an expansion.
test_that("A: advancing after exit leaves a valid inside mass for later entry", {
  source <- specify("logit", "bertrand", prices = c(2, 2.2, 2.4),
    ownerPre = c("A", "B", "C"), insideSize = 100,
    parameters = list(alpha = -1.5, meanval = c(2, 1.5, 1)),
    labels = c("a", "b", "c"))
  solved <- simulate(source, counterfactual(exit = "b"))
  next_model <- advance_market(solved,
    counterfactual(entry = entrant("e", 1, .8, 2))@steps[[1L]])
  expect_true(is.finite(next_model@shareInside))
  expect_equal(as.numeric(next_model@shares[2]), 0)
  expect_false(next_model@subset[2])
})

# A permanent: post-exit promotion must not make an ordinary sequential
# market's consumer variation undefined by summing an exited NA pre-price.
test_that("A: ordinary exit then merger has the independent logsum welfare change", {
  source <- specify("logit", "bertrand", prices = c(2, 2.2, 2.4),
    ownerPre = c("A", "B", "C"), insideSize = 100,
    parameters = list(alpha = -1.5, meanval = c(2, 1.5, 1)),
    labels = c("a", "b", "c"))
  path <- simulate(source, add_step(counterfactual(exit = "c"),
    ownership = c("A", "A", "C")))
  result <- final_result(path)
  active <- c(TRUE, TRUE, FALSE)
  mass <- function(p) 1 + sum(exp(result@slopes$meanval[active] +
    result@slopes$alpha * (p[active] - result@priceOutside)))
  expected <- result@mktSize / result@slopes$alpha *
    log(mass(result@pricePost) / mass(result@pricePre))
  expect_equal(CV(result), as.numeric(expected), tolerance = 1e-10)
  expect_true(is.finite(CV(result)))
})

# A permanent: expansion must append FALSE for a not-yet-entered product,
# then promote it to TRUE after entry; otherwise CV invents prior choices or
# recycles a stale shorter mask after a subsequent merger and exit.
test_that("A: ordinary entry then merger and exit retains aligned pre-choice sets", {
  source <- specify("logit", "bertrand", prices = c(2, 2.2, 2.4),
    ownerPre = c("A", "B", "C"), insideSize = 100,
    parameters = list(alpha = -1.5, meanval = c(2, 1.5, 1)),
    labels = c("a", "b", "c"))
  cf <- add_step(counterfactual(entry = entrant("e", 1, .8, 2)),
    ownership = c("A", "A", "C", "A"), exit = "c")
  path <- simulate(source, cf)
  for (i in 1:2) {
    result <- result_at(path, i)
    pre <- if (i == 1) c(TRUE, TRUE, TRUE, FALSE) else rep(TRUE, 4)
    post <- if (i == 1) rep(TRUE, 4) else c(TRUE, TRUE, FALSE, TRUE)
    mass <- function(p, mask) 1 + sum(exp(result@slopes$meanval[mask] +
      result@slopes$alpha * (p[mask] - result@priceOutside)))
    expected <- result@mktSize / result@slopes$alpha *
      log(mass(result@pricePost, post) / mass(result@pricePre, pre))
    expect_equal(CV(result), as.numeric(expected), tolerance = 1e-10)
  }
  expect_error(CV(final_result(path), subsetPre = c(TRUE, FALSE)), "subsetPre")
  expect_error(CV(final_result(path), subsetPre = rep(FALSE, 4)), "subsetPre")
  expect_error(CV(final_result(path), subsetPre = c(TRUE, TRUE, TRUE, NA)), "subsetPre")
})
