## A permanent economic regressions: cumulative cost wedges must not be
## reapplied to auction winning scores after promotion, and exited suppliers
## must not reappear in a later pre-state. Price FOCs alone miss the latter.
.auction_sequence_fit <- function(output = TRUE, retention = 1) {
  specify("logit", "auction2nd", prices = c(1.1, 1.2, 1.3),
    ownerPre = c("A", "A", "B"), insideSize = 100, priceOutside = 0,
    parameters = list(alpha = if (output) -6 else 6,
      meanval = log(.8 * c(40, 30, 20) / 90 / .2)),
    labels = c("a", "b", "c"), output = output, baseline = "observed",
    revenueRetentionPre = rep(retention, 3L))
}

.expect_auction_state <- function(model, value) {
  for (pre in c(TRUE, FALSE)) {
    active <- if (pre) attr(model, "antitrust_subset_pre", exact = TRUE) else model@subset
    if (is.null(active)) active <- rep(TRUE, length(model@prices))
    cost <- if (pre) model@mcPre else model@mcPost
    price <- if (pre) model@pricePre else model@pricePost
    owner <- if (pre) model@ownerPre else model@ownerPost
    eta <- value + model@slopes$alpha * cost
    weights <- exp(eta[active])
    expected_share <- weights / (1 + sum(weights))
    share <- calcShares(model, preMerger = pre)
    expect_equal(as.numeric(share[active]), as.numeric(expected_share), tolerance = 1e-10)
    expect_true(all(is.na(share[!active])))
    firm <- as.numeric(owner[active, active, drop = FALSE] %*% expected_share)
    markup <- -log1p(-firm) / (abs(model@slopes$alpha) * firm)
    actual_markup <- if (model@output) price[active] - cost[active] else cost[active] - price[active]
    expect_equal(as.numeric(actual_markup), as.numeric(markup), tolerance = 1e-8)
    expect_equal(as.numeric(calcMargins(model, preMerger = pre, level = TRUE)[active]),
                 as.numeric(markup), tolerance = 1e-8)
  }
}

test_that("A: auction cost, ownership, and no-op paths preserve values and incremental winning scores", {
  for (output in c(TRUE, FALSE)) for (retention in c(1, .8)) {
    fit <- .auction_sequence_fit(output, retention)
    value <- fit@model@slopes$meanval - fit@model@slopes$alpha * fit@model@mcPre
    first <- counterfactual(costs = c(.05, .03, .01))
    cf <- add_step(add_step(first, ownership = rep("A", 3L)))
    path <- suppressWarnings(simulate(fit, cf))
    for (i in 1:3) {
      current <- result_at(path, i)
      .expect_auction_state(current, value)
      if (i > 1L) {
        previous <- result_at(path, i - 1L)
        expect_equal(calcShares(current, TRUE), calcShares(previous, FALSE), tolerance = 1e-10)
        expect_equal(current@mcPre, previous@mcPost)
        expect_equal(current@pricePre, previous@pricePost)
      }
    }
    terminal <- final_result(path)
    expect_equal(terminal@pricePost, terminal@pricePre, tolerance = 1e-10)
    expected_cost <- fit@model@mcPre + c(.05, .03, .01) / retention
    expect_equal(as.numeric(terminal@mcPost), as.numeric(expected_cost), tolerance = 1e-10)
    direct <- suppressWarnings(simulate(fit, rep("A", 3L),
      mcDelta = c(.05, .03, .01), revenueRetentionPost = rep(retention, 3L)))
    expect_equal(terminal@pricePost, direct@pricePost, tolerance = 1e-10)
    resumed <- suppressWarnings(simulate(simulate(fit, add_step(first)),
                        counterfactual(ownership = rep("A", 3L))))
    expect_equal(final_result(resumed)@pricePost, terminal@pricePost, tolerance = 1e-10)

    advanced <- advance_market(result_at(path, 1L))
    expect_equal(calcShares(advanced, TRUE), calcShares(result_at(path, 1L), FALSE),
                 tolerance = 1e-10)
    expect_equal(advanced@mcDelta, rep(0, 3L))
    wrapped <- fit
    wrapped@model <- advanced
    replay <- suppressWarnings(simulate(wrapped, counterfactual(ownership = rep("A", 3L))))
    expect_equal(replay@pricePost, terminal@pricePost, tolerance = 1e-10)
  }
})

test_that("A: exited auction suppliers stay absent through ownership and market advancement", {
  for (output in c(TRUE, FALSE)) {
    fit <- .auction_sequence_fit(output)
    value <- fit@model@slopes$meanval - fit@model@slopes$alpha * fit@model@mcPre
    cf <- add_step(add_step(counterfactual(exit = "b"), ownership = rep("A", 3L)))
    path <- simulate(fit, cf)
    for (i in 1:3) .expect_auction_state(result_at(path, i), value)
    for (i in 2:3) {
      current <- result_at(path, i)
      expect_equal(calcShares(current, TRUE), calcShares(result_at(path, i - 1L), FALSE))
      expect_true(is.na(calcMargins(current, TRUE, level = TRUE)[2L]))
    }
    advanced <- advance_market(result_at(path, 1L))
    expect_equal(calcShares(advanced, TRUE), calcShares(result_at(path, 1L), FALSE))
    expect_true(is.na(calcMargins(advanced, TRUE, level = TRUE)[2L]))
  }
})

test_that("A: nonzero outside auction shocks cannot silently lose revenue state on promotion", {
  fit <- .auction_sequence_fit()
  fit@model@priceOutside <- .1
  solved <- suppressWarnings(simulate(fit, counterfactual(costs = rep(.05, 3L))))
  expect_error(advance_market(solved), "nonzero outside cost shocks.*revenue-state promotion")
  expect_error(suppressWarnings(simulate(fit, add_step(counterfactual(costs = rep(.05, 3L))))),
               "nonzero outside cost shocks.*revenue-state promotion")
})

test_that("A: advancing a retention-shocked auction rebases its winning index once", {
  fit <- .auction_sequence_fit()
  value <- fit@model@slopes$meanval - fit@model@slopes$alpha * fit@model@mcPre
  solved <- suppressWarnings(simulate(fit, fit@model@ownerPre,
    mcDelta = c(.05, .03, .01), revenueRetentionPost = rep(.8, 3L)))
  .expect_auction_state(solved, value)
  advanced <- advance_market(solved)
  expect_equal(calcShares(advanced, TRUE), calcShares(solved, FALSE), tolerance = 1e-10)
  expect_equal(calcShares(advanced, TRUE, revenue = TRUE),
               calcShares(solved, FALSE, revenue = TRUE), tolerance = 1e-10)
  fit@model <- advanced
  replay <- simulate(fit, counterfactual(ownership = rep("A", 3L)))
  .expect_auction_state(replay, value)
  expect_equal(calcShares(replay, FALSE), calcShares(solved, FALSE), tolerance = 1e-10)
})
