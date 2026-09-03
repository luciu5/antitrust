## Tests for sequential (multi-step) counterfactual composition: the Markov
## property (each step solves from the previous step's equilibrium, not the
## original baseline), quality shocks, entry, exit persistence, and the
## simultaneous-vs-sequential distinction between combine_counterfactuals()
## and add_step().

.seq_fit <- function() {
    calibrate(
        "logit", "bertrand", prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .2), margins = c(.4, .35, .3),
        ownerPre = c("A", "B", "C"), insideSize = 100
    )
}

test_that("sequential ownership -> ownership uses the promoted post-merger state", {
    fit <- .seq_fit()
    cf <- counterfactual(ownership = c("A", "A", "C"))
    cf <- add_step(cf, ownership = c("A", "A", "A"))
    path <- simulate(fit, cf)

    expect_s4_class(path, "CounterfactualPath")
    expect_length(path@steps, 2L)
    expect_length(path@results, 2L)

    step1 <- result_at(path, 1)
    step2 <- result_at(path, 2)

    ## Markov property: step 2 must start from step 1's solved equilibrium.
    expect_equal(step2@pricePre, step1@pricePost, tolerance = 1e-8)
    expect_equal(step2@mcPre, step1@mcPost, tolerance = 1e-8)

    ## Resetting to the originally calibrated baseline would give a
    ## different (wrong) answer for step 2 -- confirm the two diverge.
    naive <- simulate(fit, ownerPost = c("A", "A", "A"))
    expect_false(isTRUE(all.equal(naive@pricePost, step2@pricePost, tolerance = 1e-8)))

    expect_identical(final_result(path), step2)
})

test_that("sequential costs -> costs compound multiplicatively on the current cost state", {
    fit <- .seq_fit()
    cf <- counterfactual(costs = c(-.05, 0, 0))
    cf <- add_step(cf, costs = c(-.05, 0, 0))
    path <- simulate(fit, cf)

    step2 <- result_at(path, 2)
    ## Two -5% shocks compound to -9.75%, not -10%.
    reference <- simulate(fit, ownerPost = fit@model@ownerPre, mcDelta = c(-.0975, 0, 0))
    expect_equal(step2@mcPost, reference@mcPost, tolerance = 1e-8)
})

test_that("sequential ownership -> costs and costs -> ownership both use the promoted state", {
    fit <- .seq_fit()

    cf1 <- counterfactual(ownership = c("A", "A", "C"))
    cf1 <- add_step(cf1, costs = c(0, -.1, 0))
    path1 <- simulate(fit, cf1)
    step1_a <- result_at(path1, 1)
    step2_a <- result_at(path1, 2)
    expect_equal(step2_a@pricePre, step1_a@pricePost, tolerance = 1e-8)

    cf2 <- counterfactual(costs = c(0, -.1, 0))
    cf2 <- add_step(cf2, ownership = c("A", "A", "C"))
    path2 <- simulate(fit, cf2)
    step1_b <- result_at(path2, 1)
    step2_b <- result_at(path2, 2)
    expect_equal(step2_b@pricePre, step1_b@pricePost, tolerance = 1e-8)

    ## Order matters: these two orderings must not produce the same final price.
    expect_false(isTRUE(all.equal(step2_a@pricePost, step2_b@pricePost, tolerance = 1e-8)))
})

test_that("simulate() preserves the ordinary legacy result for one-step Counterfactuals", {
    fit <- .seq_fit()
    cf <- counterfactual(ownership = c("A", "A", "C"))
    result <- simulate(fit, cf)
    expect_s4_class(result, "Logit")
    expect_false(methods::is(result, "CounterfactualPath"))
})

## ---- Quality --------------------------------------------------------------

test_that("quality = 0 reproduces the baseline meanval", {
    fit <- .seq_fit()
    cf <- counterfactual(quality = c(Prod1 = 0))
    result <- simulate(fit, cf)
    expect_equal(result@slopes$meanval, fit@model@slopes$meanval, tolerance = 1e-12)
})

test_that("quality applies an exact proportional change to meanval", {
    fit <- .seq_fit()
    baseline_meanval <- fit@model@slopes$meanval

    cf <- counterfactual(quality = c(Prod1 = .10))
    result <- simulate(fit, cf)
    expect_equal(result@slopes$meanval[["Prod1"]], baseline_meanval[["Prod1"]] * 1.10, tolerance = 1e-12)

    cf_neg <- counterfactual(quality = c(Prod1 = -.05))
    result_neg <- simulate(fit, cf_neg)
    expect_equal(result_neg@slopes$meanval[["Prod1"]], baseline_meanval[["Prod1"]] * .95, tolerance = 1e-12)
})

test_that("sequential quality shocks compound: +10% then +20% gives 1.32x", {
    fit <- .seq_fit()
    baseline_meanval <- fit@model@slopes$meanval[["Prod1"]]

    cf <- counterfactual(quality = c(Prod1 = .10))
    cf <- add_step(cf, quality = c(Prod1 = .20))
    path <- simulate(fit, cf)

    final <- final_result(path)
    expect_equal(final@slopes$meanval[["Prod1"]], baseline_meanval * 1.10 * 1.20, tolerance = 1e-10)
    expect_equal(final@slopes$meanval[["Prod1"]], baseline_meanval * 1.32, tolerance = 1e-10)
})

test_that("quality persists into a later ownership change", {
    fit <- .seq_fit()
    baseline_meanval <- fit@model@slopes$meanval[["Prod1"]]

    cf <- counterfactual(quality = c(Prod1 = .10))
    cf <- add_step(cf, ownership = c("A", "A", "C"))
    path <- simulate(fit, cf)

    step2 <- result_at(path, 2)
    expect_equal(step2@slopes$meanval[["Prod1"]], baseline_meanval * 1.10, tolerance = 1e-10)
})

test_that("quality persists into a later cost change", {
    fit <- .seq_fit()
    baseline_meanval <- fit@model@slopes$meanval[["Prod1"]]

    cf <- counterfactual(quality = c(Prod1 = .10))
    cf <- add_step(cf, costs = c(0, -.1, 0))
    path <- simulate(fit, cf)

    step2 <- result_at(path, 2)
    expect_equal(step2@slopes$meanval[["Prod1"]], baseline_meanval * 1.10, tolerance = 1e-10)
})

test_that("a quality shock to an exited product errors", {
    fit <- .seq_fit()
    cf <- counterfactual(exit = "Prod1")
    cf <- add_step(cf, quality = c(Prod1 = .10))
    expect_error(simulate(fit, cf), "not active")
})

test_that("quality does not mutate the source Fit", {
    fit <- .seq_fit()
    baseline_meanval <- fit@model@slopes$meanval
    simulate(fit, counterfactual(quality = c(Prod1 = .10)))
    expect_equal(fit@model@slopes$meanval, baseline_meanval, tolerance = 0)
})

test_that("quality is rejected for unsupported model families", {
    fit <- calibrate(
        "logit_cap", "bertrand", prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .2), margins = c(.4, .35, .3),
        ownerPre = c("A", "B", "C"), insideSize = 100, mktSize = 1000,
        capacitiesPre = c(500, 400, 300)
    )
    expect_error(simulate(fit, counterfactual(quality = c(Prod1 = .1))), "does not support")
})

test_that("quality works for Logit Cournot and CES Bertrand/Cournot", {
    fit_lc <- calibrate(
        "logit", "cournot", prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .2), margins = c(.4, .35, .3),
        ownerPre = c("A", "B", "C"), insideSize = 100
    )
    result_lc <- simulate(fit_lc, counterfactual(quality = c(Prod1 = .1)))
    expect_equal(result_lc@slopes$meanval[["Prod1"]],
                 fit_lc@model@slopes$meanval[["Prod1"]] * 1.1, tolerance = 1e-10)

    fit_ces <- calibrate(
        "ces", "bertrand", prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .2), margins = c(.4, .35, .3),
        ownerPre = c("A", "B", "C"), insideSize = 100
    )
    result_ces <- simulate(fit_ces, counterfactual(quality = c(Prod1 = .1)))
    expect_equal(result_ces@slopes$meanval[["Prod1"]],
                 fit_ces@model@slopes$meanval[["Prod1"]] * 1.1, tolerance = 1e-10)

    fit_cesc <- calibrate(
        "ces", "cournot", prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .2), margins = c(.4, .35, .3),
        ownerPre = c("A", "B", "C"), insideSize = 100
    )
    result_cesc <- simulate(fit_cesc, counterfactual(quality = c(Prod1 = .1)))
    expect_equal(result_cesc@slopes$meanval[["Prod1"]],
                 fit_cesc@model@slopes$meanval[["Prod1"]] * 1.1, tolerance = 1e-10)
})

## ---- Exit ------------------------------------------------------------------

test_that("sequential exit persists and combines across steps", {
    fit <- .seq_fit()
    cf <- counterfactual(exit = "Prod1")
    cf <- add_step(cf, exit = "Prod2")
    path <- simulate(fit, cf)

    step1 <- result_at(path, 1)
    step2 <- result_at(path, 2)
    expect_equal(step1@subset, c(FALSE, TRUE, TRUE))
    ## Step 2 must retain Prod1's exit while adding Prod2's.
    expect_equal(step2@subset, c(FALSE, FALSE, TRUE))
})

test_that("simultaneous multi-exit removes both products in one step", {
    fit <- .seq_fit()
    result <- simulate(fit, counterfactual(exit = c("Prod1", "Prod2")))
    expect_equal(result@subset, c(FALSE, FALSE, TRUE))
})

test_that("a cost shock to an exited product errors, not silently ignored", {
    fit <- .seq_fit()
    cf <- counterfactual(exit = "Prod1")
    cf <- add_step(cf, costs = c(Prod1 = -.1, Prod2 = 0, Prod3 = 0))
    expect_error(simulate(fit, cf), "not active")
})

## ---- Entry -----------------------------------------------------------------

test_that("one entrant increases product and firm count by one", {
    fit <- .seq_fit()
    e1 <- entrant(label = "E1", meanval = fit@model@slopes$meanval[[1]] * .5,
                  cost = 1, priceStart = 2)
    result <- simulate(fit, counterfactual(entry = e1))

    expect_equal(length(result@labels), 4L)
    expect_true("E1" %in% result@labels)
    owner_mat <- ownerToMatrix(result, preMerger = TRUE)
    expect_equal(nrow(owner_mat), 4L)
    ## The entrant is a single-product firm: its ownership row/column has
    ## exactly one nonzero entry (itself).
    idx <- match("E1", result@labels)
    expect_equal(sum(owner_mat[idx, ] != 0), 1L)
    expect_true(is.finite(result@pricePost[idx]))
})

test_that("entrant can receive a later quality shock, cost shock, and merger", {
    fit <- .seq_fit()
    e1 <- entrant(label = "E1", meanval = fit@model@slopes$meanval[[1]] * .5,
                  cost = 1, priceStart = 2)

    cf <- counterfactual(entry = e1)
    cf <- add_step(cf, quality = c(E1 = .1))
    path <- simulate(fit, cf)
    step2 <- result_at(path, 2)
    idx <- match("E1", step2@labels)
    expect_equal(step2@slopes$meanval[[idx]],
                 result_at(path, 1)@slopes$meanval[[idx]] * 1.1, tolerance = 1e-10)

    cf2 <- counterfactual(entry = e1)
    cf2 <- add_step(cf2, costs = c(Prod1 = 0, Prod2 = 0, Prod3 = 0, E1 = -.1))
    path2 <- simulate(fit, cf2)
    step2b <- result_at(path2, 2)
    expect_true(is.finite(step2b@mcPost[match("E1", step2b@labels)]))

    cf3 <- counterfactual(entry = e1)
    ## Firm "A" (position 1, Prod1) is assigned to own the entrant (position
    ## 4, E1) post-merger; product labels are distinct from these firm IDs.
    cf3 <- add_step(cf3, ownership = c("A", "A", "C", "A"))
    path3 <- simulate(fit, cf3)
    step2c <- result_at(path3, 2)
    owner_mat <- ownerToMatrix(step2c, preMerger = FALSE)
    idxProd1 <- match("Prod1", step2c@labels)
    idxE1 <- match("E1", step2c@labels)
    expect_equal(owner_mat[idxProd1, idxE1], 1)
})

test_that("two simultaneous entrants and two sequential entrants both work", {
    fit <- .seq_fit()
    e1 <- entrant(label = "E1", meanval = fit@model@slopes$meanval[[1]] * .5, cost = 1, priceStart = 2)
    e2 <- entrant(label = "E2", meanval = fit@model@slopes$meanval[[1]] * .4, cost = 1.1, priceStart = 2)

    simultaneous <- simulate(fit, counterfactual(entry = list(e1, e2)))
    expect_equal(length(simultaneous@labels), 5L)
    expect_true(all(c("E1", "E2") %in% simultaneous@labels))

    cf <- counterfactual(entry = e1)
    cf <- add_step(cf, entry = e2)
    path <- simulate(fit, cf)
    final <- final_result(path)
    expect_equal(length(final@labels), 5L)
    expect_true(all(c("E1", "E2") %in% final@labels))
})

test_that("duplicate entrant labels error", {
    fit <- .seq_fit()
    e_dup1 <- entrant(label = "E1", meanval = .1, cost = 1, priceStart = 2)
    e_dup2 <- entrant(label = "E1", meanval = .2, cost = 1, priceStart = 2)
    expect_error(counterfactual(entry = list(e_dup1, e_dup2)), "duplicate")

    ## An entrant reusing an existing product's label must also error, but
    ## only once simulate() can check it against the model's current labels
    ## (counterfactual() alone has no fitted model to compare against).
    e_existing <- entrant(label = "Prod1", meanval = .1, cost = 1, priceStart = 2)
    expect_error(simulate(fit, counterfactual(entry = e_existing)), "duplicate")
})

test_that("entry is rejected for unsupported model families", {
    fit <- calibrate(
        "logit_cap", "bertrand", prices = c(2, 2.2, 2.5),
        shares = c(.35, .25, .2), margins = c(.4, .35, .3),
        ownerPre = c("A", "B", "C"), insideSize = 100, mktSize = 1000,
        capacitiesPre = c(500, 400, 300)
    )
    e1 <- entrant(label = "E1", meanval = .1, cost = 1, priceStart = 2)
    expect_error(simulate(fit, counterfactual(entry = e1)), "does not support")
})

test_that("entry does not mutate the source Fit", {
    fit <- .seq_fit()
    baseline_labels <- fit@model@labels
    e1 <- entrant(label = "E1", meanval = .1, cost = 1, priceStart = 2)
    simulate(fit, counterfactual(entry = e1))
    expect_equal(fit@model@labels, baseline_labels)
})

## ---- Simultaneous vs sequential -------------------------------------------

test_that("combine_counterfactuals() is simultaneous (one step); add_step() is sequential (two steps)", {
    fit <- .seq_fit()
    q <- c(Prod1 = .1)
    merger <- c("A", "A", "C")

    simultaneous <- combine_counterfactuals(
        counterfactual(quality = q),
        counterfactual(ownership = merger)
    )
    expect_length(simultaneous@steps, 1L)
    result_simultaneous <- simulate(fit, simultaneous)
    expect_false(methods::is(result_simultaneous, "CounterfactualPath"))

    sequential <- counterfactual(quality = q)
    sequential <- add_step(sequential, ownership = merger)
    expect_length(sequential@steps, 2L)
    result_sequential <- simulate(fit, sequential)
    expect_s4_class(result_sequential, "CounterfactualPath")
    expect_length(result_sequential@steps, 2L)
})

## ---- Path continuation (simulate(path, cf)) --------------------------------

test_that("simulate() can resume from a CounterfactualPath's final state", {
    fit <- .seq_fit()
    cf1 <- counterfactual(ownership = c("A", "A", "C"))
    cf1 <- add_step(cf1, ownership = c("A", "A", "A"))
    p1 <- simulate(fit, cf1)

    p2 <- simulate(p1, counterfactual(costs = c(0, -.1, 0)))
    expect_s4_class(p2, "CounterfactualPath")
    expect_length(p2@steps, 3L)
    expect_length(p2@results, 3L)

    resumed_step <- result_at(p2, 3)
    expect_equal(resumed_step@pricePre, final_result(p1)@pricePost, tolerance = 1e-8)

    ## Continuing again should compose further, not reset to any baseline.
    p3 <- simulate(p2, counterfactual(costs = c(0, -.1, 0)))
    expect_length(p3@steps, 4L)
    step4 <- result_at(p3, 4)
    expect_equal(step4@mcDelta, c(0, -.19, 0), tolerance = 1e-10)
})

## ---- Mixed sequence: entry -> quality -> merger -> exit -------------------

test_that("a mixed entry, quality, merger, and exit sequence resolves and tracks labels", {
    fit <- .seq_fit()
    e1 <- entrant(label = "E1", meanval = fit@model@slopes$meanval[[1]] * .5,
                  cost = 1, priceStart = 2)

    cf <- counterfactual(entry = e1)
    cf <- add_step(cf, quality = c(E1 = .1))
    cf <- add_step(cf, ownership = c("A", "A", "C", "A"))
    cf <- add_step(cf, exit = "Prod2")
    path <- simulate(fit, cf)

    expect_length(path@steps, 4L)
    expect_length(path@results, 4L)

    step_entry <- result_at(path, 1)
    step_quality <- result_at(path, 2)
    step_merger <- result_at(path, 3)
    step_exit <- result_at(path, 4)

    expect_equal(length(step_entry@labels), 4L)
    idxE1 <- match("E1", step_quality@labels)
    expect_equal(step_quality@slopes$meanval[[idxE1]],
                 step_entry@slopes$meanval[[idxE1]] * 1.1, tolerance = 1e-10)

    owner_mat <- ownerToMatrix(step_merger, preMerger = FALSE)
    idxProd1 <- match("Prod1", step_merger@labels)
    idxE1m <- match("E1", step_merger@labels)
    expect_equal(owner_mat[idxProd1, idxE1m], 1)

    idxProd2 <- match("Prod2", step_exit@labels)
    expect_false(step_exit@subset[idxProd2])
})
