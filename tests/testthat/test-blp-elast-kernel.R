test_that("flat LogitBLP elasticities match the independent draw loop", {
    nodes <- c(-1.5, -.4, .6, 1.8)
    weights <- c(.05, .20, .30, .45)
    prices <- c(1.3, 1.6, 1.9, 2.2)
    meanval <- c(.4, .2, -.1, -.3)
    active <- c(TRUE, FALSE, TRUE, FALSE)

    for (output in c(TRUE, FALSE)) {
        alpha <- if (output) -1.2 else 1.2
        model <- suppressWarnings(antitrust:::.blp_model(
            conduct = "bertrand", prices = prices,
            shares = c(.25, .20, .15, .15),
            margins = c(.20, .25, .30, .35),
            ownerPre = c("A", "B", "C", "D"),
            alphaMean = alpha, sigma = .35, meanval = meanval,
            draws = nodes, drawWeights = weights, s0 = .25,
            output = output, labels = LETTERS[1:4]
        ))
        model@subset <- active
        model@pricePost <- prices + c(.1, .2, .1, .3)

        draw_shares <- calcShares(model, preMerger = FALSE,
                                  aggregate = FALSE)
        draw_weights <- model@slopes$drawWeights
        alphas <- model@slopes$alphas
        expected_partial <- matrix(0, nrow = nrow(draw_shares),
                                    ncol = nrow(draw_shares))
        for (r in seq_len(ncol(draw_shares))) {
            s <- draw_shares[, r]
            active_r <- active & !is.na(s)
            s[!active_r] <- 0
            expected_partial <- expected_partial +
                draw_weights[r] * alphas[r] *
                (diag(s) - tcrossprod(s))
        }

        observed_partial <- elast(model, preMerger = FALSE, partial = TRUE)
        expect_equal(unname(observed_partial), expected_partial,
                     tolerance = 1e-13)
        expect_identical(dimnames(observed_partial),
                         list(model@labels, model@labels))

        aggregate_shares <- calcShares(model, preMerger = FALSE)
        aggregate_shares[is.na(aggregate_shares)] <- 0
        expected_elast <- expected_partial *
            outer(1 / pmax(aggregate_shares, 1e-10), model@pricePost)
        expected_elast <- expected_elast * outer(active, active)
        observed_elast <- elast(model, preMerger = FALSE)
        expect_equal(unname(observed_elast), unname(expected_elast),
                     tolerance = 1e-13)
        expect_identical(dimnames(observed_elast),
                         list(model@labels, model@labels))

        expected_market <- sum(
            aggregate_shares[active] / sum(aggregate_shares[active]) *
                rowSums(expected_elast[active, active, drop = FALSE])
        )
        expect_equal(unname(elast(model, preMerger = FALSE, market = TRUE)),
                     expected_market, tolerance = 1e-13)
    }
})
