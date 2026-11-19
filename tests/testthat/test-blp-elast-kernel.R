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

nested_kernel_model <- function(output, sigma_nest) {
    model <- suppressWarnings(antitrust:::.blp_model(
        conduct = "bertrand", prices = c(1.3, 1.6, 1.9, 2.2),
        shares = c(.25, .20, .15, .15),
        margins = c(.20, .25, .30, .35),
        ownerPre = c("A", "B", "C", "D"),
        alphaMean = if (output) -1.2 else 1.2, sigma = .35,
        meanval = c(.4, .2, -.1, -.3),
        draws = c(-1.5, -.4, .6, 1.8), drawWeights = c(.05, .20, .30, .45),
        s0 = .25, output = output, labels = LETTERS[1:4]
    ))
    model@slopes$sigmaNest <- sigma_nest
    model@subset <- c(TRUE, FALSE, TRUE, FALSE)
    model@pricePost <- model@prices + c(.1, .2, .1, .3)
    model
}

## Nested draw Jacobian written independently of the package's term2/term3
## form: d s_j / d p_k = alpha s_j [delta_jk / sigma - (1/sigma - 1) q_k - s_k],
## with q the within-nest shares.
nested_draw_loop <- function(model) {
    s_draw <- calcShares(model, preMerger = FALSE, aggregate = FALSE)
    active <- model@subset
    sig <- model@slopes$sigmaNest
    out <- matrix(0, nrow(s_draw), nrow(s_draw))
    for (r in seq_len(ncol(s_draw))) {
        s <- s_draw[, r]
        s[!active | is.na(s)] <- 0
        if (sum(s) == 0) next
        q <- s / sum(s)
        d <- model@slopes$alphas[r] * s * (diag(1 / sig, length(s)) -
            matrix((1 / sig - 1) * q + s, length(s), length(s), byrow = TRUE))
        out <- out + model@slopes$drawWeights[r] * d * outer(active, active)
    }
    out
}

test_that("nested LogitBLP elasticities match the independent draw loop", {
    ## Guards the vectorized nested kernel: a sign or nest-scaling error in the
    ## cross/own terms would pass the flat test above.
    for (output in c(TRUE, FALSE)) {
        model <- nested_kernel_model(output, .6)
        expected <- nested_draw_loop(model)
        observed <- elast(model, preMerger = FALSE, partial = TRUE)
        expect_equal(unname(observed), expected, tolerance = 1e-13)
        expect_identical(dimnames(observed), list(model@labels, model@labels))
        jac <- numDeriv::jacobian(function(p) {
            m <- model
            m@pricePost[m@subset] <- p
            calcShares(m, preMerger = FALSE)[m@subset]
        }, model@pricePost[model@subset])
        expect_equal(unname(observed[model@subset, model@subset]), jac, tolerance = 1e-7)
    }
})

test_that("a nested draw whose inside shares underflow contributes zero, not NaN", {
    ## Guards term3 = 1/s_g - 1 - (1/sigma)/s_g, which is Inf - Inf at s_g = 0.
    model <- nested_kernel_model(TRUE, .6)
    model@slopes$alphas[1] <- -900
    observed <- elast(model, preMerger = FALSE, partial = TRUE)
    expect_true(all(is.finite(observed)))
    expect_equal(unname(observed), nested_draw_loop(model), tolerance = 1e-13)
})
