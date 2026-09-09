# Independent reference-loop regression for the aggregate BargainingBLP FOC.

blp_bargaining_reference <- function(object, preMerger = TRUE) {
    active <- if (preMerger) rep(TRUE, length(object@shares)) else object@subset
    prices <- if (preMerger) object@pricePre else object@pricePost
    owner <- if (preMerger) object@ownerPre else object@ownerPost
    barg <- if (preMerger) object@bargpowerPre else object@bargpowerPost
    prices <- prices[active]
    owner <- owner[active, active, drop = FALSE]
    S <- calcShares(object, preMerger = preMerger, aggregate = FALSE)
    S <- S[active, , drop = FALSE]
    w <- object@slopes$drawWeights
    alpha <- object@slopes$alphas
    shares <- as.vector(S %*% w)

    ## Keep this draw loop independent of the vectorized implementation.
    derivative <- matrix(0, nrow = nrow(S), ncol = nrow(S))
    buyer_surplus <- numeric(nrow(S))
    for (r in seq_len(ncol(S))) {
        s <- S[, r]
        derivative <- derivative + w[r] * alpha[r] *
            (diag(s) - tcrossprod(s))
        buyer_surplus <- buyer_surplus + w[r] * log1p(-s) / alpha[r]
    }

    aggregate_elast <- derivative * outer(1 / shares, prices)
    revenue <- prices * shares
    margin_system <- t(
        diag(1 / revenue) %*%
            (t(aggregate_elast * owner) %*% diag(shares))
    )
    barg <- barg[active] / (1 - barg[active])
    own_normalized <- diag(derivative) / shares
    output <- ifelse(object@output, -1, 1)
    rhs <- own_normalized /
        (output * (own_normalized - barg * shares / buyer_surplus))
    rhs <- diag(owner) * rhs
    level <- as.vector(solve(t(margin_system), rhs))
    proportional <- level / prices
    list(level = level, proportional = proportional)
}


test_that("BargainingBLP vectorized kernel matches an independent draw loop", {
    prices <- c(1.75, 2.10, 2.45)
    shares <- c(.30, .25, .25)
    s0 <- 1 - sum(shares)
    nodes <- c(-1.5, -.25, .8, 1.75)
    weights <- c(.10, .20, .30, .40)
    alpha <- -1.35
    sigma <- .35
    delta <- antitrust:::.blp_contract(
        prices, shares, alpha, sigma, nodes, weights, s0
    )$delta
    owner <- matrix(c(
        .7, .1, 0,
        .1, .7, 0,
        0, 0, .8
    ), nrow = 3, byrow = TRUE)
    object <- suppressWarnings(antitrust:::.blp_model(
        conduct = "bargaining", prices = prices, shares = shares,
        margins = rep(.2, 3), ownerPre = owner,
        alphaMean = alpha, sigma = sigma, meanval = delta,
        draws = nodes, drawWeights = weights, s0 = s0,
        output = TRUE, bargpowerPre = c(.2, .35, .1),
        bargpowerPost = c(.25, .30, .15)
    ))
    object@ownerPost <- owner
    object@pricePost <- prices + c(.1, -.05, .08)
    object@subset <- c(TRUE, FALSE, TRUE)

    for (preMerger in c(TRUE, FALSE)) {
        expected <- blp_bargaining_reference(object, preMerger)
        observed_level <- calcMargins(object, preMerger = preMerger, level = TRUE)
        observed_proportional <- calcMargins(object, preMerger = preMerger,
                                             level = FALSE)
        active <- if (preMerger) rep(TRUE, 3) else object@subset
        expect_equal(unname(observed_level[active]), expected$level,
                     tolerance = 1e-12, info = paste("preMerger", preMerger))
        expect_equal(unname(observed_proportional[active]), expected$proportional,
                     tolerance = 1e-12, info = paste("preMerger", preMerger))
        expect_true(all(is.na(observed_level[!active])))
    }
})


test_that("adaptive multistart fallback decision is deterministic", {
    pilot <- antitrust:::.blp_multistart_decision(
        convergence_count = 1L, objective_values = 1e-4,
        artificial_boundary = FALSE, invalid_evaluations = 0L,
        strategy = "adaptive"
    )
    expect_true(pilot$fallback)
    expect_true(any(grepl("fewer than two", pilot$reasons, fixed = TRUE)))

    boundary <- antitrust:::.blp_multistart_decision(
        convergence_count = 3L, objective_values = c(1e-12, 2e-12, 3e-12),
        artificial_boundary = TRUE, invalid_evaluations = 0L,
        strategy = "adaptive"
    )
    expect_true(boundary$fallback)
    expect_true(any(grepl("artificial parameter bound", boundary$reasons,
                          fixed = TRUE)))

    exhaustive <- antitrust:::.blp_multistart_decision(
        convergence_count = 12L, objective_values = rep(1e-12, 12),
        artificial_boundary = TRUE, invalid_evaluations = 3L,
        strategy = "exhaustive"
    )
    expect_false(exhaustive$fallback)
    expect_identical(exhaustive$reasons, character())
})
