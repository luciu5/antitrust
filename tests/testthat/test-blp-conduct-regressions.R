# Test-only BLP conduct regressions.
#
# These fixtures keep the demand state fixed while comparing the new BLP
# conduct methods with the homogeneous legacy methods.  The expected
# heterogeneous margins below are assembled from draw-level primitives in
# this file; they do not call calcMargins() to create their oracle.

# Heterogeneous auction and bargaining conduct checks use independent kernels
# and belong to the extended economics gate.
qa_skip_if_not_extended()

blp_conduct_prices <- c(1.75, 2.10, 2.45)
blp_conduct_shares <- c(.30, .25, .25)
blp_conduct_owner <- c("A", "A", "B")
blp_conduct_s0 <- .20
blp_conduct_alpha <- -1.35
blp_conduct_bargpower <- rep(.40, length(blp_conduct_prices))


blp_conduct_delta <- function(alpha, prices = blp_conduct_prices,
                              shares = blp_conduct_shares,
                              s0 = blp_conduct_s0, priceOutside = 0) {
    log(shares / s0) - alpha * (prices - priceOutside)
}


blp_conduct_make_blp <- function(conduct, sigma, nodes, drawWeights,
                                 alpha = blp_conduct_alpha) {
    delta <- blp_conduct_delta(alpha)
    antitrust:::.blp_model(
        conduct = conduct,
        prices = blp_conduct_prices,
        shares = blp_conduct_shares,
        margins = c(.25, .30, .20),
        ownerPre = blp_conduct_owner,
        alphaMean = alpha,
        sigma = sigma,
        meanval = delta,
        draws = nodes,
        drawWeights = drawWeights,
        s0 = blp_conduct_s0,
        output = TRUE,
        priceOutside = 0,
        bargpowerPre = if (identical(conduct, "bargaining")) {
            blp_conduct_bargpower
        } else {
            NULL
        },
        bargpowerPost = if (identical(conduct, "bargaining")) {
            blp_conduct_bargpower
        } else {
            NULL
        }
    )
}


blp_conduct_make_legacy_auction <- function(alpha = blp_conduct_alpha) {
    result <- auction2nd.logit(
        prices = blp_conduct_prices,
        shares = blp_conduct_shares,
        margins = c(.25, .30, .20),
        ownerPre = blp_conduct_owner,
        ownerPost = blp_conduct_owner,
        output = TRUE
    )
    ## Auction2ndLogit stores the price term in meanval, whereas BLP stores
    ## the price coefficient separately.  This is the required translation.
    result@slopes$alpha <- alpha
    result@slopes$meanval <- log(blp_conduct_shares / blp_conduct_s0)
    result@pricePre <- blp_conduct_prices
    result
}


blp_conduct_make_legacy_bargaining <- function(alpha = blp_conduct_alpha) {
    result <- bargaining.logit(
        prices = blp_conduct_prices,
        shares = blp_conduct_shares,
        margins = c(.25, .30, .20),
        ownerPre = blp_conduct_owner,
        ownerPost = blp_conduct_owner,
        bargpowerPre = blp_conduct_bargpower,
        bargpowerPost = blp_conduct_bargpower,
        output = TRUE
    )
    result@slopes$alpha <- alpha
    result@slopes$meanval <- blp_conduct_delta(alpha)
    result@pricePre <- blp_conduct_prices
    result
}


blp_conduct_draw_shares <- function(object) {
    unname(calcShares(object, preMerger = TRUE, revenue = FALSE,
                      aggregate = FALSE))
}


blp_conduct_expected_auction <- function(object) {
    draw_shares <- blp_conduct_draw_shares(object)
    weights <- object@slopes$drawWeights
    alpha <- object@slopes$alphas
    firm_draw_shares <- object@ownerPre %*% draw_shares
    firm_shares <- as.vector(firm_draw_shares %*% weights)

    ## Conditional second-score margin, integrated over consumer types.
    draw_margin <- log(1 - firm_draw_shares) /
        matrix(alpha, nrow = nrow(firm_draw_shares),
               ncol = ncol(firm_draw_shares), byrow = TRUE)
    as.vector((draw_margin %*% weights) / firm_shares)
}


blp_conduct_expected_bargaining <- function(object) {
    draw_shares <- blp_conduct_draw_shares(object)
    weights <- object@slopes$drawWeights
    alpha <- object@slopes$alphas
    barg <- object@bargpowerPre / (1 - object@bargpowerPre)
    output_sign <- -1
    n <- nrow(draw_shares)
    margin_system <- matrix(0, nrow = n, ncol = n)
    right_hand_side <- numeric(n)

    ## This is the integrated bargaining kernel written from its draw-level
    ## definition.  In particular, no aggregate or effective alpha is used.
    for (r in seq_len(ncol(draw_shares))) {
        shares_r <- draw_shares[, r]
        kernel_r <- -object@ownerPre *
            matrix(rep(shares_r, times = n), nrow = n, ncol = n)
        diag(kernel_r) <- diag(object@ownerPre) + diag(kernel_r)
        margin_system <- margin_system + weights[r] * kernel_r

        diversion_odds <- shares_r / (1 - shares_r)
        term_r <- log(1 - shares_r) /
            (-1 * output_sign * alpha[r] *
                 (barg * diversion_odds - log(1 - shares_r)))
        right_hand_side <- right_hand_side + weights[r] *
            diag(object@ownerPre) * term_r
    }

    as.vector(solve(t(margin_system), right_hand_side))
}


test_that("sigma-zero Auction2ndBLP margins match homogeneous Auction2ndLogit", {
    blp <- blp_conduct_make_blp(
        conduct = "auction2nd", sigma = 0,
        nodes = c(-1, 0, 1), drawWeights = c(.2, .3, .5)
    )
    legacy <- blp_conduct_make_legacy_auction()

    expect_equal(unname(calcShares(blp, preMerger = TRUE)),
                 blp_conduct_shares, tolerance = 1e-12)
    expect_equal(unname(calcShares(legacy, preMerger = TRUE)),
                 blp_conduct_shares, tolerance = 1e-12)
    expect_equal(unname(calcMargins(blp, preMerger = TRUE, level = TRUE)),
                 unname(calcMargins(legacy, preMerger = TRUE, level = TRUE)),
                 tolerance = 1e-12)
    expect_equal(unname(calcMargins(blp, preMerger = TRUE, level = FALSE)),
                 unname(calcMargins(legacy, preMerger = TRUE, level = FALSE)),
                 tolerance = 1e-12)
})


test_that("sigma-zero BargainingBLP margins match homogeneous BargainingLogit", {
    blp <- blp_conduct_make_blp(
        conduct = "bargaining", sigma = 0,
        nodes = c(-1, 0, 1), drawWeights = c(.2, .3, .5)
    )
    legacy <- blp_conduct_make_legacy_bargaining()

    expect_equal(unname(calcShares(blp, preMerger = TRUE)),
                 blp_conduct_shares, tolerance = 1e-12)
    expect_equal(unname(calcShares(legacy, preMerger = TRUE)),
                 blp_conduct_shares, tolerance = 1e-12)
    expect_equal(unname(calcMargins(blp, preMerger = TRUE, level = FALSE)),
                 unname(calcMargins(legacy, preMerger = TRUE, level = FALSE)),
                 tolerance = 1e-12)
    expect_equal(unname(calcMargins(blp, preMerger = TRUE, level = TRUE)),
                 unname(calcMargins(legacy, preMerger = TRUE, level = TRUE)),
                 tolerance = 1e-12)
})


test_that("Auction2ndBLP integrates heterogeneous firm winning margins draw by draw", {
    blp <- blp_conduct_make_blp(
        conduct = "auction2nd", sigma = .45,
        nodes = c(-1.5, -.25, .8, 1.75),
        drawWeights = c(.10, .20, .30, .40)
    )
    expected_level <- blp_conduct_expected_auction(blp)
    observed_level <- unname(calcMargins(blp, preMerger = TRUE, level = TRUE))
    observed_proportional <- unname(calcMargins(blp, preMerger = TRUE,
                                                level = FALSE))

    expect_equal(observed_level, expected_level, tolerance = 1e-12)
    expect_equal(observed_proportional,
                 expected_level / blp_conduct_prices, tolerance = 1e-12)
})


test_that("BargainingBLP integrates its bargaining kernel draw by draw", {
    blp <- blp_conduct_make_blp(
        conduct = "bargaining", sigma = .35,
        nodes = c(-1.5, -.25, .8, 1.75),
        drawWeights = c(.10, .20, .30, .40)
    )
    expected_proportional <- blp_conduct_expected_bargaining(blp) /
        blp_conduct_prices
    observed_proportional <- unname(calcMargins(blp, preMerger = TRUE,
                                                level = FALSE))
    observed_level <- unname(calcMargins(blp, preMerger = TRUE, level = TRUE))

    expect_equal(observed_proportional, expected_proportional,
                 tolerance = 1e-12)
    expect_equal(observed_level,
                 expected_proportional * blp_conduct_prices,
                 tolerance = 1e-12)
})
