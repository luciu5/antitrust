# Test-only BLP conduct regressions.
#
# These fixtures keep the demand state fixed while comparing the new BLP
# conduct methods with the homogeneous legacy methods.  The expected
# heterogeneous margins below are assembled from independently aggregated
# demand primitives in this file; they do not call calcMargins() to create
# their oracle.

# Sigma-zero and zero-bargaining limits are deterministic fast checks.
# Heterogeneous draw-level kernels are representative extended checks.

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
    aggregate_shares <- as.vector(draw_shares %*% weights)
    derivative <- matrix(0, nrow = n, ncol = n)
    buyer_surplus <- numeric(n)
    for (r in seq_len(ncol(draw_shares))) {
        shares_r <- draw_shares[, r]
        derivative <- derivative + weights[r] * alpha[r] *
            (diag(shares_r) - tcrossprod(shares_r))
        buyer_surplus <- buyer_surplus + weights[r] *
            log1p(-shares_r) / alpha[r]
    }

    aggregate_elast <- derivative * outer(1 / aggregate_shares,
                                          object@pricePre)
    revenue <- object@pricePre * aggregate_shares
    margin_system <- t(
        diag(1 / revenue) %*%
            (t(aggregate_elast * object@ownerPre) %*%
                 diag(aggregate_shares))
    )
    own_normalized <- diag(derivative) / aggregate_shares
    right_hand_side <- own_normalized /
        (output_sign * (own_normalized - barg * aggregate_shares /
                        buyer_surplus))
    right_hand_side <- diag(object@ownerPre) * right_hand_side

    ## This is the aggregate Nash system: aggregate the demand Jacobian and
    ## buyer surplus first, then solve the ownership-adjusted FOCs.  Averaging
    ## inverses of draw-level systems is not an aggregate equilibrium.
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
    qa_skip_unless_tier("extended")
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


test_that("BargainingBLP solves the aggregate heterogeneous bargaining FOC", {
    qa_skip_unless_tier("extended")
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


test_that("zero buyer bargaining power satisfies the aggregate Bertrand FOC", {
    prices <- c(2, 2.5, 3)
    shares <- c(.2, .3, .25)
    alpha <- -1.5
    sigma <- .6
    nodes <- c(-1, 0, 1)
    weights <- rep(1 / 3, 3)
    s0 <- 1 - sum(shares)
    delta <- antitrust:::.blp_contract(
        prices, shares, alpha, sigma, nodes, weights, s0
    )$delta
    bertrand <- suppressWarnings(antitrust:::.blp_model(
        conduct = "bertrand", prices = prices, shares = shares,
        margins = rep(.2, 3), ownerPre = 1:3, alphaMean = alpha,
        sigma = sigma, meanval = delta, draws = nodes,
        drawWeights = weights, s0 = s0
    ))
    bargaining <- suppressWarnings(antitrust:::.blp_model(
        conduct = "bargaining", prices = prices, shares = shares,
        margins = rep(.2, 3), ownerPre = 1:3, alphaMean = alpha,
        sigma = sigma, meanval = delta, draws = nodes,
        drawWeights = weights, s0 = s0, bargpowerPre = rep(0, 3)
    ))

    expected <- unname(calcMargins(bertrand, level = TRUE))
    observed <- unname(calcMargins(bargaining, level = TRUE))
    ## The zero-power boundary is the documented Bertrand limit.  The
    ## aggregate bargaining system has a distinct Nash objective for positive
    ## buyer power; this regression guards the economically decisive boundary
    ## without treating the legacy Bertrand margin method as a demand FOC
    ## oracle.
    expect_equal(observed, expected, tolerance = 1e-12)

    fractional_owner <- matrix(c(
        .7, .1, 0,
        .1, .7, 0,
        0, 0, .8
    ), nrow = 3, byrow = TRUE)
    bertrand_fractional <- suppressWarnings(antitrust:::.blp_model(
        conduct = "bertrand", prices = prices, shares = shares,
        margins = rep(.2, 3), ownerPre = fractional_owner,
        alphaMean = alpha, sigma = sigma, meanval = delta,
        draws = nodes, drawWeights = weights, s0 = s0
    ))
    bargaining_fractional <- suppressWarnings(antitrust:::.blp_model(
        conduct = "bargaining", prices = prices, shares = shares,
        margins = rep(.2, 3), ownerPre = fractional_owner,
        alphaMean = alpha, sigma = sigma, meanval = delta,
        draws = nodes, drawWeights = weights, s0 = s0,
        bargpowerPre = rep(0, 3)
    ))
    expect_equal(
        unname(calcMargins(bargaining_fractional, level = TRUE)),
        unname(calcMargins(bertrand_fractional, level = TRUE)),
        tolerance = 1e-12
    )
})


test_that("sigma-zero BLP auction allocation uses values net of cost", {
    prices <- c(2, 2.5, 3)
    shares <- c(.2, .3, .1)
    s0 <- .4
    meanval <- log(shares / s0) - (-1.5) * prices
    fit <- suppressWarnings(specify(
        demand = "BLP", conduct = "auction2nd", prices = prices,
        shares = shares, ownerPre = 1:3,
        parameters = list(alpha = -1.5, sigma = 0, meanval = meanval,
                          integration = "auto"),
        bargpowerPre = rep(.5, 3)
    ))
    simulated <- suppressWarnings(simulate(fit, ownerPost = c(1, 1, 3)))
    expect_equal(unname(calcShares(simulated, preMerger = FALSE)), shares,
                 tolerance = 1e-12)
    expect_equal(unname(simulated@pricePost),
                 c(2.180384, 2.631585, 3), tolerance = 1e-6)
})
