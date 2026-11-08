# Permanent economic/numerical regression: a unit change could send Cournot's
# finite-difference calibration to a zero slope/cost bound and abort, even when
# observed quantities, margins and elasticity satisfy independent FOCs exactly.
test_that("Cournot calibration respects its finite domain across market units", {
  owner_pre <- c("a", "b", "c", "c")
  owner_post <- c("a", "a", "c", "c")
  cases <- data.frame(quantity_scale = c(1, 100, 100), price_scale = c(1, .1, 1))
  for (i in seq_len(nrow(cases))) {
    qs <- cases$quantity_scale[i]
    ps <- cases$price_scale[i]
    observed_q <- c(40, 30, 20, 10) * qs
    observed_p <- 10 * ps
    # Linear inverse demand with elasticity -1 implies slope -p/Q.
    # Each plant's observed margin is its owner's aggregate quantity share.
    margins <- c(.4, .3, .3, .3)
    model <- cournot(
      prices = observed_p, quantities = matrix(observed_q, ncol = 1),
      margins = matrix(margins, ncol = 1), mktElast = -1,
      ownerPre = owner_pre, ownerPost = owner_post
    )
    expect_equal(unname(model@intercepts / ps), 20, tolerance = 1e-4)
    expect_equal(unname(model@slopes * qs / ps), -.1, tolerance = 1e-4)
    observed_mc <- vapply(seq_along(observed_q), function(j) {
      model@mcfunPre[[j]](observed_q[j])
    }, numeric(1))
    expect_equal(observed_mc / observed_p, 1 - margins, tolerance = 1e-4)
    for (pre in c(TRUE, FALSE)) {
      q <- as.numeric(if (pre) model@quantityPre else model@quantityPost)
      owner <- if (pre) owner_pre else owner_post
      cost <- if (pre) model@mcfunPre else model@mcfunPost
      p <- model@intercepts[1] + model@slopes[1] * sum(q)
      mc <- vapply(seq_along(q), function(j) cost[[j]](q[j]), numeric(1))
      own_q <- as.numeric(outer(owner, owner, "==") %*% q)
      expect_true(all(is.finite(q) & q > 0))
      expect_lt(max(abs((p + model@slopes[1] * own_q - mc) / observed_p)), 2e-4)
    }
  }
})
