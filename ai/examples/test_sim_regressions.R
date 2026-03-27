library(antitrust)

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) {
    stop(msg, call. = FALSE)
  }
}

assert_close <- function(x, y, tol = 1e-6, msg) {
  if (!isTRUE(all.equal(x, y, tolerance = tol, check.attributes = FALSE))) {
    stop(msg, call. = FALSE)
  }
}

prices4 <- c(10, 12, 11, 9)
shares4 <- c(0.25, 0.20, 0.18, 0.17)
owner_pre4 <- c("Firm_A", "Firm_A", "Firm_B", "Firm_C")
owner_post4 <- c("Firm_A", "Firm_A", "Firm_A", "Firm_C")

blp_params <- list(
  alpha = -2,
  sigma = 0.1,
  sigmaNest = 1,
  nDraws = 80
)

set.seed(42)
blp_explicit <- sim(
  prices = prices4,
  shares = shares4,
  demand = "BLP",
  demand.param = blp_params,
  supply = "bertrand",
  ownerPre = owner_pre4,
  ownerPost = owner_post4,
  insideSize = 1000
)

set.seed(42)
blp_alias <- sim(
  prices = prices4,
  shares = shares4,
  demand = "LogitBLP",
  demand.param = blp_params,
  ownerPre = owner_pre4,
  ownerPost = owner_post4,
  insideSize = 1000
)

assert_true(is(blp_alias, "LogitBLP"), "LogitBLP alias did not create a LogitBLP object.")
assert_close(blp_alias@pricePost, blp_explicit@pricePost, msg = "LogitBLP alias does not match explicit BLP Bertrand simulation.")

set.seed(42)
cournot_explicit <- sim(
  prices = prices4,
  shares = shares4,
  demand = "BLP",
  demand.param = blp_params,
  supply = "cournot",
  ownerPre = owner_pre4,
  ownerPost = owner_post4,
  insideSize = 1000
)

set.seed(42)
cournot_alias <- sim(
  prices = prices4,
  shares = shares4,
  demand = "CournotBLP",
  demand.param = blp_params,
  ownerPre = owner_pre4,
  ownerPost = owner_post4,
  insideSize = 1000
)

assert_true(is(cournot_alias, "CournotBLP"), "CournotBLP alias did not create a CournotBLP object.")
assert_close(cournot_alias@pricePost, cournot_explicit@pricePost, msg = "CournotBLP alias does not match explicit BLP Cournot simulation.")

set.seed(7)
blp_default_sigma <- sim(
  prices = prices4,
  shares = shares4,
  demand = "BLP",
  demand.param = list(alpha = -2, sigma = 0.1, nDraws = 60),
  supply = "bertrand",
  ownerPre = owner_pre4,
  ownerPost = owner_post4,
  insideSize = 1000
)

assert_true(all(is.finite(blp_default_sigma@pricePost)), "BLP default sigmaNest path returned non-finite post-merger prices.")
assert_true(isTRUE(all.equal(unname(blp_default_sigma@slopes$sigmaNest), 1)), "BLP default sigmaNest was not set to 1.")

set.seed(11)
blp_stress <- sim(
  prices = prices4,
  shares = shares4,
  demand = "BLP",
  demand.param = list(alpha = -0.05, sigma = 1, sigmaNest = 0.9, nDraws = 100),
  supply = "bertrand",
  ownerPre = owner_pre4,
  ownerPost = owner_post4,
  insideSize = 1000
)

stress_elast <- elast(blp_stress, preMerger = TRUE)
assert_true(all(is.finite(stress_elast)), "BLP stress case still produces non-finite elasticities.")

logitcap <- sim(
  prices = c(10, 11, 9),
  margins = c(0.30, 0.25, 0.20),
  demand = "LogitCap",
  demand.param = list(alpha = -1, mktSize = 150),
  ownerPre = c("A", "B", "C"),
  ownerPost = c("A", "A", "C"),
  capacities = c(30, 25, 20)
)

assert_true(is(logitcap, "LogitCap"), "sim(LogitCap) did not create a LogitCap object.")
assert_close(logitcap@shares, c(30, 25, 20) / 150, msg = "sim(LogitCap) did not map capacities into shares correctly.")
assert_true(isTRUE(all.equal(logitcap@insideSize, 75)), "sim(LogitCap) did not set insideSize to total inside quantity.")
assert_true(all(is.finite(logitcap@pricePost)), "sim(LogitCap) returned non-finite post-merger prices.")

cat("sim()/BLP regression checks passed.\n")
