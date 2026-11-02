userlib <- Sys.getenv("R_LIBS_USER")
if (nzchar(userlib)) .libPaths(c(userlib, .libPaths()))

library(antitrust)

set.seed(123)

prices <- c(
  2.46, 2.21, 2.4, 2.29, 2.05, 2.45, 1.13, 2.39, 2.66, 2.73,
  1.83, 1.84, 2.76, 2.46, 1.1, 2.53, 2.27, 2.31, 2.72, 1.96, 2.26,
  2.63, 2.07, 2.31, 1.33, 2.67, 2.68, 2.59, 2.62, 1.85, 2.63, 1.51,
  1.8, 2.66, 1.4, 2.39, 2.35, 1.82
)

shares <- c(
  0.12499, 0.11462, 0.07748, 0.05231, 0.05168, 0.02918, 0.01622,
  0.01245, 0.00992, 0.00578, 0.00461, 0.00417, 0.00415, 0.00275,
  0.00271, 0.00214, 0.00206, 0.00184, 0.00176, 0.00155, 0.00151,
  0.00146, 0.00146, 0.00103, 0.00101, 0.00095, 8e-04, 0.00066,
  0.00062, 0.00061, 0.00057, 0.00053, 0.00053, 5e-04, 0.00049,
  0.00048, 0.00041, 0.00039
)

ownerPre <- diag(38)
ownerPost <- diag(38)
ownerPost[2, 3] <- 1
ownerPost[3, 2] <- 1

insideSize <- 4773473000
blp_params <- list(alpha = 1.0, sigma = 0.3, sigmaNest = 0.9, nDraws = 200)

is_party <- rowSums(abs(ownerPost - ownerPre)) > 0

run_sim <- function(supply, mc_delta = rep(0, length(prices))) {
  set.seed(123)
  sim(
    prices, shares,
    supply = supply, demand = "BLP",
    demand.param = blp_params,
    ownerPre = ownerPre, ownerPost = ownerPost,
    insideSize = insideSize,
    mcDelta = mc_delta,
    control.equ = list(maxit = 3000, tol = 1e-10)
  )
}

logitblp <- run_sim("bertrand")
cournotblp <- run_sim("cournot")

logit_cmcr_cost <- cmcr(logitblp, rel = "cost")
logit_cmcr_price <- cmcr(logitblp, rel = "price")
logit_cmcr_levels <- cmcr(logitblp, rel = "cost", levels = TRUE)

logit_diversion <- diversion(logitblp, preMerger = TRUE, revenue = FALSE)
diag(logit_diversion) <- -1

logit_direct_cost_full <- cmcr.bertrand(
  logitblp@pricePre,
  calcMargins(logitblp, preMerger = TRUE),
  logit_diversion,
  logitblp@ownerPre,
  ownerPost = logitblp@ownerPost,
  output = logitblp@output,
  rel = "cost",
  labels = logitblp@labels
)

logit_direct_price_full <- cmcr.bertrand(
  logitblp@pricePre,
  calcMargins(logitblp, preMerger = TRUE),
  logit_diversion,
  logitblp@ownerPre,
  ownerPost = logitblp@ownerPost,
  output = logitblp@output,
  rel = "price",
  labels = logitblp@labels
)

logit_mc_delta <- rep(0, length(prices))
logit_mc_delta[is_party] <- if (logitblp@output) -logit_cmcr_cost / 100 else logit_cmcr_cost / 100
logit_offset <- run_sim("bertrand", logit_mc_delta)

cournot_cmcr_cost <- cmcr(cournotblp, rel = "cost")
cournot_cmcr_price <- cmcr(cournotblp, rel = "price")

cournot_margin_pre <- calcMargins(cournotblp, preMerger = TRUE)
cournot_post_owner <- cournotblp
cournot_post_owner@ownerPre <- cournot_post_owner@ownerPost
cournot_margin_post <- calcMargins(cournot_post_owner, preMerger = TRUE)
cournot_manual_price <- (cournot_margin_post - cournot_margin_pre)[is_party] * 100
cournot_manual_cost <- cournot_margin_post - cournot_margin_pre
if (cournotblp@output) {
  cournot_manual_cost <- cournot_manual_cost / (1 - cournot_margin_pre)
} else {
  cournot_manual_cost <- cournot_manual_cost / (1 + cournot_margin_pre)
}
cournot_manual_cost <- cournot_manual_cost[is_party] * 100

cournot_mc_delta <- rep(0, length(prices))
cournot_mc_delta[is_party] <- if (cournotblp@output) -cournot_cmcr_cost / 100 else cournot_cmcr_cost / 100
cournot_offset <- run_sim("cournot", cournot_mc_delta)

correct_level_cost <- (logit_cmcr_cost / 100) * logitblp@mcPre[is_party]

cat("antitrust version:", as.character(packageVersion("antitrust")), "\n")
cat("PNB merging products:", paste(which(is_party), collapse = ", "), "\n\n")

cat("LogitBLP cmcr() cost %:\n")
print(logit_cmcr_cost)
cat("LogitBLP cmcr.bertrand() cost %, party subset:\n")
print(logit_direct_cost_full[is_party])
cat("Max abs diff, cost %:", max(abs(logit_cmcr_cost - logit_direct_cost_full[is_party])), "\n\n")

cat("LogitBLP cmcr() price %:\n")
print(logit_cmcr_price)
cat("LogitBLP cmcr.bertrand() price %, party subset:\n")
print(logit_direct_price_full[is_party])
cat("Max abs diff, price %:", max(abs(logit_cmcr_price - logit_direct_price_full[is_party])), "\n\n")

cat("LogitBLP levels=TRUE from cmcr():\n")
print(logit_cmcr_levels)
cat("Correct cost-level CMCR using mcPre[is_party]:\n")
print(correct_level_cost)
cat("Max abs diff, levels:", max(abs(logit_cmcr_levels - correct_level_cost)), "\n\n")

cat("LogitBLP base party post-pre prices:\n")
print(logitblp@pricePost[is_party] - logitblp@pricePre[is_party])
cat("LogitBLP with signed mcDelta=cmcr party post-pre prices:\n")
print(logit_offset@pricePost[is_party] - logit_offset@pricePre[is_party])
cat("Max abs post-pre price diff with CMCR:", max(abs(logit_offset@pricePost[is_party] - logit_offset@pricePre[is_party])), "\n\n")

cat("CournotBLP cmcr() cost %:\n")
print(cournot_cmcr_cost)
cat("CournotBLP manual post-owner FOC cost %:\n")
print(cournot_manual_cost)
cat("Max abs diff, cost %:", max(abs(cournot_cmcr_cost - cournot_manual_cost)), "\n\n")

cat("CournotBLP cmcr() party price %:\n")
print(cournot_cmcr_price)
cat("CournotBLP manual post-owner FOC price %:\n")
print(cournot_manual_price)
cat("Max abs diff, party price %:", max(abs(cournot_cmcr_price - cournot_manual_price)), "\n\n")

cat("CournotBLP base party post-pre prices:\n")
print(cournotblp@pricePost[is_party] - cournotblp@pricePre[is_party])
cat("CournotBLP with signed mcDelta=cmcr party post-pre prices:\n")
print(cournot_offset@pricePost[is_party] - cournot_offset@pricePre[is_party])
cat("Max abs post-pre price diff with CMCR:", max(abs(cournot_offset@pricePost[is_party] - cournot_offset@pricePre[is_party])), "\n")
