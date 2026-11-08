library(antitrust)
warnings_seen <- character()
checks <- list()
check <- function(name, value, detail = "") {
  checks[[length(checks) + 1L]] <<- data.frame(check = name, pass = isTRUE(value), detail = as.character(detail))
}
attempt <- function(expr) tryCatch(force(expr), error = function(e) e)
close <- function(x, y, tol = 1e-7) isTRUE(all.equal(unname(x), unname(y), tolerance = tol, check.attributes = FALSE))
foc <- function(p, mc, mv, alpha, own) {
  z <- exp(mv + alpha * p)
  s <- z / (1 + sum(z))
  D <- alpha * (diag(s) - tcrossprod(s))
  vapply(seq_along(s), function(j) (s[j] + sum(own[j, ] * (p - mc) * D[, j])) / s[j], numeric(1))
}
withCallingHandlers({
  set.seed(7261106)
  max_foc <- 0
  step_count <- 0L
  for (case in 1:24) {
    fit <- calibrate("logit", "bertrand", prices = c(2, 2.2, 2.5),
      shares = c(.35, .25, .2), margins = c(.4, .35, .3),
      ownerPre = c("A", "B", "C"), insideSize = 100)
    original <- serialize(fit, NULL)
    e <- entrant("E1", meanval = runif(1, -.4, .5), cost = runif(1, .8, 1.2), priceStart = 2)
    changes <- list(
      list(entry = e),
      list(ownership = c("A", "A", "C", "E")),
      list(ownership = c("A", "A", "C", "A")),
      list(quality = c(E1 = runif(1, -.15, .25))),
      list(costs = c(Prod1 = runif(1, -.12, -.01), E1 = runif(1, -.12, -.01))),
      list(exit = "Prod3"))
    # Change shock ordering while retaining two distinct mergers after entry.
    if (case %% 2L == 0L) changes <- changes[c(1, 6, 4, 2, 5, 3)]
    cf <- do.call(counterfactual, changes[[1]])
    for (s in changes[-1]) cf <- do.call(add_step, c(list(cf), s))
    path <- attempt(simulate(fit, cf))
    check(paste0("case_", case, "_solves"), is(path, "CounterfactualPath"), if (inherits(path, "error")) conditionMessage(path) else "")
    if (inherits(path, "error")) next
    mc <- fit@model@mcPre
    mv <- fit@model@slopes$meanval
    alpha <- fit@model@slopes$alpha
    labels <- fit@model@labels
    active <- rep(TRUE, 3)
    own <- diag(3)
    for (i in seq_along(changes)) {
      s <- changes[[i]]
      if (!is.null(s$entry)) {
        mc <- c(mc, e@cost); mv <- c(mv, e@meanval)
        labels <- c(labels, "E1"); active <- c(active, TRUE)
        next_own <- diag(length(labels)); next_own[seq_len(nrow(own)), seq_len(nrow(own))] <- own
        own <- next_own
      }
      if (!is.null(s$ownership)) own <- outer(s$ownership, s$ownership, "==") * 1
      if (!is.null(s$quality)) {
        k <- match(names(s$quality), labels); mv[k] <- mv[k] + log1p(s$quality)
      }
      if (!is.null(s$costs)) {
        k <- match(names(s$costs), labels); mc[k] <- mc[k] * (1 + s$costs)
      }
      if (!is.null(s$exit)) active[match(s$exit, labels)] <- FALSE
      m <- result_at(path, i)
      residual <- max(abs(foc(m@pricePost[active], mc[active], mv[active], alpha, own[active, active, drop = FALSE])))
      max_foc <- max(max_foc, residual)
      step_count <- step_count + 1L
      check(paste(case, i, "costs", sep = "_"), close(m@mcPost, mc))
      check(paste(case, i, "demand", sep = "_"), close(m@slopes$meanval, mv))
      check(paste(case, i, "active", sep = "_"), identical(m@subset, active))
      check(paste(case, i, "ownership", sep = "_"), close(ownerToMatrix(m, FALSE), own))
      check(paste(case, i, "independent_foc", sep = "_"), is.finite(residual) && residual < 1e-6, residual)
      if (i > 1 && is.null(s$entry)) {
        check(paste(case, i, "price_continuity", sep = "_"), close(m@pricePre, result_at(path, i - 1L)@pricePost))
      }
    }
    # Resume through the public path boundary and compare with an uninterrupted solve.
    cf_a <- do.call(counterfactual, changes[[1]])
    for (s in changes[2:3]) cf_a <- do.call(add_step, c(list(cf_a), s))
    cf_b <- do.call(counterfactual, changes[[4]])
    for (s in changes[5:6]) cf_b <- do.call(add_step, c(list(cf_b), s))
    resumed <- attempt(simulate(simulate(fit, cf_a), cf_b))
    check(paste0("case_", case, "_resume"), !inherits(resumed, "error") && close(final_result(resumed)@pricePost, final_result(path)@pricePost))
    check(paste0("case_", case, "_source_immutable"), identical(serialize(fit, NULL), original))
    # Independent cold-start solve from tracked economic primitives.
    cold <- nleqslv::nleqslv(rep(2.5, sum(active)), function(p) foc(p, mc[active], mv[active], alpha, own[active, active, drop = FALSE]))
    check(paste0("case_", case, "_cold_reference"), max(abs(cold$fvec)) < 1e-6 && close(cold$x, final_result(path)@pricePost[active], 1e-6))
  }
  cat("LONG_CHAINS", 24, "SOLVED_STEPS", step_count, "MAX_NORMALIZED_INDEPENDENT_FOC", format(max_foc, digits = 12), "\n")

  fit <- calibrate("logit", "bertrand", prices = c(2, 2.2, 2.5),
    shares = c(.35, .25, .2), margins = c(.4, .35, .3), ownerPre = c("A", "B", "C"))
  full <- counterfactual(entry = entrant("E1", meanval = .2, cost = 1, priceStart = 2))
  full <- add_step(full, ownership = c("A", "A", "C", "E"))
  full <- add_step(full, ownership = c("A", "A", "C", "A"))
  full <- add_step(full, tariff = c(0, 0, .1, .1))
  full <- add_step(full, exit = "Prod3")
  cat("ANTITRUST_EXACT_REQUEST:", conditionMessage(attempt(simulate(fit, full))), "\n")

  policy <- trade::as_trade_fit(fit, tariffPre = 0, cost_basis = "effective", margin_basis = "net_revenue")
  two <- add_step(counterfactual(ownership = c("A", "A", "C")), tariff = c(0, 0, .1))
  cat("TARIFF_GAME_PATH:", conditionMessage(attempt(simulate(policy, two))), "\n")
  cat("TARIFF_GAME_ENTRY:", conditionMessage(attempt(simulate(policy, counterfactual(entry = entrant("E1", meanval = .2, cost = 1, priceStart = 2))))), "\n")
}, warning = function(w) { warnings_seen <<- c(warnings_seen, conditionMessage(w)); invokeRestart("muffleWarning") })
out <- do.call(rbind, checks)
write.csv(out, "/tmp/iovrse-assessment-20261106/long-chain-checks.csv", row.names = FALSE)
writeLines(unique(warnings_seen), "/tmp/iovrse-assessment-20261106/long-chain-warnings.txt")
cat("CHECKS", nrow(out), "PASSED", sum(out$pass), "FAILED", sum(!out$pass), "WARNINGS", length(warnings_seen), "\n")
print(out[!out$pass, ], row.names = FALSE)
