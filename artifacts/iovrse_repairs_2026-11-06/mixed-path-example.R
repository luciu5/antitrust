library(antitrust)
source <- specify("logit", "bertrand",
  prices = c(2, 2.2, 2.4, 2.6), ownerPre = c("A", "A", "B", "C"),
  labels = c("a", "b", "c", "d"), insideSize = 100,
  parameters = list(alpha = -1.5, meanval = c(3, 2.5, 2, 1.5)))
fit <- trade::as_trade_fit(source, entrant_tariff = .12)
cf <- counterfactual(entry = entrant("e", meanval = 1, cost = .9, priceStart = 2))
cf <- add_step(cf, ownership = c("A", "A", "A", "C", "e"))
cf <- add_step(cf, ownership = c("A", "A", "A", "C", "A"))
cf <- add_step(cf, tariff = c(a = 0, b = .15, c = .15, d = .1, e = .2))
cf <- add_step(cf, exit = "d")
path <- suppressWarnings(simulate(fit, cf))
stopifnot(is(path, "CounterfactualPath"), length(path@results) == 5L,
  result_at(path, 1)@tariff_state$tariffPost[["e"]] == .12,
  final_result(path)@tariff_state$entrant_tariff == .12)
rows <- do.call(rbind, lapply(seq_along(path@results), function(i) {
  accounts <- trade::tariff_accounts(result_at(path, i))
  stopifnot(all(is.finite(accounts$q)), all(accounts$q >= 0),
    max(abs(accounts$producer_profit + accounts$government_revenue -
      (accounts$gross_revenue - accounts$physical_cost * accounts$q))) < 1e-8)
  cbind(step = i, shock = c("entry", "merger", "merger", "tariff", "exit")[[i]], accounts)
}))
print(rows[, c("step", "shock", "product", "tau", "price", "physical_cost", "q")], row.names = FALSE)
dir.create("/tmp/iovrse-fixes-20261106/final-evidence/example", recursive = TRUE, showWarnings = FALSE)
write.csv(rows, "/tmp/iovrse-fixes-20261106/final-evidence/example/mixed-path.csv", row.names = FALSE)
