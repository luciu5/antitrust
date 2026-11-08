args <- commandArgs(TRUE)
package <- args[[1]]
root <- "/tmp/iovrse-fixes-20261106"
out <- file.path(root, "final-evidence", package)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(ANTITRUST_TEST_TIER = "extended", TRADE_TEST_TIER = "extended",
  FINANCIAL_EXTENDED_TESTS = "1", ANTITRUSTBAYES_RUN_STAN_TESTS = "false",
  COMPETITIONTOOLBOX_TEST_TIER = "EXTENDED")
cat("SOURCE", file.path(root, "sources", package), "\n")
cat("LIBRARIES", paste(.libPaths(), collapse = ";"), "\n")
for (p in c("antitrust", "coordination", "trade", "vertical", "financial", "competitiontoolbox", "antitrustBayes")) {
  found <- find.package(p, quiet = TRUE)
  cat("DEPENDENCY", p, if (length(found))
    paste(as.character(packageVersion(p)), found) else "unavailable", "\n")
}
r <- testthat::test_local(file.path(root, "sources", package), reporter = "summary", stop_on_failure = FALSE)
d <- as.data.frame(r)
saveRDS(d, file.path(out, "test-results.rds"))
write.csv(d[, intersect(c("file", "test", "nb", "failed", "error", "warning", "skipped", "passed"), names(d))], file.path(out, "test-summary.csv"), row.names = FALSE)
counts <- colSums(d[, c("failed", "error", "warning", "skipped", "passed"), drop = FALSE])
cat("FINAL_COUNTS\n"); print(counts)
jsonlite::write_json(as.list(counts), file.path(out, "test-counts.json"), auto_unbox = TRUE, pretty = TRUE)
quit(status = as.integer(counts[["failed"]] > 0 || counts[["error"]] > 0))
