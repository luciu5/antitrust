#!/usr/bin/env Rscript

# Reproducibility runner for the litigation-oriented synthetic QA suite.
# Usage: Rscript ai/package_audit/run_qa.R [repository-root] [output-directory]

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else normalizePath(getwd(), mustWork = TRUE)
root <- normalizePath(root, mustWork = TRUE)
output_dir <- if (length(args) >= 2) args[[2]] else {
    file.path(root, "ai", "package_audit", "latest")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(root)

required <- c("testthat", "pkgload")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing QA packages: ", paste(missing, collapse = ", "))

git_value <- function(...) {
    value <- tryCatch(system2("git", c(...), stdout = TRUE, stderr = FALSE),
                      error = function(e) "unavailable")
    if (!length(value)) "unavailable" else paste(value, collapse = "\n")
}

fixture_files <- sort(list.files(file.path(root, "tests", "testthat", "fixtures"),
                                 full.names = TRUE, recursive = TRUE))
fixture_files <- fixture_files[file.info(fixture_files)$isdir %in% FALSE]
fixture_hash <- if (length(fixture_files)) tools::md5sum(fixture_files) else character()

build_dir <- tempfile("antitrust-qa-build-")
dir.create(build_dir)
build_log <- file.path(output_dir, "build.log")
# Build from a sanitized staging copy.  Developer checkouts often contain
# large R-devel installations or check directories; those are not package
# inputs and can otherwise make R CMD build appear to hang while scanning.
staging_dir <- file.path(build_dir, "source")
staging_archive <- file.path(build_dir, "source.tar")
stage_status <- system2("tar", c(
    "-cf", staging_archive,
    "--exclude=.git", "--exclude=.R-devel-build", "--exclude=.Rlib",
    "--exclude=.R-devel-install", "--exclude=.R-devel-src",
    "--exclude=..Rcheck", "--exclude=antitrust.Rcheck",
    "--exclude=antitrust_*.tar.gz", "--exclude=.codex",
    "-C", root, "."
), stdout = FALSE, stderr = FALSE)
dir.create(staging_dir)
if (identical(stage_status, 0L)) {
    stage_status <- system2("tar", c("-xf", staging_archive, "-C", staging_dir),
                            stdout = FALSE, stderr = FALSE)
}
oldwd <- setwd(build_dir)
build_timeout <- Sys.getenv("ANTITRUST_QA_BUILD_TIMEOUT", "30")
build_status <- tryCatch(
    system2("timeout",
            c(paste0(build_timeout, "s"), file.path(R.home("bin"), "R"),
              "CMD", "build", "--no-build-vignettes", "--no-manual",
              "--no-resave-data", "--md5", staging_dir),
            stdout = build_log, stderr = build_log),
    error = function(e) 127L
)
tarballs <- list.files(build_dir, pattern = "^antitrust_.*\\.tar\\.gz$", full.names = TRUE)
tarball <- if (length(tarballs)) tarballs[[1]] else {
    # A pre-existing package archive is still useful evidence when a local
    # filesystem prevents R CMD build from completing; build_status records
    # that distinction explicitly.
    existing <- list.files(root, pattern = "^antitrust_.*\\.tar\\.gz$", full.names = TRUE)
    if (length(existing)) existing[[which.max(file.info(existing)$mtime)]] else NA_character_
}
tarball_hash <- if (!is.na(tarball)) unname(tools::md5sum(tarball)) else NA_character_
setwd(oldwd)

warning_log <- file.path(output_dir, "warnings.csv")
if (file.exists(warning_log)) unlink(warning_log)
Sys.setenv(ANTITRUST_QA_WARNING_LOG = warning_log)

test_output <- character()
test_status <- 0L
test_results <- NULL
test_output <- tryCatch(
    capture.output({
        test_results <- testthat::test_local(root, reporter = "summary",
                                             load_package = "source")
    }),
    error = function(e) {
        test_status <<- 1L
        c(conditionMessage(e), capture.output(print(e)))
    }
)
writeLines(test_output, file.path(output_dir, "test-output.txt"))

expectation_rows <- list()
if (!is.null(test_results)) {
    for (entry in test_results) {
        if (!length(entry$results)) next
        for (condition in entry$results) {
            expectation_rows[[length(expectation_rows) + 1L]] <- data.frame(
                file = entry$file, test = entry$test,
                expectation_class = class(condition)[1],
                message = conditionMessage(condition),
                stringsAsFactors = FALSE
            )
            if (inherits(condition, c("expectation_failure", "expectation_error"))) {
                test_status <- 1L
            }
        }
    }
}
expectations <- if (length(expectation_rows)) do.call(rbind, expectation_rows) else
    data.frame(file=character(), test=character(), expectation_class=character(),
               message=character(), stringsAsFactors=FALSE)
write.csv(expectations, file.path(output_dir, "expectations.csv"), row.names = FALSE)

if (file.exists(warning_log)) {
    recorded_warnings <- tryCatch(read.csv(warning_log, stringsAsFactors = FALSE),
                                  error = function(e) data.frame())
} else {
    recorded_warnings <- data.frame()
}
if (!nrow(recorded_warnings)) {
    recorded_warnings <- data.frame(context=character(), class=character(),
                                     message=character(), stringsAsFactors=FALSE)
}

# Warnings emitted by tests that do not use qa_capture are represented by
# testthat expectation_warning results.  Merge them with helper-recorded rows.
test_warnings <- expectations[grepl("warning", expectations$expectation_class,
                                     ignore.case = TRUE), c("file", "test", "message")]
if (nrow(test_warnings)) {
    patterns <- read.csv(file.path(root, "ai", "package_audit", "warning_dispositions.csv"),
                         stringsAsFactors = FALSE)
    test_warnings$context <- paste(test_warnings$file, test_warnings$test, sep = "::")
    test_warnings$class <- vapply(test_warnings$message, function(msg) {
        hit <- which(vapply(patterns$pattern, function(p) {
            grepl(p, msg, perl = TRUE, ignore.case = TRUE)
        }, logical(1)))
        if (length(hit)) patterns$class[[hit[[1]]]] else "unclassified-review-required"
    }, character(1))
    test_warnings <- test_warnings[, c("context", "class", "message")]
    recorded_warnings <- unique(rbind(recorded_warnings, test_warnings))
}
write.csv(recorded_warnings, warning_log, row.names = FALSE)

# Keep session capture deterministic and avoid platform-specific timezone
# discovery (which can call `timedatectl` on Linux CI images).
old_tz <- Sys.getenv("TZ", unset = NA_character_)
Sys.setenv(TZ = "UTC")
session_lines <- capture.output(sessionInfo())
if (is.na(old_tz)) Sys.unsetenv("TZ") else Sys.setenv(TZ = old_tz)
writeLines(session_lines, file.path(output_dir, "session-info.txt"))
soft_version <- extSoftVersion()
soft_value <- function(name) if (name %in% names(soft_version)) soft_version[[name]] else NA_character_
metadata <- list(
    audit_timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    commit = git_value("rev-parse", "HEAD"),
    git_status = git_value("status", "--short"),
    package_tarball = if (is.na(tarball)) NA_character_ else basename(tarball),
    package_tarball_md5 = tarball_hash,
    build_status = build_status,
    R = R.version.string,
    platform = R.version$platform,
    dependency_versions = vapply(c("antitrust", "testthat", "pkgload", "BB", "numDeriv", "nleqslv", "SQUAREM"),
                                 function(pkg) if (requireNamespace(pkg, quietly=TRUE)) as.character(packageVersion(pkg)) else NA_character_,
                                 character(1)),
    BLAS = soft_value("BLAS"),
    LAPACK = soft_value("LAPACK"),
    RNGkind = RNGkind(),
    fixture_md5 = fixture_hash,
    test_status = test_status,
    test_count = nrow(expectations),
    warning_count = nrow(recorded_warnings),
    unclassified_warning_count = sum(recorded_warnings$class == "unclassified-review-required"),
    session_info = session_lines
)
saveRDS(metadata, file.path(output_dir, "metadata.rds"))

report <- c(
    "# antitrust synthetic QA run", "",
    paste0("Generated: ", metadata$audit_timestamp_utc),
    paste0("Commit: `", metadata$commit, "`"),
    paste0("R: `", metadata$R, "` on `", metadata$platform, "`"),
    paste0("BLAS: `", metadata$BLAS, "`; LAPACK: `", metadata$LAPACK, "`"),
    paste0("Package tarball: `", metadata$package_tarball, "` (MD5 `", metadata$package_tarball_md5, "`)"),
    paste0("Formal expectations: ", metadata$test_count,
           "; status: ", if (test_status == 0L) "PASS" else "FAIL"),
    paste0("Warnings captured: ", metadata$warning_count,
           "; unclassified: ", metadata$unclassified_warning_count), "",
    "## Warning dispositions", "",
    "Warnings are recorded rather than fatal under the QA policy. Every warning is retained in `warnings.csv`; an unclassified warning is a review item.", "",
    if (nrow(recorded_warnings)) paste0("- `", recorded_warnings$class, "` — ", recorded_warnings$context, ": ", recorded_warnings$message) else "- None", "",
    "## Reproduction", "",
    "Run `Rscript ai/package_audit/run_qa.R .` from a clean checkout. The runner builds a temporary tarball, runs the formal testthat suite, records fixture hashes and session details, and writes machine-readable outputs to the selected output directory."
)
writeLines(report, file.path(output_dir, "report.md"))
cat(paste(test_output, collapse = "\n"), "\n", sep = "")
cat("QA outputs: ", normalizePath(output_dir), "\n", sep = "")
if (test_status != 0L) quit(status = test_status)
