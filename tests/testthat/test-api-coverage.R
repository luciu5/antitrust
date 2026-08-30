test_that("API coverage matrix is complete and reviewable", {
    repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
    matrix_path <- file.path(repo_root, "ai", "package_audit", "api_coverage.csv")
    if (!file.exists(matrix_path)) {
        testthat::skip("The maintained API inventory is a repository-level QA artifact and is excluded from the package tarball")
    }
    coverage <- read.csv(matrix_path, stringsAsFactors = FALSE)
    testthat::expect_true(all(coverage$status %in% c("formal-test", "documented-exclusion")))
    key <- paste(coverage$kind, coverage$name, coverage$signature, sep = "::")
    testthat::expect_equal(anyDuplicated(key), 0L)

    ns <- asNamespace("antitrust")
    exports <- getNamespaceExports("antitrust")
    exports <- exports[!grepl("^\\.__[CT]__", exports)]
    listed_exports <- coverage$name[coverage$kind == "export"]
    testthat::expect_setequal(exports, listed_exports)

    classes <- getClasses(where = ns)
    classes <- setdiff(classes, "ANY")
    listed_classes <- coverage$name[coverage$kind == "S4 class"]
    testthat::expect_setequal(classes, listed_classes)

    exclusions <- coverage[coverage$status == "documented-exclusion", ]
    testthat::expect_true(all(nzchar(exclusions$exclusion_reason)))
    evidence <- unique(coverage$evidence)
    evidence_paths <- unique(unlist(strsplit(evidence, "; ", fixed = TRUE)))
    evidence_paths <- evidence_paths[grepl("^tests/|^ai/", evidence_paths)]
    testthat::expect_true(all(file.exists(file.path(repo_root, evidence_paths))))
})
