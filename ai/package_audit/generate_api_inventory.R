# Generate the maintained API coverage matrix.  The matrix is intentionally
# explicit: an API is either tied to a formal test or carries a documented
# exclusion explaining why a fixture is not appropriate.

if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
} else {
    library(antitrust)
}

ns <- asNamespace("antitrust")
exports <- sort(getNamespaceExports("antitrust"))
user_exports <- exports[!grepl("^\\.__[CT]__", exports)]
classes <- sort(getClasses(where = ns))
classes <- classes[classes != "ANY"]
generics <- sort(getGenerics(where = ns))

formal_exports <- setdiff(user_exports, c("antitrust_shiny", "plot"))
export_rows <- data.frame(
    kind = "export", name = user_exports, signature = "",
    status = ifelse(user_exports %in% formal_exports, "formal-test", "documented-exclusion"),
    evidence = ifelse(
        user_exports %in% formal_exports,
        "tests/testthat/test-model-combinations.R; tests/testthat/test-public-output-methods.R",
        "ai/package_audit/audit_report.md"
    ),
    exclusion_reason = ifelse(
        user_exports == "antitrust_shiny", "Interactive optional UI; exercised by R CMD check examples when available.",
        ifelse(user_exports == "plot", "Graphics device output is environment-dependent; documented examples remain the smoke check.", "")
    ), stringsAsFactors = FALSE
)

direct_classes <- c(
    "AIDS", "CES", "CESALM", "CESCournot", "CESCournotALM", "CESNests",
    "Cournot", "CournotBLP", "Linear", "LogLin", "Logit", "LogitALM",
    "LogitBLP", "LogitCap", "LogitCapALM", "LogitCournot", "LogitCournotALM",
    "LogitNests", "LogitNestsALM", "PCAIDS", "PCAIDSNests", "Stackelberg",
    "Auction2ndCap", "Auction2ndCES", "Auction2ndCESALM", "Auction2ndLogit",
    "Auction2ndLogitALM", "Auction2ndLogitNests", "BargainingCES",
    "BargainingCESALM", "BargainingLogit", "BargainingLogitALM",
    "Bargaining2ndCES", "Bargaining2ndLogit", "VertBargBertLogit",
    "VertBarg2ndLogit", "VertBargBertLogitNests", "VertBarg2ndLogitNests"
)
class_rows <- data.frame(
    kind = "S4 class", name = classes, signature = "",
    status = ifelse(classes %in% direct_classes, "formal-test", "documented-exclusion"),
    evidence = ifelse(classes %in% direct_classes,
                      "tests/testthat/test-model-combinations.R; tests/testthat/test-public-output-methods.R",
                      "ai/package_audit/audit_report.md"),
    exclusion_reason = ifelse(classes %in% direct_classes, "",
                              "Abstract/container class or specialized plant-level fixture; inherited behavior is covered by child classes."),
    stringsAsFactors = FALSE
)

tested_generics <- c(
    "CV", "calcBuyerExpectedCost", "calcBuyerValuation", "calcDiagnostics",
    "calcdMC", "calcExpectedLowestCost", "calcExpectedPrice", "calcMC",
    "calcMargins", "calcMarginsAG", "calcOptimalReserve", "calcPriceDelta",
    "calcPrices", "calcPricesAG", "calcProducerSurplus", "calcQuantities",
    "calcRevenues", "calcSellerCostParms", "calcShares", "calcSlopes",
    "calcVC", "cdfG", "cmcr", "diversion", "elast", "getParms", "hhi",
    "getNestsParms", "ownerToMatrix", "ownerToVec", "show", "summary", "upp"
)
method_rows <- do.call(rbind, lapply(generics, function(generic) {
    found <- findMethods(generic, where = ns)
    signatures <- names(found)
    if (!length(signatures)) return(NULL)
    formal <- generic %in% tested_generics
    data.frame(
        kind = "S4 method", name = generic, signature = signatures,
        status = if (formal) "formal-test" else "documented-exclusion",
        evidence = if (formal) "tests/testthat/test-public-output-methods.R" else "ai/package_audit/audit_report.md",
        exclusion_reason = if (formal) "" else "Specialized method is not identified by the synthetic core fixture; retain documented example as the review path.",
        stringsAsFactors = FALSE
    )
}))

inventory <- rbind(export_rows, class_rows, method_rows)
inventory <- inventory[order(inventory$kind, inventory$name, inventory$signature), ]
dir.create("ai/package_audit", showWarnings = FALSE, recursive = TRUE)
write.csv(inventory, "ai/package_audit/api_coverage.csv", row.names = FALSE)

cat("exports:", nrow(export_rows), "\n")
cat("classes:", nrow(class_rows), "\n")
cat("generics:", length(generics), "\n")
cat("methods:", nrow(method_rows), "\n")
cat("formal rows:", sum(inventory$status == "formal-test"), "\n")
cat("documented exclusions:", sum(inventory$status == "documented-exclusion"), "\n")
