# antitrust test tiers

The default `fast` tier protects public contracts, independent economics,
canonical behavioral oracles, and historical regressions.  It runs on every
push and pull request.

The additive `extended` tier runs the fast suite plus exhaustive legacy parity,
solver comparisons, BLP recovery, and compound counterfactual paths.  It runs
nightly and through manual workflow dispatch.

Run locally:

```sh
Rscript ai/package_audit/run_qa.R . qa-output
ANTITRUST_TEST_TIER=extended Rscript ai/package_audit/run_qa.R . qa-output-extended
```

| Contract | Fast canonical owner | Extended evidence |
| --- | --- | --- |
| Model specification and registry | `test-model-registry.R` | Full legacy model matrix |
| Fit API and legacy wrapper | `test-model-api-contracts.R` | `test-legacy-model-parity.R` |
| Master behavioral oracle | `test-refactor-oracle.R` | Legacy parity by model family |
| Counterfactual construction | `test-counterfactual.R` | `test-sequential-counterfactuals.R` |
| Demand translations | `test-local-demand-translations.R` | Family parity where applicable |
| BLP integration | `test-blp-integration-regressions.R` | Recovery and conduct regressions |
| Solver alternatives | Audit regressions | `test-solver-parity.R` |
| Public methods and economics | Independent/oracle/invariant tests | Specialized model parity |

Before deleting or merging an assertion, update this table or the relevant
test comment to identify the surviving canonical owner.  Historical regression
tests retain their named context even when their implementation is consolidated.
