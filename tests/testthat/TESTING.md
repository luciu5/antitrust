# antitrust test tiers

The suite has three additive tiers. `fast` is the default for local work and
push/pull-request CI; it protects public contracts, independent economics,
analytic FOCs, lifecycle invariants, registry coverage, and load-order
behavior. `extended` adds representative parameter recovery, integration,
solver, sequential, and cross-package checks. `nightly` adds broad optimizer
recovery grids, stochastic and multi-rule replications, exhaustive model and
legacy-parity combinations, and the slow CRAN-like package check.

The practical targets are two minutes for fast, five minutes for extended, and
an unrestricted scheduled window for nightly. Tier gates are block-level where
a file contains both cheap contract checks and expensive recovery cases.

Run locally:

```sh
Rscript ai/package_audit/run_qa.R . qa-output-fast
ANTITRUST_TEST_TIER=extended Rscript ai/package_audit/run_qa.R . qa-output-extended
ANTITRUST_TEST_TIER=nightly Rscript ai/package_audit/run_qa.R . qa-output-nightly
```

| Contract | Fast canonical owner | Extended evidence | Nightly evidence |
| --- | --- | --- | --- |
| Model specification and registry | `test-model-registry.R` | Representative model construction | `test-model-combinations.R` full matrix |
| Fit API and legacy wrapper | `test-model-api-contracts.R` | Representative solver checks | `test-legacy-model-parity.R` |
| Master behavioral oracle | `test-refactor-oracle.R` | Independent economics/oracles | Legacy parity by model family |
| Counterfactual construction | `test-counterfactual.R` | Representative path continuation | `test-sequential-counterfactuals.R` family matrix |
| Demand translations | `test-local-demand-translations.R` | BLP/price-leadership integration | Broad legacy translations |
| BLP integration | BLP derivative and dimension contracts | Integration reuse and representative recovery | All-conduct and multi-rule recovery |
| Solver alternatives | Independent FOC and KKT tests | Representative nleqslv/AG parity | Solver regression fixtures and grids |
| Public methods and economics | Independent/oracle/invariant tests | Specialized model parity | Exhaustive model combinations |

Before deleting or merging an assertion, update this table or the relevant
test comment to identify the surviving canonical owner. Historical regression
tests retain their named context even when their implementation is routed to a
slower tier.

## CI source provenance

The antitrust workflow records the resolved `git rev-parse HEAD` and the
`ANTITRUST_TEST_TIER` value in each job summary and QA artifact. Push and pull
request jobs run `fast`; scheduled jobs run `nightly`; manual dispatch accepts
`extended` or `nightly`. The `refactor` branch name is a mutable development
default, so a recorded commit SHA is required when the source revision itself
must be reproducible. Trade's CI accepts that antitrust SHA through its manual
`antitrust_ref` input and records the resolved core and trade SHAs together.
