# antitrust

`antitrust` calibrates structural demand and cost models from observed market
data and simulates mergers and other changes to competition. It also accepts
known structural parameters and provides screening, welfare, and synthetic-market
tools.

## Installation

Run this in a fresh R session. The shared iovrse installer installs `antitrust`
and its required R dependencies while preserving compatible cohort packages.

```r
source(url("https://raw.githubusercontent.com/luciu5/competitiontoolbox/master/tools/install_iovrse.R"))
install_iovrse(packages = "antitrust")
```

Use `install_iovrse()` to install the complete iovrse cohort. Model use does not
require the release-certification workflow or a documentation build.

## Calibrate and simulate a merger

This three-product Logit market has an outside share of 20%. Prices and gross
margins identify demand and marginal costs under Bertrand competition.

```r
fit <- antitrust::calibrate(
  "logit", "bertrand",
  prices = c(2, 2.2, 2.5), shares = c(.35, .25, .20),
  margins = c(.40, .35, .30), ownerPre = c("A", "B", "C"),
  insideSize = 100
)
merger <- antitrust::counterfactual(ownership = c("A", "A", "C"))
result <- antitrust::simulate(fit, merger)
data.frame(before = result@pricePre, after = result@pricePost)

# Apply a merger and a 10% cost reduction for product 2 together.
efficient <- antitrust::simulate(fit, antitrust::counterfactual(
  ownership = c("A", "A", "C"), costs = c(0, -.10, 0)
))
efficient@pricePost
```

Each call starts from the fitted market. To solve changes sequentially, use
`add_step()` and inspect the returned path with `result_at()` or `final_result()`.

## Choose an implemented model

```r
antitrust::supportedModels()[c("demand", "conduct", "variant",
                              "calibrate", "specify", "simulate")]
```

The registry lists complete implemented combinations. Demand families include
Logit, CES, nested and capacity variants, Linear/LogLin, AIDS/PCAIDS, and BLP.
Conduct includes Bertrand, Cournot, Stackelberg, monopolistic competition,
auctions, and bargaining where implemented. Model-specific inputs and supported
counterfactual fields are described in the [models guide](vignettes/Models.Rmd).
Vertical bargaining belongs to the sibling `vertical` package; specialized
coordinated and firm-level role games belong to `coordination`.

| Task | Public interface |
| --- | --- |
| Recover a model from observed market inputs | `calibrate()` |
| Supply known structural parameters | `specify()` |
| Solve a merger, costs, exit, entry, or another supported scenario | `counterfactual()` and `simulate()` |
| Recalibrate after revising observed inputs | `update()` |
| Translate to another supported model while retaining compatible primitives | `respecify()` and `model_transition()` |
| Construct and realize a synthetic market | `synthetic_market()`, `fake_market()`, `realize_market()` |
| Inspect model results | `calcShares()`, `calcMargins()`, `elast()`, `CV()`, `summary()` |

Read the executable [workflow guide](vignettes/Workflow.Rmd) first, then the
[economic models guide](vignettes/Models.Rmd) to choose assumptions. Package
authors can continue to [fitted objects, dispatch, and extension boundaries](vignettes/Architecture.Rmd).
Use `help("antitrust", package = "antitrust")` for the public overview. Source guides are readable without rendering the vignettes.

## Contributing and verification

The repository's routine branch is `master`. Contributor tests and the separate
extended numerical tiers are described in [the test audit](TEST_SUITE_AUDIT.md).
Exact release revisions and the canonical CI environment serve certification;
ordinary installation uses the shared installer above.
