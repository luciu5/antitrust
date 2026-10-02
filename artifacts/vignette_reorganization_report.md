# Vignette reorganization audit

## Inventory before editing

`vignettes/Reference.Rmd` has 2,395 lines, 141,909 bytes, 82 headings (35 at
levels 1–2), and several repeated explanations of simulation. The categories
below are editorial destinations, not a mechanical split. “Models” means the
new economics vignette; “Workflow” and “Architecture” name the other two.

| Old section or section group | Destination | Decision |
|---|---|---|
| Front matter and introduction | Workflow / delete | Retain the short statement of purpose and input sensitivity; remove stale address, contact, Shiny promotion, disclaimers and commented prose. |
| Separating calibration from simulation | Workflow | Rewrite for the current five operations and executable code; the old chunks are `eval=FALSE`. |
| Bertrand game and mathematical model | Models | Retain the ownership-adjusted FOC and cost-recovery logic; shorten repeated derivations. |
| Exogenous capacity constraints | Models | Retain the economic distinction for `LogitCap`; move argument details to function help. |
| Calibrating demand and costs | Models | Keep the identification logic and the useful diversion definitions; correct blanket claims about data requirements that vary by family. |
| Linear, Log-linear, LA-AIDS, nested LA-AIDS/PCAIDS, Logit, nested Logit, CES, nested CES | Models | Condense to economic forms, identified parameters, normalization and constraints. The old constructor recipes belong in function help. |
| Unknown outside share / ALM variants | Models | Keep the identification limitation and variant distinction; remove repeated constructor instructions. |
| Marginal costs | Models | Keep FOC recovery and assumptions; drop duplicate algebra. |
| Bertrand simulation, summaries, plots, efficiencies, exits, welfare | Workflow / function documentation | Use one coherent Logit market for the current API; leave method arguments and plotting options to help pages. |
| Market definition and gotchas | Function documentation / delete | `HypoMonTest()` and `diversionHypoMon()` help cover the API; old guidelines quotations and repeated caveats do not belong in these vignettes. |
| Known demand parameters | Workflow | Explain `specify()` in the lifecycle without adding a second market example. |
| Cournot game, calibration, and first-mover advantage | Models | Retain the distinct quantity and Stackelberg assumptions; delete duplicate simulation walkthrough. |
| Capacity-constrained second-price auction | Models | Keep distribution, reserve and moment-identification distinctions; leave full distribution menu to help. |
| Differentiated second-score Logit/CES auction | Models | Retain bid/valuation and normalization distinctions; remove repeated summaries and examples. |
| Nash bargaining and bargaining game | Models | Keep disagreement payoff and bargaining-power logic; remove duplicated headings and “experimental” labels that do not describe the registry. |
| Vertical supply | Architecture / historical artifact | State package boundary briefly. Migration explanation remains in existing `artifacts/refactor_architecture.md`. |
| CMCR, generalized pricing pressure, HHI | Function documentation / Architecture | Mention as public diagnostics; detailed formulas and argument lists remain in `man/`. |
| Coordinated effects and grim trigger | Models / function documentation | Keep a short economic distinction; the long repeated-game derivation and use instructions belong in dedicated help or historical Git history. |
| Under the Hood, Getting Help, extending Bertrand/auctions | Architecture / developer documentation | Replace stale constructor-step advice with the current fit/registry/dispatch contract; retain `showClass()` and `showMethods()` as discovery tools. |
| Appendix function-input table, method table, CV formula table | Function documentation / delete | Stale static catalog (omits BLP and new lifecycle). Replace support list with a table generated from `supportedModels()`; method details remain in help. |
| Class diagram | Delete | The old 1000px diagram omits `StructuralFit`, `AntitrustFit`, counterfactual/path, BLP, bargaining and newer subclasses. A compact current relationship diagram in Architecture replaces it. |

No old section needs a NEWS entry: it does not describe a new user-visible code
release. Migration reports, test architecture, parity results and rejected
designs remain in `artifacts/` rather than user-facing vignettes. Git history
retains the old manual and obsolete prose.

## Design before implementation

1. **Workflow** — one three-product Logit–Bertrand market. Fit observed
   prices/shares/margins/ownership; inspect; solve merger, efficiency, exit,
   entry and quality; combine and sequence changes; compare `specify()`,
   `update()` and registered `respecify()` transitions. State which scenario
   fields are family-specific, including capacity and bargaining/leadership.
2. **Models** — general identification and FOCs, demand families grouped by
   substitution structure, economically distinct conduct, then a compact
   support table produced directly by `supportedModels()`. Equations are kept
   only where they identify a parameter or equilibrium condition.
3. **Architecture** — fit/spec/result/path relationships, the important S4
   inheritance branches, public generics, registry dispatch and stable sibling
   extension hooks. Refer maintainers to `ARCHITECTURE.md` and help for slots
   and implementation details.

The three audiences have separate entry points. Models supplies theory to
Workflow by cross-reference; Architecture explains dispatch without repeating
the theory or the tutorial. The package's model-specific constructors remain
the behavioral source of truth.

## Implementation and verification

The three source vignettes total 26,214 bytes (roughly 82% less than the old
source). The deliberate split is:

| Vignette | Source bytes | Rendered words, including code/table | Scope |
|---|---:|---:|---|
| `Workflow.Rmd` | 8,328 | about 2,680 | One Logit–Bertrand market, repeated/simultaneous/sequential scenarios, lifecycle semantics and registered transitions. |
| `Models.Rmd` | 11,176 | about 2,560 | Identification, demand and conduct economics, registry-generated complete-model table. |
| `Architecture.Rmd` | 6,710 | about 1,450 | Fit/result/path classes, selected inheritance, public dispatch and sibling-package boundary. |

`Reference.Rmd`, its committed generated `.R` and
`.html`, the stale `ClassDiagram.png`, and `Thumbs.db` were deleted. The
bibliography was retained: all nine keys cited in `Models.Rmd` resolve in
`antitrustbib.bib`. `.Rbuildignore` no longer excludes all vignettes, and
`DESCRIPTION` declares `VignetteBuilder: knitr` so package builds actually
render them. A package-level help page is generated from
`R/AntitrustPackage.R`.

### Public-example audit

I reviewed all **19** previously generated `man/` example blocks against
their roxygen sources and checked the major constructors' execution through
the package check. The final documentation has **28** blocks: nine new
canonical lifecycle examples, 11 edited existing blocks, and eight existing
blocks unchanged. Nine of the 19 original blocks are retained specifically
as direct-constructor compatibility/convenience examples (eight shortened,
one unchanged). No complete old block was deleted; the obsolete 38-product
BLP `\dontrun{}` comparison embedded in the `sim()` block was deleted.

| Existing example block | Classification | Edit/decision |
|---|---|---|
| `Antitrust-Class`, `BertrandOther-Classes`, `BertrandRUM-Classes`, `Cournot-classes` | KEEP AS IS | Small `showClass()` discovery examples. |
| `Auction-Classes`, `AuctionCap-Methods`, `Ownership-methods` | KEEP AS IS | Small S4 method/class discovery examples. |
| `HHI-Functions` | KEEP AS IS | Partial-ownership example is useful; corrected its inaccurate firm/product comment. |
| `Auction2ndCap-Functions` | KEEP AS LEGACY/COMPATIBILITY EXAMPLE | Specialized auction constructor illustrates its distinct cost distribution; left intact. |
| `AIDS-Functions`, `Linear-Functions`, `Logit-Functions`, `CES-Functions` | KEEP AS LEGACY/COMPATIBILITY EXAMPLE | Shortened repeated output/method tours; AIDS retains three genuinely distinct calibration paths. |
| `Cournot-Functions` | KEEP AS LEGACY/COMPATIBILITY EXAMPLE | Replaced random 138-line derivation with deterministic Cournot and Stackelberg calls from existing test fixtures. |
| `Auction2ndLogit-Functions`, `BargainingLogit-Functions` | KEEP AS LEGACY/COMPATIBILITY EXAMPLE | Retained one local specialized example each; removed repetitive class tours and repeated unit experiments. |
| `Sim-Functions` | KEEP AS LEGACY/COMPATIBILITY EXAMPLE | Retained a small supplied-parameter direct call and pointed readers to `specify()`/`simulate()`; deleted the expensive `\dontrun{}` comparison. |
| `CMCRBertrand-Functions`, `CMCRCournot-Functions` | SIMPLIFY | Kept local screening calculations; removed redundant nested parameter-grid loops. |

The nine added examples are package-level orientation plus `model_spec()`,
`supportedModels()`, `calibrate()`, `specify()`, `simulate()`,
`counterfactual()`, `update()` and `respecify()`. They use small executable
markets. No examples use `\dontrun{}` or `\donttest{}` now. BLP's old large
unrun illustration is not replaced by a costly example in routine checks;
the BLP implementation and integration rules are covered by existing tiered
tests and the economics vignette.

### Documentation/API findings

* `update.AntitrustFit()` rejects changes to demand, conduct and variant.
  Its old help claimed that `...` accepted model-specification replacements;
  the help and workflow now state the actual same-specification contract.
* There are no exported `prices()` or `margins()` generics. The current public
  S4 methods are `calcPrices()` and `calcMargins()`, so Architecture names them.
* `ClassDiagram.png` predates the fit and path classes and omits current model
  families; it was removed and replaced with a selective text map.
* Package-level help previously lacked a generated page and recommended
  legacy constructors as the starting path; it now introduces the lifecycle.
* The `auction2nd.logit.alm` return-class sentence misnamed
  `auction2nd.logit`; the Logit normalization sentence said one instead of
  zero; the CES constructor described its curvature as a price coefficient.
  These documentation statements were corrected without changing economics.
* The first simplified CMCR example used a nonexistent `isOne` named
  argument. The example check caught it; changing it to `ownerPre` made the
  example execute. No production API was changed.

### QA and editorial review

`devtools::document()` regenerated help and namespace metadata. Roxygen
reported its existing deprecated `@docType "package"` style and skipped
the manually maintained `setRetention.Rd`; it completed documentation
generation. All three vignettes rendered independently. A clean staged
`R CMD build --no-manual --no-resave-data` created all three `inst/doc`
vignettes, and `R CMD check --no-manual` passed with **Status: OK**. That
check ran examples, tests, Rd cross-references, vignette checks and vignette
rebuilding. Its first pass exposed the CMCR example error above; the passing
result is from the second, corrected pass. The canonical fast QA runner also
passed 1,256 formal expectations; it captured 84 recorded
warnings, including two clipping warnings awaiting QA disposition.

`R CMD check --as-cran --no-manual` could not complete its remote CRAN
incoming-package step because DNS for CRAN/Bioconductor failed, including
when retried with network approval. The offline package check completed.
Rendered HTML was checked for raw equations, unresolved bibliography keys,
cross-vignette links, the registry table, and excessive repetition. The
workflow's expected construction warnings were suppressed in vignette
rendering so they do not obscure the example.

Editorial pass: the merger reader reaches a solved result in the first two
sections; the calibration reader gets identification and a generated support
table without class details; the extension reader gets public hooks and the
package boundary without a constructor recipe. No unsupported general
respecification or arbitrary demand-conduct composition is claimed.
The sibling-package descriptions were checked against their source and
package metadata, but their cross-package test suites were outside this
documentation check.
