# iovrse vignette implementation and practitioner review

Completed locally on 2026-11-09. This record covers documentation implementation,
deterministic example execution and independent review. It does not certify new
posterior estimates, a public website release or manual browser acceptance.

The user confirmed use of the maintained master/main branches, preserving local
changes. Three Sol implementation agents wrote the package guides; four fresh
Sol reviewers approached them as R-savvy research assistants supporting
antitrust consulting teams or competition agencies. Their actionable findings
were corrected and checked. No commits, pushes or publication were performed.

## Delivered reading paths

| Package | Guides, in suggested order |
|---|---|
| antitrust | `Workflow`, `Models`, `Architecture` |
| antitrustBayes | `pnb-merger`, `posterior-workflow`, `supply-and-demand` |
| trade | `Reference`, `quotas`, `cournot-tariffs` |
| coordination | `price-leadership`, `grim-trigger`, `role-games` |
| financial | `debt-default-merger`, `bank-liquidity-merger` |
| vertical | `linked-supply-chain-merger`, `vertical-tariffs-quotas` |
| competitiontoolbox | `gui-to-r-merger` |
| iopolicy | `reproducible-market-designs` (experimental developer infrastructure) |

There are 17 guides across the seven public installer packages and one separate
iopolicy developer guide. Existing antitrust `Workflow` and trade `Reference`
identifiers remain usable by the GUI. README/help maps, vignette metadata and
pkgdown navigation accompany the sources. Guides target short, focused tasks,
generally below ten printed pages. Actual print pagination remains unmeasured.

Ordinary lightweight teaching code executes. Bayesian fitting, bridge evidence
and interactive app commands are explicitly optional and unevaluated. The PNB
starter uses illustrative supplied parameters and requires no Stan installation;
its numbers are not presented as estimated paper results.

## Substantive corrections

- antitrust introduces the small merger before synthetic-market design, keeps
  calibration discrepancies visible, and distinguishes lifecycle operations.
- Bayes separates the supplied-state starter from optional posterior analysis;
  explains persistent quality plus recovered period demand, excluded measurement
  noise, fixed draw plans and reconstructed baselines; uses draw-wise Paasche
  harm with explicit units, invalid counts and the negative deposit-market HMT
  direction. HMT summaries also check their own returned/requested counts.
- Trade separates tariff accounting, binding quotas and quantity-setting games.
  It distinguishes physical unit cost from total cost, reports account totals,
  retains original owner labels, qualifies effective-cost equivalence, explains
  poor illustrative Cournot margin fit, and derives/checks quota shadows only
  for the example's single-product constrained firms.
- Coordination presents the PLE identification boundary honestly: a zero
  additional wedge is distinct from coordination-minus-Bertrand price changes.
  Missing timing cannot establish full repeated-game sustainability. The other
  guides cover firm portfolios and strategic versus atomistic role assumptions.
- Financial explains payment/rate units, default and survival comparisons,
  footprint-weighted bank baselines, changed-bank row alignment, tiny pooling
  effects and near-boundary deposit outcomes. Vertical explains integrated
  internal transfers, tier margin denominators, inside-quantity HHI, tariff
  accounts, quota units and KKT diagnostics.
- The GUI guide documents copying code, stored research tabs and the difference
  between inside-only Details shares and full-market script shares, with an
  executable conversion. iopolicy documents neutral designs, seeds, retries and
  rejection records without claiming an economic equilibrium.

Fifteen numbered practitioner findings across the four reviews were resolved,
along with minor source-navigation and provenance improvements. Package
equations, likelihoods, solvers and runtime economic APIs were not changed.
R-source changes for these guides are documentation comments. The installer
minimum for antitrustBayes is now 0.3.0.9002 in both identical installer copies.

## Build and distribution

The installer still disables vignette/manual building. Documentation tools are
optional build dependencies; they do not become runtime Imports. iopolicy has
not been added to the public installer.

competitiontoolbox now contains a manual documentation-review workflow and
`.ci/build_documentation.R`, described in `.ci/DOCUMENTATION.md`. Setup installs
the current source cohort once in an isolated library. The build helper records
versions, loaded paths, Git dirty state and actual source/installed file hashes;
it never installs packages. It changes into each package directory before
pkgdown initializes, because antitrust's Rd title macro reads DESCRIPTION there.

Generated sites use `.scratch/pkgdown`, preserving vertical's authored `docs/`.
The workflow has read-only repository permissions and no deployment step. Local
rendered guides and normal vignette bundles exist. Public sites were not
published; maintained source fallbacks accompany cross-guide HTML links.

## Verification actually performed

- All 18 guides rendered with their lightweight code. Small public workflows
  were independently exercised by the practitioner reviewers. Optional fitting
  and evidence code was inspected/parsed, not executed.
- All eight packages produced normal source bundles with vignettes. Installed
  vignette indexes and HTML discovery were checked. Scoped package checks used
  `R CMD check --no-manual --no-tests --no-examples`; these are documentation
  checks, not full package regression certifications. Seven packages reported
  `Status: OK`. Financial reported one warning because one existing cubature
  installation has no help database. The specific cross-reference check passes
  with another existing cubature installation containing its complete help.
- Existing Bayes price-index/decision tests passed 35 expectations. Selected
  financial accounting/default/independent FOC and vertical bargaining/welfare/
  policy tests passed. Existing vertical fixture tests emitted 29 documented
  calibration/convergence warnings, distinct from clean new guide examples.
- Installer regression/parity checks, GUI documentation-access, app clean-process
  replay and FAST model-registry tests passed. The parent's GUI run emitted 31
  existing fixture warnings and skipped one explicitly EXTENDED-only test;
  this was not a manual browser run or an extended-suite certification.
- Independent merger portfolio FOCs, draw-wise Paasche accounting, quota shadows,
  tariff accounting and seeded design/rejection/RNG checks passed. The HMT
  availability guard was exercised on a tiny supplied fit, including shortened
  and missing returns.
- pkgdown article builds/navigation passed for antitrust, antitrustBayes, trade
  and coordination. Bootstrap-5 site builds for the other packages stop at a
  missing headroom CDN asset: DNS/download failed, including an approved network
  retry using an isolated temporary cache. Their ordinary vignette HTML renders
  and bundles passed. The complete new Actions job has not run remotely.
- R/YAML parsing, article-source targets and byte-identical installer mirrors
  pass. Documentation-owned diffs pass whitespace checks. Whole-worktree checks
  still report pre-existing whitespace in trade's launcher and the toolbox's
  trade input file; those concurrent files were preserved.

All installations were explicit local setup into temporary libraries. Existing
paper libraries were read-only. No HMC, bridge estimation, posterior fitting,
cluster job, GUI browser launch or public deployment ran.

## Evidence and remaining manual acceptance

Implementation and review evidence is retained locally:

- `/tmp/iovrse-vignette-core-implementation.md`
- `/tmp/iovrse-vignette-trade-coordination-implementation.md`
- `/tmp/iovrse-vignette-financial-vertical-implementation.md`
- `/tmp/iovrse-practitioner-review-core-bayes.md`
- `/tmp/iovrse-practitioner-review-trade-coordination.md`
- `/tmp/iovrse-practitioner-review-financial-vertical.md`
- `/tmp/iovrse-practitioner-review-gui-design.md`
- `/tmp/iovrse-build-trade-coordination-20261109/report.md`
- `/tmp/iovrse-doc-final-checks-20261109/report.md`
- `/tmp/iovrse-core-final-validation/`
- `/tmp/iovrse-doc-verification-20261109/`

Sibling-repository edits were applied through approved, SHA-guarded manifests
under `/tmp/iovrse-doc-implementation`, so unrelated concurrent edits were not
replaced. Temporary reports/builds should be archived if needed beyond this
machine/session; generated HTML is not committed as authoritative source.

Before a public documentation release, run the prepared online maintainer job,
review readable-size print pagination and HTML wrapping/MathJax in a browser,
and complete the toolbox's manual browser acceptance in `.ci/READINESS.md`.
No browser executable was available here. Pages publication remains a separate
authorized release action. Keep the source links until public rendered sites
exist. A positive PLE leadership fixture remains separate economic example work;
the current guide intentionally teaches the verified identification boundary.
