# Trade MonCom handoff

This handoff accompanies the antitrust/refactor MonCom CES correction and the
trade MonCom regression additions.

The downstream audit found that the legacy `TariffMonComCES` path already
uses a constant `1/gamma` output-market margin. That is consistent with the
corrected antitrust convention: hold the CES aggregate/index fixed and use
the direct own elasticity `-gamma`, rather than the share-adjusted diagonal
`-gamma + (gamma - 1) s_j` returned by the full CES elasticity method.

The focused trade regression now checks these invariants against the current
antitrust/refactor implementation. With common `gamma = 2`, `priceOutside =
1`, revenue shares `(0.30, 0.25, 0.25)`, and a zero baseline tariff, the two
MonCom CES paths agree on baseline prices, revenue shares, quantities, and
marginal costs. The independently evaluated atomistic FOC is

\[
q_j + (p_j-c_j)(-\gamma q_j/p_j) = 0.
\]

It remains zero after a tariff shock as well. A tariff vector `(0.10, 0,
0.20)` leaves the retained baseline `mcPre` unchanged and produces
`mcPost = mcPre/(1-t)` and `pricePost = pricePre/(1-t)`, while the tariff
state remains trade-owned. Repeating the same check with different ownership
labels gives the same prices, as required by atomistic MonCom pricing.

The trade package exposes only the output-market MonCom CES route: its
parameterized `sim()` boundary requires positive `gamma` and has no `output`
orientation argument. There is therefore no trade input-market sign path to
test or silently infer. The antitrust input convention remains documented and
tested in its own package.

The noisy calibration objectives are intentionally different and are tested
as such. Core MonCom CES uses `weighted.mean(1 / margins)` (equal weights in
the focused case), whereas the mature trade calibrator minimizes squared
margin distance and yields `1 / mean(margins)`. They agree under noiseless
common margins but differ for noisy margins; the handoff does not force those
estimates to coincide or alter the mature trade ALM equations.

Trade retains its own tariff/quota adjustments and does not copy antitrust's
demand formulas into the trade implementation. If trade adds BLP MonCom
support later, it should reuse antitrust's stored integration points and
weights rather than implementing a second draw aggregator.

Trade should treat this as a shared economic invariant, not as a request to
change its mature tariff/quota equations. Any target-specific policy state
must remain trade-owned.
