# Trade MonCom handoff

This handoff accompanies the antitrust/refactor MonCom CES correction. No
trade files were edited.

The read-only downstream audit found that the legacy `TariffMonComCES` path
already uses a constant `1/gamma` output-market margin. That is consistent
with the corrected antitrust convention: hold the CES aggregate/index fixed
and use the direct own elasticity `-gamma`, rather than the share-adjusted
diagonal `-gamma + (gamma - 1) s_j` returned by the full CES elasticity
method.

When trade is next tested against the installed antitrust/refactor head:

1. regression-test tariff MonCom CES with zero and nonzero tariffs;
2. compare baseline shares, prices, marginal costs, and own-product FOC
   residuals against the corrected antitrust flat CES path when policy is
   economically null;
3. retain trade-owned tariff/quota adjustments and do not copy antitrust's
   demand formulas into trade;
4. explicitly test the input-market sign convention if a tariff/input CES
   path is exposed;
5. if trade adds BLP MonCom support, reuse antitrust's stored integration
   points and weights rather than implementing a second draw aggregator.

Trade should treat this as a shared economic invariant, not as a request to
change its mature tariff/quota equations. Any target-specific policy state
must remain trade-owned.

