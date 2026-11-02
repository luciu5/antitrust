# Price Leadership Regime Architecture

## Overview

The `calcMargins` method for `PriceLeadership` class now supports 4 distinct pricing regimes via the `regime` parameter. This provides flexible markup calculation that can be invoked either through explicit arguments or by setting ownership slots appropriately.

## The Four Regimes

### 1. Bertrand (`regime = "bertrand"`)
**Standard Bertrand competition with no coordination**

- Uses traditional ownership matrix (`control = FALSE`)
- Coalition firms do NOT coordinate
- Each firm independently maximizes profit
- Equivalent to standard Logit model

**Usage:**
```r
margins <- calcMargins(object, preMerger = TRUE, regime = "bertrand")
prices <- calcPrices(object, preMerger = TRUE, regime = "bertrand")
```

**When to use:**
- Baseline competitive scenario
- Counterfactual: "what if no coordination?"
- Grim trigger punishment phase

---

### 2. Deviation (`regime = "deviation"`)
**One firm deviates while others maintain PLE prices**

- Uses traditional ownership matrix (`control = FALSE`)
- Deviating firm optimizes independently
- Other coalition firms' prices held fixed at PLE levels
- Primarily used in IC constraint checking

**Usage:**
```r
# Typically called internally by calcDeviationProfit
margins <- calcMargins(object, preMerger = TRUE, regime = "deviation")
```

**When to use:**
- Calculating deviation profits for IC constraints
- Testing stability of coordination
- Identifying which firm has strongest incentive to deviate

---

### 3. Coordination (`regime = "coordination"`)
**Unconstrained coalition coordination - DEFAULT**

- Uses coalition control matrix (`control = TRUE`)
- Coalition products fully coordinate (act as single entity)
- NO supermarkup added yet (pure Bertrand with coordination)
- IC constraints are SLACK (non-binding)

**Usage:**
```r
# Default behavior
margins <- calcMargins(object, preMerger = TRUE)  # regime defaults to "coordination"
prices <- calcPrices(object, preMerger = TRUE)
```

**When to use:**
- Standard PLE calculation
- When supermarkup = 0 (no need for extra markup beyond coordination)
- Initial calibration before testing IC constraints

---

### 4. Constrained (`regime = "constrained"`)
**Coordination with IC binding - includes supermarkup**

- Uses coalition control matrix (`control = TRUE`)
- Coalition coordinates AND adds supermarkup to maintain IC
- At least one IC constraint BINDS
- Supermarkup determined by most constrained firm

**Usage:**
```r
margins <- calcMargins(object, preMerger = TRUE, regime = "constrained")
prices <- calcPrices(object, preMerger = TRUE, regime = "constrained")
```

**When to use:**
- When calibrated supermarkup > 0
- At least one coalition firm's IC constraint binds
- Final equilibrium with timing parameter calibrated from binding IC

---

## Implementation Details

### calcMargins Method

The method selects the appropriate ownership matrix based on regime:

```r
ownerMatrix <- switch(regime,
  "bertrand"     = ownerToMatrix(object, preMerger, control = FALSE),
  "deviation"    = ownerToMatrix(object, preMerger, control = FALSE),
  "coordination" = ownerToMatrix(object, preMerger, control = TRUE),
  "constrained"  = ownerToMatrix(object, preMerger, control = TRUE)
)
```

Then temporarily swaps ownership and calls parent (Logit) method:
```r
object@ownerPre <- ownerMatrix
result <- callNextMethod(object, preMerger, level)
object@ownerPre <- ownerOriginal
```

### calcPrices Method

Prices are calculated via FOCs where:

**Coalition products:**
```r
margin = bertrand_margin(regime) + supermarkup
```

**Fringe products:**
```r
margin = bertrand_margin("bertrand")  # Always standard Bertrand
```

The regime parameter determines which coordinated margin is used as the base.

---

## Comparison with Supreme Implementation

The architecture aligns with the Supreme (MSW 2023) implementation pattern:

**Supreme pattern:**
```r
# Step 1: Bertrand with coalition coordination
bert_rate <- bertrand_rate(alpha, s0, mrp, meanval, data_id, owner_id)

# Step 2: Add supermarkup to coalition
eval_rate <- bert_rate - supm * as.integer(data_id %in% coalition_id)

# Step 3: Solve fringe best response
ple_rate <- fix_rate(fixed_rates = eval_rate, ...)
```

**Antitrust equivalent:**
```r
# Step 1 & 2: Combined in calcMargins + FOC
margins <- calcMargins(object, preMerger, regime = "coordination")
foc <- margins - predMargin - supermarkup  # Coalition
foc <- margins - predMargin               # Fringe

# Step 3: Solved simultaneously via nleqslv/BBsolve
prices <- calcPrices(object, preMerger)
```

---

## Key Design Principles

1. **Separation of Concerns:**
   - `calcMargins`: Calculates Bertrand-level coordinated margins
   - `calcPrices`: Adds supermarkup and solves equilibrium

2. **Flexible Invocation:**
   - Via explicit `regime` argument
   - Via ownership slot settings (future enhancement)
   - Defaults match standard PLE behavior

3. **OOP Inheritance:**
   - Leverages Logit parent class FOC logic
   - Ownership matrix swap enables clean regime switching
   - No code duplication

4. **Supreme Alignment:**
   - Sequential logic (Bertrand → supermarkup → fringe response)
   - Clear distinction between coordination and supermarkup
   - Consistent with empirical implementation

---

## Usage Examples

### Example 1: Standard PLE Calibration
```r
result <- ple(
  prices = prices,
  shares = shares,
  margins = margins,
  coalitionPre = c(1, 2, 3),
  ownerPre = owners
)

# Calculates regime="coordination" by default
prices_ple <- calcPrices(result)
```

### Example 2: IC Constrained Equilibrium
```r
# If calibration finds binding IC and positive supermarkup
if(result@supermarkup > 0){
  prices_constrained <- calcPrices(result, regime = "constrained")
}
```

### Example 3: Bertrand Counterfactual
```r
# What if no coordination?
prices_bertrand <- calcPrices(result, regime = "bertrand")
welfare_gain <- calcConsumerSurplus(result) - calcConsumerSurplus_bertrand
```

### Example 4: Deviation Profit (Internal)
```r
# Called automatically in IC checking
dev_profit <- calcDeviationProfit(result, deviatingFirm = "Firm1")
# Uses regime = "deviation" internally
```

---

## Migration Notes

### Backward Compatibility

All existing code continues to work:
```r
# Old code still works - defaults to "coordination"
margins <- calcMargins(object, preMerger = TRUE)
prices <- calcPrices(object, preMerger = TRUE)
```

### New Capabilities

Code can now explicitly specify regime:
```r
# New capabilities
margins_bert <- calcMargins(object, regime = "bertrand")
margins_ple <- calcMargins(object, regime = "coordination")
margins_ic <- calcMargins(object, regime = "constrained")
```

---

## Future Enhancements

1. **Ownership Slot Setting:**
   Allow regime to be inferred from ownership slots:
   ```r
   object@regime <- "bertrand"
   prices <- calcPrices(object)  # Uses bertrand automatically
   ```

2. **Regime-Specific Methods:**
   ```r
   calcBertrandPrices(object)      # Wrapper for regime="bertrand"
   calcCoordinationPrices(object)  # Wrapper for regime="coordination"
   ```

3. **Diagnostic Tools:**
   ```r
   compareRegimes(object)  # Compare all 4 regimes side-by-side
   plotRegimes(object)     # Visualize price differences
   ```

---

## References

- Mansley, Miller, Sheu & Weinberg (2023), "A price leadership model for coordinated effects in merger analysis"
- Supreme implementation: `price_leader_functions.R` (MSW co-author code)
