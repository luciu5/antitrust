# Performance Analysis: PriceLeadershipFunctions.R

## Executive Summary

**Overall Assessment**: Code is reasonably efficient but has several significant bottlenecks, particularly in nested loops and repeated calculations during binary search iterations.

**Estimated Impact**: 
- **High Impact** improvements could reduce runtime by 50-70%
- **Medium Impact** improvements could reduce runtime by 20-40%
- **Low Impact** improvements would reduce runtime by 5-15%

---

## Critical Bottlenecks (HIGH IMPACT)

### 1. **calcSupermarkup Binary Search - Repeated Full Equilibrium Solves**
**Location**: Lines 838-950  
**Issue**: Each binary search iteration recalculates Bertrand equilibrium from scratch

```r
# CURRENT: Inside binary search loop (called ~15-20 times)
profitsBertrand <- calcProducerSurplus(object, preMerger = TRUE)  # Line 853
# ↓ This calls calcPrices which solves full nonlinear system

# Bertrand profits (line 855-856)
object@pricePre <- calcPrices(object, preMerger = TRUE, regime = "bertrand")
profitsBertrand <- calcProducerSurplus(object, preMerger = TRUE)
```

**Problem**: Bertrand profits are INDEPENDENT of supermarkup m, yet calculated in every iteration.

**Impact**: ~50% of binary search time wasted

**SOLUTION**:
```r
# Calculate ONCE before binary search loop
priceBertrand_fixed <- calcPrices(object, preMerger = preMerger, regime = "bertrand")
object@pricePre <- priceBertrand_fixed
profitsBertrand_fixed <- calcProducerSurplus(object, preMerger = preMerger)
profitsBertrand_byFirm <- sapply(coalitionFirms, function(firm){
  sum(profitsBertrand_fixed[ownerVec == firm])
})

# Then in loop, REUSE profitsBertrand_byFirm
while(m_high - m_low > tol){
  # ... calculate PL and Deviation profits only
  # Use pre-calculated profitsBertrand_byFirm
}
```

**Estimated Speedup**: 2-3x faster binary search (40-50% overall reduction)

---

### 2. **Deviation Profit Calculation in Binary Search**
**Location**: Lines 874-910 (inside binary search loop)  
**Issue**: Each coalition firm requires full equilibrium solve (calcPrices + calcProducerSurplus)

```r
profitsDeviation_byFirm <- sapply(coalitionFirms, function(firm){
  # ... setup ...
  optimalPrices <- calcPrices(object, ..., subset = subset, regime = "deviation")
  # ↑ Full nonlinear system solve
  # ...
})
```

**Problem**: For n coalition firms, each binary search iteration does:
- n full deviation equilibrium solves (calcPrices)
- n profit calculations (calcProducerSurplus)
- With 4 firms × 15 iterations = 60 equilibrium solves

**Impact**: Dominates binary search runtime

**POTENTIAL SOLUTIONS**:

**Option A - Cache Gradient Information** (Complex but powerful):
- Deviation profits change smoothly with m
- Could use previous iteration's solution as warm start
- Implement Newton-style updates

**Option B - Parallel Computation** (Simple):
```r
library(parallel)
cl <- makeCluster(detectCores() - 1)
profitsDeviation_byFirm <- parSapply(cl, coalitionFirms, function(firm){
  # ... calculation ...
})
stopCluster(cl)
```

**Option C - Reduce Solver Tolerance During Search**:
```r
# Use looser tolerance during binary search, tight only at end
control_loose <- modifyList(object@control.equ, list(tol = 1e-3))
# vs. final tol = 1e-6
```

**Estimated Speedup**: 
- Option A: 3-5x (complex implementation)
- Option B: Linear in cores (2-4x on typical machines)
- Option C: 1.5-2x (easy implementation)

---

### 3. **Ownership Matrix Conversion in ple()**
**Location**: Lines 180-207  
**Issue**: Nested loops to convert ownership vectors to matrices

```r
# CURRENT: O(n²) nested loops
if(!is.matrix(ownerPre)){
  ownerPreMat <- matrix(0, nprods, nprods)
  for(i in 1:nprods){
    for(j in 1:nprods){
      ownerPreMat[i,j] <- ifelse(ownerPre[i] == ownerPre[j], 1, 0)
    }
  }
}
```

**SOLUTION** (Vectorized):
```r
if(!is.matrix(ownerPre)){
  ownerPreMat <- outer(ownerPre, ownerPre, "==") * 1
}
```

**Estimated Speedup**: 10-50x for this operation (but small % of total time)

---

## Moderate Bottlenecks (MEDIUM IMPACT)

### 4. **Repeated ownerVec Extraction**
**Location**: Lines 473-476, 527-532, 823-828 (multiple occurrences)  
**Issue**: Same `apply()` operation repeated multiple times

```r
# Called 3+ times in calcSlopes alone
ownerVec <- if(is.vector(owner) || is.factor(owner)){
  owner
} else {
  apply(owner, 1, function(row) which(row > 0)[1])
}
```

**SOLUTION**: Cache result
```r
# At start of calcSlopes
ownerVec <- .extractOwnerVec(object@ownerPre)  # Helper function

# Helper:
.extractOwnerVec <- function(owner){
  if(is.vector(owner) || is.factor(owner)){
    return(owner)
  }
  apply(owner, 1, function(row) which(row > 0)[1])
}
```

---

### 5. **Profit Aggregation by Firm - Repeated sapply**
**Location**: Lines 539-567, 918-926  
**Issue**: Three separate `sapply` calls that iterate over same firms

```r
# CURRENT: Three separate passes
profitsPL_byFirm <- sapply(coalitionFirms, function(firm){
  sum(profitsPL[ownerVec == firm])
})
profitsDeviation_byFirm <- sapply(coalitionFirms, function(firm){
  sum(profitsDeviation[ownerVec == firm])
})
profitsBertrand_byFirm <- sapply(coalitionFirms, function(firm){
  sum(profitsBertrand[ownerVec == firm])
})
```

**SOLUTION**: Single pass with matrix operations
```r
# Pre-compute firm-product mapping matrix once
firmProductMap <- outer(ownerVec, coalitionFirms, "==")
# firmProductMap[i,f] = TRUE if product i owned by firm f

# Then vectorized aggregation
profitsPL_byFirm <- colSums(profitsPL * firmProductMap)
profitsDeviation_byFirm <- colSums(profitsDeviation * firmProductMap)
profitsBertrand_byFirm <- colSums(profitsBertrand * firmProductMap)
```

**Estimated Speedup**: 2-3x for aggregation operations

---

### 6. **Coalition Extraction - Repeated Logic**
**Location**: Lines 336-341, 632-637, 764-774  
**Issue**: Same coalition extraction code copy-pasted

```r
# Repeated pattern:
if(is.matrix(object@coalitionPre)){
  coalition <- which(rowSums(object@coalitionPre) > 1)
} else {
  coalition <- object@coalitionPre
}
```

**SOLUTION**: Helper function
```r
.getCoalitionIndices <- function(coalition_slot){
  if(is.matrix(coalition_slot)){
    return(which(rowSums(coalition_slot) > 1))
  }
  return(coalition_slot)
}

# Usage:
coalition <- .getCoalitionIndices(object@coalitionPre)
```

---

## Minor Bottlenecks (LOW IMPACT)

### 7. **Unnecessary if/else in Profit Aggregation**
**Location**: Lines 539-567  
**Issue**: Redundant conditional

```r
profitsPL_byFirm <- sapply(coalitionFirms, function(firm){
  if(is.numeric(firm)){
    sum(profitsPL[ownerVec == firm])
  } else {
    sum(profitsPL[ownerVec == firm])  # SAME CODE
  }
})
```

**SOLUTION**: Just remove condition
```r
profitsPL_byFirm <- sapply(coalitionFirms, function(firm){
  sum(profitsPL[ownerVec == firm])
})
```

---

### 8. **Coalition Product Loop Inefficiency**
**Location**: Lines 204-217  
**Issue**: Growing vector in loop (`coalitionPost <- c(coalitionPost, partner)`)

**SOLUTION**: Pre-allocate or use list
```r
# Collect new partners first
newPartners <- integer()
for(coalProd in coalitionPre){
  partners <- which(newCommonOwnership[coalProd, ])
  newPartners <- c(newPartners, setdiff(partners, coalitionPost))
}
coalitionPost <- unique(c(coalitionPost, newPartners))
```

---

### 9. **String Conversion in Named Vector Access**
**Location**: Line 591, 932  
**Issue**: Repeated `as.character(firm)` calls

```r
# CURRENT
profitsDeviation_byFirm[as.character(firm)]

# BETTER: Pre-convert names or use integer indexing
names(profitsPL_byFirm) <- as.character(coalitionFirms)
# Then access by integer position or pre-converted names
```

---

## Memory Bottlenecks

### 10. **Object Copying in Loops**
**Location**: Throughout `calcSupermarkup`  
**Issue**: R's copy-on-modify semantics may copy entire object in sapply

```r
profitsDeviation_byFirm <- sapply(coalitionFirms, function(firm){
  object@pricePre <- plPricesWithMarkup  # Potential copy
  # ... modifications ...
})
```

**SOLUTION**: Use environment to avoid copying
```r
# Create local environment
calcEnv <- new.env()
calcEnv$object <- object

# Modify within environment (no copies)
profitsDeviation_byFirm <- sapply(coalitionFirms, function(firm){
  calcEnv$object@pricePre <- plPricesWithMarkup
  # ... rest ...
})
```

---

## Recommended Implementation Priority

### Phase 1: Quick Wins (1-2 hours, 50-60% speedup)
1. ✅ **Fix #1**: Move Bertrand calculation outside binary search loop
2. ✅ **Fix #3**: Vectorize ownership matrix conversion
3. ✅ **Fix #7**: Remove redundant conditionals
4. ✅ **Fix #4**: Cache ownerVec extraction

### Phase 2: Moderate Effort (4-6 hours, additional 20-30% speedup)
5. ✅ **Fix #5**: Vectorize profit aggregation
6. ✅ **Fix #6**: Extract coalition indices helper
7. ✅ **Fix #2C**: Reduce solver tolerance during binary search

### Phase 3: Advanced Optimization (8-16 hours, additional 50-100% speedup)
8. ✅ **Fix #2B**: Parallelize deviation calculations
9. ✅ **Fix #2A**: Warm-start equilibrium solves (needs research)
10. ✅ **Fix #10**: Environment-based computation to avoid copies

---

## Profiling Recommendations

Before implementing, profile with realistic data:

```r
library(profvis)

profvis({
  result <- ple(
    prices = prices,
    shares = shares,
    margins = margins,
    coalitionPre = coalition,
    ownerPre = owners,
    ownerPost = owners
  )
})
```

Focus optimization on functions consuming >5% of runtime.

---

## Expected Overall Impact

**Conservative Estimate** (Phase 1 only):
- Runtime reduction: 50-60%
- Development time: 2-3 hours
- Risk: Very low (straightforward optimizations)

**Aggressive Estimate** (All phases):
- Runtime reduction: 75-85%
- Development time: 15-20 hours
- Risk: Medium (requires testing, potential bugs in parallel code)

---

## Code Quality Notes

**Positive**:
- ✅ Clear separation of concerns
- ✅ Proper use of S4 methods
- ✅ Good documentation
- ✅ Theoretically correct (verified separately)

**Areas for Improvement**:
- ⚠️ Excessive code duplication (coalition extraction, owner vec extraction)
- ⚠️ Mixed use of apply/sapply (inconsistent style)
- ⚠️ No caching of intermediate results
- ⚠️ No consideration of computational complexity in design

---

## Next Steps

1. **Profile current implementation** with representative data
2. **Implement Phase 1 fixes** (quick wins)
3. **Re-profile** to measure impact
4. **Decide on Phase 2/3** based on remaining bottlenecks
5. **Add performance tests** to prevent regressions
