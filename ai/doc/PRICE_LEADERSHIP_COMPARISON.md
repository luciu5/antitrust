# Comparison: Existing Grim Trigger vs. MSW Price Leadership Model

## Current Implementation: `calcProducerSurplusGrimTrigger`

### **What It Does:**
Evaluates incentive compatibility for collusion using Grim Trigger strategies in a Bertrand game.

### **Key Features:**

1. **Coalition Structure**:
   - User specifies `coalition` = vector of product indices in the collusive agreement
   - All coalition products cooperate by jointly maximizing profits
   - Full ownership matrix becomes `1` for coalition products: `ownerPre[coalition, coalition] <- 1`

2. **Three Profit States**:
   - **Cooperate** (`psCoord`): Profits when all coalition firms set joint profit-maximizing prices
   - **Defect** (`psDefect`): Profits when one firm optimally undercuts coalition prices for one period
   - **Punish** (`psPunish`): Bertrand profits forever after deviation

3. **Incentive Compatibility Check**:
   ```r
   IC = Σ(psCoord/(1-δ)) >= Σ(psDefect + psPunish*δ/(1-δ))
   ```
   - Present value of cooperation ≥ Present value of deviation + punishment

4. **Deviation Calculation**:
   - Loops over each firm in coalition
   - Temporarily sets that firm's ownership back to Bertrand
   - Recalculates prices with that firm best-responding to coalition prices
   - Stores deviation profit

5. **Outputs**:
   - Data frame with: Coalition products, Discount rates, Coord profit, Defect profit, Punish profit, IC (TRUE/FALSE)

---

## MSW (2023) Price Leadership Model

### **What It Does:**
Simulates coordinated effects where a **leader announces a supermarkup** that coalition firms adopt above Bertrand prices.

### **Key Features:**

1. **Two-Stage Game Structure**:
   - **Stage 1**: Leader announces supermarkup `m ≥ 0` (cheap talk)
   - **Stage 2**: Coalition sets `p^PL = p^B + m`; Fringe best-responds
   - Repeated infinitely with timing parameter `δ ∈ (0,1)`

2. **Coalition vs. Fringe**:
   - **Coalition** (C): Firms that follow leader's supermarkup announcement
   - **Fringe** (F): Firms that maximize static profit (like Bertrand)
   - Leader is firm 1 ∈ C

3. **Three Profit States**:
   - **π^PL_f(m)**: Profit when coalition adopts supermarkup m
   - **π^D_f(m)**: Optimal deviation profit (best-respond to p^PL(m))
   - **π^B_f**: Bertrand profit (punishment)

4. **Slack Function** (Incentive Compatibility):
   ```r
   g_f(m) = [π^PL_f(m) - π^D_f(m)] + δ/(1-δ) * [π^PL_f(m) - π^B_f]
   ```
   - Must satisfy `g_f(m) ≥ 0` for all `f ∈ C`

5. **Leader's Problem**:
   ```r
   max_{m ≥ 0} π^PL_1(m)  subject to  g_f(m) ≥ 0  ∀f ∈ C
   ```
   - Leader chooses supermarkup to maximize own profit
   - Subject to IC constraints of all coalition members

6. **Equilibrium Selection**:
   - **Constrained PLE**: At least one IC constraint binds
   - **Unconstrained PLE**: Leader sets monopoly-like markup (if δ high enough)

7. **Calibration** (with Logit demand):
   - From shares, 2 margins (1 coalition, 1 fringe), 1 diversion → identify (α, m, δ, costs)
   - Fringe margin identifies α (price coefficient)
   - Coalition margin identifies m (supermarkup)
   - Binding IC constraint identifies δ (timing parameter)

---

## Key Differences

| Feature | `calcProducerSurplusGrimTrigger` | MSW Price Leadership |
|---------|----------------------------------|----------------------|
| **Purpose** | Check if current collusion is IC | Simulate equilibrium with coordination |
| **Coordination** | Full joint profit maximization | Leader-follower with supermarkup |
| **Fringe Firms** | Not explicitly modeled | Explicitly in model, price à la Bertrand |
| **Equilibrium** | Tests one configuration | Solves for equilibrium supermarkup |
| **Calibration** | Uses existing demand/cost params | Calibrates from margins + diversion |
| **Leader** | No designated leader | Firm 1 is leader |
| **Supermarkup** | Implicit (joint pricing) | Explicit parameter `m` |
| **Output** | IC check (TRUE/FALSE) | Full merger simulation with prices |
| **Merger Analysis** | Pre/post comparison of IC | Full counterfactual with price changes |

---

## Conceptual Similarities

Both models:
1. ✅ Use **Grim Trigger** punishment (revert to Bertrand forever)
2. ✅ Calculate three profit states: Cooperation, Deviation, Punishment
3. ✅ Check incentive compatibility using discounted profit flows
4. ✅ Allow for product-level or firm-level analysis
5. ✅ Can be used for merger analysis (comparing pre/post)

---

## What MSW Price Leadership Adds

1. **Fringe Firms**: Explicitly models non-coordinating competitors
2. **Supermarkup**: Observable coordination parameter that can be calibrated
3. **Leader Selection**: Endogenous through leader's optimization
4. **Binding Constraints**: Identifies which firm limits coordination
5. **Merger Counterfactuals**: Changes in `m` due to merger
6. **Calibration Framework**: Recovers all parameters from observables

---

## Implementation Strategy

### **Option 1: Extend `calcProducerSurplusGrimTrigger`**
- Add `leader`, `fringe`, `supermarkup` parameters
- Add calibration routine
- Return full simulation object instead of data frame

**Pros**: Builds on existing code, familiar to users
**Cons**: Different design philosophy, harder to extend

### **Option 2: New `PriceLeadership` Class** ⭐ **RECOMMENDED**
- Inherit from `Logit` (or `Bertrand`)
- New slots: `@supermarkup`, `@timingParam`, `@coalition`, `@leader`
- Implement full MSW framework
- Methods: `calcPrices`, `calcSlack`, `calcSupermarkup`

**Pros**: Clean separation, follows package design, extensible to BLP
**Cons**: More code, separate from existing Grim Trigger function

### **Option 3: Hybrid**
- Keep `calcProducerSurplusGrimTrigger` for simple IC checks
- Create `PriceLeadership` class for full MSW implementation
- Both available depending on user needs

**Pros**: Flexibility for different use cases
**Cons**: Some code duplication

---

## Recommendation

**Implement Option 2** (New `PriceLeadership` class) because:

1. MSW model is fundamentally different (leader-follower vs. joint maximization)
2. Calibration requires new framework (margins → supermarkup)
3. Extensibility to LogitBLP/nested logit easier with dedicated class
4. Follows package architecture (Bertrand, Cournot, etc.)
5. Can reference existing Grim Trigger for IC calculations

The existing `calcProducerSurplusGrimTrigger` remains useful for:
- Quick IC checks without full calibration
- Scenarios without clear leader/fringe structure
- Research on different collusion mechanisms

---

## Next Steps

1. Create `PriceLeadership` class inheriting from `Logit`
2. Implement calibration (shares + 2 margins + diversion → α, m, δ, costs)
3. Implement `calcPrices` with supermarkup logic
4. Implement slack functions and IC constraint identification
5. Test with examples from MSW paper (beer industry)
6. Document differences from `calcProducerSurplusGrimTrigger`
