# WARDEN Efficiency Patterns

## Step 0 — Always Profile First

Never recommend optimisations without profiling. Understanding where time is actually spent prevents wasted effort.

```r
# Save as profile_run.R
library(WARDEN)
# [paste model setup code here]

Rprof("warden_profile.out", interval=0.005)
results <- run_sim(npats=500, n_sim=5, psa_bool=TRUE, seed=42, ...)
Rprof(NULL)

prof <- summaryRprof("warden_profile.out")
print(prof$by.self[1:20, ])   # top 20 time consumers
```

**Interpreting profiling output**:
- `self.time`: time spent IN this function (not its callees) — the primary diagnostic
- `total.time`: total time including all sub-calls
- Functions at the top of `$by.self` are the actual bottlenecks

Tell the user:
1. The top 3–5 time consumers by `self.time`
2. What each function does in the WARDEN context
3. Which pattern below addresses it

---

## Pattern E1 — Replace Cycles with `qtimecov` + `adj_val`

**Problem indicator**: profiling shows many calls to `new_event` / event scheduling; model has yearly (or sub-yearly) cycle events; large `number_events` per patient in IPD output.

**Why cycles are slow**: each cycle creates a new event entry in the C++ queue, which must be processed with a full reaction evaluation. A model with 30-year horizon and 1-year cycles generates 30+ events per patient per arm — a 1000-patient model becomes 60,000+ event evaluations instead of a handful.

**Solution**: Use `qtimecov()` to compute time-varying TTE analytically at model start; use `adj_val()` for time-varying utility/cost adjustments.

```r
# BEFORE (slow — creates 30+ events per patient):
init_event_list <- add_tte(arm=c("int","noint"), evts=c("start","cycle","death"), input={
  start <- 0; cycle <- 1; death <- qexp(luck, rate=bg_rate)
})
evt_react_list <- add_reactevt("cycle", input={
  new_event(c(cycle=curtime+1))
  modify_event(c(death=curtime+qcond_exp(luck2, rate=new_rate)))
  q_default <- u_base * age_adj[floor(age+curtime)]
}) |> ...

# AFTER (fast — 2 events per patient):
init_event_list <- add_tte(arm=c("int","noint"), evts=c("start","death"), input={
  start <- 0
  # Rate as function of time (changes annually):
  rate_fn <- function(.t) { bg_rate + covariate_effect * floor(.t) }
  death <- qtimecov(luck=luck, a_fun=rate_fn, dist="exp", dt=1)$tte
})
evt_react_list <-
  add_reactevt("start", input={
    # Utility adjustment across the whole period until death — single value:
    adj <- adj_val(curtime, next_event()$time, by=1,
                   age_adj[min(100, floor(.time + baseline_age))],
                   discount=drq)
    q_default <- u_base * adj
  }) |>
  add_reactevt("death", input={ curtime <- Inf })
```

**Expected impact**: 60–80% runtime reduction for cycle-heavy models. Undiscounted outcomes may differ slightly from the cycle approach — document this if needed.

**Tip — `vectorized_f = TRUE`**: when the expression inside `adj_val()` supports a vector `.time` argument (e.g., `age_adj[floor(.time + baseline_age)]` where `age_adj` is a named vector), pass `vectorized_f = TRUE` to avoid internal per-step looping. The expression receives a vector of time points and must return a vector of the same length. This can yield an additional 20–40% speedup for models with long time horizons.

```r
# FORWARD accumulation (unconstrained models only):
adj <- adj_val(curtime, next_event()$time, by=1,
               age_adj[pmin(100L, as.integer(.time + baseline_age))],
               discount=drq, vectorized_f=TRUE)

# BACKWARD accumulation (constrained models, or any model using accum_backwards=TRUE):
adj <- adj_val(prevtime, curtime, by=1,
               age_adj[pmin(100L, as.integer(.time + baseline_age))],
               discount=drq, vectorized_f=TRUE)
```

**Caveat — `next_event()$time` in constrained models**: `next_event()$time` is only reliable in forward accumulation (unconstrained models). In resource-constrained models (`constrained=TRUE`), the next event can be pre-empted by a resource release from another patient. Always use `accum_backwards = TRUE` with `adj_val(prevtime, curtime, ...)` in constrained models — the backward boundaries are certain.

---

## Pattern E2 — Batch `modify_event` Calls

**Problem indicator**: profiling shows `modify_event` appearing many times in `self.time`; multiple separate `modify_event` calls in the same reaction block.

**Why it's slow**: each `modify_event` call crosses the R→C++ boundary to update the priority queue. Each crossing has fixed overhead.

**Solution**: combine all modifications into a single call.

```r
# BEFORE (slow — 3 C++/R crossings):
add_reactevt("progression", input={
  modify_event(c("death"    = curtime + os_pps))
  modify_event(c("followup" = curtime + 0.5))
  modify_event(c("ae"       = curtime + 0.1))
})

# AFTER (fast — 1 crossing):
add_reactevt("progression", input={
  modify_event(c("death"    = curtime + os_pps,
                 "followup" = curtime + 0.5,
                 "ae"       = curtime + 0.1))
})
```

**Expected impact**: small but consistent (3–10%) in models with many event modifications per patient.

---

## Pattern E3 — Replace `case_when` / `if_else` with `ifelse` or `if`

**Problem indicator**: profiling shows `case_when`, `dplyr::if_else`, or similar in `self.time`; dplyr in the call stack.

**Why it's slow**: `dplyr::case_when()` and `if_else()` are designed for vectorised data frame operations and carry dplyr infrastructure overhead. Inside a per-patient, per-event loop they're called millions of times.

**Solution**: use base R `ifelse()` or `if(...) {...} else {...}`.

```r
# BEFORE (slow — dplyr overhead per event):
q_default <- case_when(
  fl.prog == 0 ~ util.pfs,
  fl.prog == 1 & fl.2l == 0 ~ util.pps1,
  TRUE ~ util.pps2
)

# AFTER (fast):
q_default <- ifelse(fl.prog == 0, util.pfs,
             ifelse(fl.2l == 0, util.pps1, util.pps2))
# OR, for readability (equally fast):
q_default <- if(fl.prog == 0) util.pfs else if(fl.2l == 0) util.pps1 else util.pps2
```

**Expected impact**: 5–20% speedup in models with complex state-dependent utility/cost expressions evaluated at every event.

---

## Pattern E4 — Pre-Allocate Random Numbers

**Problem indicator**: profiling shows `runif`, `rnorm`, `rexp`, or similar distribution functions in `self.time` inside event reactions; model uses inline `r*()` calls inside reactions.

**Why it's slow**: drawing from R distributions inside event reactions (a) disrupts the random state management WARDEN handles per patient, and (b) has per-call overhead when called millions of times in a PSA.

**Solution**: pre-draw all needed random values before simulation.

```r
# BEFORE (problematic — inline draw):
add_reactevt("ae", input={
  n_ae <- rpois(1, lambda=rate_ae * (curtime - prevtime))
  if(n_ae > 0) new_event(c(ae=curtime + rexp(1, rate_ae)))
})

# AFTER (correct — pre-drawn stream):
unique_pt_inputs <- add_item(
  rnd_ae_stream <- random_stream(200)  # pre-generate 200 uniform draws
)
add_reactevt("ae", input={
  n_ae <- qpois(rnd_ae_stream$draw_n(), lambda=rate_ae * (curtime - prevtime))
  if(n_ae > 0) {
    next_ae_time <- curtime + qexp(rnd_ae_stream$draw_n(), rate_ae)
    if(next_ae_time < get_event("death")) new_event(c(ae=next_ae_time))
  }
})
```

**Expected impact**: correctness improvement in constrained models; 5–15% in complex PSA models.

---

## Pattern E5 — Use `run_sim_parallel` for PSA

**Problem indicator**: `n_sim` ≥ 100 and each simulation takes > 2 seconds; user is running PSA with `run_sim()`.

**Why `run_sim_parallel` helps**: parallelises at the simulation level across multiple CPU cores. Each simulation is independent — no shared state between simulations — so parallelism is safe and effective.

**When NOT to use it**:
- `n_sim = 1` (deterministic): no benefit; setup cost ~2–5s outweighs gain
- `n_sim < 20` and each simulation runs < 1s: setup overhead exceeds benefit
- Still developing/debugging: parallel mode suppresses `debug=TRUE` logs
- Constrained models with very high RAM use (>2GB per simulation × n_cores may exhaust RAM)

```r
library(future)
plan(multisession, workers = parallel::detectCores() - 1)  # leave 1 core for OS

results_psa <- run_sim_parallel(
  npats    = 1000,
  n_sim    = 500,
  psa_bool = TRUE,
  seed     = 42,
  ...
)

# Reset to sequential when done:
plan(sequential)
```

**Expected impact**: 20–40% wall-clock reduction for large PSA. Not linear due to RAM duplication across sessions. A machine with 32GB RAM and 8 cores will benefit substantially; one with 8GB and 4 cores may see less gain.

---

## Pattern E6 — Move Constant Computations Out of the Event Loop

**Problem indicator**: profiling shows the same data manipulation or function call appearing repeatedly (e.g. reading from a data frame, computing survival curve values) inside event reactions.

**Why it's slow**: anything computed inside `add_reactevt` is evaluated once per event per patient per simulation. A computation run 50 times per patient × 1000 patients × 500 PSA = 25 million evaluations.

**Solution**: pre-compute any value that doesn't change during the simulation into `common_all_inputs` or `unique_pt_inputs`.

```r
# BEFORE (slow — reads data frame inside reaction):
add_reactevt("cycle", input={
  mort_rate <- life_table$rate[life_table$age == floor(age + curtime)]
  modify_event(c(death = curtime + qexp(luck, mort_rate)))
})

# AFTER (fast — vectorise lookup, pre-convert to simple vector):
common_all_inputs <- add_item(
  mort_vec <- setNames(life_table$rate, life_table$age)  # named vector indexed by age
)
add_reactevt("cycle", input={
  mort_rate <- mort_vec[as.character(floor(age + curtime))]
  modify_event(c(death = curtime + qexp(luck, mort_rate)))
})
```

---

## Pattern E7 — IPD Output Minimisation

**Problem indicator**: model runs slow and `ipd=1` is set; `input_out` vector has many variables; results object is very large.

**Why it's slow/large**: `ipd=1` returns one row per event per patient per arm per simulation — for a 1000-patient, 10-event, 2-arm, 500-sim PSA this is 10 million rows.

**Solution**: use `ipd=2` (aggregates per patient-arm) unless per-event granularity is needed; minimise `input_out` to only what's actually used.

```r
# For PSA analysis — usually only summary statistics needed, not IPD:
results <- run_sim(..., ipd=0, input_out=character())  # no IPD, fastest

# For model validation/debugging:
results <- run_sim(..., ipd=1, input_out=c("fl.prog","fl.on_tx"))  # minimal flags

# For QC of event-level data (only during development):
results <- run_sim(..., ipd=1, input_out=c("all","the","flags","needed"))
```

---

## Pattern E8 — Data.frame to Named Vector Conversion

**Problem indicator**: profiling shows data.frame indexing (`[.data.frame`, `[<-.data.frame`, `==`) or dplyr verbs (`filter`, `pull`) in `self.time`; event reactions reference `df[df$col == val, "result"]` patterns.

**Why it's slow**: data.frame row subsetting is O(n) — it scans every row to find matches. A named vector lookup is O(1) hash table access. In a 1000-patient model with 5 events each and a 50-row lookup table, that's 5000 O(50) scans vs 5000 O(1) lookups.

**Solution**: convert lookup tables to named vectors once in `common_all_inputs`.

```r
# BEFORE (slow — O(n) scan per event):
add_reactevt("progression", input={
  drug_cost <- cost_table[cost_table$state == "pps" & cost_table$arm == arm, "monthly"]
  q_default <- util_table[util_table$state == "pps", "value"]
})

# AFTER (fast — O(1) lookup per event):
common_all_inputs <- add_item(
  # Convert once at model start:
  cost_by_state_arm <- setNames(
    cost_table$monthly,
    paste0(cost_table$state, "_", cost_table$arm)
  ),
  util_by_state <- setNames(util_table$value, util_table$state)
)

add_reactevt("progression", input={
  drug_cost <- cost_by_state_arm[[paste0("pps_", arm)]]
  q_default <- util_by_state[["pps"]]
})
```

**For multi-column lookups**: create composite keys with `paste0()` at conversion time. For single-column lookups, `setNames(df$value, df$key)` is sufficient.

**Expected impact**: 20–60% speedup for lookup-heavy models (common in multi-state models with state-dependent costs/utilities from parameter tables).

---

## Pattern E9 — Named Vector Dispatch for State Logic

**Problem indicator**: multiple `if/else if/else` branches (4+) in event reactions selecting values based on a state flag; profiling shows time in reaction evaluation proportional to number of branches.

**Why it's slow**: each `if/else` is evaluated sequentially — in the worst case, all conditions are tested before reaching the matching branch. A named vector is direct hash indexing regardless of the number of states.

**Solution**: replace multi-branch logic with named vector indexing.

```r
# BEFORE (slow — sequential condition evaluation):
add_reactevt("cycle_update", input={
  if(fl.state == "pfs") {
    q_default <- 0.82
    c_default <- 1500
  } else if(fl.state == "response") {
    q_default <- 0.90
    c_default <- 2000
  } else if(fl.state == "stable") {
    q_default <- 0.75
    c_default <- 1800
  } else if(fl.state == "progressed") {
    q_default <- 0.55
    c_default <- 3500
  } else {
    q_default <- 0.40
    c_default <- 5000
  }
})

# AFTER (fast — O(1) dispatch):
common_all_inputs <- add_item(
  util_by_state <- c(pfs=0.82, response=0.90, stable=0.75,
                     progressed=0.55, bsc=0.40),
  cost_by_state <- c(pfs=1500, response=2000, stable=1800,
                     progressed=3500, bsc=5000)
)

add_reactevt("cycle_update", input={
  q_default <- util_by_state[[fl.state]]
  c_default <- cost_by_state[[fl.state]]
})
```

**When if/else is still appropriate**: 2–3 branches, complex conditions that aren't simple equality checks (e.g., `if(curtime > 5 & fl.on_tx == 1)`), or when the branch body contains more than a simple assignment.

**Expected impact**: 5–15% for models with many states (5+) and frequent state-dispatch events. The primary benefit is readability and maintainability — adding a new state is one vector element, not a new branch.

---

## Timing Template

After applying optimisations, always show before/after timing:

```r
# Time the original:
t0 <- proc.time()
results_orig <- run_sim(npats=500, n_sim=10, ...)
cat("Original:", round((proc.time()-t0)[3],2), "s\n")

# Time the optimised:
t1 <- proc.time()
results_opt  <- run_sim_parallel(npats=500, n_sim=10, ...)
cat("Optimised:", round((proc.time()-t1)[3],2), "s\n")

# Verify results are consistent:
cat("ICER original:", summary_results_det(results_orig[[1]][[1]])["ICUR","int"], "\n")
cat("ICER optimised:", summary_results_det(results_opt[[1]][[1]])["ICUR","int"], "\n")
```

Always verify that the optimised model produces results consistent with the original (within Monte Carlo noise). Report both the speedup and any caveats (e.g. slight numerical differences from `adj_val` vs cycle approach for undiscounted outcomes).

---

## Quick Reference: Optimisation Priority

| Issue found | Expected speedup | Effort | Priority |
|-------------|-----------------|--------|----------|
| Yearly cycles (30+ events/patient) | 60–80% | Medium | **Do first** |
| `run_sim` → `run_sim_parallel` for PSA | 20–40% | Low | **Do second** |
| df subsetting → named vector | 20–60% | Low | **Do third** |
| if/else chains → named vector dispatch | 5–15% | Very low | Easy win |
| Multiple `modify_event` → single call | 3–10% | Very low | Easy win |
| `case_when` → `ifelse` | 5–20% | Very low | Easy win |
| `adj_val` with `vectorized_f=TRUE` | 20–40% (of adj_val time) | Very low | Easy win |
| Inline r* draws → pre-drawn | Variable | Low | Correctness + speed |
| `ipd=1` → `ipd=0` for PSA | Memory + 5–10% | Very low | Do always for PSA |
| Move constants out of reactions | 10–30% | Medium | Model-dependent |
| `accum_backwards=TRUE` (constrained) | 0% (no speedup) | Very low | **Required** for correctness |
