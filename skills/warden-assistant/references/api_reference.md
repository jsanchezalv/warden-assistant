# WARDEN API Reference

## Package Overview

WARDEN (v2.0+) runs a nested loop: per sensitivity → per simulation → per patient → per arm → while curtime < Inf.
All inputs use delayed execution — expressions are stored as lists and evaluated at the right scope level.
The C++ event queue processes events in time order; ties broken by declaration order in `add_tte(..., evts=)`.
Default discount rates: `drc = 0.03` (costs), `drq = 0.03` (QALYs/utilities) if not set in inputs.

---

## Input Layer Functions

### `add_item(input = { ... })` or `add_item(...)`
Stores a delayed-execution expression (or named arguments) as the input block.
Can be chained with `|>`. Items defined earlier in the block are available to later ones.
Used for all four input levels: `sensitivity_inputs`, `common_all_inputs`, `common_pt_inputs`, `unique_pt_inputs`.

```r
common_all_inputs <- add_item(
  util.sick   = 0.8,
  util.sicker = 0.5,
  drc         = 0.035,
  drq         = 0.035
)
# OR with block syntax (equivalent):
common_all_inputs <- add_item(input = {
  util.sick   <- 0.8
  util.sicker <- 0.5
})
```

**Note**: Items are reset after each patient. Modifications to `common_all_inputs` items inside event reactions persist only for that patient-arm.

### `input_block(base, psa, sens=NULL, names_out, psa_indicators=NULL, sens_indicators=NULL, indicator_sens_binary=FALSE, dsa_names=NULL)`
High-level wrapper for `pick_val_v()`. Auto-detects `n_sensitivity` for DSA from block metadata.
Returns a named list of parameter values depending on whether PSA/DSA/scenario is active.

```r
common_all_inputs <- input_block(
  base      = param_df$base_value,
  psa       = pick_psa(param_df$dist, param_df$n, param_df$a, param_df$b),
  sens      = param_df,          # data frame with DSA_min, DSA_max columns
  names_out = param_df$name,
  psa_indicators = param_df$psa_flag,  # 0/1 vector, which params are PSA-active
  dsa_names = c("DSA_min", "DSA_max")
) |> add_item(other_param = 42)
```

### `pick_val_v(base, psa, sens=NULL, psa_ind, sens_ind, indicator, names_out, indicator_psa=NULL, deploy_env=TRUE)`
Lower-level input selector. Returns base, PSA draw, or sensitivity value per parameter.
- `psa_ind = psa_bool` (set by model automatically)
- `sens_ind = sensitivity_bool` (set by model automatically)
- `indicator` — vector of 0/1 for which params are active in the current sensitivity iteration

### `pick_psa(dist_vec, n_vec, a_vec, b_vec)`
Calls each distribution function (`rnorm`, `rlnorm`, `rbeta`, `rgamma`, etc.) with the given arguments.
Returns a list of one draw per parameter.

### `create_indicators(sens, n_sensitivity_total, base_indicators)`
Returns a 0/1 vector marking which parameter is active in the current sensitivity iteration.
Used internally when `sensitivity_bool = TRUE`.

---

## Event Setup Functions

### `add_tte(arm, evts, other_inp=NULL, input={ ... })`
Defines initial time-to-event for one or more arms.
- `arm`: character vector of arm names (e.g. `c("int","noint")`) or single arm
- `evts`: character vector of event names. Order determines tie-breaking priority.
- `other_inp`: character vector of intermediate objects to store (not events)
- `input`: expression evaluated per patient-arm. Events not assigned here default to `Inf`.

```r
init_event_list <- add_tte(
  arm  = c("int","noint"),
  evts = c("start","progression","death"),  # "start" processed first on tie
  input = {
    start       <- 0
    progression <- draw_tte(1, "weibull", coef1=shape, coef2=scale,
                            beta_tx = ifelse(arm=="int", HR, 1))
    death       <- draw_tte(1, "lnorm",   coef1=os_mu, coef2=os_sigma)
  }
) |> add_tte(arm="int", evts=c("ae"), input={ ae <- Inf })
```

### `add_reactevt(name_evt, input={ ... })`
Defines the reaction when a named event fires. Chained with `|>`.
Inside the block, modify items by direct assignment or call event manipulation functions.

```r
evt_react_list <-
  add_reactevt("start",       input = { q_default <- util.pfs; c_default <- cost.drug }) |>
  add_reactevt("progression", input = {
    q_default <- util.pps
    c_default <- cost.bsc
    fl.prog   <- 1
    modify_event(c("death" = min(curtime + draw_tte(1,"exp",0.5), get_event("death"))))
  }) |>
  add_reactevt("death",       input = { curtime <- Inf })
```

---

## Context Variables (available inside `add_reactevt` blocks)

| Variable | Type | Description |
|----------|------|-------------|
| `curtime` | numeric | Current event time. Set to `Inf` to terminate simulation. |
| `prevtime` | numeric | Time of previous event |
| `evt` | character | Name of current event being processed |
| `arm` | character | Current arm being iterated |
| `i` | integer | Current patient index |
| `simulation` | integer | Current simulation index (1 to n_sim) |
| `sens` | integer | Current sensitivity analysis index |
| `cur_evtlist` | XPtr | C++ event queue pointer (do not manipulate directly) |
| `npats` | integer | Total number of patients |
| `psa_bool` | logical | Whether PSA is active |
| `sensitivity_bool` | logical | Whether DSA/scenario is active |
| `sens_name_used` | character | Name of current sensitivity column |

---

## Event Manipulation Functions (use inside `add_reactevt`)

### `modify_event(events, create_if_missing=TRUE)`
Change the time of existing events. Named vector: `c("event_name" = new_time)`.
- With `create_if_missing=TRUE` (default): creates the event if it doesn't exist
- With `create_if_missing=FALSE`: errors if event doesn't exist
- Group all modifications into one call for speed

```r
modify_event(c("death" = curtime + 2, "ae" = curtime + 0.5))
```

### `new_event(named_vec)`
Add new events to the queue. Only one event per event type at a time.

```r
new_event(c("ae" = curtime + 0.25, "followup" = curtime + 1))
```

### `remove_event(event_name)`
Remove an event from the queue. Equivalent to `modify_event(c(name = Inf))`.

### `get_event(event_name)`
Returns current scheduled TTE of a named event. Returns Inf if not in queue.

### `has_event(event_name)`
Returns TRUE/FALSE whether the event exists in the queue (including Inf events).

### `next_event(n=1, ptr)`
Returns a list of the next `n` events across ALL patients in the shared queue (constrained engine). Each element has `$time`, `$event_name`, `$patient_id`. In unconstrained mode, where each patient has their own queue, this is equivalent to `next_event_pt()`.

**Usage in cost/utility calculations**: `next_event()$time` is valid in forward-accumulation (`accum_backwards=FALSE`) unconstrained models only. In constrained models, another patient's resource release can fire before the predicted next event — making `next_event()$time` unreliable as a boundary.

### `next_event_pt(n=1, ptr, patient_id)`
Returns a list of the next `n` events for the current patient only. Each element has `$time`, `$event_name`, `$patient_id`. If `patient_id` is omitted, uses `i` from the calling environment.

**Usage**: useful for inspecting what's coming next for a patient (e.g., branching logic based on event type). NOT safe as a time boundary in constrained models for the same reason as `next_event()`.

### `queue_empty()`
Returns TRUE if no events remain (Inf events count as present).

### `queue_size()`
Returns integer count of events in queue (including Inf).

---

## Time-to-Event Drawing Functions

### `draw_tte(n_chosen, dist, coef1=NULL, coef2=NULL, coef3=NULL, ..., beta_tx=1, seed=NULL)`
Draw time to event from parametric distribution, optionally applying a treatment effect.
- `dist`: `"exp"`, `"weibull"`, `"lnorm"`, `"llogis"`, `"gamma"`, `"gompertz"`, `"weibullPH"`, `"gengamma"`, `"beta"`, `"poisgamma"`
- `coef1`, `coef2`, `coef3`: distribution parameters (log-scale for flexsurv conventions — coef2 is `log(scale)` for Weibull/llogis)
- `beta_tx`: hazard ratio (default `1` = no effect). Internally applied as `log(beta_tx)` on the log-scale parameter. Pass the HR directly, not its log.
- `seed`: integer seed for reproducibility; typically `seed_vec[i]`

```r
# Weibull with HR for intervention
ttp <- draw_tte(1, "weibull",
                coef1 = shape_par,
                coef2 = log(scale_par),
                beta_tx = ifelse(arm=="int", HR, 1),
                seed = seed_vec[i])
```

### `qcond_exp(p, rate, lower_bound=0)` / `qcond_weibull()` / `qcond_lnorm()` / `qcond_gamma()` / `qcond_llogis()` / `qcond_norm()` / `qcond_gompertz()` / `qcond_weibullPH()`
Conditional quantile functions — draw TTE given patient has survived to `lower_bound`.
Use with pre-drawn uniform `luck` values for reproducibility.

```r
common_pt_inputs <- add_item(luck = runif(1))
# then in add_tte:
death <- qcond_exp(luck, rate = bg_mort_rate, lower_bound = 0)
```

### `rcond_gompertz(n, shape, rate, lower_bound=0)` / `rcond_gompertz_lu(...)`
Random draws from conditional Gompertz (useful for background mortality by age).

### `qtimecov(luck, a_fun, b_fun=NULL, dist="exp", dt=0.1, max_time=100, start_time=0)`
Compute TTE with time-varying hazard/distribution parameters. Avoids needing yearly cycles.
- `luck`: pre-drawn uniform(0,1) quantile for this patient
- `a_fun`: function of `.time` returning the first distribution parameter (e.g., rate for "exp")
- `b_fun`: function of `.time` returning the second parameter (e.g., shape for "weibull"); NULL for 1-param dists
- `dist`: `"exp"`, `"weibull"`, `"llogis"`, `"gompertz"`, `"gamma"`, `"lnorm"`, `"norm"`
- `dt`: time step for numerical integration (default 0.1; smaller = more accurate, slower)
- `max_time`: maximum horizon to search for event (default 100)
- `start_time`: offset to start integration from (default 0; use when calling mid-simulation)
- Returns list with `$tte` (time to event from start_time) and `$luck` (residual luck value for chaining)

```r
# Simple age-dependent mortality:
death <- qtimecov(luck=luck, a_fun=\(.time) base_rate * exp(0.1 * (.time + bs_age)),
                  dist="exp", dt=1)$tte

# Luck-chaining for mid-simulation recalculation (resource-constrained):
# At queue exit, advance luck through wait period, then draw remaining TTE:
result <- qtimecov(luck=luck_os, a_fun=rate_fn, dist="exp",
                   start_time=0, max_time=curtime)
remaining <- qtimecov(luck=result$luck, a_fun=new_rate_fn, dist="exp",
                      start_time=curtime, max_time=100)
death_time <- curtime + remaining$tte
```

### `luck_adj(prevsurv, cursurv, luck, condq=TRUE)`
Adjusts the luck value when distribution parameters change mid-simulation, maintaining the patient's survival rank. Use this when a hazard changes at a known time point (e.g., treatment switch, covariate update) and you need to recalculate the remaining TTE under new parameters.

- `prevsurv`: survival at the change-point under the **previous** parametrization
- `cursurv`: survival at the change-point under the **new** parametrization
- `luck`: current luck value (0-1 quantile)
- `condq`: if `TRUE` (default), uses conditional quantile approach (prevsurv/cursurv evaluated at previous and current timepoints with the SAME parametrization). If `FALSE`, unconditional approach (prevsurv/cursurv at the same timepoint with different parametrizations).
- Returns: adjusted luck value (0-1) to use with a quantile function under the new parameters

```r
# Treatment switch at time 10: Weibull rate changes from 0.02 to 0.025
new_luck <- luck_adj(prevsurv = 1 - pweibull(10, 3, 1/0.02),
                     cursurv  = 1 - pweibull(10, 3, 1/0.025),
                     luck = patient_luck, condq = FALSE)
new_tte <- qweibull(new_luck, 3, 1/0.025)
```

Note: `qtimecov()` handles continuous time-varying parameters automatically (preferred for smooth changes). Use `luck_adj()` when parameters change at discrete known time points and you want to preserve exact distributional properties between parametric families.

---

## Outcome Accumulation Functions

### `adj_val(curtime, nexttime, by, expression, discount=NULL, vectorized_f=FALSE)`
Compute a time-weighted adjustment factor for ongoing values that change with time (e.g. age utility).
Returns a scalar that, when multiplied by the base utility, gives correctly weighted results without cycles.
- `curtime`: start of interval
- `nexttime`: end of interval
- `by`: step size for evaluation (e.g., 1 for annual)
- `expression`: R expression using `.time` as the time variable (evaluated at each step)
- `discount`: discount rate (NULL = no discounting; pass `drq` for utilities, `drc` for costs)
- `vectorized_f`: if `TRUE`, evaluates expression once with `.time` as a vector (faster when expression supports vectorized input)

**Forward accumulation** (default, unconstrained models):
```r
adj_factor <- adj_val(curtime, next_event()$time, by=1,
                      u_age[min(100, floor(.time + baseline_age))],
                      discount=drq)
q_default <- u_base * adj_factor
```

**Backward accumulation** (resource-constrained models with `accum_backwards=TRUE`):
```r
# Both boundaries are known — always correct regardless of queue interrupts
adj_factor <- adj_val(prevtime, curtime, by=1,
                      u_age[min(100, floor(.time + baseline_age))],
                      discount=drq, vectorized_f=TRUE)
q_default <- u_base * adj_factor
```

**Note**: In resource-constrained models, use `accum_backwards=TRUE` in `run_sim()` and compute `adj_val(prevtime, curtime, ...)`. Forward-looking `adj_val(curtime, next_event()$time, ...)` is wrong when queue events can pre-empt the next scheduled event.

Outcome types in `run_sim()`:
- `util_ongoing_list` / `cost_ongoing_list` — accrued continuously between events (discounted with drq/drc)
- `util_instant_list` / `cost_instant_list` — accrued at point of event (one-off, discounted at event time)
- `util_cycle_list` / `cost_cycle_list` — accrued per cycle; requires companion `_cycle_l`, `_cycle_starttime`, optionally `_max_cycles` variables
- `other_ongoing_list` / `other_instant_list` — same as util but for arbitrary outputs

---

## Simulation Execution

### `run_sim(arm_list, common_all_inputs=NULL, common_pt_inputs=NULL, unique_pt_inputs=NULL, sensitivity_inputs=NULL, init_event_list=NULL, evt_react_list, util_ongoing_list=NULL, util_instant_list=NULL, util_cycle_list=NULL, cost_ongoing_list=NULL, cost_instant_list=NULL, cost_cycle_list=NULL, other_ongoing_list=NULL, other_instant_list=NULL, npats=500, n_sim=1, psa_bool=NULL, sensitivity_bool=FALSE, sensitivity_names=NULL, n_sensitivity=1, input_out=character(), ipd=1, constrained=FALSE, timed_freq=NULL, debug=FALSE, accum_backwards=FALSE, continue_on_error=FALSE, seed=NULL)`

Main simulation function. Returns a nested list: `results[[sensitivity]][[simulation]]`.

Key arguments:
- `psa_bool`: flag passed into model; does NOT automatically draw from distributions. Distributions must be wired in `common_all_inputs` with `ifelse(psa_bool, rnorm(1,...), base_val)` or `input_block()`.
- `n_sim > 1` with `psa_bool=TRUE` → PSA. Each simulation draws fresh parameter values.
- `sensitivity_names`: column names in the data frame passed to `input_block(sens=...)` or `pick_val_v(sens=...)`
- `n_sensitivity`: number of parameters in DSA (auto-detected with `input_block`)
- `input_out`: character vector of item names to export per-patient in IPD output
- `ipd=1`: full IPD; `ipd=2`: aggregated IPD (last value per patient-arm); other: no IPD
- `constrained=TRUE`: enables shared resource simulation (patients interact)
- `seed`: default 1; set explicitly for reproducibility
- `debug=TRUE`: exports log file on error
- `accum_backwards=FALSE` (default): accumulation is forward-looking (current value applied until next change)

**Protected names** — do not use as input variable names:
`arm`, `arm_list`, `categories_for_export`, `cur_evtlist`, `curtime`, `evt`, `i`, `prevtime`, `sens`, `simulation`, `sens_name_used`, `list_env`, `uc_lists`, `npats`, `ipd`

### `run_sim_parallel(...)`
Identical signature to `run_sim()`. Runs PSA simulations in parallel using `future` + `doFuture`.
- More efficient when `n_sim > 5` and each simulation takes >2s
- Expected speedup: 20–40% (not linear in core count due to RAM duplication)
- Requires `future::plan(multisession)` to be set before calling
- Do NOT use with `debug=TRUE` or when still developing/testing

```r
library(future)
plan(multisession, workers=4)
results <- run_sim_parallel(npats=1000, n_sim=500, psa_bool=TRUE, ...)
```

---

## Output Summary Functions

### `summary_results_det(sim_result)` 
Prints deterministic summary table (costs, LYs, QALYs, ICER, ICUR, INMB, discounted and undiscounted).
Input: `results[[1]][[1]]` (single simulation result).

### `summary_results_sim(results_list)`
Aggregates PSA results across simulations.

### `summary_results_sens(results)`
Aggregates across sensitivity analyses (DSA/scenarios).

### `extract_psa_result(results, outcome)`
Extracts a specific outcome across all PSA simulations as a data frame.

### `ceac_des(results, wtp_range=seq(0,100000,1000))`
Computes cost-effectiveness acceptability curve from PSA results.

### `evpi_des(results, wtp_range=seq(0,100000,1000))`
Computes Expected Value of Perfect Information from PSA results.

### `extract_elements_from_list(results, element_name)`
Utility to extract a named element from the nested results list.

---

## Resource Constraint Functions

### `resource_discrete(n, discipline="FIFO", max_queue=Inf, allow_multiple_queue=TRUE)`
Creates a discrete constrained resource with `n` units. Initialise in `common_all_inputs`.
- `discipline`: `"FIFO"` (first-in-first-out) or `"LIFO"` (within same priority level)
- `max_queue`: maximum queue length; patients are rejected (`NA`) when full
- `allow_multiple_queue`: if `FALSE`, a patient already queued is rejected on a second `seize()` attempt

**Core methods:**
- `$n_free()` — units currently available
- `$n_using()` — units currently seized
- `$size()` — total capacity
- `$utilization()` — proportion in use
- `$queue_size()` — patients currently waiting

**Patient-level methods:**
- `$had_to_queue(patient_id)` — 0L/1L: did patient ever wait?
- `$queue_wait_time(patient_id)` — final wait time (after acquisition)
- `$queue_wait_time_current(patient_id, current_time)` — growing wait time while still in queue
- `$time_in_use(patient_id, current_time)` — time since acquisition
- `$is_patient_in_queue(patient_id)` — TRUE/FALSE
- `$is_patient_using(patient_id)` — TRUE/FALSE

**Management methods:**
- `$modify_priority(patient_id, new_priority)` — change queue priority
- `$add_resource(n)` — increase capacity
- `$remove_resource(n, current_time)` — decrease capacity (evicts last-acquired patients)
- `$queue_info(n)` — data frame of first `n` queued patients (id, priority, start_time)
- `$next_patient_in_line()` — patient_id of first in queue

### `seize(resource, amount=1L)`
Attempt to seize `amount` units of a resource for the current patient. Returns `TRUE` (acquired), `FALSE` (queued), `NA` (rejected — queue full or `allow_multiple_queue=FALSE` violation).

### `release(resource, resume_event=NULL, amount=NULL)`
Release resource held by current patient. `amount=NULL` releases 1 unit (default). If `resume_event` is supplied and patients are queued, schedules that event for the next patient in line.

### `release_if_using(resource, resume_event=NULL, amount=NULL)`
Safe release: frees the resource only if the current patient is actually using it. Never errors on a non-using patient. Never removes queue entries. Use in terminal events (death) where the patient may or may not hold the resource.

### `seize_all(resources, policy="all_or_none", amounts=NULL, priorities=NULL, force_unblock=FALSE, accum_queue=TRUE)`
Atomically acquire multiple resources. Returns `TRUE` (all acquired), `FALSE` (queued on bottlenecks), `NA` (rejected — a bottleneck queue is full).
- `policy`: `"all_or_none"` — acquires all or queues for ALL bottlenecks simultaneously; `"sequential"` — acquires in order, stops at first failure
- `amounts`: integer vector of units per resource (default 1 each)
- `priorities`: integer vector of queue priorities per resource (default 0 each)
- `force_unblock`: if `TRUE` and patient is blocked by queue position (not capacity), moves patient to front of those queues and acquires. Use to resolve deadlocks.
- `accum_queue`: if `FALSE`, retries do NOT add new queue entries when the patient is already queued. Prevents phantom queue accumulation on repeated failures. Use `FALSE` for resume-event retry patterns.

### `release_all(resources, resume_event=NULL, amounts=NULL)`
Release multiple resources with all-or-none policy.
- If patient is using ALL listed resources: releases all, purges queue entries (when `amounts=NULL`), and schedules resume events for next-in-line patients.
- If patient is using SOME but not all: does NOT release the used ones (all-or-none). However, when `amounts=NULL`, it DOES purge queue entries for resources the patient is NOT using. This handles the case where a dead patient is queued on some resources but not others.
- If patient is using NONE: does nothing to usage, but `amounts=NULL` still purges all queue entries.
- `amounts=c(...)`: exact per-resource release without queue purge.
- `resume_event`: character vector of event names (one per resource). Deduplicates: if the same patient is first in queue on multiple released resources, only one resume event is scheduled.

**This is the correct death-event cleanup function** when a patient may be using some resources and queued on others. The `amounts=NULL` behavior ensures dead patients are fully removed from all queues.

### `release_all_if_using(resources, resume_event=NULL, amounts=NULL)`
Releases resources only if the patient is currently using ALL listed resources. If patient is not using all (e.g., only queued, or using some but not all), does nothing — no release, no queue purge.

Use this for **non-terminal transitions** where a patient should keep their queue position if they haven't acquired yet. Do NOT use for death events in multi-resource models where the patient might be queued — dead patients will remain in the queue.

**When to use which:**
- Death event (single resource): either works — `release_all_if_using(list(r), ...)` releases if using, does nothing if queued (but with 1 resource, a queued dead patient won't cause issues since they can't acquire after `curtime <- Inf`)
- Death event (multi-resource): use `release_all(list(r1, r2), resume_event=c("e1","e2"))` — purges queues even if patient never acquired
- Non-terminal release (keeping queue position): use `release_all_if_using()`

### `shared_input(initial_value)`
Creates a value shared across all patients within a simulation (not cloned per arm unless defined in `common_all_inputs`).

### `shared_incr(shared_obj, amount)` / `shared_decr(shared_obj, amount)`
Atomically increment/decrement a shared input.

---

## Randomness and Reproducibility

### `random_stream(n)`
Pre-generates `n` uniform random numbers for use as a stream.
Returns an R6-like object with `$draw_n()` method that pops the next value.
Use when you need a sequence of random draws inside event reactions without disrupting the random state.

```r
# In unique_pt_inputs:
rnd_ae <- random_stream(100)
# In evt_react_list:
n_ae <- qpois(rnd_ae$draw_n(), lambda=0.25*duration)
```

### Distribution helpers
- `rbeta_mse(n, mean, se)` / `qbeta_mse(p, mean, se)` — Beta parameterised by mean and SE
- `rgamma_mse(n, mean, se)` / `qgamma_mse(p, mean, se)` — Gamma parameterised by mean and SE
- `rpoisgamma(n, mean, se)` — Poisson-Gamma compound
- `rdirichlet(n, alpha)` / `cond_dirichlet(n, alpha, weights)` — Dirichlet draws
- `cond_mvn(n, mu, sigma, active)` — Conditional multivariate normal

---

## Result Output Structure

`run_sim()` returns: `list[[sensitivity_index]][[simulation_index]]`

Each element contains:
- `$total_lys`, `$total_qalys`, `$total_costs` — named numeric vectors by arm
- `$total_lys_undisc`, `$total_qalys_undisc`, `$total_costs_undisc` — undiscounted
- `$dcosts`, `$dqalys`, `$dlys` — incremental (vs first arm)
- `$ICER`, `$ICUR`, `$INMB` — cost-effectiveness metrics (NA for reference arm)
- `$sensitivity_name` — which sensitivity this is
- `$arm_list` — arm names
- IPD data frame if `ipd=1` or `ipd=2`
- Named outcome columns for any items in `util_*_list`, `cost_*_list`, `other_*_list`, `input_out`

INMB default WTP = £20,000 (UK NICE threshold). Override by setting `wtp` in `common_all_inputs`.

---

## Model Utility Functions

### `extract_from_reactions(evt_react_list)`
Returns a data frame of all assignments and event modifications in the reaction list.
Columns: `event`, `name`, `type` (item/event), `conditional_flag`, `definition`.
Useful for documentation and model review.

### `replicate_profiles(n, ...)` 
Helper to generate patient profiles.

### `sens_iterator(n_sensitivity, sensitivity_names)`
Manual DSA iteration helper (usually not needed when using `input_block`).

### `disc_cycle(value, time, dr)` / `disc_cycle_v()` / `disc_instant(value, time, dr)` / `disc_ongoing(value, from, to, dr)` etc.
Internal discounting functions. Not normally called directly by users but available.
