# CEA Patterns for WARDEN Models

## Standard Model Skeleton

Every WARDEN model follows this skeleton. Fill in the indication-specific sections.

```r
library(WARDEN)
library(dplyr)       # optional but common
library(flexsurv)    # for parametric survival models
options(scipen=999)

# ── 1. INPUTS ──────────────────────────────────────────────────────────────────

# Parameters common across all patients and arms within a simulation
# Set drc/drq here unless they need to vary across scenarios
common_all_inputs <- add_item(input = {
  drc <- 0.035       # cost discount rate (override default 0.03)
  drq <- 0.035       # QALY discount rate
  
  # Pre-draw seeds for reproducibility (one per patient per event type)
  seed_ttp_i  <- sample.int(1e6, npats, replace=FALSE)
  seed_os_i   <- runif(npats)
})

# Parameters that vary per patient but not per arm (e.g. natural death, sex)
common_pt_inputs <- add_item(input = {
  nat_death <- qcond_gompertz(seed_os_i[i], shape=0.10, rate=0.0001, lower_bound=50)
})

# Parameters that vary per patient AND per arm (starting flags, initial utility/cost)
unique_pt_inputs <- add_item(input = {
  fl.on_tx  <- 1L
  fl.prog   <- 0L
  q_default <- util.pfs
  c_default <- cost.drug.int * ifelse(arm=="int", 1, 0) + cost.drug.ref * ifelse(arm=="ref", 1, 0)
})

# ── 2. EVENTS ──────────────────────────────────────────────────────────────────

init_event_list <- add_tte(
  arm  = c("int","ref"),
  evts = c("start","progression","death"),
  input = {
    start       <- 0
    progression <- draw_tte(1, "weibull", coef1=pfs_shape, coef2=log(pfs_scale),
                            beta_tx=ifelse(arm=="int", HR_pfs, 1),
                            seed=seed_ttp_i[i])
    death       <- min(draw_tte(1, "lnorm", coef1=os_mu, coef2=os_sigma,
                                seed=as.integer(seed_os_i[i]*1e6)),
                       nat_death)
  }
)

evt_react_list <-
  add_reactevt("start",       input = {}) |>
  add_reactevt("progression", input = {
    q_default <- util.pps
    c_default <- cost.bsc
    fl.prog   <- 1L
    fl.on_tx  <- 0L
  }) |>
  add_reactevt("death",       input = { curtime <- Inf })

# ── 3. OUTCOMES ────────────────────────────────────────────────────────────────

util_ongoing <- "q_default"
cost_ongoing <- "c_default"

# ── 4. RUN ─────────────────────────────────────────────────────────────────────

results <- run_sim(
  npats             = 1000,
  n_sim             = 1,
  psa_bool          = FALSE,
  arm_list          = c("int","ref"),
  common_all_inputs = common_all_inputs,
  common_pt_inputs  = common_pt_inputs,
  unique_pt_inputs  = unique_pt_inputs,
  init_event_list   = init_event_list,
  evt_react_list    = evt_react_list,
  util_ongoing_list = util_ongoing,
  cost_ongoing_list = cost_ongoing,
  seed              = 42
)

summary_results_det(results[[1]][[1]])
```

---

## Pattern 1 — 2-State Model (On/Off Treatment)

Simple model where patients are either on treatment or not, with a single terminal event.

```
State: [On Treatment] → [Dead]
```

Key WARDEN choices:
- Two events: `start` (time 0), `death`
- `q_default` and `c_default` set at `start`, reset to 0 at `death`
- No `progression` event needed — use arm-specific cost in `unique_pt_inputs`

---

## Pattern 2 — 3-State Partitioned Survival (PFS / PPS / Dead)

Standard late-oncology structure. Two survival curves define three mutually exclusive states.

```
[PFS] → [PPS] → [Dead]
         ↗
[PFS] ─────────→ [Dead]
```

Key WARDEN choices:
- Events: `start`, `progression`, `death`
- `death = min(os_draw, nat_death)` to compete against natural mortality
- On progression: `modify_event` if death needs to be updated (e.g. post-progression survival)
- Separate OS and PFS curves; OS must be ≥ PFS at patient level (`death = max(death, progression + eps)`)
- Utility: `util.pfs` → `util.pps` at progression event
- Costs: drug cost stops at progression (set `fl.on_tx <- 0L`)

```r
evt_react_list <-
  add_reactevt("start",       input = {}) |>
  add_reactevt("progression", input = {
    q_default <- util.pps
    c_default <- cost.bsc
    fl.on_tx  <- 0L
    fl.prog   <- 1L
    # If you have PPS-specific OS, update death:
    # modify_event(c("death" = curtime + draw_tte(1,"lnorm",pps_mu,pps_sigma,seed=seed2[i])))
  }) |>
  add_reactevt("death",       input = { curtime <- Inf })
```

---

## Pattern 3 — Multi-Line Oncology (L1 → L2 → BSC)

```
[L1] → [L2] → [BSC] → [Dead]
```

Key WARDEN choices:
- Events: `start`, `prog1` (L1→L2), `prog2` (L2→BSC), `death`
- Flags: `fl.L1`, `fl.L2`, `fl.BSC`
- At `prog1`: update `c_default` to L2 costs, set `fl.L1 <- 0L`, `fl.L2 <- 1L`
- At `prog2`: update to BSC costs
- Consider: at `prog1`, draw L2 TTP and schedule `prog2` dynamically with `new_event()`
- Time on L2 relative to L1 discontinuation: use `curtime + draw_tte(...)` for `prog2`

---

## Pattern 4 — Tunnel States (Time-in-State Dependent)

When transition probabilities depend on time spent in a state (e.g. post-transplant).

Key WARDEN choices:
- Option A (simple): Use `curtime - state_entry_time` in the reaction to compute duration-dependent costs/utilities
- Option B (TTE-based): Draw multiple events at entry using conditional distributions given time already spent
- Option C (avoid cycles): Use `qtimecov()` with a time-varying hazard function

```r
# Track when patient entered current state:
evt_react_list <-
  add_reactevt("state_entry", input = {
    entry_time <- curtime  # remember when we entered
    # Draw next event conditional on time already spent if needed
  }) |>
  add_reactevt("some_later_event", input = {
    duration_in_state <- curtime - entry_time
    # Use duration to compute costs/utilities
  })
```

---

## Pattern 5 — Chronic Disease with Acute Events

```
[Stable] ←→ [Acute Event] → [Dead]
    ↓
  [Dead]
```

Key WARDEN choices:
- Background stable state: `q_default` and `c_default` set at `start`
- Acute events: use `new_event()` dynamically based on Poisson process
- `random_stream()` for drawing multiple AEs per patient
- After AE: one-off cost/disutility via `cost_instant_list` / `util_instant_list`
- Return to stable state: update back to stable utility/cost

```r
# Pre-generate AE random stream
unique_pt_inputs <- add_item(rnd_ae = random_stream(200))

add_reactevt("start", input = {
  # Schedule first AE
  n_ae_this_period <- qpois(rnd_ae$draw_n(), lambda=rate_ae * horizon)
  if(n_ae_this_period > 0) new_event(c(ae = draw_tte(1,"exp",rate_ae, seed=as.integer(rnd_ae$draw_n()*1e6))))
}) |>
add_reactevt("ae", input = {
  c_instant <- cost.ae           # one-off AE cost
  q_instant <- disutil.ae        # one-off AE disutility (if using util_instant_list)
  # Schedule next AE if still alive
  next_ae <- curtime + draw_tte(1,"exp",rate_ae, seed=as.integer(rnd_ae$draw_n()*1e6))
  if(next_ae < get_event("death")) new_event(c(ae = next_ae))
})
```

---

## Pattern 6 — Using `qtimecov` to Avoid Cycles

When hazard is time-varying (e.g. age-dependent mortality), use `qtimecov` instead of yearly cycles.
Reduces event count dramatically — can cut runtime by 60–80%.

```r
# Define rate as a function of time
common_pt_inputs <- add_item(input = {
  luck      <- runif(1)
  bs_age    <- 60
  # Rate increases with age
  mort_rate <- function(.t) { gompertz_rate(bs_age + .t) }
})

init_event_list <- add_tte(
  arm="int", evts=c("start","death"), input={
    start <- 0
    death <- qtimecov(luck=luck, a_fun=mort_rate, dist="exp", dt=0.5)$tte
  }
)
```

Use `adj_val()` for time-varying utilities (e.g. age utility multipliers):
```r
add_reactevt("start", input = {
  # Gets a scalar that gives correctly discounted outcomes without cycling:
  adj <- adj_val(curtime, next_event()$time, by=1,
                 u_age_vec[min(100, floor(.time + bs_age))],
                 discount=drq)
  q_default <- util.base * adj
})
```

**Caveat**: `next_event()$time` is only valid in forward-accumulation unconstrained models. For constrained models (or `accum_backwards=TRUE`), use `adj_val(prevtime, curtime, ...)` instead — see "Forward vs Backward Accumulation" section below.

---

## Pattern 7 — IPD (Individual Patient Data) Model

When trial IPD is available, patient-specific parameters replace population draws.

```r
# IPD stored as a data frame with columns: pat_id, TTP, OS, arm, ...
IPD <- read.csv("trial_ipd.csv")

common_pt_inputs <- add_item(input = {
  # Pull patient-specific values directly from IPD
  pt_ttp <- IPD$TTP[i]
  pt_os  <- IPD$OS[i]
  pt_arm <- IPD$arm[i]   # note: arm is still set by WARDEN
})
```

---

## Cost Structure Conventions

| Cost category | WARDEN type | Variable pattern |
|---------------|-------------|------------------|
| Annual drug cost | ongoing | `cost.drug.annual` (model auto-multiplies by time) |
| Per-administration cost | instant | `cost.admin` via `cost_instant_list` |
| Annual disease management | ongoing | `cost.dm.state` |
| One-off AE cost | instant | `cost.ae` via `cost_instant_list` |
| Per-cycle cost (e.g., annual visit) | cycle | `cost.visit` via `cost_cycle_list` + `cost.visit_cycle_l=1` |
| Terminal care cost | instant | At death event |

**Critical**: Do NOT mix ongoing and instant for the same cost — pick one and be consistent.
**Double-counting check**: If drug cost is in `unique_pt_inputs` as `c_default` AND also added at an event, it will be double-counted.

---

## Utility Structure Conventions

- EQ-5D values: typically 0–1, from published tariffs (UK EQ-5D-3L/5L, country-specific)
- `q_default` is the ongoing utility (accumulates over time) — set at `start`, update at each state transition
- Disutility of AEs: add as `util_instant_list` variable (e.g. `disutil.ae <- -0.05`)
- Age adjustment: either adjust utility values by age at model start, or use `adj_val()` for dynamic adjustment
- Utility cannot exceed 1 or be below 0 in most HTA contexts — validate parameter inputs
- Dead state: set `q_default <- 0` and `c_default <- 0` before `curtime <- Inf`

---

## Discount Rate Defaults by HTA Body

| HTA Body | Country | Costs | QALYs | Notes |
|----------|---------|-------|-------|-------|
| NICE | England | 3.5% | 3.5% | NHS+PSS perspective |
| CADTH | Canada | 1.5% | 1.5% | Healthcare payer perspective |
| IQWiG | Germany | 3.0% | 3.0% | Societal perspective available |
| HAS | France | 2.5% | 2.5% | Varies by submission type |
| SMC | Scotland | 3.5% | 3.5% | Follows NICE methodology |
| PBAC | Australia | 5.0% | 5.0% | Standard for both |
| ZIN | Netherlands | 4.0% | 1.5% | Different rates for costs/effects |
| TLV | Sweden | 3.0% | 3.0% | Societal perspective |

Set in `common_all_inputs`:
```r
drc <- 0.035   # costs
drq <- 0.035   # QALYs
```

---

## PSA Setup Pattern

```r
# In common_all_inputs with input_block():
common_all_inputs <- input_block(
  base      = param_df$base,
  psa       = pick_psa(param_df$dist, param_df$n, param_df$a, param_df$b),
  names_out = param_df$name,
  psa_indicators = param_df$psa_flag  # 0/1: which params vary in PSA
) |> add_item(drc=0.035, drq=0.035)

# Run PSA:
results_psa <- run_sim_parallel(
  npats    = 1000,
  n_sim    = 500,   # 1000+ for publication; 200 for development
  psa_bool = TRUE,
  seed     = 42,
  ...
)
```

Distribution conventions for PSA:
- Utility (bounded 0-1): `rbeta_mse(1, mean, se)` — parameterised by mean and SE
- Cost (positive): `rgamma_mse(1, mean, se)` — parameterised by mean and SE
- Log-HR: `rnorm(1, log(HR), se_log_HR)` — then `exp()` to get HR
- Probability (bounded): `rbeta_mse()` or `rdirichlet()` for multiple states
- Rate: `rgamma_mse()` or `rlnorm(1, log(rate), se_log)`

---

## DSA / Scenario Setup Pattern

```r
# param_df must have columns: name, base, DSA_low, DSA_high
common_all_inputs <- input_block(
  base      = param_df$base,
  psa       = pick_psa(param_df$dist, param_df$n, param_df$a, param_df$b),
  sens      = param_df,                    # the full data frame
  names_out = param_df$name,
  dsa_names = c("DSA_low", "DSA_high")    # column names
)

results_dsa <- run_sim(
  ...,
  sensitivity_bool = TRUE,
  sensitivity_names = c("DSA_low","DSA_high"),
  # n_sensitivity auto-detected from input_block
  seed = 42
)

summary_results_sens(results_dsa)
```

---

## Pattern 8 — Resource-Constrained Model

When patients compete for limited capacity (ICU beds, specialist appointments, infusion chairs, organs).

```
[Waiting] → [Treated/Using Resource] → [Released] → [Dead]
     ↓ (rejected/queue full)
   [Dead faster]
```

Key WARDEN choices:
- `constrained = TRUE` in `run_sim()`
- `accum_backwards = TRUE` (required — next event is unpredictable due to queuing)
- `resource_discrete()` in `common_all_inputs` (shared across patients within arm)
- `seize()` / `seize_all()` to acquire; `release()` / `release_all_if_using()` to free
- Death event must clean up: `release_all_if_using()` before `curtime <- Inf`
- Use `adj_val(prevtime, curtime, ...)` for time-varying adjustments (backward-looking, both times known)

```r
common_all_inputs <- add_item(input = {
  drc <- 0.035; drq <- 0.035
  beds <- resource_discrete(50, discipline = "FIFO", max_queue = Inf)
  seed_os_i <- runif(npats)
  seed_tx_i <- sample.int(1e6, npats, replace = FALSE)
})

common_pt_inputs <- add_item(input = {
  nat_death <- qcond_gompertz(seed_os_i[i], shape = 0.08, rate = 0.0002, lower_bound = 65)
})

unique_pt_inputs <- add_item(input = {
  fl.has_bed <- 0L
  q_default  <- util.waiting
  c_default  <- cost.waiting
})

init_event_list <- add_tte(
  arm = c("int", "ref"), evts = c("arrival", "death"), input = {
    arrival <- 0
    death   <- nat_death
  }
)

evt_react_list <-
  add_reactevt("arrival", input = {
    acquired <- seize(beds)
    if (isTRUE(acquired)) {
      fl.has_bed <- 1L
      q_default  <- util.treated
      c_default  <- cost.treated
      modify_event(c("discharge" = curtime + draw_tte(1, "exp", coef1 = rate.los,
                                                      seed = seed_tx_i[i])))
    }
    # If FALSE (queued), patient waits — resume_bed event will fire when bed freed
  }) |>
  add_reactevt("resume_bed", input = {
    acquired <- seize_all(list(beds), accum_queue = FALSE)
    if (isTRUE(acquired)) {
      fl.has_bed <- 1L
      q_default  <- util.treated
      c_default  <- cost.treated
      modify_event(c("discharge" = curtime + draw_tte(1, "exp", coef1 = rate.los,
                                                      seed = seed_tx_i[i])))
    }
  }) |>
  add_reactevt("discharge", input = {
    fl.has_bed <- 0L
    q_default  <- util.discharged
    c_default  <- cost.followup
    release(beds, resume_event = "resume_bed")
  }) |>
  add_reactevt("death", input = {
    q_default <- 0; c_default <- 0
    release_all_if_using(list(beds), resume_event = c("resume_bed"))
    curtime <- Inf
  })

results <- run_sim(
  npats = 500, n_sim = 1, psa_bool = FALSE,
  arm_list = c("int", "ref"),
  common_all_inputs = common_all_inputs,
  common_pt_inputs = common_pt_inputs,
  unique_pt_inputs = unique_pt_inputs,
  init_event_list = init_event_list,
  evt_react_list = evt_react_list,
  util_ongoing_list = "q_default",
  cost_ongoing_list = "c_default",
  constrained = TRUE,
  accum_backwards = TRUE,
  seed = 42
)
```

### How the resource lifecycle works

The constrained engine processes events from a **shared queue** across all patients. When a patient fires an event, WARDEN evaluates the reaction, then accumulates costs/utilities for the elapsed interval `(prevtime, curtime)`.

**Step-by-step lifecycle:**

1. **Arrival**: patient calls `seize(resource)`.
   - If capacity available: returns `TRUE`. Patient holds the resource. Set flags/costs/utilities for the "using" state.
   - If no capacity: returns `FALSE`. Patient is added to the resource queue. Their state (utility, cost) remains as initialized — they are "waiting". No event is scheduled; the patient is effectively frozen until another patient frees a unit.

2. **Release by another patient**: when a using patient calls `release(resource, resume_event = "X")`, WARDEN:
   - Frees 1 unit of capacity
   - Finds the next patient in queue (per discipline: FIFO = longest-waiting, LIFO = most-recent)
   - Schedules event "X" for that queued patient at the current time

3. **Resume event fires for queued patient**: the queued patient receives event "X". In this reaction, call `seize_all(list(resource), accum_queue = FALSE)`:
   - `accum_queue = FALSE` is critical: it prevents adding a duplicate queue entry. The patient is already queued from step 1; this call just checks availability and acquires without re-queuing.
   - If capacity is available (it should be, since a unit was just freed): returns `TRUE`. Set the "using" state.
   - If another patient grabbed it first (race condition in multi-resource): returns `FALSE`. Patient stays queued — will get another resume when another unit frees.

4. **Death (terminal event)**: patient may be using, queued, or neither.
   - **Single resource**: `release_all_if_using(list(resource), resume_event = c("X"))` is safe — releases if using; if only queued, the dead patient stays in queue but the engine skips events for patients with `curtime = Inf`, so it's harmless (the stale queue entry is eventually overwritten).
   - **Multi-resource** (patient might be using one and queued on another): use `release_all(list(r1, r2), resume_event = c("e1", "e2"))` — this purges queue entries for resources the patient is NOT using, preventing stale queue positions that block other patients from advancing.
   - MUST come before `curtime <- Inf`.

### Why `accum_queue = FALSE` on retry

Without it, each resume attempt that fails (because another patient grabbed the resource first) adds a NEW queue entry. After 3 failed retries, the patient has 4 queue positions. When they finally acquire, 3 phantom entries remain. These trigger extra resume events for OTHER patients, causing cascading errors.

Detection: if IPD shows the same patient receiving `resume_bed` multiple times in sequence, or if `queue_size()` grows unexpectedly, phantom entries are likely present.

### Choosing queue discipline

| Discipline | Behavior | When to use |
|-----------|----------|-------------|
| `"FIFO"` | First patient to queue is first served | Default. Fair queuing: waiting lists, transplant registers, routine appointments |
| `"LIFO"` | Last patient to queue is first served | Stack-based: emergency overflow where most recent arrival is most critical |

For priority-based queuing (e.g., sicker patients served first), use `seize_all(..., priorities = c(priority_value))`. Lower integer = higher priority (served first). Within the same priority, discipline (FIFO/LIFO) breaks ties.

### `max_queue` and rejection

When `max_queue` is finite and the queue is full:
- `seize()` returns `NA` (rejected)
- Patient is NOT queued — they never enter the wait list
- Handle with: `if (is.na(acquired)) { ... }` (e.g., divert to alternative care, accelerated decline)

Use finite `max_queue` when: physically limited waiting rooms, organ transplant lists with size caps, models where "turned away" is a distinct clinical pathway.

### Time-varying hazards after queue release (`qtimecov` restart)

When a patient waits in a queue and their baseline hazard is time-varying (e.g., age-dependent mortality), the TTE drawn at arrival becomes stale. On resume, recalculate with `qtimecov()` using `start_time = curtime` and the saved luck value:

```r
# In common_pt_inputs:
common_pt_inputs <- add_item(input = {
  luck_mort <- runif(1)
  mort_result <- qtimecov(luck = luck_mort,
                          a_fun = function(.t) mort_rate_by_age(baseline_age + .t),
                          dist = "exp", dt = 1, max_time = 60)
  nat_death <- mort_result$tte
  luck_mort_remaining <- mort_result$luck
})

# In resume event (after queue wait):
add_reactevt("resume_bed", input = {
  acquired <- seize_all(list(beds), accum_queue = FALSE)
  if (isTRUE(acquired)) {
    fl.has_bed <- 1L
    q_default <- util.treated
    # Recalculate mortality from current time with remaining luck:
    mort_result <- qtimecov(luck = luck_mort_remaining,
                            a_fun = function(.t) mort_rate_by_age(baseline_age + .t),
                            dist = "exp", dt = 1, max_time = 60,
                            start_time = curtime)
    modify_event(c("death" = mort_result$tte))
    luck_mort_remaining <- mort_result$luck
  }
})
```

The `luck` value preserves the patient's survival rank: a patient with "good luck" (high quantile) stays lucky regardless of queue duration.

### `shared_input()` for cross-patient state

When patients affect each other beyond resource competition (e.g., cumulative hospital admissions, infection prevalence, budget impact tracking):

```r
common_all_inputs <- add_item(input = {
  beds <- resource_discrete(50, discipline = "FIFO")
  total_admissions <- shared_input(0L)   # shared counter across all patients
  cumulative_cost  <- shared_input(0)    # shared cost tracker
})

# In admission event:
add_reactevt("admission", input = {
  shared_incr(total_admissions, 1L)
  # Use shared state for capacity decisions:
  if (total_admissions$value > 200L) {
    # Hospital overflow protocol...
  }
})
```

`shared_input` values are visible to all patients within the same arm-simulation. They are NOT cloned per patient — modifications by one patient are immediately visible to others.

---

## Pattern 9 — Treatment Delivery Comparison (IV vs SC)

When the decision problem is route of administration (infusion vs subcutaneous injection) with differences in time, resource use, and patient experience.

```
[Admin Event] → resource seized (chair/nurse) → [Admin Complete] → release → ... repeat
```

Key WARDEN choices:
- Infusion chairs as `resource_discrete(n_chairs)` — only for IV arm
- SC arm has no resource constraint (self-administered or brief nurse time)
- Track total admin time per patient for utility of time saved
- Costs: drug acquisition (may differ), admin cost per session, chair time cost
- Frequency: IV every N weeks vs SC every M weeks

```r
common_all_inputs <- add_item(input = {
  chairs    <- resource_discrete(10)  # infusion chairs, IV arm only
  freq.iv   <- 1/3   # every 3 weeks
  freq.sc   <- 1/2   # every 2 weeks
  dur.iv    <- 2/52  # 2 hours in years
  dur.sc    <- 0     # negligible time
  cost.iv.admin  <- 350
  cost.sc.admin  <- 50
  cost.drug.iv   <- 5000
  cost.drug.sc   <- 5200
  util.base      <- 0.82
  # Time disutility: patient loses util during infusion time
  disutil.chair_time <- 0.01  # per-session disutility for time lost
})

unique_pt_inputs <- add_item(input = {
  q_default <- util.base
  c_default <- if (arm == "iv") cost.drug.iv * freq.iv else cost.drug.sc * freq.sc
})

# IV arm: model each admin as resource event
# SC arm: no infusion events needed — costs/utils captured via ongoing defaults
evt_react_list <-
  add_reactevt("admin_iv", input = {
    if (arm == "iv") {
      acquired <- seize(chairs)
      if (isTRUE(acquired)) {
        modify_event(c("admin_complete" = curtime + dur.iv))
      }
    }
  }) |>
  add_reactevt("admin_complete", input = {
    release(chairs, resume_event = "admin_iv")
    # Schedule next admin
    next_admin <- curtime + 1/freq.iv
    if (next_admin < get_event("death")) modify_event(c("admin_iv" = next_admin))
  }) |>
  add_reactevt("death", input = {
    q_default <- 0; c_default <- 0
    if (arm == "iv") release_all_if_using(list(chairs), resume_event = c("admin_iv"))
    curtime <- Inf
  })
```

---

## Pattern 10 — Multi-Line Treatment Sequencing (L1 → L2 → L3 → BSC)

When patients progress through multiple treatment lines with conditional distributions.

```
[L1] → progression → [L2] → progression → [L3] → progression → [BSC] → [Dead]
```

Key WARDEN choices:
- Track current line with integer flag: `fl.line`
- Named vector dispatch for line-specific costs, utilities, and parameters
- At progression: draw next-line TTE with `new_event()` dynamically
- BSC is terminal treatment state (only death event remaining)
- Treatment effect may depend on prior line duration or response

```r
common_all_inputs <- add_item(input = {
  max_lines <- 3L
  util_by_line <- c("1" = 0.82, "2" = 0.70, "3" = 0.55, "bsc" = 0.40)
  cost_by_line <- c("1" = 15000, "2" = 12000, "3" = 8000, "bsc" = 3000)
  pfs_shape_by_line <- c("1" = 1.2, "2" = 1.1, "3" = 0.9)
  pfs_scale_by_line <- c("1" = log(2.5), "2" = log(1.8), "3" = log(1.0))
  seed_prog_i <- sample.int(1e6, npats, replace = FALSE)
})

unique_pt_inputs <- add_item(input = {
  fl.line   <- 1L
  q_default <- util_by_line[["1"]]
  c_default <- cost_by_line[["1"]] * ifelse(arm == "int", 1, 0.8)
})

init_event_list <- add_tte(
  arm = c("int", "ref"), evts = c("start", "progression", "death"), input = {
    start       <- 0
    progression <- draw_tte(1, "weibull",
                            coef1 = pfs_shape_by_line[["1"]],
                            coef2 = pfs_scale_by_line[["1"]],
                            beta_tx = ifelse(arm == "int", HR_L1, 1),
                            seed = seed_prog_i[i])
    death       <- min(nat_death, 40)
  }
)

evt_react_list <-
  add_reactevt("start", input = {}) |>
  add_reactevt("progression", input = {
    fl.line <- fl.line + 1L
    if (fl.line <= max_lines) {
      line_key <- as.character(fl.line)
      q_default <- util_by_line[[line_key]]
      c_default <- cost_by_line[[line_key]]
      # Draw next-line PFS
      next_prog <- curtime + draw_tte(1, "weibull",
                                      coef1 = pfs_shape_by_line[[line_key]],
                                      coef2 = pfs_scale_by_line[[line_key]],
                                      seed = seed_prog_i[i] + fl.line)
      if (next_prog < get_event("death")) {
        new_event(c("progression" = next_prog))
      }
    } else {
      # BSC — no further treatment
      q_default <- util_by_line[["bsc"]]
      c_default <- cost_by_line[["bsc"]]
    }
  }) |>
  add_reactevt("death", input = { curtime <- Inf })
```

---

## Pattern 11 — Chronic Progressive with Complications

When patients accumulate complications over time (diabetes, heart failure, obesity with comorbidities). Multiple event types as competing risks, each permanently worsening outcomes.

```
[Baseline] → MI → [Post-MI] → Stroke → [Post-MI+Stroke] → ... → [Dead]
```

Key WARDEN choices:
- Multiple complication events drawn as competing risks at baseline
- Each event modifies utility/cost AND may draw the next occurrence
- `qtimecov()` for age/covariate-dependent baseline mortality
- Utility decrements accumulate: `q_default <- util.base * (1 - sum(decrements))`
- Use `random_stream()` when number of recurring events per patient is unknown

```r
common_all_inputs <- add_item(input = {
  rate.mi    <- 0.02    # annual MI rate
  rate.stroke <- 0.015  # annual stroke rate
  disutil.mi  <- 0.05   # permanent utility decrement per MI
  disutil.stroke <- 0.10
  cost.mi.acute  <- 15000
  cost.stroke.acute <- 20000
  seed_mort_i <- runif(npats)
})

common_pt_inputs <- add_item(input = {
  bs_age <- 55
  luck_mort <- seed_mort_i[i]
  mort_rate_fn <- function(.time) { 0.001 * exp(0.08 * (.time + bs_age - 50)) }
})

unique_pt_inputs <- add_item(input = {
  rnd_stream  <- random_stream(50)
  n_mi        <- 0L
  n_stroke    <- 0L
  q_default   <- 0.85
  c_default   <- 2000  # annual disease management
})

init_event_list <- add_tte(
  arm = c("int", "ref"), evts = c("start", "mi", "stroke", "death"), input = {
    start  <- 0
    mi     <- qexp(rnd_stream$draw_n(), rate.mi * ifelse(arm == "int", HR.mi, 1))
    stroke <- qexp(rnd_stream$draw_n(), rate.stroke * ifelse(arm == "int", HR.stroke, 1))
    death  <- qtimecov(luck = luck_mort, a_fun = mort_rate_fn, dist = "exp", dt = 1)$tte
  }
)

evt_react_list <-
  add_reactevt("start", input = {}) |>
  add_reactevt("mi", input = {
    n_mi <- n_mi + 1L
    q_default <- max(0, 0.85 - n_mi * disutil.mi - n_stroke * disutil.stroke)
    c_default <- 2000 + n_mi * 500  # ongoing post-MI management
    c_instant <- cost.mi.acute
    # Schedule next MI (recurring)
    next_mi <- curtime + qexp(rnd_stream$draw_n(), rate.mi * (1 + 0.5 * n_mi))
    if (next_mi < get_event("death")) modify_event(c("mi" = next_mi))
  }) |>
  add_reactevt("stroke", input = {
    n_stroke <- n_stroke + 1L
    q_default <- max(0, 0.85 - n_mi * disutil.mi - n_stroke * disutil.stroke)
    c_default <- 2000 + n_mi * 500 + n_stroke * 1000
    c_instant <- cost.stroke.acute
    next_stroke <- curtime + qexp(rnd_stream$draw_n(), rate.stroke * (1 + 0.3 * n_stroke))
    if (next_stroke < get_event("death")) modify_event(c("stroke" = next_stroke))
  }) |>
  add_reactevt("death", input = { curtime <- Inf })
```

---

## Forward vs Backward Accumulation

WARDEN supports two accumulation modes for ongoing values (`q_default`, `c_default`):

| Mode | `run_sim()` setting | `adj_val` pattern | When to use |
|------|--------------------|--------------------|-------------|
| **Forward** (default) | `accum_backwards = FALSE` | `adj_val(curtime, next_event()$time, ...)` | Standard unconstrained models where next event time is known |
| **Backward** | `accum_backwards = TRUE` | `adj_val(prevtime, curtime, ...)` | Resource-constrained models, or any model where events can be pre-empted |

**Forward**: at each event, set `q_default` to the value that should accrue going forward until the next event. Works perfectly when you know when the next event fires (typical in unconstrained models). `adj_val(curtime, next_event()$time, ...)` is valid because no external event will interrupt.

**Backward**: at each event, set `q_default` to the value that accrued during the interval that just completed (`prevtime` → `curtime`). Both times are known with certainty. Required for constrained models because resource-release events can fire at any time, making `next_event()$time` unreliable.

**Performance**: no runtime cost difference. Choose based on correctness, not speed.

---

## Random Number Management

### Strategy 1: Pre-drawn seeds + conditional quantile (preferred for TTEs)

Pre-generate one uniform draw per patient per event type. Use `qcond_*()` or `draw_tte(seed=...)` for reproducible TTEs.

```r
common_all_inputs <- add_item(input = {
  seed_os_i  <- runif(npats)          # one per patient for OS
  seed_pfs_i <- sample.int(1e6, npats, replace = FALSE)  # integer seeds for draw_tte
})
common_pt_inputs <- add_item(input = {
  death <- qcond_gompertz(seed_os_i[i], shape = 0.08, rate = 0.0002, lower_bound = 60)
})
```

### Strategy 2: `random_stream()` for variable-count draws

When the number of random values needed per patient is unknown at input time (recurring AEs, Poisson events, dynamic decisions).

```r
unique_pt_inputs <- add_item(rnd_ae = random_stream(200))
# In reactions:
next_ae_time <- curtime + qexp(rnd_ae$draw_n(), rate = lambda_ae)
```

### Strategy 3: Never inline `r*()` in reactions

Inline `rexp()`, `rnorm()`, etc. inside `add_reactevt` blocks disrupt WARDEN's per-patient random state management. In constrained mode this causes non-reproducible results. Always use one of the above strategies instead.

### When to use each:
- **Fixed number of TTEs per patient** → Strategy 1 (pre-drawn seeds)
- **Variable number of events** (AEs, hospitalizations, recurring) → Strategy 2 (`random_stream`)
- **Constrained model** → Strategy 1 or 2 are both safe; Strategy 3 (inline) is NEVER safe
