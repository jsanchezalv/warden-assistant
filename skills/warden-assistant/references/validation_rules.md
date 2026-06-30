# WARDEN Model Validation Rules

## How to Run Validation

**Static pass**: Read the entire model code. Apply all categories below. Note each issue with: line reference, what's wrong, why it matters, fix.

**Runtime pass**: Execute the model and parse output:
```bash
Rscript model.R 2>&1
```
Check for: R warnings (not just errors), NaN in results, unexpected Inf costs/QALYs, ICER that is wildly outside clinical plausibility.

---

## Category A — API Misuse [BLOCKING]

### A1: Protected variable names used as inputs
**Check**: scan `add_item`, `add_tte`, `add_reactevt` for: `arm`, `curtime`, `prevtime`, `evt`, `i`, `simulation`, `sens`, `npats`, `ipd`, `cur_evtlist`, `list_env`, `uc_lists`, `arm_list`, `sens_name_used`, `categories_for_export`
**Why it matters**: silently overwrites WARDEN internal state, causing unpredictable simulation behaviour
**Fix**: rename any clashing variables

### A2: `modify_event` called on non-existent or already-fired events
**Check**: any `modify_event(c("X" = ...))` where event "X" is not declared in `add_tte(evts=...)` or created via `new_event()` before this point
**Why it matters**: runtime error that terminates the simulation
**Fix**: use `modify_event(..., create_if_missing=TRUE)` or `new_event()` first; check event exists with `has_event()`

### A3: `curtime` not set to `Inf` in terminal event
**Check**: the terminal event (usually "death" or "os") must contain `curtime <- Inf`
**Why it matters**: simulation runs forever (or until some very large time), producing incorrect results
**Fix**: add `curtime <- Inf` to the terminal event reaction

### A4: Cost/utility variables of length > 1 passed to outcome lists
**Check**: variables in `util_ongoing_list`, `cost_ongoing_list`, etc. should always be length 1 at the time of accumulation
**Why it matters**: WARDEN expands data per element, creating multiple rows per event — results are then wrong without post-processing
**Fix**: extract the single relevant value before assignment (e.g., `q_default <- util_vec[state]` not `q_default <- util_vec`)

### A5: `add_tte` used for the same event name in conflicting arm calls
**Check**: if `add_tte(arm=c("int","ref"), ...)` and then a separate `add_tte(arm="int", ...)` both define the same event, the second overwrites the first for that arm
**Why it matters**: silent event time overwrite with no error
**Fix**: consolidate into one call or ensure second call uses `other_inp` only

### A6: `draw_tte` coef2 not on log scale for Weibull/lnorm/llogis
**Check**: `draw_tte(dist="weibull", coef2=scale_val)` — coef2 should be the log of the scale parameter for Weibull (as in flexsurv parameterisation), not the scale directly. For lnorm, coef2 is the log-SD.
**Why it matters**: TTEs will be systematically wrong
**Fix**: use `coef2=log(scale)` for Weibull; for lnorm `coef1=log(mean_t)` and `coef2=log_sd`

### A7: `n_sim > 1` with `psa_bool = FALSE`
**Check**: if `n_sim > 1` and `psa_bool = FALSE`, no warning is produced but results vary across simulations only due to patient sampling
**Why it matters**: user may interpret this as PSA when it's not; parameter uncertainty is not captured
**Fix**: set `psa_bool = TRUE` for PSA, or `n_sim = 1` for deterministic

### A8: Cycle-list variable missing `_cycle_l` or `_cycle_starttime` companion
**Check**: any variable in `cost_cycle_list` or `util_cycle_list` must have corresponding `{name}_cycle_l` (cycle length) and `{name}_cycle_starttime` variables accessible in the model environment
**Why it matters**: runtime error or silent wrong discounting
**Fix**: declare `my_cost_cycle_l <- 1` and `my_cost_cycle_starttime <- 0` (or the appropriate start time) in inputs

---

## Category B — DES Logic Errors [BLOCKING]

### B1: Events defined in `evts` but timing depends on events not yet scheduled
**Check**: in `add_tte` input block, if event B's time depends on event A's time, but A isn't computed yet in the expression
**Why it matters**: A will be `Inf` when referenced; B gets wrong time
**Fix**: compute A before B in the `input` block (order within the block matters)

### B2: Terminal event may never fire
**Check**: is there a code path where `curtime` is never set to `Inf`? E.g. a conditional `if(condition){ curtime <- Inf }` where condition might never be TRUE
**Why it matters**: simulation runs indefinitely or until a very large number
**Fix**: ensure terminal event is unconditionally terminal; add a fallback horizon event

### B3: `new_event` scheduling creates circular event chains
**Check**: event A schedules event B, event B schedules event A — without a termination condition
**Why it matters**: infinite loop, simulation never completes
**Fix**: add a guard condition before scheduling the next event (check `curtime < death_time`, flag, etc.)

### B4: Competing events not properly handled
**Check**: if patient can die from cause A OR cause B, is `death = min(cause_A, cause_B)` computed? Or are two separate death events in the queue?
**Why it matters**: if two competing events both set `curtime <- Inf`, whichever fires first wins, but the second may fire if queue handling is wrong
**Fix**: use `min()` in `add_tte` input to pick the minimum, or use a single death event that is updated as causes change

### B5: `get_event` called on event that may not exist
**Check**: `get_event("X")` inside a reaction where event "X" may have been removed or never scheduled for some patients
**Why it matters**: runtime error
**Fix**: wrap with `if(has_event("X")) { ... }` first

### B6: Random draws in event reactions that should be pre-drawn
**Check**: calls to `rexp()`, `rnorm()`, `runif()`, `rbeta()`, etc. inside `add_reactevt` blocks (without `random_stream`)
**Why it matters**: in unconstrained mode, random state is managed per patient-arm — inline random calls within reactions can be inconsistent across PSA runs. In constrained mode, this is a serious reproducibility problem.
**Fix**: pre-draw in `common_pt_inputs` or `unique_pt_inputs` using `runif(1)` seeds, then use `qcond_*()` functions; or use `random_stream()` object

### B7: `modify_event` called after the event has already fired
**Check**: inside a reaction for event "X", attempting to modify event "X" itself
**Why it matters**: runtime error — you cannot modify an event currently being processed
**Fix**: modify OTHER events from within "X"; never self-reference

---

## Category C — Health Economics Errors [HIGH]

### C1: Double-counted costs
**Check**: is the same cost appearing in both `c_default` (ongoing) AND `cost_instant_list`? Or in `unique_pt_inputs` AND modified in an event reaction?
**Why it matters**: inflated cost estimates — could double total drug costs
**Fix**: choose one accumulation mechanism per cost category; document in comments

### C2: Wrong discount rates
**Check**: are `drc` and `drq` defined in the model? If not, the default is 3.0% for both.
**Check**: are rates appropriate for the target HTA body? (see cea_patterns.md for body-specific rates)
**Why it matters**: incorrect discounting leads to wrong ICERs

### C3: Costs accumulated after terminal event
**Check**: are `c_default` and `q_default` set to 0 BEFORE `curtime <- Inf` in the terminal event?
**Why it matters**: if set after, the model may have already closed accumulation. Convention: set to 0, then set `curtime <- Inf`
**Fix**: `q_default <- 0; c_default <- 0; curtime <- Inf`

### C4: Time horizon insufficient for condition
**Check**: does the model run until patient death (natural mortality included)? Or is it artificially capped?
**Check**: for chronic conditions, lifetime horizon expected by most HTA bodies
**Why it matters**: truncation of costs and outcomes biases ICER

### C5: Drug costs not stopping at discontinuation
**Check**: when a patient discontinues treatment (`fl.on_tx <- 0L`), is `c_default` updated to remove the drug cost?
**Why it matters**: ongoing drug cost continues accumulating after patient is off treatment

### C6: Utility not accounting for dead time
**Check**: at death, `q_default` should be set to 0 (it usually is via the terminal event reaction, but verify)
**Why it matters**: if patient dies mid-period and utility isn't zeroed, time after death contributes to QALYs (very small effect but technically wrong)

### C7: Natural mortality not competing with disease-specific death
**Check**: is `nat_death` drawn in `common_pt_inputs` and taken as `min()` with disease-specific OS?
**Why it matters**: without competing natural mortality, patients never die of non-disease causes — inflates LYs/QALYs
**Fix**: `death <- min(disease_os, nat_death)` in `add_tte`

### C8: Perspective mismatch
**Check**: if NHS perspective is stated, are indirect costs, productivity losses, or patient travel costs excluded?
**Check**: if societal perspective, are productivity losses included?
**Why it matters**: including costs from wrong perspective invalidates submission

---

## Category D — PSA / DSA Issues [HIGH]

### D1: `psa_bool = TRUE` but distributions not wired in inputs
**Check**: does `common_all_inputs` actually draw from distributions when `psa_bool=TRUE`? Or do all parameters use fixed values regardless?
**Why it matters**: PSA runs but all simulations return identical results
**Fix**: use `input_block()`, or `ifelse(psa_bool, rnorm(1,...), base_value)` for each parameter

### D2: Beta distribution parameterisation wrong
**Check**: `rbeta_mse(n, mean, se)` expects mean and SE (standard error), not alpha/beta
**Check**: if using `rbeta(n, alpha, beta)` directly, verify alpha/beta are computed from moments correctly: `alpha = mean*(mean*(1-mean)/var - 1)`, `beta = (1-mean)*(mean*(1-mean)/var - 1)`
**Why it matters**: wrong distribution → wrong PSA uncertainty
**Fix**: use `rbeta_mse()` with mean and SE instead of direct parameterisation

### D3: LogNormal parameterisation on wrong scale
**Check**: `rlnorm(n, meanlog, sdlog)` — both arguments must be ON THE LOG SCALE. If user provides natural-scale mean and SD, must convert.
**Check**: `draw_tte(dist="lnorm", coef1=X, coef2=Y)` — coef1 is meanlog, coef2 is sdlog (both log-scale)
**Fix**: `meanlog = log(natural_mean) - 0.5*log_var`; or use `rlnorm` and ensure inputs are already on log scale

### D4: HR PSA using normal distribution on log scale
**Check**: HR uncertainty should be `log(HR) ~ N(log(HR_est), SE_log_HR)`, not `HR ~ N(HR_est, SE_HR)`
**Why it matters**: normal on natural scale can produce negative HRs; log-normal correct for positive ratios
**Fix**: `hr_psa <- exp(rnorm(1, log(HR_base), se_log_hr))`

### D5: Correlated parameters treated as independent in PSA
**Check**: if two parameters are estimated from the same regression/model (e.g., Weibull shape and scale), they are correlated. Independent draws will understate or overstate joint uncertainty.
**Why it matters**: biased PSA uncertainty bounds
**Fix**: use `MASS::mvrnorm()` with the variance-covariance matrix, or `cond_mvn()`

### D6: Seed not set in `run_sim()`
**Check**: `seed` argument should be set explicitly
**Why it matters**: results not reproducible between runs
**Fix**: `run_sim(..., seed=42)`

### D7: DSA range not symmetric or implausible
**Check**: are DSA ranges (95% CI or ± 25% or similar) appropriate for each parameter?
**Check**: do DSA ranges respect parameter constraints (e.g. utility cannot exceed 1, HR cannot be negative)?
**Why it matters**: DSA with implausible ranges produces misleading tornado charts

---

## Category E — Code Quality [MEDIUM]

### E1: Multiple separate `modify_event` calls that could be grouped
**Check**: two or more `modify_event(c("X" = ...))` calls in the same event reaction block
**Why it matters**: each call has overhead from interacting with the C++ queue; a single call with multiple entries is meaningfully faster in large models
**Fix**: `modify_event(c("X" = val1, "Y" = val2))` — single call, named vector

### E2: `case_when()` or `dplyr::if_else()` used inside reactions
**Check**: any use of `dplyr::case_when()`, `dplyr::if_else()`, or similar tidyverse functions inside `add_reactevt` blocks
**Why it matters**: these load dplyr overhead in a tight loop; `ifelse()` or standard `if`/`else` is faster
**Fix**: replace with `ifelse()` or `if(cond){x}else{y}` 

### E3: Named vectors assigned to `q_default` or `c_default`
**Check**: `q_default <- c(pfs=0.8)` instead of `q_default <- 0.8`
**Why it matters**: named vector of length 1 causes data expansion (category A4 variant)
**Fix**: `q_default <- unname(util_vec[state])` or `q_default <- util_vec[[state]]`

### E4: Magic numbers without named constants
**Check**: bare numbers in event reactions (e.g. `cost.drug <- 85000`, threshold values) without explanation
**Why it matters**: maintainability — someone reviewing the model can't understand assumptions
**Fix**: define named constants in `common_all_inputs` with comments citing source

### E5: Events created dynamically that could be in `init_event_list`
**Check**: `new_event()` calls in reactions for events that fire for every patient
**Why it matters**: pre-declaring in `add_tte` is faster; dynamic creation via `new_event` has overhead
**Fix**: if the event fires for all (or nearly all) patients, put it in `add_tte` with `Inf` default and update via `modify_event`

### E6: `assign()` or `get()` with dynamic names inside reactions
**Check**: patterns like `assign(paste0("fl_", state), 1)` inside event reactions
**Why it matters**: WARDEN's debug logger cannot capture dynamically named assignments; harder to maintain
**Fix**: use static explicit assignments; pre-compute the variable name outside if needed

### E7: `tibble`/`data.frame` operations inside event reactions
**Check**: any call to `tibble()`, `data.frame()`, `dplyr::filter()`, `dplyr::mutate()`, `dplyr::select()`, `subset()`, or `merge()` inside `add_reactevt` blocks
**Why it matters**: data.frame construction and subsetting is expensive (memory allocation, class dispatch). In a DES, event reactions fire thousands to millions of times — even microsecond overhead compounds into minutes.
**Fix**: pre-compute all lookups as named vectors or lists in `common_all_inputs`. Replace `df[df$x == val, "y"]` with `setNames(df$y, df$x)` then `vec[[as.character(val)]]` inside reactions.

### E8: Data.frame row subsetting for single-value lookup
**Check**: patterns like `df[df$col == value, "result"]` or `df %>% filter(x == val) %>% pull(y)` inside event reactions for extracting a single value
**Why it matters**: linear O(n) scan on every lookup. Named vector indexing is O(1) hash lookup, 20-60x faster for tables with 10+ rows. Over millions of event firings this dominates runtime.
**Fix**: in `common_all_inputs`, convert once: `lookup_vec <- setNames(df$result, df$col)`. In reactions: `value <- lookup_vec[[as.character(key)]]`

---

## Category F — Constrained Model Specific [HIGH — if `constrained=TRUE`]

### F1: "r" distribution functions called inside event reactions
**Check**: `rexp()`, `rnorm()`, `rpois()`, etc. directly in `add_reactevt` blocks when `constrained=TRUE`
**Why it matters**: in constrained mode, the event loop processes patients from the shared queue; random state is NOT patient-specific, so inline draws produce non-reproducible results that may differ between `constrained=FALSE` and `constrained=TRUE` runs
**Fix**: pre-draw all random values in `unique_pt_inputs`; use `random_stream()` for streams needed during simulation

### F2: Resource not initialised in `common_all_inputs`
**Check**: `resource_discrete()` should be in `common_all_inputs` (shared within arm, cloned across arms)
**Why it matters**: if in `unique_pt_inputs`, each patient gets their own resource — not shared
**Fix**: move resource creation to `common_all_inputs`

### F3: `release()` called without corresponding `seize()`
**Check**: if a patient can reach the death event without ever having acquired the resource, `release()` will error
**Fix**: use `release_all_if_using(resource)` in the terminal event, which is safe to call even if the patient never acquired

### F4: Dead patients left in resource queues
**Check**: terminal event reaction in a constrained model — does it call `release_all_if_using()` or `release_all()` for every resource the patient could be queued on?
**Why it matters**: dead patients remain in resource queues and get rescheduled when another patient releases, causing phantom events for dead patients and incorrect resource allocation
**Fix**: call `release_all_if_using(list(resource1, resource2, ...), resume_event = c("resume_evt1", "resume_evt2"))` BEFORE `curtime <- Inf` in the terminal event. This purges both usage and queue entries.

### F5: Forward `adj_val` in constrained models
**Check**: any call to `adj_val` using `next_event()$time` or a predicted future time as the second argument when `constrained = TRUE`
**Why it matters**: in resource-constrained models, the next event can be pre-empted by another patient's resource release. The time boundary used at accumulation may not match actual elapsed time, producing wrong utility/cost calculations.
**Fix**: use `accum_backwards = TRUE` in `run_sim()` and `adj_val(prevtime, curtime, ...)` in event reactions. Both boundaries are known with certainty at the time of computation.

### F6: `seize_all()` retry without `accum_queue = FALSE`
**Check**: if a patient can fail `seize_all()` and retry later (e.g., on a resume event triggered by another patient's release), is `accum_queue = FALSE` set on the retry call?
**Why it matters**: without it, each failed retry adds a new queue entry for the same patient on each bottleneck resource, leading to phantom queue positions. When the patient finally acquires, previous extra entries remain and may cause double-scheduling of resume events for other patients.
**Fix**: use `seize_all(resources, accum_queue = FALSE)` on retry attempts (the first call can use the default `TRUE`)

---

## Category G — HTA Submission Risks [CONTEXTUAL]

### G1: No natural mortality
**Flag if**: no general population mortality (Gompertz, life tables) competing with disease-specific OS
**HTA bodies**: NICE, CADTH, IQWiG all expect general population mortality to be included
**Caveat**: acceptable if disease OS already incorporates all-cause mortality from trial data

### G2: Time horizon < lifetime for chronic condition
**Flag if**: time horizon < 40+ years for conditions in patients under 60
**HTA bodies**: most expect lifetime horizon for conditions with long-term mortality impact

### G3: No scenario analysis for structural assumptions
**Flag if**: key structural assumptions (e.g., choice of survival distribution, treatment waning, extrapolation method) not tested in scenarios
**HTA bodies**: NICE expects sensitivity to extrapolation assumptions for oncology submissions

### G4: Single parametric distribution without justification
**Flag if**: one survival distribution used without comparison to alternatives (AIC/BIC, clinical plausibility check)
**HTA bodies**: NICE TSD 14 guidance requires considering multiple distributions

### G5: Utility source not UK EQ-5D for NICE submissions
**Flag if**: utilities from non-EQ-5D instrument, or non-UK value set
**HTA bodies**: NICE requires EQ-5D with UK preference weights (3L or 5L)

### G6: Extrapolation beyond trial data without flagging
**Flag if**: model time horizon substantially exceeds trial follow-up and no sensitivity analysis on post-trial extrapolation assumptions
**Why it matters**: key area of scrutiny by ERGs/HTA committees
