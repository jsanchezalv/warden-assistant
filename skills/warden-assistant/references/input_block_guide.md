# `input_block()` Complete Guide

## What it does

`input_block()` is a high-level wrapper that generates the code needed to automatically select the correct parameter value depending on whether the model is running:
- **Deterministic** (base case only)
- **PSA** (draw from distributions)
- **DSA** (vary one parameter/group at a time through min/max)
- **Probabilistic DSA** (PSA + DSA simultaneously)
- **Scenario** (apply all alternative values from a named column at once)

It replaces the manual `pick_val_v()` + `create_indicators()` + `sens_iterator()` pattern with a single call that auto-detects `n_sensitivity` and handles offsets across multiple blocks.

---

## The data structure

`input_block()` expects parameters organized as a **named list of lists** (or equivalently, a list-of-columns structure like a data.frame):

```r
l_inputs <- list(
  parameter_name = list("util.pfs", "util.pps", "cost.drug", "HR"),
  base_value     = list(0.82, 0.55, 1500, 0.75),
  PSA_dist       = list("rbeta_mse", "rbeta_mse", "rgamma_mse", "rlnorm"),
  a              = list(0.82, 0.55, 1500, log(0.75)),
  b              = list(0.04, 0.03, 300, 0.15),
  n              = list(1, 1, 1, 1),
  DSA_min        = list(0.70, 0.40, 1000, 0.50),
  DSA_max        = list(0.90, 0.70, 2000, 0.95),
  scenario_opt   = list(0.85, 0.60, 1200, 0.65),
  scenario_pess  = list(0.75, 0.45, 2500, 0.90),
  psa_indicators = list(1, 1, 1, 1)
)
```

Column meanings:
- `parameter_name`: output variable names (what gets assigned in the model environment)
- `base_value`: deterministic base case values
- `PSA_dist`: distribution function name (e.g., `"rnorm"`, `"rbeta_mse"`, `"rgamma_mse"`, `"rlnorm"`, `"mvrnorm"`)
- `a`, `b`, `n`: arguments passed to the distribution function via `pick_psa(dist, n, a, b)`
- `DSA_min`, `DSA_max`: column names for DSA (one parameter varied at a time)
- `scenario_*`: column names for scenario analyses (all params take their alternative value simultaneously)
- `psa_indicators`: 0/1 per parameter — which ones draw from PSA (0 = always use base)

---

## Basic usage: PSA only (no DSA or scenarios)

When you only need deterministic + PSA, omit `sens` and `dsa_names`:

```r
common_all_inputs <- input_block(
  base           = l_inputs[["base_value"]],
  psa            = pick_psa(l_inputs[["PSA_dist"]], l_inputs[["n"]],
                            l_inputs[["a"]], l_inputs[["b"]]),
  names_out      = l_inputs[["parameter_name"]],
  psa_indicators = l_inputs[["psa_indicators"]]
) |> add_item(drc = 0.035, drq = 0.035)

# Deterministic:
run_sim(..., psa_bool = FALSE)
# PSA:
run_sim(..., psa_bool = TRUE, n_sim = 500)
```

---

## DSA: one-at-a-time sensitivity

For DSA, add `sens` (the full data object) and `dsa_names` (which columns are DSA directions):

```r
common_all_inputs <- input_block(
  base           = l_inputs[["base_value"]],
  psa            = pick_psa(l_inputs[["PSA_dist"]], l_inputs[["n"]],
                            l_inputs[["a"]], l_inputs[["b"]]),
  sens           = l_inputs,
  names_out      = l_inputs[["parameter_name"]],
  psa_indicators = l_inputs[["psa_indicators"]],
  dsa_names      = c("DSA_min", "DSA_max")
)

# DSA (deterministic, one param at a time):
results_dsa <- run_sim(
  ...,
  psa_bool         = FALSE,
  sensitivity_bool = TRUE
  # n_sensitivity and sensitivity_names auto-detected from input_block metadata
)

# Probabilistic DSA (PSA + one param at a time):
results_prob_dsa <- run_sim(
  ...,
  psa_bool         = TRUE,
  n_sim            = 100,
  sensitivity_bool = TRUE
)
```

**Key point**: `dsa_names` tells the engine which `sensitivity_names` iterate one-at-a-time through parameters. Everything else is treated as a scenario.

---

## Scenarios: all parameters change simultaneously

For scenario analyses, add columns to the data structure with scenario values, but do NOT include those column names in `dsa_names`:

```r
# l_inputs has columns: scenario_opt, scenario_pess (in addition to DSA_min, DSA_max)

common_all_inputs <- input_block(
  base           = l_inputs[["base_value"]],
  psa            = pick_psa(l_inputs[["PSA_dist"]], l_inputs[["n"]],
                            l_inputs[["a"]], l_inputs[["b"]]),
  sens           = l_inputs,
  names_out      = l_inputs[["parameter_name"]],
  psa_indicators = l_inputs[["psa_indicators"]],
  dsa_names      = c("DSA_min", "DSA_max")  # only these iterate one-at-a-time
)

# Scenario analysis:
results_scen <- run_sim(
  ...,
  psa_bool          = FALSE,
  sensitivity_bool  = TRUE,
  sensitivity_names = c("scenario_opt", "scenario_pess")
  # These are NOT in dsa_names, so ALL parameters take their scenario value at once
)
```

**The critical distinction**:
- Column names listed in `dsa_names` → one-parameter-at-a-time iteration (DSA)
- Column names NOT in `dsa_names` but passed via `sensitivity_names` in `run_sim()` → all parameters use their alternative value simultaneously (scenario)

---

## Covarying parameters in DSA (grouped mode)

By default (`indicator_sens_binary = TRUE` or omitted), each parameter is varied independently in DSA. To vary multiple parameters together (e.g., utilities together, costs together), use **grouped mode**:

```r
l_inputs <- list(
  parameter_name = list("util.pfs", "util.pps", "cost.drug", "cost.admin", "HR"),
  base_value     = list(0.82, 0.55, 1500, 200, 0.75),
  # ... PSA_dist, a, b, n, DSA_min, DSA_max as before ...
  psa_indicators = list(1, 1, 1, 1, 1),
  dsa_indicators = list(1, 1, 2, 2, 3)
  #                      ^--^     ^--^   ^
  #                      group1   group2  group3
  # group 1: both utilities varied together
  # group 2: both costs varied together
  # group 3: HR varied alone
)

common_all_inputs <- input_block(
  base                  = l_inputs[["base_value"]],
  psa                   = pick_psa(l_inputs[["PSA_dist"]], l_inputs[["n"]],
                                   l_inputs[["a"]], l_inputs[["b"]]),
  sens                  = l_inputs,
  names_out             = l_inputs[["parameter_name"]],
  psa_indicators        = l_inputs[["psa_indicators"]],
  dsa_names             = c("DSA_min", "DSA_max"),
  sens_indicators       = l_inputs[["dsa_indicators"]],
  indicator_sens_binary = FALSE,
  distributions         = l_inputs[["PSA_dist"]],
  covariances           = l_inputs[["b"]]
)
```

**How `sens_indicators` works in grouped mode:**
- `0` = never varied (always uses base/PSA value)
- Same integer = varied together in the same DSA iteration
- The total number of DSA iterations = number of unique non-zero groups (here: 3)

**Required extra arguments in grouped mode:**
- `distributions`: list of distribution names (needed to know which params are vectors)
- `covariances`: list of variance/covariance values (scalars for univariate, matrices for `mvrnorm`)

---

## Vector-valued parameters (e.g., `mvrnorm`)

Some parameters are vectors (e.g., transition probabilities drawn from a multivariate normal or Dirichlet). Handle these by making the relevant entries vectors:

```r
l_inputs_pat <- list(
  parameter_name = list("age", "sex", "v_trans_probs"),
  base_value     = list(60, 1, c(0.3, 0.5)),          # v_trans_probs is length 2
  PSA_dist       = list("rnorm", "rbinom", "mvrnorm"),
  a              = list(60, 1, c(0.3, 0.5)),           # mean vector for mvrnorm
  b              = list(10, 0.5, matrix(c(0.01, 0.005, 0.005, 0.02), 2, 2)),  # covariance matrix
  n              = list(1, 1, 1),
  DSA_min        = list(30, 0, c(0.2, 0.3)),
  DSA_max        = list(80, 1, c(0.4, 0.7)),
  psa_indicators = list(1, 1, c(1, 0)),                # vary first element in PSA, not second
  dsa_indicators = list(1, 1, c(2, 2))                 # both elements of v_trans_probs covaried
)

i_pat <- input_block(
  base                  = l_inputs_pat[["base_value"]],
  psa                   = pick_psa(l_inputs_pat[["PSA_dist"]], l_inputs_pat[["n"]],
                                   l_inputs_pat[["a"]], l_inputs_pat[["b"]]),
  sens                  = l_inputs_pat,
  names_out             = l_inputs_pat[["parameter_name"]],
  psa_indicators        = l_inputs_pat[["psa_indicators"]],
  dsa_names             = c("DSA_min", "DSA_max"),
  sens_indicators       = l_inputs_pat[["dsa_indicators"]],
  indicator_sens_binary = FALSE,
  distributions         = l_inputs_pat[["PSA_dist"]],
  covariances           = l_inputs_pat[["b"]]
)
```

**Key rules for vector parameters:**
- `psa_indicators` can be a vector per element (e.g., `c(1, 0)` means first component varies in PSA, second doesn't)
- `dsa_indicators` uses same-integer grouping for sub-elements
- The output variable (`v_trans_probs`) will be a vector in the model environment

---

## Multiple blocks at different levels

Parameters often live at different levels:
- **Simulation-level** (`common_all_inputs`): fixed per simulation (costs, utilities, HRs)
- **Patient-level** (`common_pt_inputs`): drawn per patient (age, sex, comorbidities)
- **Arm-level** (`unique_pt_inputs`): different per arm (starting treatment, dose)

Use `input_block()` at each level. The engine automatically computes DSA offsets across all blocks:

```r
# Simulation-level (7 params):
i_sim <- input_block(
  base      = l_inputs[["base_value"]],
  psa       = pick_psa(l_inputs[["PSA_dist"]], l_inputs[["n"]],
                       l_inputs[["a"]], l_inputs[["b"]]),
  sens      = l_inputs,
  names_out = l_inputs[["parameter_name"]],
  psa_indicators = l_inputs[["psa_indicators"]],
  dsa_names = c("DSA_min", "DSA_max")
)

# Patient-level (2 params):
i_pat <- input_block(
  base      = l_inputs_pat[["base_value"]],
  psa       = pick_psa(l_inputs_pat[["PSA_dist"]], l_inputs_pat[["n"]],
                       l_inputs_pat[["a"]], l_inputs_pat[["b"]]),
  sens      = l_inputs_pat,
  names_out = l_inputs_pat[["parameter_name"]],
  psa_indicators = l_inputs_pat[["psa_indicators"]],
  dsa_names = c("DSA_min", "DSA_max")
)

results <- run_sim(
  ...,
  common_all_inputs = i_sim,
  common_pt_inputs  = i_pat,
  sensitivity_bool  = TRUE
  # n_sensitivity auto-detected: 7 + 2 = 9 total DSA iterations
  # Offsets computed automatically: patient-level block starts at index 8
)
```

**No manual `n_sensitivity` or `n_sens_before` needed** — the engine accumulates across all `input_block()` calls in the model.

---

## Mixing `input_block()` with `add_item()`

`input_block()` returns a `{}` expression, same as `add_item()`. They compose with `|>`:

```r
common_all_inputs <- input_block(
  base      = l_inputs[["base_value"]],
  psa       = pick_psa(...),
  sens      = l_inputs,
  names_out = l_inputs[["parameter_name"]],
  dsa_names = c("DSA_min", "DSA_max")
) |> add_item(
  # Additional derived params not in the input table:
  drc = 0.035,
  drq = 0.035,
  mort_vec = setNames(life_table$rate, life_table$age)
)
```

Order matters: `input_block()` output variables are available in subsequent `add_item()` blocks.

---

## Common mistakes

### 1. Forgetting `dsa_names` — everything becomes a scenario
```r
# WRONG: without dsa_names, "DSA_min" and "DSA_max" are treated as scenarios
input_block(base=..., psa=..., sens=l_inputs, names_out=...)

# RIGHT: explicitly declare which columns are DSA (one-at-a-time)
input_block(base=..., psa=..., sens=l_inputs, names_out=...,
            dsa_names = c("DSA_min", "DSA_max"))
```

### 2. DSA column values that violate parameter constraints
```r
# WRONG: utility DSA_max > 1 (impossible for a utility)
DSA_max = list(1.1, ...)
# WRONG: HR DSA_min is negative
DSA_min = list(-0.2, ...)
```

### 3. Using grouped mode without `distributions` and `covariances`
```r
# WRONG: grouped mode requires distributions and covariances
input_block(..., sens_indicators = list(1,1,2,2), indicator_sens_binary = FALSE)
# Error: distributions cannot be NULL in grouped mode

# RIGHT:
input_block(..., sens_indicators = list(1,1,2,2), indicator_sens_binary = FALSE,
            distributions = list("rnorm","rnorm","rgamma_mse","rgamma_mse"),
            covariances = list(0.04, 0.03, 300, 200))
```

### 4. Passing `sensitivity_names` for DSA when using `input_block()`
```r
# UNNECESSARY: input_block auto-detects sensitivity_names for DSA
run_sim(..., sensitivity_bool = TRUE,
        sensitivity_names = c("DSA_min", "DSA_max"))  # redundant but harmless

# REQUIRED: for scenarios, you MUST pass sensitivity_names in run_sim()
run_sim(..., sensitivity_bool = TRUE,
        sensitivity_names = c("scenario_opt", "scenario_pess"))
```

### 5. Mismatched `psa_indicators` with `sens_indicators` in grouped mode
```r
# The psa_indicators control PSA draws, sens_indicators control DSA grouping
# They are independent — a parameter can be excluded from PSA but included in DSA
psa_indicators = list(0, 0, 1, 1)     # first two excluded from PSA
dsa_indicators = list(1, 1, 2, 2)     # but still varied in DSA (groups 1 and 2)
```

### 6. Setting `sens_indicators` to 0 expecting the parameter to use base in scenarios
```r
# 0 in sens_indicators means the parameter is NEVER varied in DSA or scenarios
# If you want it varied in scenarios but not DSA, use indicator_sens_binary = TRUE
# and don't include it in dsa_names columns — or leave it as 1 in sens_indicators
```

---

## Decision tree: which mode to use

```
Do you have sensitivity analysis at all?
├─ No → omit sens, dsa_names, sens_indicators
│       (just base + PSA)
│
└─ Yes
   ├─ Only scenarios (all params change together)?
   │  → set sens = l_inputs, omit dsa_names
   │    In run_sim(): sensitivity_names = c("scenario_1", "scenario_2")
   │
   ├─ Only DSA (one param at a time)?
   │  → set sens = l_inputs, dsa_names = c("DSA_min", "DSA_max")
   │    In run_sim(): sensitivity_bool = TRUE (auto-detected)
   │
   ├─ Both DSA and scenarios?
   │  → set sens = l_inputs, dsa_names = c("DSA_min", "DSA_max")
   │    For DSA: run_sim(sensitivity_bool = TRUE) — auto-detected
   │    For scenarios: run_sim(sensitivity_bool = TRUE,
   │                           sensitivity_names = c("scenario_opt", ...))
   │
   └─ DSA with covaried parameters?
      → add: sens_indicators = list(1,1,2,2,3,...),
             indicator_sens_binary = FALSE,
             distributions = list(...),
             covariances = list(...)
```

---

## How `run_sim()` interacts with `input_block()`

| `run_sim()` argument | With `input_block()` | Without (manual `pick_val_v`) |
|---------------------|---------------------|-------------------------------|
| `n_sensitivity` | Auto-detected from block metadata | Must set manually |
| `sensitivity_names` (DSA) | Auto-detected from `dsa_names` | Must set manually |
| `sensitivity_names` (scenarios) | Must set manually | Must set manually |
| `sensitivity_inputs` | Not needed | Must build with `sens_iterator()` + `create_indicators()` |
| `sensitivity_bool` | Set to `TRUE` to activate | Set to `TRUE` to activate |
| `psa_bool` | Set to `TRUE` for PSA | Set to `TRUE` for PSA |

---

## Validation checklist

When reviewing a model's `input_block()` setup:

1. **Column consistency**: all list entries in `l_inputs` must have the same length (number of parameters)
2. **PSA distribution args**: `a` and `b` must match what the distribution expects:
   - `rnorm(n, mean, sd)` → a=mean, b=sd
   - `rbeta_mse(n, mean, se)` → a=mean, b=se
   - `rgamma_mse(n, mean, se)` → a=mean, b=se
   - `rlnorm(n, meanlog, sdlog)` → a=meanlog, b=sdlog (LOG scale)
   - `mvrnorm(n, mu, Sigma)` → a=mu (vector), b=Sigma (matrix)
3. **`dsa_names` columns exist** in `l_inputs` (column names must match exactly)
4. **Scenario columns exist** in `l_inputs` (and match what's passed to `sensitivity_names`)
5. **Parameter constraints respected** in DSA/scenario columns (utilities ≤ 1, costs ≥ 0, HR > 0)
6. **`sens_indicators` groups** are sequential integers starting from 1 when covarying
7. **Vector parameters** have matching vector lengths in base, DSA, and scenario columns
