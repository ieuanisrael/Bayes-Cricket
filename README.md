# Bayes Cricket

Hierarchical Bayesian model for **runs per innings** using **Stan** (via **rstan**): player-specific log-ability (`theta`), contextual predictors (`X β`), and a **negative binomial** likelihood. The goal is partially pooled player rankings and interpretable context effects.

---

## Requirements

| Item | Notes |
|------|--------|
| **R** | 4.x recommended; Docker image uses `rocker/r-ver:4.4.2`. |
| **rstan** | `install.packages("rstan", dependencies = TRUE)`. First run compiles the Stan model to C++ (one-time per model change). |
| **ggplot2** | For diagnostic and context-effect plots. |
| **C++ toolchain** | On **macOS**, Apple Command Line Tools are required. The repo includes **`.R/Makevars`** so the compiler sees the system SDK. `fit_runs_ranking()` sets `R_MAKEVARS_USER` automatically when that file is present. |

---

## Data format

Your CSV must include at least:

- **`player_id`** — any character/id per batter.
- **`runs`** — non-negative integer, runs in that innings.
- **`season`** (recommended) — any label that identifies the competition season for that row; required for **season-level** age curves / FPCA and for several actionable plots. The example file uses short codes such as `P001_01` (player + block index).

Optional **context** columns (passed via `--context` or `context_cols` in R):

- **Numeric** — centered and scaled inside the model.
- **Character / factor** — expanded to one-hot (full rank) **treatment** contrasts; the **first** level (alphabetical by default) is the **reference** level.

The example file is **`data/example_innings.csv`**.

---

## Quick start (local R)

From the **repository root** (so `stan/` and `R/` resolve correctly). Recommended: set `BAYES_CRICKET_ROOT` to the project path if you run scripts from elsewhere.

```bash
export BAYES_CRICKET_ROOT="$(pwd)"   # optional, from repo root
Rscript R/run_example.R
```

This fits the model with default context columns, then prints a player ranking and P(best) summaries.

### Reproducible CLI (recommended for experiments)

Fits the model, saves **`fitted.rds`**, a **`context_betas.csv`**, **`run_manifest.txt`**, and **plots** of each context coefficient:

```bash
Rscript R/run_reproducible_model.R \
  --data data/example_innings.csv \
  --context is_home,opp_strength,format \
  --out outputs/my_run
```

**Useful options**

| Option | Default | Description |
|--------|---------|-------------|
| `--data` | `data/example_innings.csv` | Path to CSV (absolute or relative to project). |
| `--context` | `is_home,opp_strength,format` | Comma-separated context column names. |
| `--out` | `outputs/repro` | Output directory. |
| `--chains` | `2` (or `CHAINS` env) | HMC chains. |
| `--iter` / `--warmup` / `--seed` | 800 / 400 / 42 | Sampling controls. |
| `--stan` | `stan/runs_nb_hierarchical.stan` | Alternative Stan file. |
| `--adapt-delta` | `0.9` | HMC; increase if you see divergences. |
| `--fpca-age-k` | `0` | If `K > 0`, aggregates to **player–season**, builds curves of prior-season `(mean_age, mean_runs)` on an age grid, PCA → `fpc_age_*`, merged to innings. Needs **`age`** + **`season`** (see `--fpca-season-col`). |
| `--fpca-season-col` | `season` | Season id for aggregation. |
| `--fpca-season-fallback` | off | If `TRUE`, use `floor(age)` as season when missing (**demo only**). |
| `--fpca-statistic` | `mean` | `mean` or `median` within-season runs for the curve. |

Environment variables **`DATA`**, **`CONTEXT`**, **`MC_CORES`**, **`ITER`**, etc., override defaults when a flag is omitted. Run **`Rscript R/run_reproducible_model.R --help`** for the full list.

**Season-level age curves:** Curves are built from **seasonal** summaries (mean/median runs per innings in the season vs mean age in the season), using only **previous seasons** for that player, then discretized **PCA**; scores are attached to every innings in that season. Stan still models **innings-level** outcomes.

**Outputs** (under `--out`):

- **`fitted.rds`** — R object from `fit_runs_ranking()` (reload with `readRDS()`).
- **`context_betas.csv`** — posterior summaries for each `β` component.
- **`run_manifest.txt`** — R version, data path, context list, seed, and sampling settings.
- **`plots/`** — forest plots (log and `exp(β)` scale) and per-coefficient density plots.

**Interpretation:** `η = log(μ) = θⱼ + Xβ` for the negative binomial mean. Coefficients `β` are on the **log-mean** scale; **`exp(β)`** is a **multiplier** on expected runs, holding the player and other columns fixed. Numeric predictors are per **one SD**; factors are **vs the reference** level.

### Actionable views (after a fit)

From a saved **`fitted.rds`** and the **same** innings CSV (and design: re-apply **`--fpca-age-k`** if the fit used it), run:

```bash
Rscript R/run_actionable_views.R \
  --fitted outputs/my_run/fitted.rds \
  --data data/example_innings.csv \
  --out outputs/actionable_1 \
  --fpca-age-k 0
```

This writes **`plots/`**: **season-level** mean age vs mean residual (`season_age_vs_residual.png`), shrinkage (two panels), season observed vs posterior-mean expected totals (uses **`season`**, else **`match_date`** year, else **`floor(age)`**), plus **`next_season_threshold.csv`**: Monte Carlo **`P(sum of n_future innings runs > threshold)`** per player. If **`fpc_age_*`** exist, **`trajectory_clusters.*`** runs k-means on **one row per player–season** (not averaged over players).

---

## Other R scripts

| Script | Purpose |
|--------|--------|
| **`R/fit_runs_model.R`** | Source this for **`fit_runs_ranking()`**, ranking summaries, **`compare_recent_vs_posterior_predictive()`** (last *k* innings vs posterior predictive). |
| **`R/context_effects.R`** | Source with **`fit_runs_model.R`** for **`summarize_context_betas()`** and **`plot_context_effects()`** (used by the reproducible runner). |
| **`R/fpca_age_predictors.R`** | **`append_historical_age_curve_pc_scores()`** — age–run functional inputs for the hierarchical model (see `--fpca-age-k`). |
| **`R/actionable_views.R`** | **`run_all_actionable_views()`** — residuals vs age, shrinkage plots, trajectory clustering (if `fpc_age_*` present), next-season exceedance table. |
| **`R/run_actionable_views.R`** | CLI wrapper: **`--fitted`**, **`--data`**, **`--out`**, optional **`--fpca-age-k`**. |
| **`R/plot_priors_data_posteriors.R`** | Priors, data plots, a sample fit, and top-player posteriors under defaults. |
| **`R/investigate_age_effect.R`** | Age / age² in context; writes plots to **`outputs/plots/`** (expects `age` in the CSV). |

Settings like **`MC_CORES`**, **`CHAINS`**, **`ITER`**, **`WARMUP`** are read where documented in each script.

---

## Docker

Build and run the containerized example (same as **`R/run_example.R`** on the R side):

```bash
docker build -t bayes-cricket .
docker run --rm -t bayes-cricket
# More CPU for parallel chains:
docker run --rm -e MC_CORES=4 -t bayes-cricket
```

With **docker compose**:

```bash
docker compose build
docker compose run --rm model
```

To mount your own data directory (see **`Dockerfile`**), point **`--data`** at the file inside the mount. The default image entrypoint runs **`R/run_example.R`** only. For the reproducible CLI, run a shell in the container and call **`Rscript R/run_reproducible_model.R ...`**, or extend **`docker/entrypoint.sh`**.

---

## Model (short)

- **Likelihood:** `Runs ~ NegBinomial₂_log(μ, φ)` with `log(μ) = θ[player] + X β`.
- **Hierarchy:** `θⱼ = μ_global + σ_player * zⱼ`, with priors in **`stan/runs_nb_hierarchical.stan`**.
- **Posterior predictive replicates** `y_rep` are generated in Stan for model checking.

For serious inference, check **R̂**, **ESS**, and **divergences** in the Stan output; increase **`iter`**, **`warmup`**, or **`adapt_delta`** if needed.

---

## Troubleshooting

- **Stan C++ build fails on macOS** (e.g. missing `<new>`): install/fix **Xcode Command Line Tools** (`xcode-select --install`). Ensure the repo’s **`.R/Makevars`** is used (`R_MAKEVARS_USER` is set by **`fit_runs_ranking()`** when possible).
- **Divergences or high R̂:** raise **`--adapt-delta`**, more iterations, or simplify predictors.
- **“Recent innings”** logic assumes row or time order; pass a time column in **`compare_recent_vs_posterior_predictive()`** when you have dates.
