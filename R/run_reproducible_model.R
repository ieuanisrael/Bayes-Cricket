#!/usr/bin/env Rscript
# Reproducible CLI: fit the runs model and write outputs + context-effect plots.
#
# Examples (from repository root):
#   Rscript R/run_reproducible_model.R \
#     --data data/example_innings.csv \
#     --context is_home,opp_strength,format \
#     --out outputs/run1
#
#   Rscript R/run_reproducible_model.R  # uses defaults below
#
# Environment (optional):
#   BAYES_CRICKET_ROOT, MC_CORES, CHAINS, ITER, WARMUP, SEED, ADAPT_DELTA

parse_cli <- function(args = NULL) {
  if (is.null(args)) {
    args <- commandArgs(trailingOnly = TRUE)
  }
  opts <- list()
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (grepl("^--", a)) {
      a <- sub("^--", "", a)
      if (grepl("=", a, fixed = TRUE)) {
        kv <- strsplit(a, "=", fixed = TRUE)[[1L]]
        key <- kv[[1L]]
        val <- paste(kv[-1L], collapse = "=")
        opts[[key]] <- val
        i <- i + 1L
      } else {
        if (i < length(args) && !grepl("^--", args[[i + 1L]])) {
          opts[[a]] <- args[[i + 1L]]
          i <- i + 2L
        } else {
          opts[[a]] <- "TRUE"
          i <- i + 1L
        }
      }
    } else {
      i <- i + 1L
    }
  }
  opts
}

log_msg <- function(...) {
  cat(paste0(..., collapse = ""), "\n", sep = "", file = stderr())
  flush(stderr())
}

`%||%` <- function(x, y) {
  if (!is.null(x) && nzchar(as.character(x))) {
    return(x)
  }
  y
}

find_project_root <- function() {
  env <- Sys.getenv("BAYES_CRICKET_ROOT", unset = "")
  if (nzchar(env) && file.exists(file.path(env, "stan", "runs_nb_hierarchical.stan"))) {
    return(normalizePath(env, winslash = "/"))
  }
  args <- commandArgs()
  f <- sub("^--file=", "", args[grepl("^--file=", args)][1L])
  if (is.na(f) || !nzchar(f)) {
    f <- NA_character_
  } else {
    f <- normalizePath(f, winslash = "/", mustWork = TRUE)
  }
  from_script <- if (!is.na(f)) {
    p <- dirname(f)
    p <- if (grepl("(^|/)R$", p)) {
      normalizePath(file.path(p, ".."), winslash = "/")
    } else {
      p
    }
    p
  } else {
    getwd()
  }
  for (cand in c(from_script, getwd(), normalizePath(file.path(getwd(), ".."), mustWork = FALSE))) {
    if (file.exists(file.path(cand, "stan", "runs_nb_hierarchical.stan"))) {
      return(normalizePath(cand, winslash = "/"))
    }
  }
  from_script
}

proj <- find_project_root()
setwd(proj)

raw_args <- commandArgs(trailingOnly = TRUE)
opt <- parse_cli(raw_args)
if (any(raw_args %in% c("-h", "--help"))) {
  opt$help <- "TRUE"
}

# Defaults: env overrides, then static defaults
data_path <- opt$data %||% Sys.getenv("DATA", unset = "data/example_innings.csv")
# relative paths resolve from project root
if (!grepl("^(/|~)", data_path)) {
  data_path <- file.path(proj, data_path)
}

context_str <- opt$context %||% Sys.getenv("CONTEXT", unset = "is_home,opp_strength,format")
context_str <- gsub(" ", "", context_str, fixed = FALSE)
if (!nzchar(context_str)) {
  context_cols <- character()
} else {
  context_cols <- strsplit(context_str, ",", fixed = TRUE)[[1L]]
  context_cols <- context_cols[nzchar(context_cols)]
}

out_arg <- opt$out %||% file.path("outputs", "repro")
if (!grepl("^(/|~)", out_arg)) {
  out_dir <- file.path(proj, out_arg)
} else {
  out_dir <- path.expand(out_arg)
}
plot_subdir <- file.path(out_dir, "plots")
dir.create(plot_subdir, recursive = TRUE, showWarnings = FALSE)

makevars <- file.path(proj, ".R", "Makevars")
if (file.exists(makevars) && !nzchar(Sys.getenv("R_MAKEVARS_USER", ""))) {
  Sys.setenv(R_MAKEVARS_USER = normalizePath(makevars, winslash = "/"))
}

mc <- Sys.getenv("MC_CORES", unset = "")
if (nzchar(mc)) {
  v <- suppressWarnings(as.integer(mc))
  if (!is.na(v) && v > 0L) {
    options(mc.cores = v)
  }
}

source(file.path(proj, "R", "fit_runs_model.R"))
source(file.path(proj, "R", "context_effects.R"))

stan_file <- opt$stan %||% file.path(proj, "stan", "runs_nb_hierarchical.stan")
if (!is.null(opt$help) || !is.null(opt$h)) {
  cat("
Bayes-Cricket: fit hierarchical runs model (reproducible)

  --data PATH     CSV with player_id, runs, plus context columns (default: data/example_innings.csv)
  --context A,B   Comma-separated context column names (default: is_home,opp_strength,format)
                  Use empty or omit for no context (not recommended: weak beta meaning).
  --out DIR       Output root under project, or absolute path (default: outputs/repro)
  --chains N      (default: env CHAINS or 2)
  --iter N
  --warmup N
  --seed N
  --stan PATH     Stan file (default: stan/runs_nb_hierarchical.stan)
  --adapt-delta   Passed to rstan (default: env ADAPT_DELTA or 0.9)
  --fpca-age-k K  If K>0, append K historical age-curve PC scores (see R/fpca_age_predictors.R)
                  and add them to --context automatically. Requires numeric column `age`.
  --fpca-time-col NAME  Optional column for time ordering within player (else row order).
  -h, --help

Writes: <out>/fitted.rds, <out>/context_betas.csv, <out>/run_manifest.txt, <out>/plots/*.png
")
  quit(status = 0L)
}

chains <- as.integer(opt$chains %||% Sys.getenv("CHAINS", unset = "2"))
iter <- as.integer(opt$iter %||% Sys.getenv("ITER", unset = "800"))
warmup <- as.integer(opt$warmup %||% Sys.getenv("WARMUP", unset = "400"))
seed <- as.integer(opt$seed %||% Sys.getenv("SEED", unset = "42"))
ad_delta <- as.numeric(opt$`adapt-delta` %||% Sys.getenv("ADAPT_DELTA", unset = "0.9"))
if (anyNA(chains) || anyNA(iter) || anyNA(warmup) || anyNA(seed)) {
  stop("Invalid --chains, --iter, --warmup, or --seed", call. = FALSE)
}
if (is.na(ad_delta)) {
  ad_delta <- 0.9
}

fpca_k <- suppressWarnings(as.integer(opt[["fpca-age-k"]] %||% Sys.getenv("FPCA_AGE_K", unset = "0")))
if (is.na(fpca_k)) {
  fpca_k <- 0L
}

log_msg("Project root: ", normalizePath(proj, winslash = "/"))
log_msg("Data: ", data_path)
log_msg("Context columns: ", if (length(context_cols)) {
  paste(context_cols, collapse = ", ")
} else {
  "(none — Stan uses zero column; prior-only beta if J=1)"
})
log_msg("Output: ", out_dir)

if (!file.exists(data_path)) {
  stop("Data file not found: ", data_path, call. = FALSE)
}

innings <- read.csv(data_path, stringsAsFactors = FALSE)

if (fpca_k > 0L) {
  source(file.path(proj, "R", "fpca_age_predictors.R"))
  tcol <- opt[["fpca-time-col"]]
  if (is.null(tcol) || !nzchar(as.character(tcol)) || tcol == "TRUE") {
    tcol <- NULL
  }
  innings <- append_historical_age_curve_pc_scores(
    innings,
    time_order_col = tcol,
    K = fpca_k
  )
  fc_names <- grep("^fpc_age_[0-9]+$", names(innings), value = TRUE)
  context_cols <- unique(c(context_cols, fc_names))
  log_msg("Historical age-curve PC columns: ", paste(fc_names, collapse = ", "))
}

log_msg("Rows: ", nrow(innings), "; fitting ...")

t0 <- Sys.time()
set.seed(seed)
fitted <- fit_runs_ranking(
  innings,
  context_cols = context_cols,
  chains = chains,
  iter = iter,
  warmup = warmup,
  seed = seed,
  stan_file = stan_file,
  control = list(adapt_delta = ad_delta),
  refresh = 50L
)
log_msg("Sampling finished in ", round(difftime(Sys.time(), t0, units = "mins"), 2), " min")

if (is.null(fitted$fit)) {
  stop("fit_runs_ranking returned no Stan fit.", call. = FALSE)
}

saveRDS(fitted, file = file.path(out_dir, "fitted.rds"))
log_msg("Saved: ", file.path(out_dir, "fitted.rds"))

beta_sum <- summarize_context_betas(fitted, prob = 0.9)
write.csv(
  beta_sum,
  file = file.path(out_dir, "context_betas.csv"),
  row.names = FALSE
)
log_msg("Saved: ", file.path(out_dir, "context_betas.csv"))

plot_context_effects(fitted, plot_subdir, prob = 0.9)
if (nrow(beta_sum)) {
  log_msg("Plots in: ", normalizePath(plot_subdir, winslash = "/"))
}

manifest <- c(
  paste0("r_version: ", R.version.string),
  paste0("working_dir: ", normalizePath(proj, winslash = "/")),
  paste0("date_utc: ", format(as.POSIXct(Sys.time(), tz = "UTC"), usetz = TRUE)),
  paste0("data_file: ", data_path),
  paste0("context_columns: ", paste(context_cols, collapse = ", ")),
  paste0(
    "fpca_age_k: ",
    if (!is.na(fpca_k) && fpca_k > 0L) {
      as.character(fpca_k)
    } else {
      "0"
    }
  ),
  paste0("stan_file: ", normalizePath(stan_file, winslash = "/")),
  paste0("chains: ", chains, "  iter: ", iter, "  warmup: ", warmup, "  seed: ", seed, "  adapt_delta: ", ad_delta)
)
writeLines(manifest, con = file.path(out_dir, "run_manifest.txt"))
log_msg("Saved: ", file.path(out_dir, "run_manifest.txt"))

log_msg("Done.")
quit(status = 0L)
