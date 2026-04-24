#!/usr/bin/env Rscript
# Build actionable plots/tables from a saved fit and the same innings CSV used for fitting.
#
#   Rscript R/run_actionable_views.R \
#     --fitted outputs/run_with_fpca_age/fitted.rds \
#     --data data/example_innings.csv \
#     --out outputs/actionable_run1
#
# If the fit used --fpca-age-k, pass the same CSV after recomputing FPC columns, or use the
# CSV that was saved alongside the run (you can add a save step). For a one-shot pipeline,
# run run_reproducible_model.R then this script with the same --data file and re-apply FPCA
# flags below to match the fit.

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
        opts[[kv[[1L]]]] <- paste(kv[-1L], collapse = "=")
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
    if (grepl("(^|/)R$", p)) {
      normalizePath(file.path(p, ".."), winslash = "/")
    } else {
      p
    }
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

opt <- parse_cli(commandArgs(trailingOnly = TRUE))
if (any(commandArgs(trailingOnly = TRUE) %in% c("-h", "--help")) || !is.null(opt$help)) {
  cat("
Usage:
  Rscript R/run_actionable_views.R --fitted PATH --data PATH --out DIR

Options:
  --fitted PATH   fitted.rds from fit_runs_ranking()
  --data PATH     Same innings CSV (and row order) as used for the fit
  --out DIR       Output directory for plots/ and CSVs
  --fpca-age-k K  If >0, append historical FPC columns before views (must match fit)
  --threshold N   For next-season table (default 350)
  --n-future N    Innings to sum in simulation (default 20)
  --cluster-k K   k-means clusters (default 4)
")
  quit(status = 0L)
}

fitted_path <- opt[["fitted"]]
data_path <- opt[["data"]]
if (is.null(fitted_path) || !nzchar(as.character(fitted_path))) {
  stop("--fitted is required", call. = FALSE)
}
if (is.null(data_path) || !nzchar(as.character(data_path))) {
  stop("--data is required", call. = FALSE)
}
out_dir <- opt[["out"]] %||% file.path("outputs", "actionable")
if (!grepl("^(/|~)", out_dir)) {
  out_dir <- file.path(proj, out_dir)
}
if (!grepl("^(/|~)", data_path)) {
  data_path <- file.path(proj, data_path)
}

fpca_k <- suppressWarnings(as.integer(opt[["fpca-age-k"]] %||% "0"))
if (is.na(fpca_k)) {
  fpca_k <- 0L
}
threshold <- as.numeric(opt$threshold %||% "350")
n_future <- as.integer(opt[["n-future"]] %||% "20")
cluster_k <- as.integer(opt[["cluster-k"]] %||% "4")

source(file.path(proj, "R", "fit_runs_model.R"))
source(file.path(proj, "R", "actionable_views.R"))

fitted <- readRDS(fitted_path)
innings <- read.csv(data_path, stringsAsFactors = FALSE)

if (fpca_k > 0L) {
  source(file.path(proj, "R", "fpca_age_predictors.R"))
  tcol <- opt[["fpca-time-col"]]
  if (is.null(tcol) || !nzchar(as.character(tcol)) || tcol == "TRUE") {
    tcol <- NULL
  }
  innings <- append_historical_age_curve_pc_scores(innings, time_order_col = tcol, K = fpca_k)
}

run_all_actionable_views(
  fitted,
  innings,
  out_dir = out_dir,
  cluster_k = cluster_k,
  n_future = n_future,
  threshold = threshold
)

cat("Wrote actionable views to:", normalizePath(out_dir, winslash = "/"), "\n")
quit(status = 0L)
