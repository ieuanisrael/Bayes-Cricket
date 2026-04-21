# Example: set working directory to the Bayes Cricket project root, then:
#   Rscript R/run_example.R

log_msg <- function(...) {
  cat(paste0(..., collapse = ""), "\n", sep = "", file = stderr())
  flush(stderr())
}

proj <- Sys.getenv("BAYES_CRICKET_ROOT", unset = "")
if (!nzchar(proj)) {
  proj <- normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
  if (!file.exists(file.path(proj, "stan", "runs_nb_hierarchical.stan"))) {
    proj <- getwd()
  }
}
setwd(proj)

log_msg("Bayes Cricket: project root = ", normalizePath(proj, winslash = "/", mustWork = TRUE))

mc <- Sys.getenv("MC_CORES", unset = "")
if (nzchar(mc)) {
  v <- suppressWarnings(as.integer(mc))
  if (!is.na(v) && v > 0L) {
    options(mc.cores = v)
  }
}

makevars <- file.path(proj, ".R", "Makevars")
if (file.exists(makevars)) {
  Sys.setenv(R_MAKEVARS_USER = normalizePath(makevars))
}

source("R/fit_runs_model.R")

log_msg("Loading data/example_innings.csv ...")
innings <- read.csv("data/example_innings.csv", stringsAsFactors = FALSE)
log_msg(
  "Rows: ", nrow(innings),
  "; compiling/fitting Stan model (first run compiles C++; may take a few minutes) ..."
)

fitted <- fit_runs_ranking(
  innings,
  context_cols = c("is_home", "opp_strength", "format"),
  chains = 2L,
  iter = 800L,
  warmup = 400L,
  seed = 42L,
  stan_file = file.path(getwd(), "stan", "runs_nb_hierarchical.stan"),
  control = list(adapt_delta = 0.9),
  verbose = TRUE,
  refresh = 50L
)

log_msg("Done sampling. Posterior ranking (by theta_mean):")
print(summarize_player_ranking(fitted, prob = 0.9))
flush(stdout())
log_msg("P(best player):")
print(best_player_probabilities(fitted))
flush(stdout())
