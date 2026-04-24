proj <- Sys.getenv("BAYES_CRICKET_ROOT", unset = "")
if (!nzchar(proj)) {
  proj <- normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
  if (!file.exists(file.path(proj, "stan", "runs_nb_hierarchical.stan"))) {
    proj <- getwd()
  }
}
setwd(proj)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Please install ggplot2.", call. = FALSE)
}

mc <- Sys.getenv("MC_CORES", unset = "")
if (nzchar(mc)) {
  v <- suppressWarnings(as.integer(mc))
  if (!is.na(v) && v > 0L) options(mc.cores = v)
}

source("R/fit_runs_model.R")
source("R/fpca_age_predictors.R")

# Innings-level `age` / `age_sq` in the Stan model below are separate from the
# **season-aggregated** age–runs summaries used for exploratory curves / FPCA.

innings <- read.csv("data/example_innings.csv", stringsAsFactors = FALSE)
if (!("age" %in% names(innings))) {
  stop("Column `age` not found in data/example_innings.csv", call. = FALSE)
}
innings$age <- as.numeric(innings$age)
if (anyNA(innings$age)) {
  stop("Column `age` must be numeric without NA.", call. = FALSE)
}
innings$age_sq <- innings$age^2

out_dir <- file.path("outputs", "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chains <- as.integer(Sys.getenv("CHAINS", unset = "2"))
iter <- as.integer(Sys.getenv("ITER", unset = "700"))
warmup <- as.integer(Sys.getenv("WARMUP", unset = "350"))

fit <- fit_runs_ranking(
  innings = innings,
  context_cols = c("is_home", "opp_strength", "format", "age", "age_sq"),
  chains = chains,
  iter = iter,
  warmup = warmup,
  seed = 42,
  control = list(adapt_delta = 0.9),
  refresh = 50L
)

beta_post <- rstan::extract(fit$fit, pars = "beta", permuted = TRUE)$beta
beta_names <- fit$context_info$beta_names
idx_age <- which(beta_names == "age")
idx_age_sq <- which(beta_names == "age_sq")

if (length(idx_age) != 1L || length(idx_age_sq) != 1L) {
  stop("Could not find both age and age_sq coefficients in posterior beta.", call. = FALSE)
}

# Undo standardization done in fit_runs_model.R
age_center <- mean(innings$age)
age_scale <- stats::sd(innings$age)
age_sq_center <- mean(innings$age_sq)
age_sq_scale <- stats::sd(innings$age_sq)

age_grid <- seq(min(innings$age), max(innings$age), length.out = 120)

age_effect_draws <- sapply(age_grid, function(a) {
  z_age <- (a - age_center) / age_scale
  z_age_sq <- (a^2 - age_sq_center) / age_sq_scale
  beta_post[, idx_age] * z_age + beta_post[, idx_age_sq] * z_age_sq
})

age_mu <- apply(age_effect_draws, 2, mean)
age_lo <- apply(age_effect_draws, 2, stats::quantile, probs = 0.05)
age_hi <- apply(age_effect_draws, 2, stats::quantile, probs = 0.95)

age_df <- data.frame(
  age = age_grid,
  delta_log_runs_mean = age_mu,
  delta_log_runs_lo = age_lo,
  delta_log_runs_hi = age_hi,
  runs_multiplier_mean = exp(age_mu),
  runs_multiplier_lo = exp(age_lo),
  runs_multiplier_hi = exp(age_hi)
)

p_age_multiplier <- ggplot2::ggplot(age_df, ggplot2::aes(x = age, y = runs_multiplier_mean)) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = runs_multiplier_lo, ymax = runs_multiplier_hi),
    fill = "#93c5fd",
    alpha = 0.35
  ) +
  ggplot2::geom_line(color = "#1d4ed8", linewidth = 1.1) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#111827") +
  ggplot2::labs(
    title = "Age effect on expected runs (posterior)",
    x = "Age",
    y = "Expected runs multiplier\n(relative to baseline on covariate scale)"
  ) +
  ggplot2::theme_minimal()

ggplot2::ggsave(
  file.path(out_dir, "posterior_age_runs_multiplier.png"),
  p_age_multiplier,
  width = 8,
  height = 5,
  dpi = 120
)

p_data_age <- ggplot2::ggplot(innings, ggplot2::aes(x = age, y = runs)) +
  ggplot2::geom_point(alpha = 0.08, color = "#374151") +
  ggplot2::geom_smooth(method = "loess", se = TRUE, color = "#dc2626") +
  ggplot2::labs(
    title = "Observed runs vs age",
    x = "Age",
    y = "Runs"
  ) +
  ggplot2::theme_minimal()

ggplot2::ggsave(
  file.path(out_dir, "data_runs_vs_age.png"),
  p_data_age,
  width = 8,
  height = 5,
  dpi = 120
)

# Posterior distribution of age coefficients
coef_df <- data.frame(
  value = c(beta_post[, idx_age], beta_post[, idx_age_sq]),
  term = c(
    rep("beta_age (standardized age)", nrow(beta_post)),
    rep("beta_age_sq (standardized age^2)", nrow(beta_post))
  )
)

p_coef <- ggplot2::ggplot(coef_df, ggplot2::aes(x = value, fill = term)) +
  ggplot2::geom_density(alpha = 0.35) +
  ggplot2::facet_wrap(~term, scales = "free", ncol = 1) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
  ggplot2::guides(fill = "none") +
  ggplot2::labs(
    title = "Posterior coefficients for age terms",
    x = "Coefficient value",
    y = "Density"
  ) +
  ggplot2::theme_minimal()

ggplot2::ggsave(
  file.path(out_dir, "posterior_age_coefficients.png"),
  p_coef,
  width = 8,
  height = 7,
  dpi = 120
)

peak_idx <- which.max(age_df$runs_multiplier_mean)
summary_txt <- c(
  paste0("Age range in data: ", round(min(innings$age), 2), " to ", round(max(innings$age), 2)),
  paste0("Posterior mean peak expected runs at age ~", round(age_df$age[peak_idx], 2)),
  paste0(
    "Runs multiplier at peak age (mean, 90% CI): ",
    round(age_df$runs_multiplier_mean[peak_idx], 3), " [",
    round(age_df$runs_multiplier_lo[peak_idx], 3), ", ",
    round(age_df$runs_multiplier_hi[peak_idx], 3), "]"
  )
)
writeLines(summary_txt, con = file.path(out_dir, "age_effect_summary.txt"))

if ("season" %in% names(innings)) {
  ps <- player_season_summary(
    innings,
    season_col = "season",
    age_col = "age",
    runs_col = "runs",
    player_col = "player_id",
    statistic = "mean"
  )
  p_season <- ggplot2::ggplot(ps, ggplot2::aes(x = .data$mean_age, y = .data$mean_runs)) +
    ggplot2::geom_point(alpha = 0.2, color = "#4338ca", size = 1) +
    ggplot2::geom_smooth(method = "loess", se = TRUE, color = "#b91c1c", linewidth = 0.9) +
    ggplot2::labs(
      title = "Player–season: mean age vs mean runs per innings",
      subtitle = "Each point is one player–season aggregate (not innings-level).",
      x = "Mean age in season",
      y = "Mean runs per innings in season"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(
    file.path(out_dir, "season_mean_age_vs_mean_runs.png"),
    p_season,
    width = 8,
    height = 5,
    dpi = 120
  )
}

cat("Saved age investigation plots to:", normalizePath(out_dir), "\n")
