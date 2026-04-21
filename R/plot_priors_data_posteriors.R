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
  if (!is.na(v) && v > 0L) {
    options(mc.cores = v)
  }
}

source("R/fit_runs_model.R")

out_dir <- file.path("outputs", "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

innings <- read.csv("data/example_innings.csv", stringsAsFactors = FALSE)

chains <- as.integer(Sys.getenv("CHAINS", unset = "2"))
iter <- as.integer(Sys.getenv("ITER", unset = "600"))
warmup <- as.integer(Sys.getenv("WARMUP", unset = "300"))

set.seed(42)

# Priors from the Stan model specification:
# mu_global ~ normal(log(25), 1.5)
# sigma_player ~ exponential(1)
# phi_inv ~ exponential(0.5), phi = 1/phi_inv
n_prior <- 5000
prior_df <- data.frame(
  mu_global = rnorm(n_prior, mean = log(25), sd = 1.5),
  sigma_player = rexp(n_prior, rate = 1),
  phi = 1 / rexp(n_prior, rate = 0.5)
)

p_mu <- ggplot2::ggplot(prior_df, ggplot2::aes(x = mu_global)) +
  ggplot2::geom_density(fill = "#3b82f6", alpha = 0.35) +
  ggplot2::labs(
    title = "Prior: mu_global",
    x = "mu_global (log expected runs)",
    y = "Density"
  ) +
  ggplot2::theme_minimal()

p_sigma <- ggplot2::ggplot(prior_df, ggplot2::aes(x = sigma_player)) +
  ggplot2::geom_density(fill = "#10b981", alpha = 0.35) +
  ggplot2::labs(
    title = "Prior: sigma_player",
    x = "sigma_player",
    y = "Density"
  ) +
  ggplot2::theme_minimal()

p_phi <- ggplot2::ggplot(prior_df, ggplot2::aes(x = phi)) +
  ggplot2::geom_density(fill = "#f59e0b", alpha = 0.35) +
  ggplot2::coord_cartesian(xlim = c(0, stats::quantile(prior_df$phi, 0.99))) +
  ggplot2::labs(
    title = "Prior: phi (dispersion)",
    x = "phi",
    y = "Density"
  ) +
  ggplot2::theme_minimal()

ggplot2::ggsave(file.path(out_dir, "prior_mu_global.png"), p_mu, width = 8, height = 5, dpi = 120)
ggplot2::ggsave(file.path(out_dir, "prior_sigma_player.png"), p_sigma, width = 8, height = 5, dpi = 120)
ggplot2::ggsave(file.path(out_dir, "prior_phi.png"), p_phi, width = 8, height = 5, dpi = 120)

# Data plots
p_runs <- ggplot2::ggplot(innings, ggplot2::aes(x = runs)) +
  ggplot2::geom_histogram(bins = 60, fill = "#6366f1", alpha = 0.8) +
  ggplot2::labs(title = "Observed runs distribution", x = "Runs", y = "Count") +
  ggplot2::theme_minimal()

p_format <- ggplot2::ggplot(innings, ggplot2::aes(x = format, y = runs, fill = format)) +
  ggplot2::geom_boxplot(alpha = 0.7, outlier.alpha = 0.15) +
  ggplot2::guides(fill = "none") +
  ggplot2::labs(title = "Runs by format", x = "Format", y = "Runs") +
  ggplot2::theme_minimal()

p_home <- ggplot2::ggplot(innings, ggplot2::aes(x = factor(is_home), y = runs, fill = factor(is_home))) +
  ggplot2::geom_boxplot(alpha = 0.7, outlier.alpha = 0.15) +
  ggplot2::guides(fill = "none") +
  ggplot2::labs(title = "Runs by home/away", x = "is_home", y = "Runs") +
  ggplot2::theme_minimal()

ggplot2::ggsave(file.path(out_dir, "data_runs_hist.png"), p_runs, width = 8, height = 5, dpi = 120)
ggplot2::ggsave(file.path(out_dir, "data_runs_by_format.png"), p_format, width = 8, height = 5, dpi = 120)
ggplot2::ggsave(file.path(out_dir, "data_runs_by_home.png"), p_home, width = 8, height = 5, dpi = 120)

# Fit model and posterior plots
fit <- fit_runs_ranking(
  innings = innings,
  context_cols = c("is_home", "opp_strength", "format"),
  chains = chains,
  iter = iter,
  warmup = warmup,
  seed = 42,
  control = list(adapt_delta = 0.9),
  refresh = 50L
)

mu_post <- rstan::extract(fit$fit, pars = "mu_global", permuted = TRUE)$mu_global
sig_post <- rstan::extract(fit$fit, pars = "sigma_player", permuted = TRUE)$sigma_player
phi_post <- rstan::extract(fit$fit, pars = "phi", permuted = TRUE)$phi

post_hyper <- data.frame(
  mu_global = mu_post,
  sigma_player = sig_post,
  phi = phi_post
)

p_post_mu <- ggplot2::ggplot(post_hyper, ggplot2::aes(x = mu_global)) +
  ggplot2::geom_density(fill = "#2563eb", alpha = 0.4) +
  ggplot2::labs(title = "Posterior: mu_global", x = "mu_global", y = "Density") +
  ggplot2::theme_minimal()

p_post_sigma <- ggplot2::ggplot(post_hyper, ggplot2::aes(x = sigma_player)) +
  ggplot2::geom_density(fill = "#059669", alpha = 0.4) +
  ggplot2::labs(title = "Posterior: sigma_player", x = "sigma_player", y = "Density") +
  ggplot2::theme_minimal()

p_post_phi <- ggplot2::ggplot(post_hyper, ggplot2::aes(x = phi)) +
  ggplot2::geom_density(fill = "#d97706", alpha = 0.4) +
  ggplot2::coord_cartesian(xlim = c(0, stats::quantile(post_hyper$phi, 0.99))) +
  ggplot2::labs(title = "Posterior: phi", x = "phi", y = "Density") +
  ggplot2::theme_minimal()

ggplot2::ggsave(file.path(out_dir, "posterior_mu_global.png"), p_post_mu, width = 8, height = 5, dpi = 120)
ggplot2::ggsave(file.path(out_dir, "posterior_sigma_player.png"), p_post_sigma, width = 8, height = 5, dpi = 120)
ggplot2::ggsave(file.path(out_dir, "posterior_phi.png"), p_post_phi, width = 8, height = 5, dpi = 120)

# Posterior home-field advantage effect
beta_post <- rstan::extract(fit$fit, pars = "beta", permuted = TRUE)$beta
if (!is.null(dim(beta_post))) {
  beta_names <- fit$context_info$beta_names
  idx_home <- which(beta_names == "is_home")
  if (length(idx_home) == 1L) {
    sd_home <- stats::sd(innings$is_home)
    if (is.finite(sd_home) && sd_home > 0) {
      # is_home is standardized in fit_runs_model.R, so convert back to 1-vs-0 scale.
      home_delta_log <- beta_post[, idx_home] * (1 / sd_home)
      home_ratio <- exp(home_delta_log)
      home_df <- data.frame(
        delta_log_runs = home_delta_log,
        run_ratio_home_vs_away = home_ratio
      )

      p_home_adv <- ggplot2::ggplot(home_df, ggplot2::aes(x = run_ratio_home_vs_away)) +
        ggplot2::geom_density(fill = "#7c3aed", alpha = 0.4) +
        ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "#111827") +
        ggplot2::labs(
          title = "Posterior home-field advantage",
          x = "Expected runs ratio (home / away)",
          y = "Density"
        ) +
        ggplot2::theme_minimal()

      ggplot2::ggsave(
        file.path(out_dir, "posterior_homefield_advantage.png"),
        p_home_adv,
        width = 8,
        height = 5,
        dpi = 120
      )
    }
  }
}

rank_df <- summarize_player_ranking(fit, prob = 0.9)
top_n <- min(30, nrow(rank_df))
top_df <- rank_df[seq_len(top_n), ]
top_df$player_id <- factor(top_df$player_id, levels = rev(top_df$player_id))

p_rank <- ggplot2::ggplot(top_df, ggplot2::aes(x = player_id, y = theta_mean)) +
  ggplot2::geom_point(color = "#111827", size = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = theta_qlo, ymax = theta_qhi),
    width = 0.2,
    color = "#6b7280"
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Posterior player strength (top 30 by theta_mean)",
    x = "Player",
    y = "theta"
  ) +
  ggplot2::theme_minimal()

ggplot2::ggsave(file.path(out_dir, "posterior_top30_players.png"), p_rank, width = 9, height = 8, dpi = 120)

cat("Saved plots to:", normalizePath(out_dir), "\n")
