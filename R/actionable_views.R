#' Actionable diagnostics: residuals vs age, shrinkage, trajectory clusters, next-season risk.
#'
#' Requires the same `innings` row order and columns used in [fit_runs_ranking].
#' Plots need ggplot2.

validate_fitted_innings <- function(fitted, innings) {
  if (is.null(fitted$stan_data) || is.null(fitted$player_key)) {
    stop("fitted must be from fit_runs_ranking().", call. = FALSE)
  }
  if (nrow(innings) != fitted$stan_data$N) {
    stop(
      "innings must have N = ", fitted$stan_data$N, " rows matching the fit.",
      call. = FALSE
    )
  }
  cols <- fitted$context_info$columns
  des <- stan_design_from_innings(innings, cols)
  X0 <- as.matrix(fitted$stan_data$X)
  X1 <- as.matrix(des$X)
  if (!identical(dim(X0), dim(X1))) {
    stop("Rebuilt design matrix X has wrong dimensions vs fitted$stan_data.", call. = FALSE)
  }
  diff <- max(abs(X0 - X1))
  if (!is.finite(diff) || diff > 1e-3) {
    stop(
      "Innings table does not reproduce the same scaled design matrix X as the fit ",
      "(max abs diff = ", diff, "). ",
      "Use the exact same CSV and context columns (e.g. re-run --fpca-age-k to match).",
      call. = FALSE
    )
  }
}

#' Posterior draws of linear predictor eta (log mean runs) per row, matrix S x N.
posterior_eta_draws <- function(fitted) {
  theta <- rstan::extract(fitted$fit, pars = "theta", permuted = TRUE)$theta
  beta <- rstan::extract(fitted$fit, pars = "beta", permuted = TRUE)$beta
  X <- fitted$stan_data$X
  pid <- fitted$stan_data$player_idx
  S <- nrow(theta)
  N <- nrow(X)
  eta <- matrix(NA_real_, nrow = S, ncol = N)
  lin <- beta %*% t(X)
  for (j in seq_len(N)) {
    eta[, j] <- theta[, pid[j]] + lin[, j]
  }
  eta
}

posterior_mean_expected_runs <- function(fitted) {
  eta <- posterior_eta_draws(fitted)
  colMeans(exp(eta))
}

#' Log of model-implied typical mean runs per player (theta + average X' beta for that player).
player_log_typical_mean_runs <- function(fitted) {
  theta <- rstan::extract(fitted$fit, pars = "theta", permuted = TRUE)$theta
  beta <- rstan::extract(fitted$fit, pars = "beta", permuted = TRUE)$beta
  X <- fitted$stan_data$X
  pid <- fitted$stan_data$player_idx
  P <- ncol(theta)
  log_typ <- numeric(P)
  for (j in seq_len(P)) {
    rows <- which(pid == j)
    if (!length(rows)) {
      log_typ[j] <- NA_real_
      next
    }
    xb <- colMeans(X[rows, , drop = FALSE])
    nu <- theta[, j] + as.numeric(beta %*% xb)
    log_typ[j] <- log(mean(exp(pmin(pmax(nu, -6), 12))))
  }
  names(log_typ) <- fitted$player_key$player_id
  log_typ
}

#' 1) Age vs residual (runs minus posterior mean expected runs).
plot_age_vs_residual <- function(innings, fitted, age_col = "age", out_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2.", call. = FALSE)
  }
  validate_fitted_innings(fitted, innings)
  if (!age_col %in% names(innings)) {
    stop("Column ", age_col, " not found (needed for age vs residual).", call. = FALSE)
  }
  mu_hat <- posterior_mean_expected_runs(fitted)
  df <- data.frame(
    age = as.numeric(innings[[age_col]]),
    runs = as.integer(innings$runs),
    expected = mu_hat,
    residual = innings$runs - mu_hat,
    player_id = innings$player_id
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$age, y = .data$residual)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#9ca3af") +
    ggplot2::geom_point(alpha = 0.12, color = "#1d4ed8", size = 0.8) +
    ggplot2::geom_smooth(method = "loess", se = TRUE, color = "#dc2626", linewidth = 1) +
    ggplot2::labs(
      title = "Innings residual vs age (model mean vs observed runs)",
      subtitle = "Residual = observed − posterior mean E[runs | model]; LOESS trend.",
      x = age_col,
      y = "Residual (runs)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(out_path, p, width = 9, height = 5.5, dpi = 120)
  invisible(out_path)
}

#' 2) Shrinkage: crude log(mean runs) vs model log(typical mean); distance vs n innings.
#' Writes two PNGs: `{out_prefix}_crude_vs_model.png` and `{out_prefix}_offset_vs_n.png`.
plot_shrinkage_dashboard <- function(innings, fitted, out_prefix) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2.", call. = FALSE)
  }
  validate_fitted_innings(fitted, innings)
  emp <- stats::aggregate(
    runs ~ player_id,
    data = data.frame(player_id = innings$player_id, runs = innings$runs),
    FUN = mean
  )
  n_inn <- stats::aggregate(
    runs ~ player_id,
    data = data.frame(player_id = innings$player_id, runs = innings$runs),
    FUN = length
  )
  names(n_inn)[2] <- "n_innings"
  emp$log_crude <- log(pmax(emp$runs, 0.5))
  mod <- player_log_typical_mean_runs(fitted)
  emp$log_model <- as.numeric(mod[as.character(emp$player_id)])
  emp <- merge(emp, n_inn, by = "player_id")
  emp <- emp[is.finite(emp$log_model) & is.finite(emp$log_crude), , drop = FALSE]

  p1 <- ggplot2::ggplot(emp, ggplot2::aes(x = .data$log_crude, y = .data$log_model)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#9ca3af") +
    ggplot2::geom_point(ggplot2::aes(size = .data$n_innings), alpha = 0.45, color = "#2563eb") +
    ggplot2::scale_size_area(max_size = 6) +
    ggplot2::labs(
      title = "Shrinkage: crude vs model-typical log mean runs",
      subtitle = "Each point = one player; diagonal = no shift from crude to hierarchical model.",
      x = "log(empirical mean runs per innings)",
      y = "log(model typical mean at player-mean covariates)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  emp$shrink_offset <- emp$log_model - emp$log_crude
  p2 <- ggplot2::ggplot(emp, ggplot2::aes(x = .data$n_innings, y = .data$shrink_offset)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#9ca3af") +
    ggplot2::geom_point(alpha = 0.4, color = "#059669") +
    ggplot2::geom_smooth(method = "loess", se = TRUE, color = "#111827", linewidth = 0.9) +
    ggplot2::labs(
      title = "Shrinkage offset vs sample size",
      subtitle = "log(model) − log(crude); sparse careers often move more toward the group.",
      x = "Innings in data",
      y = "Shrinkage offset (log scale)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  p1_path <- paste0(out_prefix, "_crude_vs_model.png")
  p2_path <- paste0(out_prefix, "_offset_vs_n.png")
  ggplot2::ggsave(p1_path, p1, width = 8, height = 6, dpi = 120)
  ggplot2::ggsave(p2_path, p2, width = 8, height = 6, dpi = 120)
  invisible(c(p1_path, p2_path))
}

#' 3) k-means on mean fPC scores per player; scatter + save assignments.
plot_trajectory_clusters <- function(innings, fitted, k = 4L, out_png, out_csv) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2.", call. = FALSE)
  }
  validate_fitted_innings(fitted, innings)
  fcols <- grep("^fpc_age_[0-9]+$", names(innings), value = TRUE)
  if (length(fcols) < 2L) {
    message("Skipping trajectory clusters: need at least two fpc_age_* columns on innings.")
    return(invisible(NULL))
  }
  ag <- stats::aggregate(innings[, fcols, drop = FALSE], by = list(player_id = innings$player_id), FUN = mean)
  set.seed(42L)
  km <- stats::kmeans(ag[, fcols, drop = FALSE], centers = as.integer(k), nstart = 25L)
  ag$cluster <- km$cluster
  write.csv(ag, out_csv, row.names = FALSE)

  pc1 <- fcols[[1L]]
  pc2 <- if (length(fcols) >= 2L) fcols[[2L]] else fcols[[1L]]
  p <- ggplot2::ggplot(ag, ggplot2::aes(x = .data[[pc1]], y = .data[[pc2]], color = factor(.data$cluster))) +
    ggplot2::geom_point(size = 2, alpha = 0.75) +
    ggplot2::labs(
      title = "Trajectory clusters (k-means on mean fPC scores per player)",
      subtitle = paste0("k = ", k, "; requires fpc_age_* columns from append_historical_age_curve_pc_scores()."),
      x = pc1,
      y = pc2,
      color = "Cluster"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(out_png, p, width = 8.5, height = 6, dpi = 120)
  invisible(ag)
}

#' 4) Monte Carlo P(sum of n_future innings runs > threshold) per player.
#'
#' Uses player-average covariate row `xb` and [MASS::rnegbin] matching Stan's NegBin2
#' (mean = exp(eta), variance = mu + mu^2/phi).
next_season_threshold_table <- function(
    fitted,
    innings,
    n_future = 20L,
    threshold = 350,
    n_draws = 1500L,
    seed = 42L
) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Install MASS (usually installed with R) for next-season simulation.", call. = FALSE)
  }
  validate_fitted_innings(fitted, innings)
  set.seed(seed)
  phi <- as.vector(rstan::extract(fitted$fit, pars = "phi", permuted = TRUE)$phi)
  theta <- rstan::extract(fitted$fit, pars = "theta", permuted = TRUE)$theta
  beta <- rstan::extract(fitted$fit, pars = "beta", permuted = TRUE)$beta
  X <- fitted$stan_data$X
  pid <- fitted$stan_data$player_idx
  S <- nrow(theta)
  P <- ncol(theta)
  players <- fitted$player_key$player_id
  idx_draw <- sample.int(S, size = min(as.integer(n_draws), S), replace = TRUE)

  rows_out <- vector("list", P)
  for (j in seq_len(P)) {
    rows <- which(pid == j)
    n_hist <- length(rows)
    if (!n_hist) {
      next
    }
    xb <- colMeans(X[rows, , drop = FALSE])
    totals <- numeric(length(idx_draw))
    for (ii in seq_along(idx_draw)) {
      s <- idx_draw[ii]
      mu <- exp(min(12, max(-6, theta[s, j] + sum(beta[s, ] * xb))))
      yf <- MASS::rnegbin(n_future, mu = mu, theta = phi[s])
      totals[ii] <- sum(yf)
    }
    rows_out[[j]] <- data.frame(
      player_id = players[j],
      n_innings_hist = n_hist,
      n_future_innings = as.integer(n_future),
      threshold = as.numeric(threshold),
      P_total_gt_threshold = mean(totals > threshold),
      median_sim_total = stats::median(totals),
      sim_q10 = as.numeric(stats::quantile(totals, 0.1)),
      sim_q90 = as.numeric(stats::quantile(totals, 0.9)),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows_out)
  rownames(out) <- NULL
  out
}

#' Ensure a season-like grouping column exists.
ensure_season_column <- function(innings, season_col = "season") {
  if (season_col %in% names(innings)) {
    innings[[season_col]] <- as.character(innings[[season_col]])
    return(list(innings = innings, col = season_col, note = NULL))
  }
  if ("match_date" %in% names(innings)) {
    d <- as.Date(innings$match_date)
    innings$season <- format(d, "%Y")
    return(list(innings = innings, col = "season", note = NULL))
  }
  if (!"age" %in% names(innings)) {
    stop("Need column `season`, `match_date`, or `age` to build season grouping.", call. = FALSE)
  }
  innings$season <- paste0("age_", floor(as.numeric(innings$age)))
  list(
    innings = innings,
    col = "season",
    note = "Synthetic season from floor(age); add a real `season` column for production."
  )
}

#' Observed vs posterior-mean expected total runs per player-season.
plot_season_observed_vs_projected <- function(innings, fitted, out_path, season_col = "season") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2.", call. = FALSE)
  }
  validate_fitted_innings(fitted, innings)
  es <- ensure_season_column(innings, season_col)
  inn <- es$innings
  sc <- es$col
  mu_hat <- posterior_mean_expected_runs(fitted)
  df <- data.frame(
    player_id = inn$player_id,
    season = inn[[sc]],
    runs = inn$runs,
    expected = mu_hat
  )
  agg <- stats::aggregate(
    cbind(runs, expected) ~ player_id + season,
    data = df,
    FUN = sum
  )
  n_players <- length(unique(agg$player_id))
  if (n_players > 40L) {
    top_players <- names(sort(table(inn$player_id), decreasing = TRUE)[seq_len(40L)])
    agg <- agg[agg$player_id %in% top_players, , drop = FALSE]
  }
  agg$player_id <- factor(agg$player_id)
  p <- ggplot2::ggplot(agg, ggplot2::aes(x = .data$season, y = .data$runs, group = .data$player_id)) +
    ggplot2::geom_line(color = "#94a3b8", linewidth = 0.3, alpha = 0.5) +
    ggplot2::geom_point(color = "#1d4ed8", size = 0.8, alpha = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = .data$expected), color = "#dc2626", linewidth = 0.35, alpha = 0.6) +
    ggplot2::geom_point(ggplot2::aes(y = .data$expected), color = "#dc2626", size = 0.7, alpha = 0.35) +
    ggplot2::facet_wrap(~ player_id, scales = "free_y", ncol = 5) +
    ggplot2::labs(
      title = "Season total runs: observed (blue) vs model expected (red)",
      subtitle = if (!is.null(es$note) && nzchar(es$note)) {
        es$note
      } else {
        paste0("Grouping column: ", sc)
      },
      x = NULL,
      y = "Total runs (sum within player-season)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 6)
    )
  ggplot2::ggsave(out_path, p, width = 14, height = 10, dpi = 120)
  invisible(out_path)
}

`%||%` <- function(x, y) {
  if (is.null(x) || (is.character(x) && !nzchar(x))) {
    return(y)
  }
  x
}

#' Run all actionable exports (best-effort if optional columns missing).
run_all_actionable_views <- function(
    fitted,
    innings,
    out_dir,
    age_col = "age",
    season_col = "season",
    cluster_k = 4L,
    n_future = 20L,
    threshold = 350,
    n_draws = 1500L
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2.", call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plot_dir <- file.path(out_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  validate_fitted_innings(fitted, innings)

  if (age_col %in% names(innings)) {
    plot_age_vs_residual(
      innings,
      fitted,
      age_col = age_col,
      out_path = file.path(plot_dir, "age_vs_residual.png")
    )
  } else {
    message("Skipping age_vs_residual: column `", age_col, "` not found.")
  }

  plot_shrinkage_dashboard(
    innings,
    fitted,
    out_prefix = file.path(plot_dir, "shrinkage")
  )

  traj <- plot_trajectory_clusters(
    innings,
    fitted,
    k = cluster_k,
    out_png = file.path(plot_dir, "trajectory_clusters.png"),
    out_csv = file.path(out_dir, "trajectory_clusters.csv")
  )

  ns <- next_season_threshold_table(
    fitted,
    innings,
    n_future = n_future,
    threshold = threshold,
    n_draws = n_draws
  )
  write.csv(ns, file.path(out_dir, "next_season_threshold.csv"), row.names = FALSE)

  es <- ensure_season_column(innings, season_col)
  plot_season_observed_vs_projected(
    es$innings,
    fitted,
    out_path = file.path(plot_dir, "season_observed_vs_projected.png"),
    season_col = es$col
  )

  invisible(list(
    plots_dir = normalizePath(plot_dir),
    trajectory = traj,
    next_season = ns
  ))
}
