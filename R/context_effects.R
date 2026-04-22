#' Summarise and plot posterior `beta` (context / covariate effects on log mean runs).
#'
#' The model is `neg_binomial_2_log(eta, phi)` with `eta = theta[player] + X %*% beta`.
#' Coefficients apply to **scaled numeric** predictors and **treatment-contrast**
#' factor columns (see [fit_runs_ranking]).
#'
#' @param fitted Output of [fit_runs_ranking].
#' @param prob Central credible mass for `beta` and for `exp(beta)`.
#' @return `data.frame` with one row per `beta` component.
#' @noRd
summarize_context_betas <- function(fitted, prob = 0.9) {
  if (is.null(fitted$context_info$beta_names)) {
    stop("fitted$context_info$beta_names missing.", call. = FALSE)
  }
  bn <- fitted$context_info$beta_names
  if (identical(bn, "_intercept_padding") || !length(bn)) {
    return(data.frame(
      term = character(),
      beta_mean = numeric(),
      beta_lo = numeric(),
      beta_hi = numeric(),
      exp_beta_mean = numeric(),
      exp_beta_lo = numeric(),
      exp_beta_hi = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  b <- rstan::extract(fitted$fit, pars = "beta", permuted = TRUE)$beta
  if (is.null(b) || !length(dim(b)) || dim(b)[2] != length(bn)) {
    stop("beta posterior dimension mismatch to beta_names.", call. = FALSE)
  }
  lo <- (1 - prob) / 2
  hi <- 1 - lo
  rows <- lapply(seq_along(bn), function(j) {
    v <- b[, j]
    ev <- exp(v)
    data.frame(
      term = bn[[j]],
      beta_mean = mean(v),
      beta_lo = stats::quantile(v, lo),
      beta_hi = stats::quantile(v, hi),
      exp_beta_mean = mean(ev),
      exp_beta_lo = stats::quantile(ev, lo),
      exp_beta_hi = stats::quantile(ev, hi),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Write plots for each context coefficient: forest (log scale) and density, plus exp scale.
#' @param fitted From [fit_runs_ranking].
#' @param out_dir Directory for PNG files (created if needed).
#' @param prob Credible level for error bars and ribbons.
#' @return Invisibly, path to `out_dir`.
#' @noRd
plot_context_effects <- function(fitted, out_dir, prob = 0.9) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 for context effect plots: install.packages(\"ggplot2\")", call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  sm <- summarize_context_betas(fitted, prob = prob)
  if (!nrow(sm)) {
    message("No context columns in fit (or only padding); skipping context effect plots.")
    return(invisible(out_dir))
  }

  sm$term <- factor(sm$term, levels = sm$term)

  p_forest <- ggplot2::ggplot(
    sm,
    ggplot2::aes(
      x = .data$term,
      y = .data$beta_mean,
      ymin = .data$beta_lo,
      ymax = .data$beta_hi
    )
  ) +
    ggplot2::geom_pointrange(color = "#1d4ed8", linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#9ca3af") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Context effects on log expected runs (posterior mean & interval)",
      subtitle = "Numeric predictors: one SD. Factors: level vs reference.",
      x = NULL,
      y = "beta (log scale, additive in linear predictor)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  p_exp <- ggplot2::ggplot(
    sm,
    ggplot2::aes(
      x = .data$term,
      y = .data$exp_beta_mean,
      ymin = .data$exp_beta_lo,
      ymax = .data$exp_beta_hi
    )
  ) +
    ggplot2::geom_pointrange(color = "#059669", linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#9ca3af") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Multiplicative effect on expected runs (exp(beta))",
      y = "Multiplier (1 = no effect on mean)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  b <- rstan::extract(fitted$fit, pars = "beta", permuted = TRUE)$beta
  bn <- fitted$context_info$beta_names
  long <- data.frame(
    value = as.vector(b),
    term = rep(factor(bn, levels = bn), each = nrow(b))
  )
  p_density <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value, fill = .data$term)) +
    ggplot2::geom_density(alpha = 0.35) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#6b7280") +
    ggplot2::facet_wrap(~ term, scales = "free", ncol = 2) +
    ggplot2::labs(
      title = "Posterior of each beta (context effects, log link)",
      x = "beta",
      y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )

  ggplot2::ggsave(
    file.path(out_dir, "context_effects_forest_log.png"),
    p_forest,
    width = 9,
    height = max(4, 0.3 * nrow(sm) + 2),
    dpi = 120
  )
  ggplot2::ggsave(
    file.path(out_dir, "context_effects_forest_exp.png"),
    p_exp,
    width = 9,
    height = max(4, 0.3 * nrow(sm) + 2),
    dpi = 120
  )
  ggplot2::ggsave(
    file.path(out_dir, "context_effects_densities.png"),
    p_density,
    width = 10,
    height = min(3 + 2.5 * ceiling(nrow(sm) / 2), 24),
    dpi = 120
  )

  invisible(normalizePath(out_dir))
}
