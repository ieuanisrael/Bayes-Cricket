#' Fit hierarchical negative-binomial model for runs per innings (rstan).
#'
#' @param innings data.frame with at least: player_id, runs (integer >= 0).
#'        Optional columns passed to `context_cols` (numeric or factor).
#' @param context_cols character vector of column names for contextual predictors.
#'        Factors are converted to treatment contrasts (first level = reference).
#' @param stan_file path to `runs_nb_hierarchical.stan`. Default: `./stan/runs_nb_hierarchical.stan` under `getwd()`.
#' @param chains, iter, warmup passed to `rstan::stan`.
#' @param ... additional arguments to `rstan::stan`.
#' @return list with `fit` (stanfit), `player_key` (id -> name), `context_info`.
#' @export
fit_runs_ranking <- function(
    innings,
    context_cols = character(),
    stan_file = NULL,
    chains = 4L,
    iter = 2000L,
    warmup = 1000L,
    ...
) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Install rstan: install.packages(\"rstan\")", call. = FALSE)
  }

  req <- c("player_id", "runs")
  miss <- setdiff(req, names(innings))
  if (length(miss)) {
    stop("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  innings <- as.data.frame(innings)
  innings$player_id <- as.character(innings$player_id)
  players <- sort(unique(innings$player_id))
  n_players <- length(players)
  pid_to_idx <- setNames(seq_len(n_players), players)

  y <- as.integer(innings$runs)
  if (any(y < 0 | is.na(y))) {
    stop("runs must be non-negative integers without NA.", call. = FALSE)
  }

  player_idx <- as.integer(pid_to_idx[innings$player_id])

  J <- 1L
  X <- matrix(0, nrow = nrow(innings), ncol = 1L)
  beta_names <- "_intercept_padding"

  if (length(context_cols)) {
    unknown <- setdiff(context_cols, names(innings))
    if (length(unknown)) {
      stop("Unknown context_cols: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    X_list <- list()
    for (nm in context_cols) {
      v <- innings[[nm]]
      if (is.factor(v) || is.character(v)) {
        v <- factor(v)
        mm <- model.matrix(~ 0 + v)
        colnames(mm) <- paste0(nm, colnames(mm))
        X_list[[length(X_list) + 1L]] <- mm
      } else {
        x <- as.numeric(v)
        if (anyNA(x)) {
          stop("Numeric context column has NA: ", nm, call. = FALSE)
        }
        x <- as.matrix(scale(x, center = TRUE, scale = TRUE))
        colnames(x) <- nm
        X_list[[length(X_list) + 1L]] <- x
      }
    }
    X <- do.call(cbind, X_list)
    beta_names <- colnames(X)
    J <- ncol(X)
  }

  stan_data <- list(
    N = nrow(innings),
    y = y,
    n_players = n_players,
    player_idx = player_idx,
    J = J,
    X = X
  )

  if (is.null(stan_file)) {
    stan_file <- file.path(getwd(), "stan", "runs_nb_hierarchical.stan")
  }
  if (!file.exists(stan_file)) {
    stop("Stan file not found: ", stan_file, call. = FALSE)
  }

  fit <- rstan::stan(
    file = stan_file,
    data = stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    ...
  )

  player_key <- data.frame(
    player_id = players,
    player_index = seq_len(n_players),
    stringsAsFactors = FALSE
  )

  list(
    fit = fit,
    player_key = player_key,
    context_info = list(
      columns = context_cols,
      beta_names = beta_names,
      J = J
    ),
    stan_data = stan_data
  )
}

#' Posterior summaries for player log-ability `theta` (higher = stronger batting).
#'
#' @param fitted output of `fit_runs_ranking`.
#' @param prob mass for central credible interval.
#' @return data.frame sorted by posterior mean `theta_mean` descending.
summarize_player_ranking <- function(fitted, prob = 0.9) {
  theta <- rstan::extract(fitted$fit, pars = "theta", permuted = TRUE)$theta
  lo <- (1 - prob) / 2
  hi <- 1 - lo
  mu <- colMeans(theta)
  qlo <- apply(theta, 2, stats::quantile, probs = lo)
  qhi <- apply(theta, 2, stats::quantile, probs = hi)

  out <- data.frame(
    player_id = fitted$player_key$player_id,
    theta_mean = mu,
    theta_qlo = qlo,
    theta_qhi = qhi,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$theta_mean, decreasing = TRUE), ]
  rownames(out) <- NULL
  out
}

#' Probability each player has the highest theta (pairwise MC comparison).
best_player_probabilities <- function(fitted) {
  ex <- rstan::extract(fitted$fit, pars = "theta", permuted = TRUE)$theta
  # iterations x n_players
  best_counts <- apply(ex, 1, function(row) which.max(row))
  tab <- tabulate(best_counts, nbins = ncol(ex))
  out <- data.frame(
    player_id = fitted$player_key$player_id,
    P_best = as.numeric(tab) / nrow(ex),
    stringsAsFactors = FALSE
  )
  out[order(out$P_best, decreasing = TRUE), ]
}
