#' Fit hierarchical negative-binomial model for runs per innings (rstan).
#'
#' @param innings data.frame with at least: player_id, runs (integer >= 0).
#'        Optional columns passed to `context_cols` (numeric or factor).
#' @param context_cols character vector of column names for contextual predictors.
#'        Factors are converted to treatment contrasts (first level = reference).
#' @param stan_file path to `runs_nb_hierarchical.stan`. Default: `./stan/runs_nb_hierarchical.stan` under `getwd()`.
#' @param chains, iter, warmup passed to `rstan::stan`.
#' @param ... additional arguments to `rstan::stan`.
#' Build Stan design matrices from an innings table (same logic as [fit_runs_ranking]).
#'
#' @param innings `data.frame` with `player_id`, `runs`, and any `context_cols`.
#' @param context_cols character vector of predictor column names.
#' @return List with `innings`, `players`, `player_idx`, `y`, `X`, `J`, `beta_names`, `n_players`.
#' @export
stan_design_from_innings <- function(innings, context_cols = character()) {
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

  list(
    innings = innings,
    players = players,
    n_players = n_players,
    player_idx = player_idx,
    y = y,
    X = X,
    J = J,
    beta_names = beta_names
  )
}

#' @return list with `fit` (stanfit), `player_key` (id -> name), `context_info`.
#' @export
fit_runs_ranking <- function(
    innings,
    context_cols = character(),
    stan_file = NULL,
    chains = 1L,
    iter = 200L,
    warmup = 100L,
    ...
) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Install rstan: install.packages(\"rstan\")", call. = FALSE)
  }

  des <- stan_design_from_innings(innings, context_cols)
  innings <- des$innings
  players <- des$players
  n_players <- des$n_players
  player_idx <- des$player_idx
  y <- des$y
  X <- des$X
  J <- des$J
  beta_names <- des$beta_names

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

  # macOS: rstan needs the Command Line Tools SDK on the compiler line; the repo
  # ships `.R/Makevars` for that. If unset, apply it when present (same as R/run_example.R).
  if (!nzchar(Sys.getenv("R_MAKEVARS_USER", ""))) {
    proj_root <- dirname(dirname(normalizePath(stan_file, winslash = "/", mustWork = TRUE)))
    makevars_proj <- file.path(proj_root, ".R", "Makevars")
    if (file.exists(makevars_proj)) {
      Sys.setenv(R_MAKEVARS_USER = normalizePath(makevars_proj, winslash = "/"))
    }
  }

  fit <- rstan::stan(
    file = stan_file,
    data = stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    ...
  )

  ex <- suppressWarnings(tryCatch(
    rstan::extract(fit, pars = "mu_global", permuted = TRUE)$mu_global,
    error = function(e) NULL
  ))
  if (is.null(ex) || !length(ex)) {
    stop(
      "Stan produced no samples (sampling failed or was interrupted). ",
      "Check printed errors (e.g. divergences, RNG overflow in generated quantities). ",
      "Try more chains, higher adapt_delta, tighter priors on beta, or fewer predictors.",
      call. = FALSE
    )
  }

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

#' Row indices for a player's most recent `k` innings.
#'
#' Uses `time_col` ascending so "recent" means latest in time. If `time_col` is
#' missing, uses the order rows appear in `innings` for that player (assumed
#' chronological).
#'
#' @param innings data.frame with `player_id`.
#' @param player_id character, one player.
#' @param k integer, number of innings.
#' @param time_col optional column name for ordering (numeric, Date, or POSIXct).
#' @return integer vector of row indices in `innings` (1..nrow), length <= k.
recent_innings_row_indices <- function(innings, player_id, k, time_col = NULL) {
  rows <- which(innings$player_id == player_id)
  if (!length(rows)) {
    return(integer())
  }
  if (!is.null(time_col)) {
    if (!time_col %in% names(innings)) {
      stop("time_col not found: ", time_col, call. = FALSE)
    }
    t <- innings[[time_col]][rows]
    if (anyNA(t)) {
      stop("time_col has NA for some rows.", call. = FALSE)
    }
    rows <- rows[order(t)]
  }
  n_take <- min(as.integer(k), length(rows))
  rows[(length(rows) - n_take + 1L):length(rows)]
}

#' Compare recent innings to in-sample posterior predictive (`y_rep` from Stan).
#'
#' Requires the same `innings` data.frame (same row order) used in
#' [fit_runs_ranking]. For each posterior draw, replicates `y_rep` are
#' independent across innings conditional on parameters; the function forms the
#' posterior predictive distribution of the **sum** or **mean** of runs over
#' the selected innings and compares it to the observed statistic.
#'
#' @param fitted output of [fit_runs_ranking].
#' @param innings same `data.frame` passed to [fit_runs_ranking].
#' @param k number of recent innings per player (default 10).
#' @param time_col passed to [recent_innings_row_indices]; NULL uses row order.
#' @param stat `"sum"` or `"mean"` over the k innings.
#' @param prob central interval mass for predictive quantiles (e.g. 0.9).
#' @return data.frame with one row per player: observed stat, predictive quantiles,
#'   and `ppc_p_upper` = posterior probability that replicate stat >= observed
#'   (one-sided predictive p-value for "more extreme high" streaks).
#' @export
compare_recent_vs_posterior_predictive <- function(
    fitted,
    innings,
    k = 10L,
    time_col = NULL,
    stat = c("sum", "mean"),
    prob = 0.9
) {
  stat <- match.arg(stat)
  if (!inherits(fitted, "list") || is.null(fitted$fit) || is.null(fitted$stan_data)) {
    stop("fitted must be the return value of fit_runs_ranking().", call. = FALSE)
  }
  N <- fitted$stan_data$N
  if (nrow(innings) != N) {
    stop(
      "innings must have the same number of rows as in the fit (",
      N, "); got ", nrow(innings), ".",
      call. = FALSE
    )
  }

  y_rep <- rstan::extract(fitted$fit, pars = "y_rep", permuted = TRUE)$y_rep
  if (is.null(y_rep) || ncol(y_rep) != N) {
    stop("Could not extract y_rep with ncol equal to N.", call. = FALSE)
  }

  players <- fitted$player_key$player_id
  lo <- (1 - prob) / 2
  hi <- 1 - lo

  out <- lapply(players, function(pid) {
    idx <- recent_innings_row_indices(innings, pid, k, time_col)
    if (!length(idx)) {
      return(NULL)
    }
    obs_runs <- innings$runs[idx]
    if (stat == "sum") {
      obs <- sum(obs_runs)
      pred <- rowSums(y_rep[, idx, drop = FALSE])
    } else {
      obs <- mean(obs_runs)
      pred <- rowMeans(y_rep[, idx, drop = FALSE])
    }
    ppc_p_upper <- mean(pred >= obs)
    data.frame(
      player_id = pid,
      n_innings_used = length(idx),
      stat = stat,
      observed = obs,
      pred_mean = mean(pred),
      pred_qlo = stats::quantile(pred, probs = lo),
      pred_qhi = stats::quantile(pred, probs = hi),
      ppc_p_upper = ppc_p_upper,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}
