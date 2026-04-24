#' Historical age–run curves and PCA scores (functional predictors for Stan).
#'
#' For each innings row, builds a **retrospective** smooth of that batter’s
#' **past** `(age, runs)` pairs only, evaluated on a common age grid, then runs
#' **PCA** on the stacked curves. The first `K` principal component scores are
#' appended as numeric columns `fpc_age_1`, …, `fpc_age_K` suitable for
#' [fit_runs_ranking] `context_cols`.
#'
#' This is **discretized functional PCA** (PCA on a grid of function values), a
#' standard approximation; for irregular sparse grids you can swap the inner PCA
#' for `fdapace::FPCA` and merge scores by row.
#'
#' **Ordering:** Rows must be in **chronological order within each player** (or
#' pass `time_order_col`). If order is wrong, “history” is wrong.
#'
#' **Leakage:** Scores for row `r` use only innings with index `< r` for that
#' player, so they are usable as **prospective** inputs for that innings’ outcome.
#'
#' @param innings `data.frame` with `player_id`, `runs`, and `age` (numeric).
#' @param age_col,runs_col,player_col column names.
#' @param time_order_col optional numeric or Date; within player, smaller = earlier.
#'        If `NULL`, uses current row order after sorting by `player_col`.
#' @param K number of PC scores to retain (default 3).
#' @param n_age_grid number of equally spaced evaluation points between observed
#'        min and max age (clamped to sensible range).
#' @return `innings` with added columns `fpc_age_1` … `fpc_age_K` (same row order as input).
append_historical_age_curve_pc_scores <- function(
    innings,
    age_col = "age",
    runs_col = "runs",
    player_col = "player_id",
    time_order_col = NULL,
    K = 3L,
    n_age_grid = 40L
) {
  if (!age_col %in% names(innings) || !runs_col %in% names(innings) ||
    !player_col %in% names(innings)) {
    stop("innings must contain age_col, runs_col, player_col.", call. = FALSE)
  }
  df <- as.data.frame(innings)
  df[[player_col]] <- as.character(df[[player_col]])
  df[[age_col]] <- as.numeric(df[[age_col]])
  df[[runs_col]] <- as.numeric(df[[runs_col]])
  if (anyNA(df[[age_col]]) || anyNA(df[[runs_col]])) {
    stop("age and runs must be non-missing.", call. = FALSE)
  }

  if (is.null(time_order_col)) {
    o <- order(df[[player_col]], seq_len(nrow(df)))
  } else {
    if (!time_order_col %in% names(df)) {
      stop("time_order_col not found: ", time_order_col, call. = FALSE)
    }
    o <- order(df[[player_col]], df[[time_order_col]], seq_len(nrow(df)))
  }
  df <- df[o, ]
  rownames(df) <- NULL

  age_min <- max(stats::quantile(df[[age_col]], 0.01), min(df[[age_col]]))
  age_max <- min(stats::quantile(df[[age_col]], 0.99), max(df[[age_col]]))
  if (age_max <= age_min) {
    stop("Insufficient age variation for a grid.", call. = FALSE)
  }
  t_grid <- seq(age_min, age_max, length.out = as.integer(n_age_grid))

  # Fallback curve when a player has no history: pooled (age, runs) smooth on grid.
  global_y <- stats::approx(
    df[[age_col]],
    df[[runs_col]],
    xout = t_grid,
    rule = 2
  )$y

  players <- split(df, df[[player_col]])
  curve_list <- vector("list", length(players))
  names(curve_list) <- names(players)

  for (nm in names(players)) {
    P <- players[[nm]]
    n <- nrow(P)
    M <- matrix(NA_real_, nrow = n, ncol = length(t_grid))
    for (r in seq_len(n)) {
      if (r == 1L) {
        M[r, ] <- global_y
      } else {
        ag <- P[[age_col]][seq_len(r - 1L)]
        rs <- P[[runs_col]][seq_len(r - 1L)]
        ok <- is.finite(ag) & is.finite(rs)
        ag <- ag[ok]
        rs <- rs[ok]
        ord <- order(ag)
        ag <- ag[ord]
        rs <- rs[ord]
        dfu <- data.frame(ag = ag, rs = rs)
        dfu <- stats::aggregate(rs ~ ag, data = dfu, FUN = mean)
        ag <- dfu$ag
        rs <- dfu$rs
        # duplicate ages / too few distinct support points → constant extension
        if (length(ag) < 2L || length(unique(ag)) < 2L) {
          mu0 <- if (length(rs)) mean(rs) else NA_real_
          M[r, ] <- rep(if (is.finite(mu0)) mu0 else mean(global_y, na.rm = TRUE), length(t_grid))
        } else {
          M[r, ] <- stats::approx(ag, rs, xout = t_grid, rule = 2)$y
        }
      }
    }
    curve_list[[nm]] <- M
  }

  # Reassemble in same row order as df
  curve_mat <- matrix(NA_real_, nrow = nrow(df), ncol = length(t_grid))
  idx <- 0L
  for (nm in names(players)) {
    P <- players[[nm]]
    nr <- nrow(P)
    curve_mat[(idx + 1L):(idx + nr), ] <- curve_list[[nm]]
    idx <- idx + nr
  }

  K <- as.integer(K)
  if (K < 1L) {
    stop("K must be >= 1.", call. = FALSE)
  }
  K <- min(K, ncol(curve_mat), nrow(curve_mat))

  pc <- stats::prcomp(curve_mat, center = TRUE, scale. = FALSE)
  n_comp <- min(K, ncol(pc$x), nrow(pc$x))
  scores <- pc$x[, seq_len(n_comp), drop = FALSE]
  cn <- paste0("fpc_age_", seq_len(n_comp))
  colnames(scores) <- cn

  out <- df
  for (j in seq_len(n_comp)) {
    out[[cn[j]]] <- scores[, j]
  }
  # Z-score each score column so magnitudes match other scaled numeric predictors in Stan.
  for (nm in cn) {
    v <- out[[nm]]
    s <- stats::sd(v, na.rm = TRUE)
    out[[nm]] <- if (!is.finite(s) || s < 1e-12) {
      rep(0, length(v))
    } else {
      as.numeric(scale(v))
    }
  }

  # Restore original row order of `innings`
  inv <- integer(length(o))
  inv[o] <- seq_along(o)
  out[inv, , drop = FALSE]
}
