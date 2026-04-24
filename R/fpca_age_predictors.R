#' Season-level age curves and PCA scores (functional predictors for Stan).
#'
#' Innings are aggregated to **player–season** summaries (`mean_runs`, `mean_age`,
#' etc.). For each season, a retrospective curve uses only **prior seasons** for
#' that player: points `(mean_age, statistic)` are interpolated on a common age
#' grid, then **PCA** across all player–season rows produces `fpc_age_*` scores.
#' Those scores are merged back onto **every innings row** in that player–season.
#'
#' Requires a **`season`** (or custom) column on `innings`. If unavailable, set
#' `synthetic_season_from_age = TRUE` to use `paste0("Y", floor(age))` as a demo
#' label (not a real competition season).
#'
#' @param innings `data.frame` with `player_id`, `runs`, `age`, and `season_col`.
#' @param season_col column name for season identifier (character or factor).
#' @param statistic `"mean"` or `"median"` innings runs within the season for the curve.
#' @param age_col,runs_col,player_col column names.
#' @param synthetic_season_from_age if `TRUE` and `season_col` missing, creates
#'        `innings$season` from `floor(age)`.
#' @param K number of PC scores.
#' @param n_age_grid evaluation points on age axis.
#' @return `innings` with added `fpc_age_1` … (same row order as input).
append_historical_age_curve_pc_scores <- function(
    innings,
    season_col = "season",
    age_col = "age",
    runs_col = "runs",
    player_col = "player_id",
    statistic = c("mean", "median"),
    synthetic_season_from_age = FALSE,
    K = 3L,
    n_age_grid = 40L
) {
  statistic <- match.arg(statistic)
  if (!age_col %in% names(innings) || !runs_col %in% names(innings) ||
    !player_col %in% names(innings)) {
    stop("innings must contain age_col, runs_col, player_col.", call. = FALSE)
  }
  if (!season_col %in% names(innings)) {
    if (isTRUE(synthetic_season_from_age) && age_col %in% names(innings)) {
      innings <- as.data.frame(innings)
      innings[[season_col]] <- paste0("Y", floor(as.numeric(innings[[age_col]])))
    } else {
      stop(
        "Missing column `", season_col, "`. Add a real season id, or pass ",
        "synthetic_season_from_age = TRUE for a floor(age) fallback (demo only).",
        call. = FALSE
      )
    }
  }

  df <- as.data.frame(innings)
  df[[player_col]] <- as.character(df[[player_col]])
  df[[season_col]] <- as.character(df[[season_col]])
  df[[age_col]] <- as.numeric(df[[age_col]])
  df[[runs_col]] <- as.numeric(df[[runs_col]])
  if (anyNA(df[[age_col]]) || anyNA(df[[runs_col]])) {
    stop("age and runs must be non-missing.", call. = FALSE)
  }

  o_inn <- order(df[[player_col]], df[[season_col]], seq_len(nrow(df)))
  inv_inn <- integer(length(o_inn))
  inv_inn[o_inn] <- seq_along(o_inn)

  ps <- player_season_summary(
    df,
    season_col = season_col,
    age_col = age_col,
    runs_col = runs_col,
    player_col = player_col,
    statistic = statistic
  )
  if (nrow(ps) < 2L) {
    stop("Need at least two player–season rows after aggregation.", call. = FALSE)
  }

  season_levels <- sort_season_levels(unique(as.character(ps$season)))
  ps$season_rank <- match(ps$season, season_levels)
  if (anyNA(ps$season_rank)) {
    stop("Season values could not be ranked; check season_col format.", call. = FALSE)
  }
  ps <- ps[order(ps[[player_col]], ps$season_rank), , drop = FALSE]
  rownames(ps) <- NULL
  ps$i_ps <- seq_len(nrow(ps))

  age_min <- max(stats::quantile(ps$mean_age, 0.01, na.rm = TRUE), min(ps$mean_age))
  age_max <- min(stats::quantile(ps$mean_age, 0.99, na.rm = TRUE), max(ps$mean_age))
  if (age_max <= age_min) {
    stop("Insufficient age variation across seasons.", call. = FALSE)
  }
  t_grid <- seq(age_min, age_max, length.out = as.integer(n_age_grid))

  y_stat <- if (statistic == "mean") {
    ps$mean_runs
  } else {
    ps$median_runs
  }
  global_y <- stats::approx(
    ps$mean_age,
    y_stat,
    xout = t_grid,
    rule = 2
  )$y

  curve_mat_ps <- matrix(NA_real_, nrow = nrow(ps), ncol = length(t_grid))
  players_ps <- split(ps, ps[[player_col]])
  for (nm in names(players_ps)) {
    P <- players_ps[[nm]]
    P <- P[order(P$season_rank), , drop = FALSE]
    n <- nrow(P)
    M <- matrix(NA_real_, nrow = n, ncol = length(t_grid))
    for (r in seq_len(n)) {
      if (r == 1L) {
        M[r, ] <- global_y
      } else {
        past <- P[seq_len(r - 1L), , drop = FALSE]
        ag <- past$mean_age
        rs <- if (statistic == "mean") past$mean_runs else past$median_runs
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
        if (length(ag) < 2L || length(unique(ag)) < 2L) {
          mu0 <- if (length(rs)) mean(rs) else NA_real_
          M[r, ] <- rep(if (is.finite(mu0)) mu0 else mean(global_y, na.rm = TRUE), length(t_grid))
        } else {
          M[r, ] <- stats::approx(ag, rs, xout = t_grid, rule = 2)$y
        }
      }
    }
    idx_rows <- P$i_ps
    curve_mat_ps[idx_rows, ] <- M
  }

  K <- as.integer(K)
  if (K < 1L) {
    stop("K must be >= 1.", call. = FALSE)
  }
  K <- min(K, ncol(curve_mat_ps), nrow(curve_mat_ps))
  pc <- stats::prcomp(curve_mat_ps, center = TRUE, scale. = FALSE)
  n_comp <- min(K, ncol(pc$x), nrow(pc$x))
  scores_ps <- pc$x[, seq_len(n_comp), drop = FALSE]
  cn <- paste0("fpc_age_", seq_len(n_comp))
  colnames(scores_ps) <- cn
  ps_scores <- cbind(ps[, c(player_col, season_col), drop = FALSE], scores_ps)

  df$.row_original <- seq_len(nrow(df))
  mrg <- merge(
    df,
    ps_scores,
    by.x = c(player_col, season_col),
    by.y = c(player_col, season_col),
    all.x = TRUE,
    sort = FALSE
  )
  mrg <- mrg[order(mrg$.row_original), , drop = FALSE]
  mrg$.row_original <- NULL
  mrg$i_ps <- NULL
  for (nm in cn) {
    if (anyNA(mrg[[nm]])) {
      mrg[[nm]][is.na(mrg[[nm]])] <- 0
    }
  }
  for (nm in cn) {
    v <- mrg[[nm]]
    s <- stats::sd(v, na.rm = TRUE)
    mrg[[nm]] <- if (!is.finite(s) || s < 1e-12) {
      rep(0, nrow(mrg))
    } else {
      as.numeric(scale(v))
    }
  }

  mrg[inv_inn, , drop = FALSE]
}

#' One row per player–season with summary stats for curve building.
#'
#' @param statistic Which within-season innings statistic traces the age curve.
#' @keywords internal
player_season_summary <- function(
    innings,
    season_col,
    age_col = "age",
    runs_col = "runs",
    player_col = "player_id",
    statistic = c("mean", "median")
) {
  statistic <- match.arg(statistic)
  inn <- as.data.frame(innings)
  key <- interaction(inn[[player_col]], inn[[season_col]], drop = TRUE, sep = "|||")
  spl <- split(seq_len(nrow(inn)), key)
  rows <- lapply(spl, function(ii) {
    chunk <- inn[ii, , drop = FALSE]
    rs <- chunk[[runs_col]]
    data.frame(
      player_id = as.character(chunk[[player_col]][1L]),
      season = as.character(chunk[[season_col]][1L]),
      mean_runs = mean(rs),
      median_runs = stats::median(rs),
      sd_runs = if (length(rs) > 1L) stats::sd(rs) else 0,
      n_innings = length(rs),
      mean_age = mean(chunk[[age_col]]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Chronological order of season labels (global calendar).
sort_season_levels <- function(u) {
  u <- unique(as.character(u))
  nu <- suppressWarnings(as.integer(u))
  if (length(u) == length(nu) && all(!is.na(nu))) {
    return(u[order(nu)])
  }
  if (all(grepl("^[0-9]{4}$", u))) {
    return(u[order(as.integer(u))])
  }
  sort(u)
}
