# CSurf initializers.
#
# CSurf needs an initial score f_hat0 that orders observations along the change
# surface. Each initializer regresses the signal y (DR pseudo-outcome or a
# first-stage CATE estimate) on X and converts the fitted values into a
# rank-smoothed score: the ordering of the fit is kept, but its magnitudes are
# replaced by a linear function of rank. The attribute "csurf_f0_source"
# records which initializer (or fallback) produced the score.

.resolve_csurf_f_hat0 <- function(X, y, f_hat0) {
  if (is.null(f_hat0)) {
    return(.csurf_mark_f0(y, "pseudo_outcome"))
  }
  if (is.character(f_hat0)) {
    init <- switch(f_hat0,
                   SuperLearner = .csurf_sl_f0,
                   xgb_smooth   = .csurf_xgb_f0,
                   gam_smooth   = .csurf_gam_f0,
                   change_plane = .csurf_change_plane_f0,
                   stop("Unknown CSurf initializer: ", f_hat0, call. = FALSE))
    return(init(X, y))
  }
  if (length(f_hat0) != length(y)) {
    stop("f_hat0 must have the same length as y.", call. = FALSE)
  }
  .csurf_mark_f0(f_hat0, "provided")
}

.csurf_mark_f0 <- function(f0, source) {
  f0 <- as.numeric(f0)
  attr(f0, "csurf_f0_source") <- source
  f0
}

.csurf_rank_smooth <- function(score, source) {
  score <- as.numeric(score)
  n <- length(score)
  if (n <= 2L || any(!is.finite(score)) ||
      stats::var(score) <= .Machine$double.eps) {
    return(NULL)
  }

  ord <- order(score)
  f0 <- numeric(n)
  f0[ord] <- stats::lm.fit(cbind(1, seq_len(n)), score[ord])$fitted.values / 2

  if (!all(is.finite(f0)) || stats::var(f0) <= .Machine$double.eps) {
    return(NULL)
  }
  .csurf_mark_f0(f0, source)
}

# SuperLearner ensemble of shallow XGBoost and an additive smoother (SL.gam,
# or SL.glm if the gam package is unavailable). Primary feasible initializer.
.csurf_sl_f0 <- function(X, y) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  ok <- complete.cases(X) & is.finite(y)
  if (sum(ok) <= 10L) {
    return(.csurf_linear_f0_fallback(X, y))
  }

  pred <- tryCatch({
    X_df <- as.data.frame(X[ok, , drop = FALSE])
    sl_lib <- if (requireNamespace("gam", quietly = TRUE)) {
      c("SL.xgboost", "SL.gam")
    } else {
      c("SL.xgboost", "SL.glm")
    }
    fit <- SuperLearner::SuperLearner(Y = y[ok], X = X_df, family = gaussian(),
                                      SL.library = sl_lib,
                                      cvControl = list(V = 5L))
    if (isTRUE(sum(fit$coef, na.rm = TRUE) == 0)) fit$coef[1] <- 1

    X_new <- as.data.frame(X)
    colnames(X_new) <- colnames(X_df)
    as.numeric(stats::predict(fit, newdata = X_new, onlySL = TRUE)$pred)
  }, error = function(e) NULL)

  f0 <- if (is.null(pred) || !all(is.finite(pred))) NULL else
    .csurf_rank_smooth(pred, "SuperLearner")
  if (is.null(f0)) .csurf_xgb_f0(X, y) else f0
}

# Shallow boosted trees (depth 2); the number of rounds is chosen by 5-fold CV.
.csurf_xgb_f0 <- function(X, y) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  ok <- complete.cases(X) & is.finite(y)
  if (sum(ok) <= 10L) {
    return(.csurf_linear_f0_fallback(X, y))
  }

  params <- list(objective = "reg:squarederror", max_depth = 2L, eta = 0.1,
                 subsample = 0.8, colsample_bytree = 0.8, nthread = 1L)
  valid_n <- function(x) !is.null(x) && length(x) > 0L && is.finite(x)

  pred <- tryCatch({
    dtrain <- xgboost::xgb.DMatrix(data = X[ok, , drop = FALSE], label = y[ok])
    cv_fit <- xgboost::xgb.cv(data = dtrain, params = params, nrounds = 200L,
                              nfold = 5L, early_stopping_rounds = 20L,
                              verbose = 0L)
    # The CV result fields differ across xgboost versions.
    best_n <- cv_fit$best_iteration
    if (!valid_n(best_n)) best_n <- cv_fit$niter
    if (!valid_n(best_n)) {
      eval_log <- cv_fit$evaluation_log
      if (!is.null(eval_log) && "test_rmse_mean" %in% names(eval_log)) {
        best_n <- eval_log$iter[which.min(eval_log$test_rmse_mean)]
      }
    }
    if (!valid_n(best_n)) best_n <- 200L
    best_n <- max(1L, as.integer(best_n[1L]))

    fit <- xgboost::xgb.train(data = dtrain, params = params,
                              nrounds = best_n, verbose = 0L)
    as.numeric(stats::predict(fit, xgboost::xgb.DMatrix(data = X)))
  }, error = function(e) NULL)

  f0 <- if (is.null(pred) || !all(is.finite(pred))) NULL else
    .csurf_rank_smooth(pred, "xgb_smooth")
  if (is.null(f0)) .csurf_linear_f0_fallback(X, y) else f0
}

# Additive model with one penalized smooth per covariate (mgcv, select = TRUE).
.csurf_gam_f0 <- function(X, y) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  p <- ncol(X)
  ok <- complete.cases(X) & is.finite(y)
  if (sum(ok) <= 10L) {
    return(.csurf_linear_f0_fallback(X, y))
  }

  pred <- tryCatch({
    df <- as.data.frame(X[ok, , drop = FALSE])
    colnames(df) <- paste0("V", seq_len(p))
    df$y <- y[ok]
    k <- min(5L, max(3L, floor(sum(ok) / (3L * p))))
    fml <- stats::as.formula(paste(
      "y ~", paste0("s(V", seq_len(p), ", k=", k, ")", collapse = " + ")))
    fit <- mgcv::gam(fml, data = df, select = TRUE, method = "REML")

    X_new <- as.data.frame(X)
    colnames(X_new) <- paste0("V", seq_len(p))
    as.numeric(stats::predict(fit, newdata = X_new))
  }, error = function(e) NULL)

  f0 <- if (is.null(pred) || !all(is.finite(pred))) NULL else
    .csurf_rank_smooth(pred, "gam_smooth")
  if (is.null(f0)) .csurf_xgb_f0(X, y) else f0
}

# Stabilized change-plane initializer. CSurf::f0_Fun is fitted to the 1, 2, 3,
# 5 or 8 covariates most correlated with the (winsorized, standardized) signal,
# using an evenly spaced subsample of at most `csurf_max_n` rows. If all fits
# fail, a native change-plane search, f0_Fun on leading principal components,
# and finally the SuperLearner initializer are tried in turn.
.csurf_change_plane_f0 <- function(X, y, max_p = 8L, csurf_max_n = 120L) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- length(y)

  fallback <- function() {
    f0 <- .csurf_sl_f0(X, y)
    attr(f0, "csurf_f0_source") <- paste0("change_plane_fallback_",
                                          attr(f0, "csurf_f0_source"))
    f0
  }

  if (n <= 10L || nrow(X) != n || ncol(X) < 1L) {
    return(fallback())
  }
  ok_y <- is.finite(y)
  if (sum(ok_y) <= 10L || stats::var(y[ok_y]) <= .Machine$double.eps) {
    return(fallback())
  }

  # Standardize non-constant columns, imputing non-finite values by the median.
  keep_col <- vapply(seq_len(ncol(X)), function(j) {
    z <- X[, j]
    finite_z <- is.finite(z)
    any(finite_z) && stats::var(z[finite_z]) > .Machine$double.eps
  }, logical(1L))
  if (!any(keep_col)) {
    return(fallback())
  }

  X_use <- X[, keep_col, drop = FALSE]
  for (j in seq_len(ncol(X_use))) {
    z <- X_use[, j]
    finite_z <- is.finite(z)
    fill <- stats::median(z[finite_z], na.rm = TRUE)
    if (!is.finite(fill)) fill <- 0
    z[!finite_z] <- fill
    z_sd <- stats::sd(z)
    if (!is.finite(z_sd) || z_sd <= .Machine$double.eps) z_sd <- 1
    X_use[, j] <- (z - mean(z)) / z_sd
  }

  # Winsorize the signal at its 1st/99th percentiles and standardize it.
  y_fill <- stats::median(y[ok_y], na.rm = TRUE)
  if (!is.finite(y_fill)) y_fill <- 0
  y_work <- y
  y_work[!ok_y] <- y_fill
  y_q <- stats::quantile(y_work[ok_y], probs = c(0.01, 0.99),
                         na.rm = TRUE, names = FALSE)
  if (all(is.finite(y_q)) && y_q[1L] < y_q[2L]) {
    y_work <- pmin(pmax(y_work, y_q[1L]), y_q[2L])
  }
  y_sd <- stats::sd(y_work[ok_y])
  if (!is.finite(y_sd) || y_sd <= .Machine$double.eps) {
    return(fallback())
  }
  y_work <- (y_work - mean(y_work[ok_y])) / y_sd

  cor_score <- vapply(seq_len(ncol(X_use)), function(j) {
    value <- suppressWarnings(stats::cor(X_use[ok_y, j], y_work[ok_y]))
    if (is.finite(value)) abs(value) else 0
  }, numeric(1L))
  ord <- order(cor_score, decreasing = TRUE)
  p_use <- min(as.integer(max_p), ncol(X_use))
  subset_sizes <- unique(pmin(c(1L, 2L, 3L, 5L, 8L, p_use), p_use))

  init_n <- min(n, as.integer(csurf_max_n))
  init_n <- max(40L, 10L * floor(init_n / 10L))
  init_n <- min(init_n, n)
  init_idx <- if (n > init_n) {
    unique(as.integer(round(seq(1L, n, length.out = init_n))))
  } else {
    seq_len(n)
  }
  if (length(init_idx) < min(n, init_n)) {
    init_idx <- sort(unique(c(init_idx, seq_len(n))))
    init_idx <- init_idx[seq_len(min(length(init_idx), init_n))]
  }

  # Project all n rows onto the fitted direction when f0_Fun returns one.
  extract_full_f0 <- function(fit, Z) {
    if (is.list(fit) && !is.null(fit$omeg.hat)) {
      omega <- as.numeric(fit$omeg.hat)
      if (length(omega) == ncol(Z) && all(is.finite(omega)) &&
          sum(omega^2) > .Machine$double.eps) {
        return(as.numeric(Z %*% omega))
      }
    }
    if (nrow(Z) == n && is.list(fit) && !is.null(fit$f0)) {
      return(as.numeric(fit$f0))
    }
    if (nrow(Z) == n) {
      return(as.numeric(fit))
    }
    NULL
  }

  try_f0_fun <- function(Z, label) {
    Z_init <- Z[init_idx, , drop = FALSE]
    # f0_Fun's k-means sieve needs more distinct rows than centers.
    if (nrow(unique(as.data.frame(Z_init))) <= ceiling(nrow(Z_init) / 10)) {
      return(NULL)
    }
    fit <- tryCatch(
      suppressWarnings({
        value <- NULL
        invisible(utils::capture.output(
          value <- CSurf::f0_Fun(as.matrix(Z_init), y_work[init_idx]),
          type = "output"
        ))
        value
      }),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)

    f0_raw <- tryCatch(extract_full_f0(fit, Z), error = function(e) NULL)
    if (is.null(f0_raw) || length(f0_raw) != n || !all(is.finite(f0_raw)) ||
        stats::var(f0_raw) <= .Machine$double.eps) {
      return(NULL)
    }
    .csurf_rank_smooth(f0_raw, paste0("change_plane_CSurf_", label))
  }

  for (k in subset_sizes) {
    f0 <- try_f0_fun(X_use[, ord[seq_len(k)], drop = FALSE], paste0("top", k))
    if (!is.null(f0)) return(f0)
  }

  idx <- ord[seq_len(p_use)]
  f0 <- .csurf_native_change_plane_f0(
    X_use[, idx, drop = FALSE], y_work,
    source = paste0("change_plane_native_top", p_use)
  )
  if (!is.null(f0)) return(f0)

  f0 <- tryCatch({
    pc <- stats::prcomp(X_use[, idx, drop = FALSE], center = FALSE, scale. = FALSE)
    npc <- min(3L, ncol(pc$x))
    try_f0_fun(pc$x[, seq_len(npc), drop = FALSE], paste0("pc", npc))
  }, error = function(e) NULL)
  if (!is.null(f0)) return(f0)

  fallback()
}

# Native change-plane search: candidate directions (signed coordinate axes,
# ridge, mean difference, sliced inverse regression, principal components) are
# scored by the BIC of a two-regime ridge fit over a grid of thresholds.
.csurf_native_change_plane_f0 <- function(X, y, source) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- length(y)
  p <- ncol(X)

  if (n <= 30L || p < 1L || nrow(X) != n || any(!is.finite(X)) ||
      any(!is.finite(y)) || stats::var(y) <= .Machine$double.eps) {
    return(NULL)
  }

  y_centered <- y - mean(y)
  candidates <- list()
  add_candidate <- function(name, omega) {
    omega <- as.numeric(omega)
    if (length(omega) != p || !all(is.finite(omega))) return(invisible(NULL))
    norm_omega <- sqrt(sum(omega^2))
    if (!is.finite(norm_omega) || norm_omega <= .Machine$double.eps) {
      return(invisible(NULL))
    }
    candidates[[name]] <<- omega / norm_omega
    invisible(NULL)
  }

  cor_score <- vapply(seq_len(p), function(j) {
    value <- suppressWarnings(stats::cor(X[, j], y))
    if (is.finite(value)) value else 0
  }, numeric(1L))
  for (j in seq_len(p)) {
    axis_j <- numeric(p)
    axis_j[j] <- ifelse(cor_score[j] >= 0, 1, -1)
    add_candidate(paste0("axis", j), axis_j)
  }

  ridge_beta <- tryCatch({
    lambda <- 1e-4 * mean(diag(crossprod(X)))
    if (!is.finite(lambda) || lambda <= 0) lambda <- 1e-4
    as.numeric(solve(crossprod(X) + diag(lambda, p), crossprod(X, y_centered)))
  }, error = function(e) NULL)
  if (!is.null(ridge_beta)) add_candidate("ridge", ridge_beta)

  mean_diff <- tryCatch({
    lo <- y <= stats::quantile(y, 0.35, na.rm = TRUE)
    hi <- y >= stats::quantile(y, 0.65, na.rm = TRUE)
    colMeans(X[hi, , drop = FALSE]) - colMeans(X[lo, , drop = FALSE])
  }, error = function(e) NULL)
  if (!is.null(mean_diff)) add_candidate("mean_diff", mean_diff)

  sir_dirs <- tryCatch({
    h <- min(5L, max(2L, floor(n / 40L)))
    slice <- cut(rank(y, ties.method = "average"),
                 breaks = unique(round(seq(0, n, length.out = h + 1L))),
                 include.lowest = TRUE, labels = FALSE)
    xbar <- colMeans(X)
    sirmat <- matrix(0, p, p)
    for (s in stats::na.omit(unique(slice))) {
      id <- which(slice == s)
      if (length(id) <= 2L) next
      diff <- matrix(colMeans(X[id, , drop = FALSE]) - xbar, ncol = 1L)
      sirmat <- sirmat + length(id) / n * (diff %*% t(diff))
    }
    sigma <- stats::cov(X)
    lambda <- 1e-3 * mean(diag(sigma))
    if (!is.finite(lambda) || lambda <= 0) lambda <- 1e-3
    eig <- eigen(solve(sigma + diag(lambda, p), sirmat))
    Re(eig$vectors[, seq_len(min(3L, ncol(eig$vectors))), drop = FALSE])
  }, error = function(e) NULL)
  if (!is.null(sir_dirs)) {
    for (j in seq_len(ncol(sir_dirs))) {
      add_candidate(paste0("sir", j), sir_dirs[, j])
    }
  }

  pc_dirs <- tryCatch({
    pc <- stats::prcomp(X, center = FALSE, scale. = FALSE)
    pc$rotation[, seq_len(min(3L, ncol(pc$rotation))), drop = FALSE]
  }, error = function(e) NULL)
  if (!is.null(pc_dirs)) {
    for (j in seq_len(ncol(pc_dirs))) {
      direction <- pc_dirs[, j]
      sign_j <- suppressWarnings(stats::cor(as.numeric(X %*% direction), y))
      if (is.finite(sign_j) && sign_j < 0) direction <- -direction
      add_candidate(paste0("pc", j), direction)
    }
  }

  if (length(candidates) == 0L) return(NULL)

  # BIC of the ridge fit y ~ X + X * 1(z >= gamma).
  two_regime_bic <- function(Xreg) {
    p_reg <- ncol(Xreg)
    lambda <- 1e-4 * mean(diag(crossprod(Xreg)))
    if (!is.finite(lambda) || lambda <= 0) lambda <- 1e-4
    beta <- tryCatch(
      as.numeric(solve(crossprod(Xreg) + diag(lambda, p_reg),
                       crossprod(Xreg, y_centered))),
      error = function(e) NULL
    )
    if (is.null(beta) || !all(is.finite(beta))) return(NULL)
    mse <- mean(as.numeric(y_centered - Xreg %*% beta)^2)
    if (!is.finite(mse) || mse <= 0) return(NULL)
    n * log(mse) + sum(abs(beta) > 1e-6) * log(n)
  }

  best <- NULL
  min_group <- max(10L, ceiling(0.10 * n))
  for (name in names(candidates)) {
    z <- as.numeric(X %*% candidates[[name]])
    if (!all(is.finite(z)) || stats::var(z) <= .Machine$double.eps) next

    grid <- unique(as.numeric(stats::quantile(
      z, probs = seq(0.15, 0.85, length.out = 15L), na.rm = TRUE, names = FALSE
    )))
    for (gamma in grid) {
      ind <- as.numeric(z >= gamma)
      if (sum(ind) < min_group || sum(1 - ind) < min_group) next
      bic <- two_regime_bic(cbind(X, X * ind))
      if (!is.null(bic) && (is.null(best) || bic < best$bic)) {
        best <- list(bic = bic, name = name)
      }
    }
  }

  if (is.null(best)) return(NULL)
  .csurf_rank_smooth(X %*% candidates[[best$name]],
                     paste0(source, "_", best$name))
}

# Near-unpenalized ridge regression on standardized X (last-resort fallback).
.csurf_linear_f0_fallback <- function(X, y) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  ok <- complete.cases(X) & is.finite(y)
  if (sum(ok) <= 2L) {
    return(.csurf_mark_f0(y, "pseudo_outcome"))
  }

  center <- colMeans(X[ok, , drop = FALSE])
  scale <- apply(X[ok, , drop = FALSE], 2L, stats::sd)
  scale[!is.finite(scale) | scale <= 0] <- 1
  Xs <- sweep(sweep(X, 2L, center, "-"), 2L, scale, "/")
  Xs_fit <- Xs[ok, , drop = FALSE]
  y_fit <- y[ok] - mean(y[ok])

  p <- ncol(Xs_fit)
  lambda <- 1e-6 * mean(diag(crossprod(Xs_fit)))
  if (!is.finite(lambda) || lambda <= 0) lambda <- 1e-6
  beta <- tryCatch(
    as.numeric(solve(crossprod(Xs_fit) + diag(lambda, p),
                     crossprod(Xs_fit, y_fit))),
    error = function(e) rep(0, p)
  )

  f0 <- as.numeric(Xs %*% beta)
  if (all(is.finite(f0))) {
    .csurf_mark_f0(f0, "ridge_linear")
  } else {
    .csurf_mark_f0(y, "pseudo_outcome")
  }
}
