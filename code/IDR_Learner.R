#' IDR-Learner: Interpretable Doubly Robust Learner
#'
#' Fits CSurf directly to cross-fitted doubly robust (AIPW) pseudo-outcomes.
#'
#' @param X Covariate matrix or data frame.
#' @param W Binary treatment indicator (0/1).
#' @param Y Outcome vector.
#' @param learner_ps SuperLearner library for the propensity score.
#' @param learner_y SuperLearner library for the outcome regressions.
#' @param n_folds Number of cross-fitting folds.
#' @param f_hat0 CSurf initializer: NULL (the pseudo-outcome itself), a numeric
#'   score (e.g. the oracle score), or one of "SuperLearner", "xgb_smooth",
#'   "gam_smooth", "change_plane" (see CSurf_init.R).
#'
#' @return A CSurf fit with the cross-fitted pseudo-outcomes (`dr_pseudo`),
#'   fold labels (`fold_ids`) and the initializer used (`f_hat0_source`).
IDR_Learner <- function(X, W, Y, learner_ps = c("SL.glm", "SL.xgboost"),
                        learner_y = c("SL.glm", "SL.xgboost"), n_folds = 5,
                        f_hat0 = NULL) {
  X <- as.data.frame(X)
  dr <- .dr_pseudo_outcome(X, W, Y, learner_ps, learner_y, n_folds)
  fit <- .fit_csurf(dr$pseudo, X, f_hat0)
  fit$dr_pseudo <- dr$pseudo
  fit$fold_ids <- dr$fold_ids
  fit
}

# SuperLearner without internal cross-validation; if every learner receives
# zero weight, the first learner is used.
.fit_sl <- function(Y, X, family, SL.library, ...) {
  fit <- SuperLearner::SuperLearner(Y = Y, X = X, family = family,
                                    SL.library = SL.library,
                                    cvControl = list(V = 0), ...)
  if (isTRUE(sum(fit$coef, na.rm = TRUE) == 0)) fit$coef[1] <- 1
  fit
}

.sl_predict <- function(fit, newX) {
  as.vector(predict(fit, newX, onlySL = TRUE)$pred)
}

# Out-of-fold propensity scores and arm-specific outcome predictions.
.cross_fit_nuisance <- function(X, W, Y, learner_ps, learner_y, n_folds) {
  n <- length(Y)
  fold_ids <- sample(rep(1:n_folds, length.out = n))
  ps <- mu0 <- mu1 <- rep(NA_real_, n)

  for (k in 1:n_folds) {
    te <- fold_ids == k
    X_tr <- X[!te, ]
    W_tr <- W[!te]
    Y_tr <- Y[!te]
    ps[te] <- .sl_predict(.fit_sl(W_tr, X_tr, binomial(), learner_ps), X[te, ])
    mu0[te] <- .sl_predict(.fit_sl(Y_tr[W_tr == 0], X_tr[W_tr == 0, ],
                                   gaussian(), learner_y), X[te, ])
    mu1[te] <- .sl_predict(.fit_sl(Y_tr[W_tr == 1], X_tr[W_tr == 1, ],
                                   gaussian(), learner_y), X[te, ])
  }

  list(fold_ids = fold_ids, ps = pmin(pmax(ps, 1e-6), 1 - 1e-6),
       mu0 = mu0, mu1 = mu1)
}

# Cross-fitted AIPW pseudo-outcome.
.dr_pseudo_outcome <- function(X, W, Y, learner_ps, learner_y, n_folds) {
  nu <- .cross_fit_nuisance(X, W, Y, learner_ps, learner_y, n_folds)
  apo_1 <- nu$mu1 + W * (Y - nu$mu1) / nu$ps
  apo_0 <- nu$mu0 + (1 - W) * (Y - nu$mu0) / (1 - nu$ps)
  list(pseudo = apo_1 - apo_0, fold_ids = nu$fold_ids)
}

.fit_csurf <- function(y, X, f_hat0) {
  f_hat0 <- .resolve_csurf_f_hat0(X, y, f_hat0)
  fit <- CSurf::CSurf(y = y, X = as.matrix(X), f_hat0 = f_hat0)
  fit$f_hat0_source <- attr(f_hat0, "csurf_f0_source", exact = TRUE)
  fit
}
