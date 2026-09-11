# Comparator estimators.
#
# Post hoc pipelines estimate CATEs with a DR-, R- or X-learner or a causal
# forest and then summarize them with CSurf, CART or CRE. The causal tree (CT)
# is a standalone reference; iCF is in ICF.R.

# ---- First-stage CATE estimators ---------------------------------------------

.drl_cate <- function(X, W, Y, learner_ps, learner_y, n_folds) {
  pseudo <- .dr_pseudo_outcome(X, W, Y, learner_ps, learner_y, n_folds)$pseudo
  .sl_predict(.fit_sl(pseudo, X, gaussian(), learner_y), X)
}

.rl_cate <- function(X, W, Y, learner_ps, learner_y, n_folds) {
  nu <- .cross_fit_nuisance(X, W, Y, learner_ps, learner_y, n_folds)
  m_hat <- nu$ps * nu$mu1 + (1 - nu$ps) * nu$mu0
  w_resid <- W - nu$ps
  pseudo <- (Y - m_hat) / w_resid
  fit <- .fit_sl(pseudo, X, gaussian(), learner_y, obsWeights = w_resid^2)
  .sl_predict(fit, X)
}

.xl_cate <- function(X, W, Y, learner_ps, learner_y) {
  ps <- .sl_predict(.fit_sl(W, X, binomial(), learner_ps), X)
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  mu0 <- .fit_sl(Y[W == 0], X[W == 0, ], gaussian(), learner_y)
  mu1 <- .fit_sl(Y[W == 1], X[W == 1, ], gaussian(), learner_y)

  d_1 <- Y[W == 1] - .sl_predict(mu0, X[W == 1, ])
  d_0 <- .sl_predict(mu1, X[W == 0, ]) - Y[W == 0]
  tau_1 <- .fit_sl(d_1, X[W == 1, ], gaussian(), learner_y)
  tau_0 <- .fit_sl(d_0, X[W == 0, ], gaussian(), learner_y)

  (1 - ps) * .sl_predict(tau_1, X) + ps * .sl_predict(tau_0, X)
}

.cf_cate <- function(X, W, Y, num.trees, min.node.size, honesty, seed) {
  if (!is.null(seed)) set.seed(seed)
  cf <- grf::causal_forest(X = as.matrix(X), Y = Y, W = W,
                           num.trees = num.trees,
                           min.node.size = min.node.size,
                           honesty = honesty)
  as.vector(predict(cf)$predictions)
}

# Regression tree on the first-stage CATE, pruned at the CV-optimal cp.
.fit_cart <- function(pseudo_y, X, minsize, maxdepth, cp, prune) {
  tree_dat <- data.frame(pseudo_y = pseudo_y, X)
  tree <- rpart::rpart(
    pseudo_y ~ .,
    data = tree_dat,
    method = "anova",
    control = rpart::rpart.control(minbucket = minsize, maxdepth = maxdepth,
                                   cp = cp)
  )
  if (prune && !is.null(tree$cptable) && nrow(tree$cptable) > 0) {
    best_cp <- tree$cptable[which.min(tree$cptable[, "xerror"]), "CP"]
    tree <- rpart::prune(tree, cp = best_cp)
  }
  tree
}

# ---- DR-learner ------------------------------------------------------------------

DRL_CSurf <- function(X, W, Y, learner_ps = c("SL.glm", "SL.xgboost"),
                      learner_y = c("SL.glm", "SL.xgboost"), n_folds = 5,
                      f_hat0 = NULL) {
  X <- as.data.frame(X)
  ite <- .drl_cate(X, W, Y, learner_ps, learner_y, n_folds)
  .fit_csurf(ite, X, f_hat0)
}

DRL_CART <- function(X, W, Y, learner_ps = c("SL.glm", "SL.xgboost"),
                     learner_y = c("SL.glm", "SL.xgboost"), n_folds = 5,
                     minsize = 40, maxdepth = 5, cp = 0.001, prune = TRUE) {
  X <- as.data.frame(X)
  ite <- .drl_cate(X, W, Y, learner_ps, learner_y, n_folds)
  .fit_cart(ite, X, minsize, maxdepth, cp, prune)
}

DRL_CRE <- function(X, W, Y, learner_ps = c("SL.glm", "SL.xgboost"),
                    learner_y = c("SL.glm", "SL.xgboost"), n_folds = 5) {
  X <- as.data.frame(X)
  ite <- .drl_cate(X, W, Y, learner_ps, learner_y, n_folds)
  CRE::cre(Y, W, X, ite = ite)
}

# ---- R-learner -------------------------------------------------------------------

RL_CSurf <- function(X, W, Y, learner_ps = c("SL.glm", "SL.xgboost"),
                     learner_y = c("SL.glm", "SL.xgboost"), n_folds = 5,
                     f_hat0 = NULL) {
  X <- as.data.frame(X)
  ite <- .rl_cate(X, W, Y, learner_ps, learner_y, n_folds)
  .fit_csurf(ite, X, f_hat0)
}

RL_CART <- function(X, W, Y, learner_ps = c("SL.glm", "SL.xgboost"),
                    learner_y = c("SL.glm", "SL.xgboost"), n_folds = 5,
                    minsize = 40, maxdepth = 5, cp = 0.001, prune = TRUE) {
  X <- as.data.frame(X)
  ite <- .rl_cate(X, W, Y, learner_ps, learner_y, n_folds)
  .fit_cart(ite, X, minsize, maxdepth, cp, prune)
}

RL_CRE <- function(X, W, Y, learner_ps = c("SL.glm", "SL.xgboost"),
                   learner_y = c("SL.glm", "SL.xgboost"), n_folds = 5) {
  X <- as.data.frame(X)
  ite <- .rl_cate(X, W, Y, learner_ps, learner_y, n_folds)
  CRE::cre(Y, W, X, ite = ite)
}

# ---- X-learner (full-sample nuisance fits, no cross-fitting) ---------------------

XL_CSurf <- function(X, W, Y, learner_ps = c("SL.glm", "SL.xgboost"),
                     learner_y = c("SL.glm", "SL.xgboost"), f_hat0 = NULL) {
  X <- as.data.frame(X)
  ite <- .xl_cate(X, W, Y, learner_ps, learner_y)
  .fit_csurf(ite, X, f_hat0)
}

XL_CART <- function(X, W, Y, learner_ps = c("SL.glm", "SL.xgboost"),
                    learner_y = c("SL.glm", "SL.xgboost"),
                    minsize = 40, maxdepth = 5, cp = 0.001, prune = TRUE) {
  X <- as.data.frame(X)
  ite <- .xl_cate(X, W, Y, learner_ps, learner_y)
  .fit_cart(ite, X, minsize, maxdepth, cp, prune)
}

XL_CRE <- function(X, W, Y, learner_ps = c("SL.glm", "SL.xgboost"),
                   learner_y = c("SL.glm", "SL.xgboost")) {
  X <- as.data.frame(X)
  ite <- .xl_cate(X, W, Y, learner_ps, learner_y)
  CRE::cre(Y, W, X, ite = ite)
}

# ---- Causal forest (out-of-bag CATE predictions) -----------------------------------

CF_CSurf <- function(X, W, Y, num.trees = 2000, min.node.size = 5,
                     honesty = TRUE, seed = 1234, f_hat0 = NULL) {
  ite <- .cf_cate(X, W, Y, num.trees, min.node.size, honesty, seed)
  .fit_csurf(ite, X, f_hat0)
}

CF_CART <- function(X, W, Y, num.trees = 2000, min.node.size = 5,
                    honesty = TRUE, minsize = 40, maxdepth = 5, cp = 0.001,
                    prune = TRUE, seed = NULL) {
  ite <- .cf_cate(X, W, Y, num.trees, min.node.size, honesty, seed)
  .fit_cart(ite, X, minsize, maxdepth, cp, prune)
}

CF_CRE <- function(X, W, Y, num.trees = 2000, min.node.size = 5,
                   honesty = TRUE, seed = NULL) {
  ite <- .cf_cate(X, W, Y, num.trees, min.node.size, honesty, seed)
  CRE::cre(Y, W, X, ite = ite)
}

# ---- Causal tree -----------------------------------------------------------------

CT <- function(X, W, Y, split.Rule = "CT", cv.option = "CT",
               split.Honest = TRUE, cv.Honest = TRUE, minsize = 40,
               cp = 0.001, prune = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dat <- data.frame(Y = Y, as.data.frame(X))

  ct_fit <- causalTree::causalTree(
    formula = Y ~ .,
    data = dat,
    treatment = W,
    split.Rule = split.Rule,
    cv.option = cv.option,
    split.Honest = split.Honest,
    cv.Honest = cv.Honest,
    minsize = minsize,
    cp = cp
  )
  if (prune && !is.null(ct_fit$cptable) && nrow(ct_fit$cptable) > 0) {
    best_cp <- ct_fit$cptable[which.min(ct_fit$cptable[, "xerror"]), "CP"]
    ct_fit <- prune(ct_fit, cp = best_cp)
  }
  ct_fit
}
