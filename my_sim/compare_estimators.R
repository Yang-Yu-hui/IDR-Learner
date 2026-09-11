# The 15 estimators of the simulation study. Each takes training and test data
# from dgp_Fun() and returns test-set CATE predictions (`yhat`) and subgroup
# labels (`group_hat`). CSurf-based methods use the initializer in
# data$f_hat0_init. The list order is the order in which they are fitted.

.make_estimator <- function(fit_fun, predict_fun, csurf = FALSE) {
  force(fit_fun)
  force(predict_fun)
  function(data, data.test) {
    fit <- if (csurf) {
      fit_fun(data$X, data$W, data$Y, f_hat0 = data$f_hat0_init)
    } else {
      fit_fun(data$X, data$W, data$Y)
    }
    predict_fun(fit, data.test$X)
  }
}

estimators <- list(
  estimator_IDR_Learner = .make_estimator(IDR_Learner, predict_CSurf, csurf = TRUE),
  estimator_DRL_CRE     = .make_estimator(DRL_CRE, predict_CRE),
  estimator_DRL_CSurf   = .make_estimator(DRL_CSurf, predict_CSurf, csurf = TRUE),
  estimator_DRL_CART    = .make_estimator(DRL_CART, predict_CART),
  estimator_CF_CRE      = .make_estimator(CF_CRE, predict_CRE),
  estimator_CF_CSurf    = .make_estimator(CF_CSurf, predict_CSurf, csurf = TRUE),
  estimator_CF_CART     = .make_estimator(CF_CART, predict_CART),
  estimator_XL_CRE      = .make_estimator(XL_CRE, predict_CRE),
  estimator_XL_CSurf    = .make_estimator(XL_CSurf, predict_CSurf, csurf = TRUE),
  estimator_XL_CART     = .make_estimator(XL_CART, predict_CART),
  estimator_RL_CRE      = .make_estimator(RL_CRE, predict_CRE),
  estimator_RL_CSurf    = .make_estimator(RL_CSurf, predict_CSurf, csurf = TRUE),
  estimator_RL_CART     = .make_estimator(RL_CART, predict_CART),
  estimator_CT          = .make_estimator(CT, predict_CT),
  estimator_ICF         = function(data, data.test) {
    fit <- ICF(data$X, data$W, data$Y, treeNo = 200, iterationNo = 50,
               min.split.var = 4, split_val_round_posi = 1,
               icf_path = get0("icf_dir", ifnotfound = "iCF"))
    predict_ICF(fit, data.test$X)
  }
)
