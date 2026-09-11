# One Monte Carlo replicate: simulate a test set (seed i) and a training set
# (seed i + 2025), fit every estimator in `estimators`, and record the test-set
# MSE, ARI, NMI and runtime. A failing estimator yields NA instead of stopping
# the replicate.

myfun <- function(i, n, example_id, pi, overlap, f_hat0_type, save_path = NULL) {
  set.seed(i)
  data.test <- dgp_Fun(n = 2000, example_id = example_id,
                       rho = 0, pi = pi, overlap = overlap)
  set.seed(i + 2025)
  data <- dgp_Fun(n = n, example_id = example_id,
                  rho = 0, pi = pi, overlap = overlap)

  # "oracle" keeps the true-score initializer returned by dgp_Fun().
  if (f_hat0_type != "oracle") {
    if (!f_hat0_type %in% c("SuperLearner", "xgb_smooth", "gam_smooth",
                            "change_plane")) {
      stop("Unknown f_hat0_type: ", f_hat0_type, call. = FALSE)
    }
    data$f_hat0_init <- f_hat0_type
  }

  true.cate <- data.test$cate
  gr.true <- data.test$gr.true
  n_test <- length(true.cate)

  rows <- lapply(names(estimators), function(name) {
    res <- NULL
    elapsed <- system.time(
      # Console output of the estimators is discarded.
      invisible(capture.output(
        res <- tryCatch(
          suppressMessages(estimators[[name]](data, data.test)),
          error = function(e) {
            warning("[rep ", i, "] ", name, " failed: ", conditionMessage(e))
            list(yhat = rep(NA_real_, n_test),
                 group_hat = rep(NA_integer_, n_test))
          }
        ),
        type = "output"
      ))
    )["elapsed"]

    g <- res$group_hat
    valid <- which(!is.na(g) & !is.na(gr.true))
    data.frame(
      rep            = i,
      estimator.name = name,
      true.cate.var  = var(true.cate, na.rm = TRUE),
      mse            = mean((res$yhat - true.cate)^2, na.rm = TRUE),
      ari            = if (length(valid) > 1) aricode::ARI(g[valid], gr.true[valid]) else NA_real_,
      nmi            = if (length(valid) > 1) aricode::NMI(g[valid], gr.true[valid]) else NA_real_,
      time_sec       = as.numeric(elapsed),
      n              = n,
      example_id     = example_id,
      pi             = pi,
      overlap        = overlap,
      f_hat0_type    = f_hat0_type,
      f_hat0_source  = if (is.null(res$f_hat0_source)) NA_character_ else res$f_hat0_source,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)

  if (!is.null(save_path)) {
    dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
    saveRDS(df, file.path(save_path, sprintf(
      "result_ex%d_n%d_pi%.2f_%s_%s_rep%04d.rds",
      example_id, n, pi, overlap, f_hat0_type, i)))
  }
  df
}
