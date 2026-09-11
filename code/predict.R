# Prediction on new covariates. Each function returns the CATE prediction
# (`yhat`) and the subgroup label (`group_hat`); predict_ICF() is in ICF.R.

predict_CSurf <- function(fit, X_new, rule = 2) {
  X_new <- as.matrix(X_new)
  n <- nrow(X_new)
  p <- ncol(X_new)
  tau <- fit$tau.hat
  alp <- fit$alp.hat
  f_hat0_source <- if (is.null(fit$f_hat0_source)) NA_character_ else
    fit$f_hat0_source

  # No change surface was selected: a single linear effect model.
  if (is.null(fit$Xord)) {
    return(list(yhat = as.vector(X_new %*% alp[seq_len(p)]),
                group_hat = rep(1L, n), f_hat = rep(0, n),
                f_hat0_source = f_hat0_source))
  }
  stopifnot(ncol(fit$Xord) == p)

  # Interpolate the fitted additive score components at the new points.
  f_all <- matrix(NA_real_, nrow = n, ncol = p)
  for (j in seq_len(p)) {
    f_all[, j] <- approx(x = fit$Xord[, j], y = fit$f_all.hat[, j],
                         xout = X_new[, j], rule = rule)$y
  }
  f_hat <- rowSums(f_all)
  keep <- complete.cases(f_all)
  X_keep <- X_new[keep, , drop = FALSE]

  # Baseline effect plus one increment per threshold crossed.
  yhat <- rep(NA_real_, n)
  yhat[keep] <- as.vector(X_keep %*% alp[seq_len(p)])
  for (k in seq_along(tau)) {
    above <- as.numeric(f_hat[keep] > tau[k])
    yhat[keep] <- yhat[keep] + as.vector(X_keep %*% alp[k * p + seq_len(p)]) * above
  }

  group_hat <- rep(NA_integer_, n)
  group_hat[keep] <- if (length(tau) >= 1) {
    1L + rowSums(outer(f_hat[keep], tau, `>`))
  } else {
    1L
  }

  list(yhat = yhat, group_hat = group_hat, f_hat = f_hat,
       f_hat0_source = f_hat0_source)
}

# Subgroups are the leaves of the fitted tree, numbered 1..K.
predict_CART <- function(fit, X_new) {
  X_df <- as.data.frame(X_new)
  node <- as.integer(predict(partykit::as.party(fit), newdata = X_df,
                             type = "node"))
  list(yhat = as.vector(predict(fit, newdata = X_df)),
       group_hat = as.integer(factor(node, levels = sort(unique(fit$where)))))
}

predict_CT <- function(fit, X_new) {
  X_df <- as.data.frame(X_new)
  node <- as.integer(predict(partykit::as.party(fit), newdata = X_df,
                             type = "node"))
  list(yhat = as.vector(predict(fit, newdata = X_df)),
       group_hat = as.integer(factor(node, levels = sort(unique(node)))))
}

# Each unit is assigned to the first CRE rule it satisfies (0 if none).
predict_CRE <- function(fit, X_new) {
  X_df <- as.data.frame(X_new)
  yhat <- predict(fit, X = X_df)

  env <- list2env(list(X = X_df), parent = parent.frame())
  group_hat <- rep(0L, nrow(X_df))
  for (i in seq_along(fit$rules)) {
    hit <- eval(parse(text = fit$rules[i]), envir = env) & group_hat == 0L
    group_hat[hit] <- i
  }

  list(yhat = yhat, group_hat = group_hat)
}
