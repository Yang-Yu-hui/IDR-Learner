# Data-generating processes (DGPs 1-5) of the simulation study.
#
# Covariates come from a Gaussian copula with Uniform(0, 1) margins: three
# confounders X_C, twelve nuisance predictors X_N and p effect modifiers Z
# (p = 20 in DGP 2, otherwise 5). The outcome is
#   Y = -0.6 + X_C'(0.6, 0.4, 0.2) + X_N'eta_N + W * tau(Z) + N(0, 1),
# with tau(Z) = Z'beta_g in subgroup g. Subgroups are defined by thresholds on
# a standardized additive score f(Z) (DGPs 1-4) or by a 2 x 2 rectangular
# partition of (Z_1, Z_2) (DGP 5). Treatment follows a logistic model in X_C;
# its intercept is calibrated to the prevalence `pi` and its slope is set by
# `overlap`.
#
# Rows are sorted by f. Besides the data, the true CATE (`cate`), subgroups
# (`gr.true`) and the oracle CSurf initializer (`f_hat0_init`) are returned.

dgp_Fun <- function(n = 500, example_id = 4, rho = 0,
                    pi = 0.5, overlap = "moderate") {
  if (!(example_id %in% 1:5))
    stop("example_id must be one of 1, 2, 3, 4, 5")
  if (!(overlap %in% c("weak", "moderate", "strong")))
    stop("overlap must be 'weak', 'moderate', or 'strong'")
  if (pi <= 0 || pi >= 1)
    stop("pi must be strictly between 0 and 1")

  p <- if (example_id == 2L) 20L else 5L
  pcov <- 15L
  s_ps <- switch(overlap, strong = 1.0, moderate = 1.5, weak = 2.0)

  # Columns: [X_C (3) | X_N (12) | Z (p)]; rho is the common correlation.
  total_cols <- pcov + p
  Sigma_all <- matrix(rho, total_cols, total_cols) + diag(total_cols) * (1 - rho)
  pvars_all <- pnorm(MASS::mvrnorm(n = n, mu = rep(0, total_cols), Sigma = Sigma_all))
  XC <- pvars_all[, 1:3, drop = FALSE]
  nuisance_term <- as.vector(pvars_all[, 4:pcov, drop = FALSE] %*%
                               c(rep(0.22, 4L), rep(0, pcov - 7L)))

  # Effect modifiers: truncated N(0, 1) on [-1, 1] (DGP 1) or Uniform(-1, 1).
  U <- pvars_all[, (pcov + 1L):(pcov + p), drop = FALSE]
  if (example_id == 1) {
    Fa <- pnorm(-1)
    Fb <- pnorm(1)
    X <- qnorm(Fa + U * (Fb - Fa))
  } else {
    X <- 2 * U - 1
  }

  # f = sum_j f_j(Z_j) / sd, with each component centered.
  additive_score <- function(comps) {
    comps <- lapply(comps, function(v) v - mean(v))
    sd_f <- sqrt(Reduce(`+`, lapply(comps, var)))
    Reduce(`+`, lapply(comps, function(v) v / sd_f))
  }

  if (example_id == 1) {
    # Smooth nonlinear score, 3 subgroups.
    f <- additive_score(list(X[, 1], 0.5 * (X[, 2] - 1)^2, 3 * pnorm(X[, 3])))
    tau.true <- c(-0.5, 0.2)
    beta_list <- list(c(2, 2, 2, 2, 0), c(2, -2, -1, 0, 0), c(2, 1, 2, 0, 2))

  } else if (example_id == 2) {
    # High-dimensional nonlinear score (8 of 20 modifiers), 3 subgroups.
    f_funs <- list(function(x) 2 * sin(1.5 * x),
                   function(x) 2.5 * x^2,
                   function(x) -2 * x,
                   function(x) exp(-x) - 0.5 * sinh(2))
    f <- additive_score(lapply(1:8, function(j) f_funs[[(j - 1) %% 4 + 1]](X[, j])))
    tau.true <- c(-0.5, 0.2)
    beta1 <- c(rep(c(1, -1), length.out = 9), rep(0, p - 10L), 1) * 1.5
    beta_list <- list(beta1, -1 * beta1, 2 * beta1)

  } else if (example_id == 3) {
    # Linear score, 2 subgroups.
    f <- additive_score(list(X[, 1], -2 * X[, 2], 1.5 * X[, 3]))
    tau.true <- 0
    beta_list <- list(c(2, 2, 2, 2, 0), c(2, -1, -2, 0, 2))

  } else if (example_id == 4) {
    # Step-function score, 2 subgroups; f takes three values (about -1.41, 0,
    # 1.41), and the threshold separates the highest value from the rest.
    f <- additive_score(list(ifelse(X[, 1] > 0, 1, 0), ifelse(X[, 2] > -0.5, 1, 0)))
    tau.true <- 0.7
    beta_list <- list(c(2, 2, 2, 2, 0), c(2, -1, -2, 0, 3))

  } else {
    # Rectangular partition G = 1 + 1(Z_1 > 0) + 2 * 1(Z_2 > 0), 4 subgroups;
    # not representable by an additive change score.
    f_raw <- as.integer(X[, 1] > 0) + 2 * as.integer(X[, 2] > 0)
    f <- (f_raw - mean(f_raw)) / sd(f_raw)
    f_vals <- sort(unique(f_raw))
    tau_raw <- (f_vals[-length(f_vals)] + f_vals[-1]) / 2
    tau.true <- (tau_raw - mean(f_raw)) / sd(f_raw)
    beta_list <- list(c(2, 2, 1, 0, 0), c(2, -1, -1, 0, 0),
                      c(-1, 2, 0, 1, 0), c(-1, -1, 0, 0, 2))
  }

  gr.true <- if (example_id == 5) {
    1L + as.integer(X[, 1] > 0) + 2L * as.integer(X[, 2] > 0)
  } else {
    1L + as.integer(rowSums(outer(f, tau.true, ">")))
  }
  cate <- rowSums(X * do.call(rbind, beta_list)[gr.true, , drop = FALSE])

  # Treatment: logistic in X_C, intercept calibrated so that mean(P(W = 1)) = pi.
  lp_raw <- s_ps * (XC[, 1] - XC[, 2] + XC[, 3])
  gamma_0 <- uniroot(function(g) mean(plogis(g + lp_raw)) - pi,
                     interval = c(-20, 20))$root
  W <- rbinom(n, 1L, prob = plogis(gamma_0 + lp_raw))

  Y <- -0.6 + as.vector(XC %*% c(0.6, 0.4, 0.2)) +
    cate * W + nuisance_term + rnorm(n, 0, 1)

  ord <- order(f)
  pvars <- pvars_all[ord, 1:pcov, drop = FALSE]
  colnames(pvars) <- paste0("cov", 1:pcov)

  list(
    X           = cbind(pvars, X[ord, , drop = FALSE]),
    W           = W[ord],
    Y           = Y[ord],
    f           = f[ord],
    tau.true    = tau.true,
    gr.true     = gr.true[ord],
    cate        = cate[ord],
    # Oracle initializer: rows are sorted by f, so seq_len(n) is the rank of f.
    f_hat0_init = lm(f[ord] ~ seq_len(n))$fitted.values / 2
  )
}
