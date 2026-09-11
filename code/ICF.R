# Iterative causal forest (iCF) comparator.
#
# The iCF scripts in iCF/ are sourced into a private environment, and a few of
# their functions are patched so that they run unattended inside
# cross-validation folds.

.load_icf_env <- local({
  cache_env <- NULL
  cache_path <- NULL

  function(path = "iCF") {
    if (!is.null(cache_env) && identical(cache_path, path)) {
      return(cache_env)
    }

    scripts <- c("best_tree_MSegar.R", "iCF_CV.R", "iCF_TREE_build.R",
                 "iCF_PARENT_node.R", "iCF_PRE_majority.R",
                 "iCF_MAJORITY_VOTE.R", "iCF_SUBGROUP_DECISION.R",
                 "iCF_SUBGROUP_PIPELINE.R", "iCF_SUBGROUP_ANALYSIS.R",
                 "iCF_SUBGROUP_MODEL.R", "iCF_GG_toolbox.R", "sim_Truth_tree.R")
    r_files <- file.path(path, scripts)
    missing <- r_files[!file.exists(r_files)]
    if (length(missing) > 0) {
      stop("Missing iCF scripts: ", paste(basename(missing), collapse = ", "))
    }

    icf_env <- new.env(parent = globalenv())
    invisible(lapply(r_files, source, local = icf_env))

    # The scripts call rlang's global_env() and caller_env() unqualified.
    icf_env$global_env <- function() rlang::global_env()
    icf_env$caller_env <- function(n = 1L) rlang::caller_env(n = n)

    icf_env$prop.func <- function(x, trt) {
      propens.model <- glmnet::cv.glmnet(y = trt, x = x, family = "binomial")
      predict(propens.model, s = "lambda.min", newx = x, type = "response")[, 1]
    }

    # CF() reads Y, W, X, selected_cf.idx and Train from the iCF environment;
    # refresh them with the values of the current fold.
    orig_subgroup_pipeline <- get("SUBGROUP_PIPELINE", envir = icf_env, inherits = FALSE)
    icf_env$SUBGROUP_PIPELINE <- function(...) {
      dots <- list(...)
      if (!is.null(dots[["X"]]) && length(dots[["selected_cf.idx"]]) == 0L) {
        dots[["selected_cf.idx"]] <- seq_len(ncol(as.data.frame(dots[["X"]])))
      }
      for (key in c("Y", "W", "Y.hat", "W.hat", "selected_cf.idx")) {
        if (!is.null(dots[[key]])) assign(key, dots[[key]], envir = icf_env)
      }
      if (!is.null(dots[["X"]])) {
        assign("X", dots[["X"]], envir = icf_env)
        # Factor levels can disappear within a fold.
        assign("vars_catover2", icf_env$find_level_over2(as.data.frame(dots[["X"]])),
               envir = icf_env)
      }
      if (!is.null(dots[["X"]]) && !is.null(dots[["Y"]]) && !is.null(dots[["W"]])) {
        assign("Train", cbind.data.frame(
          Y = as.numeric(dots[["Y"]]),
          W = as.numeric(dots[["W"]]),
          as.data.frame(dots[["X"]])
        ), envir = icf_env)
      }
      do.call(orig_subgroup_pipeline, dots)
    }

    # Drop contrasts for factors that have a single level within a fold.
    if (exists("SGMODEL_DATA", envir = icf_env, inherits = FALSE)) {
      orig_sgmodel_data <- get("SGMODEL_DATA", envir = icf_env, inherits = FALSE)
      icf_env$SGMODEL_DATA <- function(dat, outcome_type_at, W.hat) {
        out <- orig_sgmodel_data(dat = dat, outcome_type_at = outcome_type_at, W.hat = W.hat)
        if (is.list(out) && !is.null(out$dat_ID_SG_df) && length(names(out$contr)) > 0) {
          keep <- vapply(names(out$contr), function(vn) {
            if (!vn %in% names(out$dat_ID_SG_df)) return(FALSE)
            x <- out$dat_ID_SG_df[[vn]]
            !is.factor(x) || nlevels(x) >= 2
          }, logical(1))
          out$contr <- out$contr[keep]
        }
        out
      }
    }

    # Return a model for every subgroup term, even when a G_k term is absent.
    if (exists("CF_GROUP_DECISION", envir = icf_env, inherits = FALSE)) {
      patched_cf_group_decision <- function(HTE_P_cf.raw, dat, method_CF,
                                            outcome_type_at, P_threshold, W.hat) {
        if (round(HTE_P_cf.raw, 1) > P_threshold) {
          return(list(model.basic = "NA", model.g2 = "NA", model.g3 = "NA",
                      model.g4 = "NA", model.g5 = "NA"))
        }

        SGmodel_dat <- SGMODEL_DATA(dat, outcome_type_at, W.hat)
        dat_ID_SG_df <- SGmodel_dat$dat_ID_SG_df
        contr <- SGmodel_dat$contr
        if (is.null(contr)) contr <- list()
        model.g2 <- model.g3 <- model.g4 <- model.g5 <- NULL

        y_n_unique <- length(unique(dat_ID_SG_df$Y))
        if (y_n_unique >= 8) {
          model.basic <- stats::lm(formula_basic, data = dat_ID_SG_df)
          if ("G2" %in% names(contr)) model.g2 <- stats::lm(formula_g2, data = dat_ID_SG_df)
          if ("G3" %in% names(contr)) model.g3 <- stats::lm(formula_g3, data = dat_ID_SG_df)
          if ("G4" %in% names(contr)) model.g4 <- stats::lm(formula_g4, data = dat_ID_SG_df)
          if ("G5" %in% names(contr)) model.g5 <- stats::lm(formula_g5, data = dat_ID_SG_df)
        } else if (y_n_unique == 2) {
          model.basic <- stats::glm(formula_basic, data = dat_ID_SG_df, family = stats::binomial())
          if ("G2" %in% names(contr)) model.g2 <- stats::glm(formula_g2, data = dat_ID_SG_df, family = stats::binomial())
          if ("G3" %in% names(contr)) model.g3 <- stats::glm(formula_g3, data = dat_ID_SG_df, family = stats::binomial())
          if ("G4" %in% names(contr)) model.g4 <- stats::glm(formula_g4, data = dat_ID_SG_df, family = stats::binomial())
          if ("G5" %in% names(contr)) model.g5 <- stats::glm(formula_g5, data = dat_ID_SG_df, family = stats::binomial())
        } else {
          model.basic <- stats::lm(formula_basic, data = dat_ID_SG_df)
        }

        if (is.null(model.g2)) model.g2 <- model.basic
        if (is.null(model.g3)) model.g3 <- model.basic
        if (is.null(model.g4)) model.g4 <- model.basic
        if (is.null(model.g5)) model.g5 <- model.basic
        list(model.basic = model.basic, model.g2 = model.g2, model.g3 = model.g3,
             model.g4 = model.g4, model.g5 = model.g5)
      }
      environment(patched_cf_group_decision) <- icf_env
      icf_env$CF_GROUP_DECISION <- patched_cf_group_decision
    }

    cache_path <<- path
    cache_env <<- icf_env
    icf_env
  }
})

# Minimum leaf size tuning for tree depths 2-5.
.run_mls_tuning <- function(icf_env, train_df, treeNo, iterationNo,
                            split_val_round_posi) {
  depths <- c("D2", "D3", "D4", "D5")
  denominators <- c(25, 45, 65, 85)
  res <- setNames(lapply(seq_along(depths), function(i) {
    icf_env$MinLeafSizeTune(
      dat                  = train_df,
      denominator          = denominators[i],
      treeNo               = treeNo,
      iterationNo          = iterationNo,
      split_val_round_posi = split_val_round_posi,
      depth                = depths[i]
    )
  }), depths)
  res$leafsize <- setNames(lapply(depths, function(d) res[[d]]$denominator), depths)
  res
}

# Larger leaf-size denominators to retry with; the iCF documentation suggests
# this for "parent_sign" errors.
.leafsize_retry_grid <- function(leafsize) {
  base <- lapply(leafsize, as.numeric)
  lapply(c(0, 10, 20, 30), function(add) {
    list(D5 = base$D5 + add, D4 = base$D4 + add,
         D3 = base$D3 + add, D2 = base$D2 + add)
  })
}

ICF <- function(X, W, Y,
                K = 5,
                treeNo = 200,
                iterationNo = 50,
                min.split.var = 4,
                split_val_round_posi = 2,
                P_threshold = 1,
                variable_type = "non-hd",
                hdPctTop = 0.95,
                tune_leafsize = TRUE,
                icf_path = "iCF",
                seed = 1234) {
  suppressMessages(suppressWarnings({
    set.seed(seed)

    X_df <- as.data.frame(X)
    W <- as.numeric(W)
    Y <- as.numeric(Y)
    if (length(W) != nrow(X_df) || length(Y) != nrow(X_df)) {
      stop("X, W and Y must have the same number of observations.")
    }
    if (!all(W %in% c(0, 1))) stop("W must be binary (0/1).")
    if (anyNA(X_df) || anyNA(W) || anyNA(Y)) stop("iCF does not accept missing values.")
    variable_type <- match.arg(tolower(variable_type), c("non-hd", "hd", "hdshrink"))

    train_df <- cbind.data.frame(data.frame(Y = Y, W = W, check.names = FALSE), X_df)
    icf_env <- .load_icf_env(icf_path)

    leafsize <- list(D5 = 100, D4 = 80, D3 = 60, D2 = 40)
    icf_env$vars_forest <- colnames(X_df)
    icf_env$X <- X_df
    icf_env$Y <- Y
    icf_env$W <- W
    icf_env$leafsize <- leafsize
    icf_env$vars_catover2 <- icf_env$find_level_over2(X_df)

    # Step 1: raw causal forest (variable screening and HTE test).
    raw_cf <- icf_env$CF_RAW_key(
      Train_cf = train_df,
      min.split.var = min.split.var,
      variable_type = variable_type,
      hdPctTop = hdPctTop
    )
    if (length(raw_cf$selected_cf.idx) == 0L) {
      raw_cf$selected_cf.idx <- seq_len(ncol(X_df))
    }
    icf_env$selected_cf.idx <- raw_cf$selected_cf.idx
    icf_env$Y.hat <- raw_cf$Y.hat
    icf_env$W.hat <- raw_cf$W.hat
    icf_env$Train <- train_df
    icf_env$HTE_P_cf.raw <- raw_cf$HTE_P_cf.raw

    # Step 2 (optional): tune minimum leaf sizes; keep defaults on failure.
    if (isTRUE(tune_leafsize)) {
      tuned <- tryCatch(
        .run_mls_tuning(icf_env, train_df, treeNo, iterationNo, split_val_round_posi),
        error = function(e) NULL
      )
      if (!is.null(tuned$leafsize)) leafsize <- tuned$leafsize
    }

    # Step 3: cross-validated iCF, retried with larger leaf sizes on known
    # iCF failures.
    for (leafsize_try in .leafsize_retry_grid(leafsize)) {
      icf_env$leafsize <- leafsize_try
      fit <- tryCatch(
        icf_env$iCFCV(
          dat = train_df,
          K = K,
          treeNo = treeNo,
          iterationNo = iterationNo,
          min.split.var = min.split.var,
          split_val_round_posi = split_val_round_posi,
          P_threshold = P_threshold,
          variable_type = variable_type,
          hdPctTop = hdPctTop,
          HTE_P_cf.raw = raw_cf$HTE_P_cf.raw
        ),
        error = function(e) e
      )
      if (!inherits(fit, "error")) return(fit)
      if (!grepl("parent_sign|doesn't exist|no rows to aggregate",
                 conditionMessage(fit), ignore.case = TRUE)) break
    }
    stop("iCF failed: ", conditionMessage(fit), call. = FALSE)
  }))
}

# ---- Prediction ----------------------------------------------------------------

# iCF reports subgroup effects as strings such as "1.23(0.45,2.01)".
.parse_cate_iptw <- function(x) {
  as.numeric(sub("\\(.*$", "", as.character(x)))
}

# Assign units to iCF subgroups by evaluating the subgroup definitions; units
# matched by no definition go to the subgroup with the smallest ID.
.apply_subgroup_rules <- function(X_df, selected_sg) {
  if (is.null(selected_sg) || nrow(selected_sg) == 0) {
    return(rep(1L, nrow(X_df)))
  }
  nm <- colnames(selected_sg)
  if (all(c("subgroupID", "subgroup") %in% nm)) {
    ids <- selected_sg$subgroupID
    rules <- selected_sg$subgroup
  } else if (all(c("SubgroupID", "Definition") %in% nm)) {
    ids <- selected_sg$SubgroupID
    rules <- selected_sg$Definition
  } else {
    return(rep(1L, nrow(X_df)))
  }
  ids <- suppressWarnings(as.integer(ids))

  group_hat <- rep(NA_integer_, nrow(X_df))
  for (i in seq_along(ids)) {
    rule <- as.character(rules[i])
    if (!is.finite(ids[i]) || is.na(rule) || nchar(rule) == 0) next
    idx <- tryCatch(which(with(X_df, eval(parse(text = rule)))),
                    error = function(e) integer(0))
    group_hat[idx] <- ids[i]
  }

  valid_ids <- ids[is.finite(ids)]
  group_hat[is.na(group_hat)] <- if (length(valid_ids) > 0) min(valid_ids) else 1L
  as.integer(group_hat)
}

.extract_ate <- function(fit) {
  tbl <- fit$ATE_table
  ate <- if (is.data.frame(tbl) && nrow(tbl) > 0 && "CATE_iptw" %in% colnames(tbl)) {
    .parse_cate_iptw(tbl$CATE_iptw[1])
  } else {
    NA
  }
  if (is.finite(ate)) ate else 0
}

predict_ICF <- function(fit, X_new) {
  X_df <- as.data.frame(X_new)
  if (nrow(X_df) == 0) {
    return(list(yhat = numeric(0), group_hat = integer(0)))
  }

  selected_sg <- fit$selectedSG_ori
  if (is.list(selected_sg) && !is.data.frame(selected_sg)) {
    selected_sg <- selected_sg$majority
  }
  if (!is.data.frame(selected_sg)) selected_sg <- NULL
  group_hat <- .apply_subgroup_rules(X_df, selected_sg)

  yhat <- rep(NA_real_, length(group_hat))
  cate_tbl <- fit$CATE_t2_ori
  if (is.data.frame(cate_tbl) && all(c("SubgroupID", "CATE_iptw") %in% colnames(cate_tbl))) {
    cate_map <- unique(as.data.frame(cate_tbl)[, c("SubgroupID", "CATE_iptw")])
    cate_ids <- suppressWarnings(as.integer(cate_map$SubgroupID))
    yhat <- .parse_cate_iptw(cate_map$CATE_iptw)[match(group_hat, cate_ids)]
  }

  # Fall back to the overall ATE (or 0) where no subgroup effect is available.
  if (anyNA(yhat)) {
    if (all(is.na(yhat))) group_hat <- rep(1L, length(group_hat))
    yhat[is.na(yhat)] <- .extract_ate(fit)
  }

  list(yhat = as.numeric(yhat), group_hat = as.integer(group_hat))
}
