# Monte Carlo simulation study (14 conditions x 100 replicates).
#
#   Rscript my_sim/run_sim.R                          # all conditions
#   Rscript my_sim/run_sim.R <scenario>               # one row of sim_grid
#   Rscript my_sim/run_sim.R <scenario> <from> <to>   # a range of replicates
#
# Each replicate is saved to results/<condition>/ and skipped if it already
# exists, so interrupted runs can be resumed. Replicates run in parallel on
# N_CORES workers (environment variable; default: all cores but one).

# ---- Paths ------------------------------------------------------------------
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
my_sim_dir <- normalizePath(
  if (length(file_arg) > 0) dirname(sub("^--file=", "", file_arg[1])) else getwd(),
  winslash = "/"
)
if (basename(my_sim_dir) != "my_sim" && dir.exists(file.path(my_sim_dir, "my_sim"))) {
  my_sim_dir <- file.path(my_sim_dir, "my_sim")
}
project_dir <- dirname(my_sim_dir)
icf_dir <- file.path(project_dir, "iCF")
results_dir <- file.path(project_dir, "results")
dir.create(results_dir, showWarnings = FALSE)

# ---- Packages -----------------------------------------------------------------
# Packages are attached, not only loaded: SuperLearner looks up its learners
# (SL.glm, SL.xgboost, SL.gam) by name, and the iCF scripts call several
# packages without namespace prefixes. The same list is attached on workers.
pkgs <- c("CSurf", "CRE", "SuperLearner", "xgboost", "aricode", "causalTree",
          "glmnet", "gam", "mgcv", "grf", "partykit", "rpart", "MASS", "dplyr",
          "rlang", "rlist", "tidyverse", "caret", "knitr", "cowplot", "ggridges")
all_pkgs <- c("doFuture", pkgs)
missing_pkgs <- all_pkgs[!vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Please install: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
for (pkg in all_pkgs) library(pkg, character.only = TRUE)

# ---- Code ---------------------------------------------------------------------
for (f in list.files(file.path(project_dir, "code"), pattern = "\\.R$",
                     full.names = TRUE)) {
  source(f)
}
source(file.path(my_sim_dir, "myfun.R"))
source(file.path(my_sim_dir, "compare_estimators.R"))

# ---- Simulation conditions -----------------------------------------------------
# Reference condition: DGP 1, n = 2000, moderate overlap, pi = 0.5.
cond <- function(example_id = 1L, n = 2000L, pi = 0.50, overlap = "moderate",
                 f_hat0_type = "oracle") {
  data.frame(example_id, n, pi, overlap, f_hat0_type, stringsAsFactors = FALSE)
}
sim_grid <- rbind(
  # Initialization study in the reference condition
  cond(f_hat0_type = "gam_smooth"), cond(f_hat0_type = "xgb_smooth"),
  cond(f_hat0_type = "SuperLearner"), cond(f_hat0_type = "change_plane"),
  cond(),
  # One factor at a time, oracle initialization
  cond(overlap = "strong"), cond(n = 500L),
  cond(example_id = 2L), cond(example_id = 3L),
  cond(example_id = 4L), cond(example_id = 5L),
  cond(n = 10000L), cond(pi = 0.75), cond(overlap = "weak")
)
n_reps <- 100L

args <- as.integer(commandArgs(trailingOnly = TRUE))
scenarios <- if (length(args) >= 1L) args[1L] else seq_len(nrow(sim_grid))
if (any(is.na(scenarios) | scenarios < 1L | scenarios > nrow(sim_grid))) {
  stop("scenario must be between 1 and ", nrow(sim_grid), call. = FALSE)
}
reps <- seq(if (length(args) >= 2L) args[2L] else 1L,
            if (length(args) >= 3L) args[3L] else n_reps)

n_cores <- as.integer(Sys.getenv(
  "N_CORES", unset = max(1L, parallel::detectCores() - 1L, na.rm = TRUE)))

# Objects shipped to the parallel workers.
export_list <- c(Filter(function(x) is.function(get(x)), ls(all.names = TRUE)),
                 "estimators", "icf_dir")

# ---- Run ------------------------------------------------------------------------
run_scenario <- function(sc) {
  tag <- sprintf("ex%d_n%d_pi%.2f_%s_%s",
                 sc$example_id, sc$n, sc$pi, sc$overlap, sc$f_hat0_type)
  sc_dir <- file.path(results_dir, tag)
  dir.create(sc_dir, showWarnings = FALSE)

  done <- as.integer(sub(".*_rep(\\d{4})\\.rds$", "\\1",
                         list.files(sc_dir, pattern = "_rep\\d{4}\\.rds$")))
  todo <- setdiff(reps, done)
  message(tag, ": ", length(todo), " replicate(s) to run")
  if (length(todo) == 0L) return(invisible(NULL))

  run_rep <- function(i) {
    myfun(i, n = sc$n, example_id = sc$example_id, pi = sc$pi,
          overlap = sc$overlap, f_hat0_type = sc$f_hat0_type,
          save_path = sc_dir)
  }

  if (n_cores > 1L && length(todo) > 1L) {
    registerDoFuture()
    plan(multisession, workers = min(n_cores, length(todo)))
    on.exit(plan(sequential), add = TRUE)
    options(doFuture.rng.onMisuse = "ignore")
    foreach(i = todo, .packages = pkgs, .export = export_list) %dopar% run_rep(i)
  } else {
    lapply(todo, run_rep)
  }
  invisible(NULL)
}

for (s in scenarios) run_scenario(sim_grid[s, ])
