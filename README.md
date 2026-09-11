# IDR-Learner

Simulation code for

> Y. Yang et al., "IDR-Learner: Interpretable Doubly Robust Subgroup Discovery in Observational Biomedical Data."

## Repository layout

```
code/
  IDR_Learner.R         IDR-Learner and the shared cross-fitting helpers
  CSurf_init.R          CSurf initializers: SuperLearner, XGBoost, GAM, change-plane
  comparators.R         DR-/R-/X-learner and causal forest + CSurf / CART / CRE; causal tree
  ICF.R                 iterative causal forest (iCF) wrapper and predict_ICF()
  predict.R             CATE and subgroup predictions for new data
  dgp_Fun.R             data-generating processes (DGPs 1-5)
my_sim/
  run_sim.R             runs the simulation study
  myfun.R               one Monte Carlo replicate
  compare_estimators.R  the 15 compared estimators
```

## Requirements

R (>= 4.3) with

- `CSurf`
- `SuperLearner`, `xgboost`, `gam`, `mgcv`, `grf`, `rpart`, `partykit`, `CRE`, `MASS`, `aricode`, `doFuture`
- `causalTree`: `remotes::install_github("susanathey/causalTree")`
- for iCF: `glmnet`, `dplyr`, `rlang`, `rlist`, `tidyverse`, `caret`, `knitr`, `cowplot`, `ggridges`

### iCF scripts

The iCF comparator uses the scripts of Wang et al. (2024), which are not included in this repository. Download them from <https://github.com/tianshengwang/iCF> and copy the following files from its `R/` folder into a folder named `iCF/` at the top level of this repository:

```
best_tree_MSegar.R        iCF_SUBGROUP_ANALYSIS.R
iCF_CV.R                  iCF_SUBGROUP_DECISION.R
iCF_GG_toolbox.R          iCF_SUBGROUP_MODEL.R
iCF_MAJORITY_VOTE.R       iCF_SUBGROUP_PIPELINE.R
iCF_PARENT_node.R         iCF_TREE_build.R
iCF_PRE_majority.R        sim_Truth_tree.R
```

These scripts are only needed for the iCF comparator; IDR-Learner and the other estimators run without them.

## Using IDR-Learner

```r
library(SuperLearner)   # the learners SL.glm / SL.xgboost are looked up by name
for (f in list.files("code", pattern = "\\.R$", full.names = TRUE)) source(f)

set.seed(1)
train <- dgp_Fun(n = 2000, example_id = 1)
test  <- dgp_Fun(n = 2000, example_id = 1)

fit  <- IDR_Learner(train$X, train$W, train$Y, f_hat0 = "SuperLearner")
pred <- predict_CSurf(fit, test$X)
head(pred$yhat)       # CATE estimates
table(pred$group_hat) # estimated subgroups
```

`f_hat0` selects the CSurf initializer: `"SuperLearner"` (primary), `"xgb_smooth"`, `"gam_smooth"`, `"change_plane"`, or a numeric score (for example the oracle score `train$f_hat0_init`).

## Reproducing the simulation study

The study has 14 conditions with 100 replicates each:

| Scenario | DGP | n | Overlap | Prevalence | Initializer |
|---|---|---|---|---|---|
| 1-4 | 1 | 2000 | moderate | 0.50 | GAM, XGBoost, SuperLearner, change-plane |
| 5 | 1 | 2000 | moderate | 0.50 | oracle (reference condition) |
| 6 | 1 | 2000 | strong | 0.50 | oracle |
| 7 | 1 | 500 | moderate | 0.50 | oracle |
| 8-11 | 2-5 | 2000 | moderate | 0.50 | oracle |
| 12 | 1 | 10000 | moderate | 0.50 | oracle |
| 13 | 1 | 2000 | moderate | 0.75 | oracle |
| 14 | 1 | 2000 | weak | 0.50 | oracle |

```bash
Rscript my_sim/run_sim.R            # all scenarios
Rscript my_sim/run_sim.R 5          # one scenario
Rscript my_sim/run_sim.R 5 1 50     # replicates 1-50 of scenario 5
```

The number of parallel workers is set with the `N_CORES` environment variable. Each replicate is saved as `results/<condition>/result_..._repXXXX.rds`, and existing files are skipped, so an interrupted run can be resumed. To combine the results:

```r
files   <- list.files("results", pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
results <- do.call(rbind, lapply(files, readRDS))
results$relmse <- results$mse / results$true.cate.var
```
