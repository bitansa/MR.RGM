###############################################################################
## pilot_competitors.R  --  one-off check of the competitor wrappers (run before the simulations)
##
## Generates one replicate of the "roles" design with a strong, known structure (4 exposures -> 4 outcomes) and runs
## each competitor once with short settings. It prints, for each method, the raw output shapes and the estimated
## [to, from] matrix next to the truth, so that the orientation used in the wrappers of sim_functions.R
## can be confirmed by inspection: the large entries must appear where the truth has edges.
##
## Usage (from this folder): source("pilot_competitors.R")      # 1-5 minutes
###############################################################################

source("sim_functions.R")
set <- SIM_SETTINGS
set$mr2g  <- modifyList(set$mr2g, list(num_pert = 10, num_pert_M = 10))
set$mr2   <- modifyList(set$mr2, list(niter = 1500, burnin = 500, thin = 2))
set$mrdag <- modifyList(set$mrdag, list(niter = 1500, burnin = 500, thin = 2))
set$mcmc  <- list(nIter = 3000, nBurnin = 1000, Thin = 2)

for (pk in c("MR2G", "MR2", "MrDAG")) message(pk, ": ", if (has_pkg(pk)) paste("installed, version", packageVersion(pk)) else "NOT installed")

set.seed(1); truth0 <- gen_truth("roles", 8, set)
Xi <- truth0$roles$X_idx; Yi <- truth0$roles$Y_idx
n <- 5000
set.seed(2); dat0 <- gen_data(truth0, n, set)
std <- standardize_data(dat0, truth0, set); X <- std$X; Y <- std$Y; truth <- std$truth; D <- truth$D
ms <- marginal_stats(X, Y); attr(ms, "X") <- X
cat(sprintf("\nTRUE A [to, from] (rows = outcomes %s receive edges from exposures %s and earlier outcomes):\n", paste(Yi, collapse = ","), paste(Xi, collapse = ","))); print(round(truth$A, 2))
cat("\nTrue confounding |Sigma*| (standardized):\n"); print(round(truth$Sigma_star, 2))
show <- function(name, est, score = NULL) {
  cat("\n==== ", name, " ====\n"); cat("estimate [to, from]:\n"); print(round(est, 2))
  if (!is.null(score)) { cat("score (PIP or 1 - p):\n"); print(round(score, 2)) }
}

## MR.RGM
set.seed(11); f <- fit_rgm(X, Y, D, "SSSL", set, selectB = FALSE); show("MR.RGM", f$AEst, f$GammaEst); cat("ZEst:\n"); print(round(f$ZEst, 2))
## MVMR
set.seed(12); r <- fit_mvmr(ms, D); show("MVMR", r$est, 1 - r$pval)
## MR2G (wrapper; ScreenAug igraph >= 2.0 compatible version, see sim_functions.R)
if (has_pkg("MR2G")) {
  set.seed(13); r <- tryCatch(fit_mr2g(ms, Y, D, set), error = function(e) { message("fit_mr2g failed: ", conditionMessage(e)); NULL })
  if (!is.null(r)) show("MR2G (wrapper)", r$est, 1 - r$pval)
}
cat("\nTrue exposure -> outcome block A[outcomes, exposures] (for MR2 and MrDAG):\n"); print(round(truth$A[Yi, Xi], 2))
## MR2
if (has_pkg("MR2")) {
  set.seed(14); r <- tryCatch(fit_mr2(ms, D, truth$roles, set), error = function(e) { message("fit_mr2 failed: ", conditionMessage(e)); NULL })
  if (!is.null(r)) { show("MR2 (wrapper; rows = outcomes, columns = exposures should hold the effects)", r$est, r$score)
    cat("residual correlation among outcomes:\n"); print(round(r$resid_cor, 2)) }
}
## MrDAG
if (has_pkg("MrDAG")) {
  set.seed(15); r <- tryCatch(fit_mrdag(ms, D, truth$roles, set), error = function(e) { message("fit_mrdag failed: ", conditionMessage(e)); NULL })
  if (!is.null(r)) show("MrDAG (wrapper; exposure -> outcome effects, PEIP among outcomes and among exposures)", r$est, r$score)
  cat("The true outcome -> outcome edges are:\n"); A <- truth$A; yy <- which(A[Yi, Yi] != 0, arr.ind = TRUE); print(data.frame(from = Yi[yy[, 2]], to = Yi[yy[, 1]]))
}
cat("\nPilot finished.\n")
