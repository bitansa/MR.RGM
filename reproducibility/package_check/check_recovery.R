###############################################################################
## check_recovery.R  --  check of the MR.RGM samplers before the full runs
##
## Fits all 24 sampler variants (2 priors x 3 error models x instruments yes/no x
## covariates yes/no) to one simulated data set with a known cyclic network and
## correlated errors, each with the known instrument map (MR.RGM, selectB = FALSE) and,
## for the variants with instruments, also with D = 1 and selection (MR.RGM+), and prints:
##   AUC     : area under the ROC curve of GammaEst for the true graph
##   MAD     : mean absolute deviation |AEst - A| over off-diagonal entries
##   confAUC : AUC of ZEst for the true error-covariance pattern (SSSL only; meaningful only
##             with covariates = TRUE, because the omitted covariates induce correlation
##             between all pairs of responses)
## A correct sampler gives AUC close to 1 and small MAD at n = 5000.
##
## Usage: source("check_recovery.R")      (about 5-10 minutes)
###############################################################################

library(MR.RGM)
library(pROC)
message("MR.RGM version: ", packageVersion("MR.RGM"))
if (packageVersion("MR.RGM") < "0.1.3" || !"gibbsB" %in% names(formals(MR.RGM::RGM))) stop("MR.RGM >= 0.1.3 is required (argument gibbsB of RGM())")

set.seed(2026)
n <- 5000; p <- 5; ivs <- 3; k <- p * ivs; l <- 2

## Cyclic network with a 2-cycle (1 <-> 2) and a 3-cycle (3 -> 4 -> 5 -> 3); A[to, from]
A <- matrix(0, p, p)
A[2, 1] <- 0.3; A[1, 2] <- -0.3; A[4, 3] <- 0.3; A[5, 4] <- 0.3; A[3, 5] <- -0.3
D <- matrix(0, p, k); for (i in 1:p) D[i, ((i - 1) * ivs + 1):(i * ivs)] <- 1

## Errors with two correlated pairs: (1, 3) and (2, 5)
Sigma <- diag(p); Sigma[1, 3] <- Sigma[3, 1] <- 0.5; Sigma[2, 5] <- Sigma[5, 2] <- -0.5
X <- matrix(rnorm(n * k), n, k); U <- matrix(rnorm(n * l), n, l)
E <- matrix(rnorm(n * p), n, p) %*% chol(Sigma)          # rows ~ N(0, Sigma)
Y <- (X %*% t(D) + U %*% matrix(c(0.5, -0.5), l, p) + E) %*% t(solve(diag(p) - A))
X <- scale(X); Y <- scale(Y); U <- scale(U)

true_graph <- (A != 0)[row(A) != col(A)] * 1
true_conf  <- (Sigma != 0)[upper.tri(Sigma)] * 1
auc <- function(truth, score) as.numeric(pROC::auc(pROC::roc(truth, score, levels = c(0, 1), direction = "<", quiet = TRUE)))

grid <- expand.grid(prior = c("Spike and Slab", "Threshold"), Sigma = c("SSSL", "IW", "diagonal"),
                    instruments = c(TRUE, FALSE), covariates = c(TRUE, FALSE), stringsAsFactors = FALSE)
res <- vector("list", nrow(grid))
for (r in seq_len(nrow(grid))) {
  g <- grid[r, ]
  args <- list(Y = Y, D = D, prior = g$prior, SigmaStarModel = g$Sigma, nIter = 6000, nBurnin = 2000, Thin = 2, selectB = FALSE)
  if (g$instruments) args$X <- X
  if (g$covariates)  args$U <- U
  set.seed(100 + r)
  t0  <- proc.time()[["elapsed"]]
  ## warnings of the package are recorded per fit (column warnings) instead of printed at the end
  fit_w <- function(args) { w <- character(0)
    f <- withCallingHandlers(do.call(RGM, args), warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
    f$warn <- paste(unique(sub("[:(].*", "", w)), collapse = " | "); f }
  fit <- fit_w(args)
  mkrow <- function(fit, version, t0) data.frame(g, version = version,
                         AUC     = round(auc(true_graph, fit$GammaEst[row(A) != col(A)]), 3),
                         MAD     = round(mean(abs(fit$AEst - A)[row(A) != col(A)]), 3),
                         confAUC = if (g$Sigma == "SSSL") round(auc(true_conf, fit$ZEst[upper.tri(Sigma)]), 3) else NA,
                         ivAUC   = if (version == "MR.RGM+") round(auc(as.vector(D), as.vector(fit$PhiEst)), 3) else NA,
                         accA = round(fit$AccptA, 1), halves = round(if (is.null(fit$diagnostics$halvesMaxDiff)) NA else fit$diagnostics$halvesMaxDiff, 2),
                         seconds = round(proc.time()[["elapsed"]] - t0, 1), warnings = fit$warn)
  res[[r]] <- mkrow(fit, "MR.RGM", t0)
  if (g$instruments) {                                   # MR.RGM+: D = 1, instrument effects selected
    args$D <- matrix(1, p, k); args$selectB <- TRUE
    set.seed(200 + r); t0 <- proc.time()[["elapsed"]]
    fit <- fit_w(args)
    res[[r]] <- rbind(res[[r]], mkrow(fit, "MR.RGM+", t0))
  }
  print(res[[r]], row.names = FALSE)
}
res <- do.call(rbind, res)
options(width = 250)
print(res[, setdiff(names(res), "warnings")], row.names = FALSE)
write.csv(res, "check_recovery_results.csv", row.names = FALSE)     # every column, incl. the warnings of each fit
message("Saved check_recovery_results.csv")

## Without instruments the network is identified only through the error structure, so lower
## AUC is expected in the instruments = FALSE rows; with instruments every row should have AUC
## near 1 and MAD well below 0.1; the MR.RGM+ rows should also have ivAUC near 1 (PhiEst recovers D).
