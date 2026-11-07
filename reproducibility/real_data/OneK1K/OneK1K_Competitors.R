###############################################################################
## OneK1K_Competitors.R  —  every competing method of the simulation study on the OneK1K data of the MR.RGM fit
##
## Prerequisite: OneK1K_ModelFitting.R (MR.RGM), which saves the analysed data to
##   <ONEK_OUT_DIR>/output/OneK1K_analysis_data.rds  (standardized expression Y, genotypes X, covariates U, own-gene map D)
## so that every method sees exactly the same genes, SNPs and donors.
##
## Methods (wrappers and package settings of the simulation study, ../../simulation/sim_functions.R):
##   MR.RGM_NoConf          MR.RGM with a diagonal error covariance (no confounding), covariates as in the main fit;
##                          four chains pooled when converged, as MR.RGM (fit_chains() in ../real_data_functions.R)
##   MVMR                   multivariable IVW, one gene at a time as outcome, all other genes as exposures
##   MR2G                   network of total effects deconvolved to direct effects (package defaults, ScreenAug)
##   MrDAG                  one gene at a time as outcome, the other genes as exposures (by_outcome())
##   GrIVET                 causal DAG with hidden confounding from instruments (all SNPs candidates for every gene);
##                          the tuning grid of the simulation study, fitted to the covariate-adjusted UNSTANDARDIZED data
##                          (its tuning is not scale free; fit_grivet() in sim_functions.R), effects reported on the
##                          scale of the other summary-level methods
##   MR-IVW, MR-SimpleMedian, MR-WeightedMedian (MendelianRandomization), OneSampleMR (TSRI), mrbayes   every ordered pair
##   (MR2 needs an exposure/outcome split with at least two responses and does not apply to a network of genes.)
## Covariates: the summary-level competitors have no covariate argument, so expression and genotypes are adjusted for
## the covariates (residuals of the regression on U, then standardized; Frisch-Waugh), which gives the same
## covariate-adjusted associations that MR.RGM uses through its covariate model.
## MCMC methods use the chain length of the real-data MR.RGM fit (100,000 iterations, 20,000 burn-in, thinning 10;
## mrbayes: 3 chains of 20,000 burn-in + 80,000 iterations); everything else is the package default.
## Edge rules as in the simulations: p-value < 0.05 (no multiplicity correction), posterior inclusion probability > 0.5.
##
## Output: <ONEK_OUT_DIR>/output/OneK1K_competitors.rds  (per method: effects and scores [to, from], edge matrix, rule,
##         seconds, failure message); the comparison tables are written by OneK1K_Results.R.
## Usage (this folder, after OneK1K_ModelFitting.R):  source("OneK1K_Competitors.R")
## Option: ONEK_COMP_METHODS (e.g. "GrIVET"): fit only these methods and keep the other methods' results from the
##   existing output/OneK1K_competitors.rds (to refit a subset of the methods). GrIVET needs the scale of the data before
##   standardization: sd_Y / sd_X in OneK1K_analysis_data.rds, or output/OneK1K_scales.rds (written by
##   OneK1K_ModelFitting.R with ONEK_STRENGTH_ONLY <- TRUE).
###############################################################################

suppressPackageStartupMessages(library(MR.RGM))
source("../real_data_functions.R")
OUT_DIR <- if (exists("ONEK_OUT_DIR")) ONEK_OUT_DIR else "."
dat_file <- file.path(OUT_DIR, "output", "OneK1K_analysis_data.rds")
if (!file.exists(dat_file)) stop("Run OneK1K_ModelFitting.R (MR.RGM) first: ", dat_file, " not found.")
dat <- readRDS(dat_file)
Y <- dat$Y; X <- dat$X; U <- dat$U; D <- dat$D; labels <- dat$labels; p <- ncol(Y)
stopifnot(identical(colnames(Y), dat$symbols), identical(dim(D), c(ncol(Y), ncol(X))))

## The simulation study's wrappers and settings, in reproducibility/simulation (the second path is an alternative layout)
sim_file <- c(file.path("..", "..", "simulation", "sim_functions.R"), file.path("..", "..", "..", "Simulation Studies", "sim_functions.R"))
sim_file <- sim_file[file.exists(sim_file)][1]
if (is.na(sim_file)) stop("sim_functions.R of the simulation study not found")
sf <- new.env(); sys.source(sim_file, envir = sf, toplevel.env = sf)

MCMC <- dat$settings$mcmc                                     # chain length of the MR.RGM fit
set  <- sf$SIM_SETTINGS
set$mrdag   <- modifyList(set$mrdag, list(niter = MCMC$nIter, burnin = MCMC$nBurnin, thin = MCMC$Thin))
set$mrbayes <- list(n.chains = 3, n.burn = MCMC$nBurnin, n.iter = MCMC$nIter - MCMC$nBurnin)
SEED <- dat$settings$seed

## covariate-adjusted, standardized data for the summary-level and pairwise methods
resid_on <- function(M) M - U %*% solve(crossprod(U), crossprod(U, M))
Yr <- standardize_cols(resid_on(Y), "Expression (covariate-adjusted)")
Xr <- standardize_cols(resid_on(X), "Genotypes (covariate-adjusted)")
stopifnot(ncol(Yr) == p, ncol(Xr) == ncol(X))
ms <- sf$marginal_stats(Xr, Yr); attr(ms, "X") <- Xr

run <- function(name, expr_fun) {
  set.seed(SEED + sf$METHOD_SEED[[name]])
  t0 <- proc.time()[["elapsed"]]
  out <- tryCatch(expr_fun(), error = function(e) structure(list(error = conditionMessage(e)), class = "failed"))
  secs <- proc.time()[["elapsed"]] - t0
  if (inherits(out, "failed")) message(sprintf("%s failed: %s", name, out$error)) else message(sprintf("%s: %.0f s", name, secs))
  attr(out, "seconds") <- secs
  out
}
## results [to, from]; score = 1 - p-value or the inclusion probability (NA: pair not estimated by the method);
## seconds and the number of failed pairs / outcomes of the fit are kept
pval_result <- function(r, est, pval, rule = "$p < 0.05$")
  list(est = est, score = 1 - pval, pval = pval, edge = !is.na(pval) & pval < set$alpha, rule = rule,
       seconds = attr(r, "seconds"), failed_pairs = if (is.null(r$failed)) 0 else r$failed)
pip_result <- function(r, est, pip, rule = "PIP $> 0.50$")
  list(est = est, score = pip, edge = !is.na(pip) & pip > set$pip_thr, rule = rule,
       seconds = attr(r, "seconds"), failed_pairs = if (is.null(r$failed)) 0 else r$failed)

RUN <- if (exists("ONEK_COMP_METHODS")) ONEK_COMP_METHODS else
  c("MR.RGM_NoConf", "MVMR", "MR2G", "MrDAG", "MendelianRandomization", "OneSampleMR", "mrbayes", "GrIVET")
out_file <- file.path(OUT_DIR, "output", "OneK1K_competitors.rds")
res <- if (exists("ONEK_COMP_METHODS") && file.exists(out_file)) readRDS(out_file)$methods else list()
S <- summary_stats(Y, X, U)
noconf_chain <- function(seed) {
  set.seed(seed); w <- character(0)
  f <- withCallingHandlers(RGM(Syy = S$Syy, Syx = S$Syx, Sxx = S$Sxx, Syu = S$Syu, Sxu = S$Sxu, Suu = S$Suu, D = D, n = S$n,
                               selectB = FALSE, nIter = MCMC$nIter, nBurnin = MCMC$nBurnin, Thin = MCMC$Thin,
                               prior = "Spike and Slab", SigmaStarModel = "diagonal"),
                           warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  f$warnings <- w; f }
seeds_nc <- (if (is.null(dat$settings$seeds)) SEED else dat$settings$seeds) + sf$METHOD_SEED[["MR.RGM_NoConf"]]
if ("MR.RGM_NoConf" %in% RUN) {
  r <- run("MR.RGM_NoConf", function() fit_chains(noconf_chain, seeds_nc))
  res$MR.RGM_NoConf <- if (inherits(r, "failed")) r else {
    x <- pip_result(r, r$fit$AEst, r$fit$GammaEst); x$convergence <- r$convergence; x$seconds <- attr(r, "seconds"); x }
}
if ("MVMR" %in% RUN)  { r <- run("MVMR", function() sf$fit_mvmr(ms, D));                         res$MVMR  <- if (inherits(r, "failed")) r else pval_result(r, r$est, r$pval) }
if ("MR2G" %in% RUN)  { r <- run("MR2G", function() suppressWarnings(sf$fit_mr2g(ms, Yr, D, set))); res$MR2G <- if (inherits(r, "failed")) r else pval_result(r, r$est, r$pval) }
if ("MrDAG" %in% RUN) { r <- run("MrDAG", function() sf$by_outcome(sf$fit_mrdag, ms, D, set));      res$MrDAG <- if (inherits(r, "failed")) r else pip_result(r, r$est, r$score) }
if ("MendelianRandomization" %in% RUN) {
  r <- run("MendelianRandomization", function() sf$fit_mr_package(Xr, Yr, D, set))
  if (inherits(r, "failed")) { for (m in c("MR-IVW", "MR-SimpleMedian", "MR-WeightedMedian")) res[[m]] <- r } else
    for (m in c("IVW", "SimpleMedian", "WeightedMedian")) res[[paste0("MR-", m)]] <- pval_result(r, t(r$est[[m]]), t(r$pval[[m]]))
}
if ("OneSampleMR" %in% RUN) { r <- run("OneSampleMR", function() sf$fit_onesamplemr(Xr, Yr, D)); res$OneSampleMR <- if (inherits(r, "failed")) r else pval_result(r, t(r$est), t(r$pval)) }
if ("mrbayes" %in% RUN)     { r <- run("mrbayes", function() sf$fit_mrbayes(Xr, Yr, D, set));   res$mrbayes <- if (inherits(r, "failed")) r else pval_result(r, t(r$est), t(r$pval), "two-sided posterior tail probability $< 0.05$") }
## GrIVET: the covariate-adjusted data on their original scale (centered expression and genotype dosages), all SNPs
## candidates for every gene; its effects are converted to the scale of Yr (the covariate-adjusted, standardized
## expression of the other summary-level methods) as in the simulation study
if ("GrIVET" %in% RUN) {
  scl <- if (!is.null(dat$sd_Y)) dat[c("sd_Y", "sd_X")] else {
    sf_ <- file.path(OUT_DIR, "output", "OneK1K_scales.rds")
    if (!file.exists(sf_)) stop("GrIVET needs the scale of the data before standardization: run OneK1K_ModelFitting.R with ONEK_STRENGTH_ONLY <- TRUE (writes ", sf_, ")")
    readRDS(sf_) }
  stopifnot(identical(names(scl$sd_Y), colnames(Y)), identical(names(scl$sd_X), colnames(X)))
  Y_orig <- sweep(Y, 2, scl$sd_Y, `*`); X_orig <- sweep(X, 2, scl$sd_X, `*`)
  Yr_orig <- resid_on(Y_orig); Xr_orig <- resid_on(X_orig)
  r <- run("GrIVET", function() sf$fit_grivet(Xr_orig, Yr_orig, apply(Yr_orig, 2, sd), set))
  res$GrIVET <- if (inherits(r, "failed")) r else
    list(est = r$est, score = r$score, edge = r$edge & !is.na(r$edge), rule = "selected by \\texttt{GrIVET} (effect $\\neq 0$)",
         seconds = attr(r, "seconds"), failed_pairs = 0)
}

ord <- c("MR.RGM_NoConf", "MVMR", "MR2G", "MrDAG", "GrIVET", "MR-IVW", "MR-SimpleMedian", "MR-WeightedMedian", "OneSampleMR", "mrbayes")
res <- res[c(intersect(ord, names(res)), setdiff(names(res), ord))]   # the order of the simulation study
for (m in names(res)) if (!inherits(res[[m]], "failed")) {
  for (k in intersect(c("est", "score", "pval", "edge"), names(res[[m]]))) dimnames(res[[m]][[k]]) <- list(labels, labels)
  diag(res[[m]]$edge) <- FALSE
}
out <- list(methods = res, labels = labels, n = nrow(Y), mcmc = MCMC, covariates = colnames(U),
            settings = set[c("alpha", "pip_thr", "mr_boot", "mrbayes", "mr2g", "mrdag", "grivet")],
            package_versions = vapply(c("MR.RGM", "MendelianRandomization", "OneSampleMR", "mrbayes", "rjags", "MR2G", "MrDAG", "grivet"),
                                      function(pk) if (requireNamespace(pk, quietly = TRUE)) as.character(packageVersion(pk)) else NA_character_, character(1)),
            session_info = sessionInfo())
of <- out_file
saveRDS(out, of)
message("Saved ", of, "; edges per method: ",
        paste(sprintf("%s %s", names(res), vapply(res, function(x) if (inherits(x, "failed")) "failed" else as.character(sum(x$edge)), "")), collapse = ", "))
