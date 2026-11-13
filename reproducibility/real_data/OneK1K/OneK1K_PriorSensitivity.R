###############################################################################
## OneK1K_PriorSensitivity.R  --  MR.RGM on the OneK1K data under a uniform prior on the causal inclusion probabilities
##
## The main analysis uses the package default Beta(3, 1) (prior inclusion probability 0.75 for a causal edge). This
## script refits MR.RGM with Beta(1, 1) (prior inclusion probability 0.5) on exactly the same data, with the chain
## length, seeds and convergence check of the main fit, and compares the two fits:
##   - convergence of the four chains (fit_chains() in ../real_data_functions.R);
##   - the causal PIPs of every ordered pair under both priors, their Spearman correlation, and the edges at
##     PIP > 0.50 and > 0.75 under each prior;
##   - the Bayes factors (BF) implied by each fit (posterior odds / prior odds), which should agree if the data, not the
##     prior, drive the ranking.
## Output: <ONEK_OUT_DIR>/output/OneK1K_fit_MRRGM_Beta11.rds, <ONEK_OUT_DIR>/tables/onek_prior_sensitivity.csv and
##         onek_prior_sensitivity.tex (every edge with PIP > 0.5 under either prior).
## Prerequisite: OneK1K_ModelFitting.R (MR.RGM), which writes output/OneK1K_analysis_data.rds and OneK1K_fit_MRRGM.rds.
## Usage (this folder):  ONEK_OUT_DIR <- "path/to/Results/OneK1K"; source("OneK1K_PriorSensitivity.R")
###############################################################################

suppressPackageStartupMessages(library(MR.RGM))
source("../real_data_functions.R")
OUT_DIR <- if (exists("ONEK_OUT_DIR")) ONEK_OUT_DIR else "."
f_dat <- file.path(OUT_DIR, "output", "OneK1K_analysis_data.rds"); f_main <- file.path(OUT_DIR, "output", "OneK1K_fit_MRRGM.rds")
if (!file.exists(f_dat) || !file.exists(f_main)) stop("Run OneK1K_ModelFitting.R (MR.RGM) first: ", f_dat, " / ", f_main)
dat  <- readRDS(f_dat); main <- readRDS(f_main)
S    <- summary_stats(dat$Y, dat$X, dat$U)
MCMC <- dat$settings$mcmc; SEEDS <- dat$settings$seeds
A_RHO <- 1; B_RHO <- 1                                                    # uniform prior on the inclusion probability

run_chain <- function(seed) {
  set.seed(seed); w <- character(0)
  f <- withCallingHandlers(
    RGM(Syy = S$Syy, Syx = S$Syx, Sxx = S$Sxx, Syu = S$Syu, Sxu = S$Sxu, Suu = S$Suu, D = dat$D, n = S$n, selectB = FALSE,
        nIter = MCMC$nIter, nBurnin = MCMC$nBurnin, Thin = MCMC$Thin, prior = "Spike and Slab", SigmaStarModel = "SSSL",
        aRho = A_RHO, bRho = B_RHO),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  f$warnings <- w; f
}
message(sprintf("MR.RGM, Beta(%g, %g) prior on the causal inclusion probabilities: %d chains of %s iterations",
                A_RHO, B_RHO, length(SEEDS), format(MCMC$nIter, big.mark = ",")))
fc <- fit_chains(run_chain, SEEDS)
saveRDS(list(fit = fc$fit, convergence = fc$convergence, prior = c(aRho = A_RHO, bRho = B_RHO), labels = dat$labels,
             package_version = as.character(packageVersion("MR.RGM"))),
        file.path(OUT_DIR, "output", "OneK1K_fit_MRRGM_Beta11.rds"))

## ---- comparison with the main fit (Beta(3, 1)) --------------------------------------------------------------------
lab <- dat$labels; p <- length(lab); off <- which(row(diag(p)) != col(diag(p)))
G1 <- main$fit$GammaEst; G2 <- fc$fit$GammaEst; A1 <- main$fit$AEst; A2 <- fc$fit$AEst
prior1 <- 3 / 4; prior2 <- A_RHO / (A_RHO + B_RHO)
bf <- function(pip, prior) (pip / (1 - pip)) / (prior / (1 - prior))
idx <- arrayInd(off, c(p, p))
cmp <- data.frame(edge = sprintf("%s -> %s", lab[idx[, 2]], lab[idx[, 1]]),
                  PIP_Beta31 = G1[off], PIP_Beta11 = G2[off], BF_Beta31 = bf(G1[off], prior1), BF_Beta11 = bf(G2[off], prior2),
                  effect_Beta31 = A1[off], effect_Beta11 = A2[off])
cmp <- cmp[order(-cmp$PIP_Beta31), ]
dir.create(file.path(OUT_DIR, "tables"), showWarnings = FALSE)
write.csv(cmp, file.path(OUT_DIR, "tables", "onek_prior_sensitivity.csv"), row.names = FALSE)

rho <- suppressWarnings(cor(cmp$PIP_Beta31, cmp$PIP_Beta11, method = "spearman"))
rho_bf <- suppressWarnings(cor(log(cmp$BF_Beta31), log(cmp$BF_Beta11)))
top19 <- head(cmp$edge, sum(cmp$PIP_Beta31 > 0.5))
rank2 <- cmp$edge[order(-cmp$PIP_Beta11)]
message(sprintf("Convergence under Beta(1,1): %s (split-R-hat %.3f / %.3f, smallest ESS %.0f, chains' PIP correlation %.3f)",
                if (fc$convergence$converged) "yes" else "NO", fc$convergence$rhat_loglik, fc$convergence$rhat_edges_max,
                fc$convergence$ess_bulk_edges_min, fc$convergence$pip_pairwise_min_cor))
message(sprintf("Edges at PIP > 0.50: %d under Beta(3,1), %d under Beta(1,1); at PIP > 0.75: %d and %d",
                sum(cmp$PIP_Beta31 > 0.5), sum(cmp$PIP_Beta11 > 0.5), sum(cmp$PIP_Beta31 > 0.75), sum(cmp$PIP_Beta11 > 0.75)))
message(sprintf("Spearman correlation of the causal PIPs: %.3f; correlation of the log Bayes factors: %.3f", rho, rho_bf))
message(sprintf("Of the %d edges at PIP > 0.50 under Beta(3,1), %d are among the top %d under Beta(1,1)",
                length(top19), sum(top19 %in% head(rank2, length(top19))), length(top19)))
print(head(transform(cmp, PIP_Beta31 = round(PIP_Beta31, 3), PIP_Beta11 = round(PIP_Beta11, 3), BF_Beta31 = round(BF_Beta31, 2),
                     BF_Beta11 = round(BF_Beta11, 2), effect_Beta31 = round(effect_Beta31, 3), effect_Beta11 = round(effect_Beta11, 3)),
           max(20, sum(cmp$PIP_Beta11 > 0.5))), row.names = FALSE)

## LaTeX table: every edge with PIP > 0.5 under either prior
tt <- cmp[cmp$PIP_Beta31 > 0.5 | cmp$PIP_Beta11 > 0.5, ]
L <- c("\\begin{table}[htb]", "\\centering", "\\small",
       sprintf("\\caption{OneK$1$K B cells, \\texttt{MR.RGM}: sensitivity to the prior on the causal inclusion probabilities. Edges with posterior inclusion probability (PIP) $> 0.50$ under the default $\\mathrm{Beta}(3,1)$ prior (prior inclusion probability $0.75$) or under $\\mathrm{Beta}(1,1)$ ($0.50$), with the PIP and the Bayes factor (BF, posterior odds divided by prior odds) under each, and the posterior mean effects. Same data, chain length and seeds. Four chains, %s under $\\mathrm{Beta}(1,1)$. Spearman correlation of the PIPs of all %d ordered pairs: %.2f.}",
               if (fc$convergence$converged) "converged" else "not converged", length(off), rho),
       "\\label{tab:onek_prior_sensitivity}", "\\begin{tabular}{|l|c|c|c|c|c|c|}", "\\hline",
       "\\textbf{Edge} & \\textbf{PIP, Beta(3,1)} & \\textbf{BF} & \\textbf{Effect} & \\textbf{PIP, Beta(1,1)} & \\textbf{BF} & \\textbf{Effect} \\\\", "\\hline",
       sprintf("%s & %.3f & %.2f & %+.3f & %.3f & %.2f & %+.3f \\\\", gsub(" -> ", " $\\\\rightarrow$ ", tt$edge),
               tt$PIP_Beta31, tt$BF_Beta31, tt$effect_Beta31, tt$PIP_Beta11, tt$BF_Beta11, tt$effect_Beta11),
       "\\hline", "\\end{tabular}", "\\end{table}")
writeLines(display_text(L), file.path(OUT_DIR, "tables", "onek_prior_sensitivity.tex"))   # displayed names (real_data_functions.R)
message("Saved output/OneK1K_fit_MRRGM_Beta11.rds and tables/onek_prior_sensitivity.{csv,tex} in ", OUT_DIR)
