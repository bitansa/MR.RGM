###############################################################################
## OneK1K_Bootstrap.R  --  how reliable is the fitted OneK1K network at this sample size and instrument strength?
##
## Parametric bootstrap with the real genotypes and covariates: the fitted MR.RGM network is taken as the truth and
## new expression data are generated from it, R times, for the same 891 donors; MR.RGM is refitted to every data set.
##   Truth   A0 = posterior mean of A on the edges with PIP > 0.5 (zero elsewhere), B0, C0 = posterior means of B and C,
##           Sigma0 = posterior mean of Sigma* (all from the pooled four-chain fit, OneK1K_fit_MRRGM.rds)
##   Data    the analysed (standardized) genotypes X and covariates U of OneK1K_analysis_data.rds;
##           Y = (B0 X + C0 U + E)(I - A0)^{-T} row-wise, E ~ N(0, Sigma0); every column of Y centered and scaled,
##           as in the analysis; the truth is expressed on the scale of the standardized Y of each replicate
##   Fit     MR.RGM with the own-gene map, the package defaults and the MCMC settings of the simulation study
##           (one chain of 50,000 iterations, 10,000 burn-in, thinning 10)
## Reported (mean and SD over the replicates):
##   network AUC of the causal PIPs; edges selected at PIP > 0.5; TPR and false discovery proportion at PIP > 0.5;
##   false discovery proportion at PIP > 0.85; MAD of the estimated effects on the true edges; for every true 2-cycle,
##   the share of replicates in which both directions are selected; selected 2-cycles that are not true 2-cycles
##   (all, and those with effects of opposite signs); share of fits with each warning of the package (weak
##   instruments, error variances, spectral radius, chain halves, acceptance rate) and with any warning. Per edge of the fitted network: mean PIP (SD), share of replicates in which it is selected,
##   and share in which it is selected with the sign of the fitted effect; false edges selected in at least 20% of
##   the replicates are listed as well.
## The bootstrap is conditional on the fitted network being the truth; the effects of selected edges are, if anything,
## inflated by the selection, so the error rates are, if anything, optimistic.
##
## Prerequisite: OneK1K_ModelFitting.R (MR.RGM), which writes output/OneK1K_analysis_data.rds and OneK1K_fit_MRRGM.rds.
## Output: <ONEK_OUT_DIR>/output/OneK1K_bootstrap.rds, <ONEK_OUT_DIR>/tables/onek_bootstrap.{csv,tex} (summary) and
##         onek_bootstrap_edges.{csv,tex} (per edge).
## Usage (this folder):  ONEK_OUT_DIR <- "path/to/Results/OneK1K"; source("OneK1K_Bootstrap.R")
##   Options before sourcing: BOOT_N_REP (default 100), N_CORES (default: all cores but one).
##   Time: about one minute per replicate and core (p = 10, k = 48).
###############################################################################

suppressPackageStartupMessages({ library(MR.RGM); library(pROC); library(parallel) })
source("../real_data_functions.R")
OUT_DIR <- if (exists("ONEK_OUT_DIR")) ONEK_OUT_DIR else "."
f_dat <- file.path(OUT_DIR, "output", "OneK1K_analysis_data.rds"); f_main <- file.path(OUT_DIR, "output", "OneK1K_fit_MRRGM.rds")
if (!file.exists(f_dat) || !file.exists(f_main)) stop("Run OneK1K_ModelFitting.R (MR.RGM) first: ", f_dat, " / ", f_main)
dat  <- readRDS(f_dat); main <- readRDS(f_main); fit0 <- main$fit
N_REP   <- if (exists("BOOT_N_REP")) BOOT_N_REP else 100
N_CORES <- if (exists("N_CORES")) N_CORES else max(1, parallel::detectCores() - 1)
MCMC    <- list(nIter = 50000, nBurnin = 10000, Thin = 10)                # as in the simulation study
THR     <- 0.5; THR_HIGH <- 0.85
labels  <- dat$labels; p <- ncol(dat$Y); n <- nrow(dat$Y)

## ---- Truth: the fitted network ------------------------------------------------------
A0 <- fit0$AEst * (fit0$GammaEst > THR); diag(A0) <- 0
B0 <- fit0$BEst; C0 <- fit0$CEst; S0 <- fit0$SigmaEst
stopifnot(identical(dim(A0), c(p, p)), identical(dim(B0), c(p, ncol(dat$X))), identical(dim(C0), c(p, ncol(dat$U))),
          identical(dim(S0), c(p, p)))
S0 <- (S0 + t(S0)) / 2; R0 <- chol(S0)
M0inv <- solve(diag(p) - A0)
truth <- (A0 != 0) * 1; off <- row(A0) != col(A0)
cyc0  <- which(truth == 1 & t(truth) == 1 & upper.tri(truth), arr.ind = TRUE)          # true 2-cycles [to, from]
message(sprintf("Truth: %d edges (fitted network, PIP > %.1f), %d 2-cycles; spectral radius %.2f; n = %d, k = %d, l = %d; %d replicates on %d cores",
                sum(truth), THR, nrow(cyc0), max(Mod(eigen(A0)$values)), n, ncol(dat$X), ncol(dat$U), N_REP, N_CORES))

## ---- One replicate --------------------------------------------------------------------
boot_fit <- function(rep) {
  set.seed(20261103 + rep)
  E <- matrix(rnorm(n * p), n, p) %*% R0
  Yraw <- (dat$X %*% t(B0) + dat$U %*% t(C0) + E) %*% t(M0inv)
  s <- apply(Yraw, 2, sd); Y <- scale(Yraw); attr(Y, "scaled:center") <- NULL; attr(Y, "scaled:scale") <- NULL
  colnames(Y) <- colnames(dat$Y)
  At <- A0 * outer(1 / s, s)                                              # truth on the scale of the standardized Y
  S <- summary_stats(Y, dat$X, dat$U); w <- character(0)
  fit <- withCallingHandlers(
    RGM(Syy = S$Syy, Syx = S$Syx, Sxx = S$Sxx, Syu = S$Syu, Sxu = S$Sxu, Suu = S$Suu, D = dat$D, n = S$n, selectB = FALSE,
        nIter = MCMC$nIter, nBurnin = MCMC$nBurnin, Thin = MCMC$Thin, prior = "Spike and Slab", SigmaStarModel = "SSSL"),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  G <- fit$GammaEst; Ah <- fit$AEst * (G > THR); sel <- off & G > THR; selh <- off & G > THR_HIGH
  cyc <- sel & t(sel) & upper.tri(G)
  spur <- cyc & !(truth == 1 & t(truth) == 1)
  list(G = G, sign = sign(Ah),
       stats = c(auc = as.numeric(auc(roc(truth[off], G[off], levels = c(0, 1), direction = "<", quiet = TRUE))),
                 n_sel = sum(sel), tpr = sum(sel & truth == 1) / sum(truth),
                 fdp = if (sum(sel)) mean(truth[sel] == 0) else NA, fdp_high = if (sum(selh)) mean(truth[selh] == 0) else NA,
                 mad_edges = mean(abs(Ah[truth == 1] - At[truth == 1])),
                 spurious_cycles = sum(spur), spurious_cycles_opposite = sum(spur & sign(Ah) != sign(t(Ah))),
                 warn_weakIV = as.numeric(any(grepl("Weak instruments", w))), warn_errvar = as.numeric(any(grepl("error variances exceed", w))),
                 warn_spectral = as.numeric(any(grepl("spectral radius", w))), warn_halves = as.numeric(any(grepl("two halves", w))),
                 warn_accept = as.numeric(any(grepl("acceptance rate", w))), warn_any = as.numeric(length(w) > 0)),
       cycles = if (nrow(cyc0)) setNames(as.numeric(sel[cyc0] & t(sel)[cyc0]), apply(cyc0, 1, paste, collapse = "_")) else numeric(0))
}

t0 <- Sys.time()
rl <- mclapply(seq_len(N_REP), function(r) tryCatch(boot_fit(r), error = function(e) structure(conditionMessage(e), class = "try-error")),
               mc.cores = N_CORES, mc.preschedule = FALSE)
bad <- vapply(rl, function(x) inherits(x, "try-error") || is.null(x), logical(1))
if (any(bad)) warning(sprintf("%d of %d replicates failed: %s", sum(bad), N_REP, paste(unique(unlist(rl[bad])), collapse = "; ")))
rl <- rl[!bad]; R <- length(rl)
message(sprintf("%d replicates in %.1f minutes", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))

## ---- Summaries --------------------------------------------------------------------------
st <- do.call(rbind, lapply(rl, `[[`, "stats"))
cy <- if (nrow(cyc0)) do.call(rbind, lapply(rl, `[[`, "cycles")) else NULL
Gs <- simplify2array(lapply(rl, `[[`, "G")); Sg <- simplify2array(lapply(rl, `[[`, "sign"))
ms <- function(x) sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
cyc_lab <- if (nrow(cyc0)) sprintf("%s $\\leftrightarrows$ %s", labels[cyc0[, 2]], labels[cyc0[, 1]]) else character(0)
summ <- data.frame(
  quantity = c("Network AUC", sprintf("Edges selected (PIP $> %.1f$)", THR), sprintf("TPR (PIP $> %.1f$)", THR),
               sprintf("False discovery proportion (PIP $> %.1f$)", THR), sprintf("False discovery proportion (PIP $> %.2f$)", THR_HIGH),
               "MAD of the effects on the true edges",
               if (nrow(cyc0)) sprintf("2-cycle %s recovered", cyc_lab),
               "Selected 2-cycles that are not true 2-cycles", "\\quad of which with effects of opposite signs",
               "Warning: weak instruments (first-stage $F < 10$)", "Warning: error variances", "Warning: spectral radius",
               "Warning: chain halves", "Warning: acceptance rate", "Any warning"),
  value = c(ms(st[, "auc"]), ms(st[, "n_sel"]), ms(st[, "tpr"]), ms(st[, "fdp"]), ms(st[, "fdp_high"]), ms(st[, "mad_edges"]),
            if (nrow(cyc0)) sprintf("%.0f\\%%", 100 * colMeans(cy)),
            ms(st[, "spurious_cycles"]), ms(st[, "spurious_cycles_opposite"]),
            vapply(c("warn_weakIV", "warn_errvar", "warn_spectral", "warn_halves", "warn_accept", "warn_any"),
                   function(v) sprintf("%.0f\\%%", 100 * mean(st[, v])), character(1))),
  stringsAsFactors = FALSE)
if (any(is.na(st[, "fdp_high"]))) summ$value[5] <- paste0(summ$value[5], sprintf(" [%d]", sum(!is.na(st[, "fdp_high"]))))

## per edge: the edges of the fitted network, then false edges selected in at least 20% of the replicates
sel_rate <- apply(Gs > THR, c(1, 2), mean); pip_mean <- apply(Gs, c(1, 2), mean); pip_sd <- apply(Gs, c(1, 2), sd)
sign_ok  <- apply(sweep(Sg, c(1, 2), sign(A0), "==") & Gs > THR, c(1, 2), mean)
idx_true  <- which(truth == 1, arr.ind = TRUE); idx_true <- idx_true[order(-fit0$GammaEst[idx_true]), , drop = FALSE]
idx_false <- which(truth == 0 & off & sel_rate >= 0.2, arr.ind = TRUE); idx_false <- idx_false[order(-sel_rate[idx_false]), , drop = FALSE]
edge_row <- function(ix, true_edge) data.frame(
  edge = sprintf("%s $\\rightarrow$ %s", labels[ix[, 2]], labels[ix[, 1]]), in_fit = true_edge,
  fitted_pip = round(fit0$GammaEst[ix], 3), fitted_effect = round(fit0$AEst[ix], 3),
  boot_pip = sprintf("%.2f (%.2f)", pip_mean[ix], pip_sd[ix]), selected = round(100 * sel_rate[ix]),
  selected_sign = if (true_edge) round(100 * sign_ok[ix]) else NA, stringsAsFactors = FALSE)
edges <- rbind(edge_row(idx_true, TRUE), if (nrow(idx_false)) edge_row(idx_false, FALSE))

## ---- Save -------------------------------------------------------------------------------
dir.create(file.path(OUT_DIR, "tables"), showWarnings = FALSE)
saveRDS(list(stats = st, cycles = cy, pip = Gs, sign = Sg, truth = list(A = A0, B = B0, C = C0, Sigma = S0), labels = labels,
             n_rep = R, failed = sum(bad), mcmc = MCMC, threshold = c(THR, THR_HIGH),
             package_version = as.character(packageVersion("MR.RGM")), session_info = sessionInfo()),
        file.path(OUT_DIR, "output", "OneK1K_bootstrap.rds"))
write.csv(summ, file.path(OUT_DIR, "tables", "onek_bootstrap.csv"), row.names = FALSE)
write.csv(edges, file.path(OUT_DIR, "tables", "onek_bootstrap_edges.csv"), row.names = FALSE)

cap_common <- sprintf("the fitted \\texttt{MR.RGM} network ($%d$ edges with PIP $> %.1f$, posterior-mean effects) taken as the truth, expression generated %d times for the %d donors from their genotypes and covariates, and \\texttt{MR.RGM} refitted with the simulation-study MCMC settings", sum(truth), THR, R, n)
L1 <- c("\\begin{table}[htb]", "\\centering", "\\small",
        sprintf("\\caption{OneK$1$K B cells, parametric bootstrap: %s. Mean (SD) over the replicates%s.}", cap_common,
                if (any(is.na(st[, "fdp_high"]))) sprintf(". The false discovery proportion at PIP $> %.2f$ is averaged over the replicates that select at least one such edge (their number in brackets)", THR_HIGH) else ""),
        "\\label{tab:onek_bootstrap}", "\\begin{tabular}{|l|c|}", "\\hline", "\\textbf{Quantity} & \\textbf{Value} \\\\", "\\hline",
        sprintf("%s & %s \\\\", summ$quantity, summ$value), "\\hline", "\\end{tabular}", "\\end{table}")
writeLines(display_text(L1), file.path(OUT_DIR, "tables", "onek_bootstrap.tex"))   # displayed names (real_data_functions.R)
fmt_row <- function(i) with(edges[i, ], sprintf("%s & %s & %.2f & %+.3f & %s & %d & %s \\\\", edge, if (in_fit) "yes" else "no",
                                               fitted_pip, fitted_effect, boot_pip, selected, if (is.na(selected_sign)) "--" else as.character(selected_sign)))
L2 <- c("\\begin{table}[htb]", "\\centering", "\\footnotesize",
        sprintf("\\caption{OneK$1$K B cells, parametric bootstrap (%s): for every edge of the fitted network, its PIP and effect in the fit, its mean (SD) PIP over the replicates, and the percentage of replicates in which it is selected (PIP $> %.1f$) and selected with the sign of the fitted effect. Below are the false edges selected in at least $20\\%%$ of the replicates.}", cap_common, THR),
        "\\label{tab:onek_bootstrap_edges}", "\\begin{tabular}{|l|c|c|c|c|c|c|}", "\\hline",
        "\\textbf{Edge} & \\textbf{In fit} & \\textbf{PIP} & \\textbf{Effect} & \\textbf{Bootstrap PIP} & \\textbf{Selected \\%} & \\textbf{Same sign \\%} \\\\", "\\hline",
        vapply(seq_len(nrow(edges)), fmt_row, character(1)), "\\hline", "\\end{tabular}", "\\end{table}")
writeLines(display_text(L2), file.path(OUT_DIR, "tables", "onek_bootstrap_edges.tex"))
print(summ, row.names = FALSE); print(edges, row.names = FALSE)
message("Saved output/OneK1K_bootstrap.rds and tables/onek_bootstrap{,_edges}.{csv,tex} in ", OUT_DIR)
