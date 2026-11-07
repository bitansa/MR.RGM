###############################################################################
## check_design_matched.R  --  can MR.RGM recover a network with the instrument strengths observed in the real data?
##
## Giving every gene the same instrument strength is unrealistic: in the real data the strength varies
## a lot between genes, so here each simulated gene gets the partial R^2 and the number of instruments of the
## corresponding real gene, at the real n. The design is read from the candidate genes of the real data
## (OneK1K_candidates.rds, written by OneK1K_ModelFitting.R with ONEK_STRENGTH_ONLY <- TRUE: every gene whose best
## cis-SNP passes the Bonferroni check, with its partial R^2 and number of instruments), or from a saved fit.
## Gene sets of increasing size (strongest genes first, by mu^2) are compared, and the largest supported gene set is
## the largest size such that it and every smaller size tested meet every criterion below (pre-specified; the size
## used for the real-data fit, ONEK_N_GENES):
##   auc_mean >= 0.80        the network is recovered (AUC of the causal PIPs)
##   halves_mean >= 0.90     the chain is stable (correlation of the causal PIPs from the two halves of the chain)
##   warn_errvar_mean <= 0.10  at most 10% of the fits trigger the package's error-variance (weak-identification) warning
## Truth: random graph with up to 2p directed edges (duplicates and self-loops removed) of size +/-0.2, 3 confounders,
## every column standardized;
## MR.RGM (own-gene instrument map D, no selection of the instrument effects; selectB = FALSE).
##
## Reported per gene set: mean and SD over replicates (columns *_mean, *_sd):
##   sigma_med     median estimated error variance (truth ~ 1 - R^2 on the standardized data)
##   auc           AUC of the causal PIPs for the true graph
##   auc_strong    AUC restricted to edges out of genes with mu2 >= 50
##   auc_weak      AUC restricted to edges out of genes with mu2 <  50
##   rho_A         spectral radius of AEst
##   n_rep50       edges with PIP > 0.5 (the reporting rule of the real-data analysis); fdr50_strong / fdr50_weak:
##                 share of them that are false, for edges out of strong (mu2 >= 50) / weak genes (NA when no such
##                 edge is reported in a replicate; averaged over the replicates that report at least one)
##   n_rep85, fdr85_strong, fdr85_weak   the same at PIP > 0.85
##   recip_opposite  reported 2-cycles whose two effects have opposite signs (weak-identification artefact)
##   halves        correlation of causal PIPs from the two halves of the chain
##   warn_*        share of fits in which the package's weak-instrument / spectral-radius / chain-halves / error-variance
##                 warning fired
## Settings as in the simulation study: 100 replicates; MCMC 50,000 iterations, 10,000 burn-in, thinning 10.
## Usage: source("check_design_matched.R")   (about 3-5 hours on 8-10 cores; on a SLURM cluster:
##   sbatch -p <partition> check_design_matched.sbatch, minutes on 100 cores; replicates in forked processes with a
##   time limit, run_parallel() of ../simulation/sim_functions.R)
##   Options before sourcing: DESIGN_FITS (named character vector of candidate or fit files; default the OneK1K
##   candidates), DESIGN_SIZES (gene-set sizes; default 5, 8, 10, 12, 15, 20 and all candidates), DESIGN_N_REP, N_CORES.
###############################################################################

library(MR.RGM); library(pROC); library(parallel)
if (packageVersion("MR.RGM") < "0.1.3" || !"gibbsB" %in% names(formals(MR.RGM::RGM))) stop("MR.RGM >= 0.1.3 is required (arguments selectB and gibbsB of RGM())")
## run_parallel() (replicates in forked processes with a time limit), without loading the simulation settings
## and display_text() (displayed method names and wording of the LaTeX table)
sim_env <- local({ e <- new.env(); sys.source(file.path("..", "simulation", "sim_functions.R"), envir = e,
                   toplevel.env = e) ; e })
run_parallel <- sim_env$run_parallel; display_text <- sim_env$display_text
if (!exists("N_CORES")) N_CORES <- max(1, parallel::detectCores() - 1)
if (!exists("DESIGN_FITS")) DESIGN_FITS <- c(OneK1K = file.path(if (exists("ONEK_OUT_DIR")) ONEK_OUT_DIR else file.path("..", "real_data", "OneK1K"),
                                                   "output", "OneK1K_candidates.rds"))
N_REP <- if (exists("DESIGN_N_REP")) DESIGN_N_REP else 100   # as in the simulation study

sim_fit <- function(n, r2, ivs, seed, nIter = 50000, nBurnin = 10000, Thin = 10) {   # as in the simulation study
  set.seed(seed)
  p <- length(r2); A <- matrix(0, p, p)
  edges <- cbind(sample(p, 2 * p, TRUE), sample(p, 2 * p, TRUE)); edges <- unique(edges[edges[, 1] != edges[, 2], , drop = FALSE])
  A[edges] <- sample(c(-0.2, 0.2), nrow(edges), TRUE)          # A[j, i] != 0: edge i -> j
  k <- sum(ivs); B <- matrix(0, p, k); D <- B; start <- c(0, cumsum(ivs))
  for (j in 1:p) { cols <- (start[j] + 1):start[j + 1]; D[j, cols] <- 1; B[j, cols] <- sqrt(r2[j] / ivs[j]) }
  X <- matrix(rnorm(n * k), n, k)
  U <- matrix(rnorm(n * 3), n, 3); C <- matrix(sample(c(-0.3, 0, 0.3), p * 3, TRUE), p, 3)
  E <- matrix(rnorm(n * p), n, p) %*% diag(sqrt(1 - r2), p)
  Y <- (X %*% t(B) + U %*% t(C) + E) %*% t(solve(diag(p) - A))
  X <- scale(X); Y <- scale(Y)
  w <- character(0)
  fit <- withCallingHandlers(
    RGM(Syy = crossprod(Y) / n, Syx = crossprod(Y, X) / n, Sxx = crossprod(X) / n, D = D, n = n,
        nIter = nIter, nBurnin = nBurnin, Thin = Thin, prior = "Spike and Slab", SigmaStarModel = "SSSL", selectB = FALSE),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  G <- fit$GammaEst; GP <- fit$GammaPst; m <- dim(GP)[3]; h <- m %/% 2
  off <- row(A) != col(A); truth <- (A != 0) * 1
  strong_out <- matrix(n * r2 / (1 - r2) >= 50, p, p, byrow = TRUE)   # column i = exposure gene i
  auc_on <- function(sel) { t <- truth[sel]; if (length(unique(t)) < 2) return(NA)
    as.numeric(auc(roc(t, G[sel], levels = c(0, 1), direction = "<", quiet = TRUE))) }
  rep50 <- off & G > 0.50                                             # edges reported in the real-data analysis
  rep85 <- off & G > 0.85
  fdr_on <- function(sel) if (sum(sel)) mean(truth[sel] == 0) else NA
  recip <- sum(rep50 & t(rep50) & sign(fit$AEst) != sign(t(fit$AEst)) & upper.tri(G))   # opposite-sign 2-cycles
  c(sigma_med = median(diag(fit$SigmaEst)), rho_A = max(Mod(eigen(fit$AEst)$values)), auc = auc_on(off),
    auc_strong = auc_on(off & strong_out), auc_weak = auc_on(off & !strong_out),
    n_rep50 = sum(rep50), fdr50_strong = fdr_on(rep50 & strong_out), fdr50_weak = fdr_on(rep50 & !strong_out),
    n_rep85 = sum(rep85), fdr85_strong = fdr_on(rep85 & strong_out), fdr85_weak = fdr_on(rep85 & !strong_out),
    recip_opposite = recip,
    halves = cor(apply(GP[, , 1:h], c(1, 2), mean)[off], apply(GP[, , (h + 1):m], c(1, 2), mean)[off]),
    warn_weakIV = as.numeric(any(grepl("Weak instruments", w))), warn_spectral = as.numeric(any(grepl("spectral radius", w))),
    warn_halves = as.numeric(any(grepl("two halves", w))), warn_errvar = as.numeric(any(grepl("error variances exceed", w))))
}

## Observed instrument strengths of the candidate (or analysed) genes, strongest first
designs <- lapply(names(DESIGN_FITS), function(d) {
  f <- DESIGN_FITS[[d]]
  if (!file.exists(f)) stop("File not found: ", f, " (run OneK1K_ModelFitting.R with ONEK_STRENGTH_ONLY <- TRUE first, or set DESIGN_FITS)")
  obj <- readRDS(f); st <- obj$instrument_strength
  st <- st[st$gene %in% obj$symbols, ]; st <- st[order(-st$mu2), ]
  p <- nrow(st)
  list(n = obj$n, gene = st$gene, r2 = if (!is.null(st$partial_R2_exact)) st$partial_R2_exact else st$partial_R2, ivs = st$instruments,
       sizes = if (exists("DESIGN_SIZES")) DESIGN_SIZES[DESIGN_SIZES <= p] else sort(unique(c(c(5, 8, 10, 12, 15, 20)[c(5, 8, 10, 12, 15, 20) <= p], p))),
       source = basename(f), package_version = obj$package_version)
})
names(designs) <- names(DESIGN_FITS)
for (d in names(designs)) { ds <- designs[[d]]
  message(sprintf("%s (%s, fitted with MR.RGM %s): n = %d, %d genes; R2 %s; instruments %s", d, ds$source, ds$package_version, ds$n, length(ds$gene),
                  paste(sprintf("%.3f", ds$r2), collapse = " "), paste(ds$ivs, collapse = " "))) }

res <- do.call(rbind, lapply(names(designs), function(d) {
  ds <- designs[[d]]
  do.call(rbind, lapply(ds$sizes, function(p) {
    message(sprintf("[%s] %s: top %d genes (n = %d)", format(Sys.time(), "%H:%M:%S"), d, p, ds$n))
    rl <- run_parallel(seq_len(N_REP), function(rep) sim_fit(ds$n, ds$r2[1:p], ds$ivs[1:p], seed = 1000 * p + rep),
                       cores = N_CORES, timeout_min = 60, label = sprintf("design-matched %s p=%d", d, p))
    bad <- vapply(rl, function(x) inherits(x, "try-error") || is.null(x), logical(1))
    if (any(bad)) warning(sprintf("%d replicate(s) failed for %s, p = %d", sum(bad), d, p))
    r <- simplify2array(rl[!bad])
    mn <- apply(r, 1, mean, na.rm = TRUE); sdv <- apply(r, 1, sd, na.rm = TRUE)
    out <- data.frame(data = d, n = ds$n, p = p, last_gene = ds$gene[p],
                      min_mu2 = round(min(ds$n * ds$r2[1:p] / (1 - ds$r2[1:p]))), reps = sum(!bad),
                      t(round(setNames(mn, paste0(names(mn), "_mean")), 3)), t(round(setNames(sdv, paste0(names(sdv), "_sd")), 3)))
    print(out, row.names = FALSE); out
  }))
}))
print(res, row.names = FALSE)
## Largest supported gene set per data set (criteria in the header)
res$supported <- (res$auc_mean >= 0.80 & res$halves_mean >= 0.90 & res$warn_errvar_mean <= 0.10) %in% TRUE
for (d in unique(res$data)) { r <- res[res$data == d, ]; r <- r[order(r$p), ]
  ok <- cumprod(r$supported) == 1                               # this size and every smaller one supported
  m <- if (any(ok)) max(r$p[ok]) else NA
  message(sprintf("%s: largest supported gene set = %s (sizes tested: %s; supported: %s)", d, m,
                  paste(r$p, collapse = ", "), paste(r$p[r$supported], collapse = ", "))) }
attr(res, "designs") <- designs
write.csv(res, sprintf("check_design_matched_%s.csv", paste(names(designs), collapse = "_")), row.names = FALSE)
saveRDS(list(results = res, designs = designs, n_rep = N_REP, package_version = as.character(packageVersion("MR.RGM"))),
        sprintf("check_design_matched_%s.rds", paste(names(designs), collapse = "_")))

## LaTeX table of the results (Online Appendix): one table per data set, mean (SD) over the replicates
design_matched_tex <- function(r, file, data_label, n, n_rep) {
  r <- r[order(r$p), ]
  f2 <- function(m, s) sprintf("%.2f (%.2f)", m, s)
  rows <- sprintf("%d & %s & %d & %s & %s & %s & %s & %.0f & %.1f & %.2f / %.2f & %.2f / %.2f & %.1f & %s \\\\",
                  r$p, r$last_gene, r$min_mu2, f2(r$auc_mean, r$auc_sd), f2(r$halves_mean, r$halves_sd),
                  f2(r$sigma_med_mean, r$sigma_med_sd), f2(r$rho_A_mean, r$rho_A_sd), 100 * r$warn_errvar_mean,
                  r$n_rep50_mean, r$fdr50_strong_mean, r$fdr50_weak_mean, r$fdr85_strong_mean, r$fdr85_weak_mean,
                  r$recip_opposite_mean, ifelse(cumprod(r$supported) == 1, "yes", "no"))
  L <- c(sprintf("%% Generated by check_design_matched.R (%s).", basename(file)), "\\begin{table}[htb]", "\\centering", "\\small",
         sprintf("\\caption{Design-matched simulation for %s ($n=%d$; $%d$ replicates per gene set; \\texttt{MR.RGM} with the own-gene map and the simulation-study MCMC settings). Gene sets of increasing size $m$, strongest genes first; last gene added and the smallest $\\mu^2$ in the set; mean (SD) over replicates of the network AUC, the chain-half correlation of the causal PIPs, the median fitted error variance and the spectral radius of $\\widehat{\\mathbf{A}}$; percentage of fits with the error-variance warning; mean number of edges with PIP $>0.5$; false discovery proportion (FDP) at PIP $>0.5$ and at PIP $>0.85$ for edges out of genes with $\\mu^2\\ge50$ / $\\mu^2<50$; reported opposite-sign 2-cycles per fit. Supported: AUC $\\ge0.80$, chain halves $\\ge0.90$ and warning rate $\\le10\\%%$ for the set and every smaller set.}", data_label, n, n_rep),
         "\\label{tab:design_matched}", "\\begin{adjustbox}{max width=\\textwidth}", "\\begin{tabular}{|c|l|r|c|c|c|c|r|r|c|c|r|c|}", "\\hline",
         "$m$ & \\textbf{Last gene} & \\textbf{min $\\mu^2$} & \\textbf{AUC} & \\textbf{Halves} & \\textbf{Error var.} & \\textbf{Spectral radius} & \\textbf{Warn.\\ \\%} & \\textbf{Edges} & \\textbf{FDP $0.5$} & \\textbf{FDP $0.85$} & \\textbf{Opp.\\ 2-cycles} & \\textbf{Supported} \\\\",
         "\\hline", rows, "\\hline", "\\end{tabular}", "\\end{adjustbox}", "\\end{table}")
  writeLines(display_text(L), file)
}
for (d in unique(res$data))
  design_matched_tex(res[res$data == d, ], sprintf("design_matched_%s.tex", d), if (d == "OneK1K") "OneK$1$K" else d,
                     res$n[res$data == d][1], N_REP)
