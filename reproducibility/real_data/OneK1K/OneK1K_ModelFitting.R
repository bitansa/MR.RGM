###############################################################################
## OneK1K_ModelFitting.R  —  OneK1K B cells, MR.RGM (main) / MR.RGM+ (sensitivity) with covariates
##
## Prerequisite: OneK1K_PreProcessing.R (same R session)
##
## Steps
##   1) Candidate genes: KEGG pathway hsa04662, B cell receptor signaling (OneK1K_nodes.R)
##   2) Instruments: for each gene, the cis-SNPs (transcription start site +/- CIS_KB kb, default 1000) with
##      F > F_MIN (10) in log-expression ~ SNP + sex + age (decade bins) + PC1 + PC2 + PC3 (the covariates of the
##      model), LD-clumped in R (greedy by |t|: a SNP is kept only if its r^2 with every SNP already kept for that
##      gene is <= LD_R2, default 0.5), at most TOP_PER_GENE (10) per gene; a SNP selected for two genes is assigned
##      to the gene where its |t| is larger; instruments of different genes in LD (r^2 > LD_R2) are resolved in the
##      same way (the SNP with the smaller |t| is dropped), so that the instruments of different genes are not
##      near-collinear; genes without any instrument are removed. For every selected SNP the t-statistics on all
##      other analysed genes are recorded (cross-gene associations, reported with the instruments).
##   3) Preprocessing: log1p expression (from preprocessing), genotypes and covariates [SEX (male = 1),
##      AGE (decade bins 1-7), PC1-PC3]; every column centered and scaled
##   4) Fit from summary statistics: MR.RGM (default; each SNP affects only its own gene) or,
##      with ONEK_MODEL <- "MR.RGM+", MR.RGM+ (every SNP may affect every gene, D = 1); the package's warnings and
##      diagnostics are recorded; four chains with different seeds (same starting point) are checked for convergence
##      (split-R-hat, effective sample size, between-chain PIPs) and pooled
##   5) Save the fit and everything needed to interpret it to <ONEK_OUT_DIR>/output/OneK1K_fit_MRRGM.rds
##      (or OneK1K_fit_MRRGMplus.rds), and the analysed data to output/OneK1K_analysis_data.rds (used by
##      OneK1K_Competitors.R); ONEK_OUT_DIR defaults to this folder
##   ONEK_STRENGTH_ONLY <- TRUE: steps 1-3 only; the instrument strengths of every candidate gene that passes the
##      Bonferroni check are saved to output/OneK1K_candidates.rds, the design of ../../package_check/check_design_matched.R,
##      which gives the number of genes of the analysis (ONEK_N_GENES)
##
## Usage (this folder):
##   data_dir <- "path/to/zenodo/folder"; source("OneK1K_PreProcessing.R")
##   source("OneK1K_ModelFitting.R")      # slow: many SNP-wise regressions + MCMC
##   source("OneK1K_Results.R")
###############################################################################

suppressPackageStartupMessages({ library(MR.RGM); library(GenomicRanges); library(dplyr) })
source("../real_data_functions.R")

needed <- c("donor", "gene.ref", "vcf.ref", "gene.promoter.ref", "RNA.log", "genotype.mat", "Subsetted_SNP_Names")
missing <- needed[!vapply(needed, exists, logical(1))]
if (length(missing)) stop("Run OneK1K_PreProcessing.R first. Missing: ", paste(missing, collapse = ", "))

## Settings (reported in the manuscript)
TOP_PER_GENE <- if (exists("ONEK_TOP_PER_GENE")) ONEK_TOP_PER_GENE else 10   # instruments per gene (override: ONEK_TOP_PER_GENE)
F_MIN        <- if (exists("ONEK_F_MIN")) ONEK_F_MIN else 10                 # instrument strength F = t^2 (override: ONEK_F_MIN)
LD_R2        <- if (exists("ONEK_LD_R2")) ONEK_LD_R2 else 0.5                 # LD clumping: max r^2 between kept SNPs (override: ONEK_LD_R2)
## Chain length per model (four chains each, 8,000 kept draws per chain): MR.RGM 100,000 iterations, 20,000 burn-in,
## thinning 10 (longer than in the simulation study because the real-data fit is reported in detail); MR.RGM+, which
## estimates every SNP-gene effect (about ten times as many parameters), five times longer (500,000 / 100,000 / 50),
## so that a lack of convergence cannot be attributed to the length of the chains
MCMC_BY_MODEL <- list(`MR.RGM`  = list(nIter = 100000, nBurnin = 20000,  Thin = 10),
                      `MR.RGM+` = list(nIter = 500000, nBurnin = 100000, Thin = 50))
SEEDS  <- c(12345, 54321, 1001, 2002)                     # four chains (convergence diagnostics; pooled results)
SEED   <- SEEDS[1]
OUT_DIR <- if (exists("ONEK_OUT_DIR")) ONEK_OUT_DIR else "."  # outputs: <OUT_DIR>/output (fits), tables, figures

## ---- 1) Genes ------------------------------------------------------------------
source("OneK1K_nodes.R")                                    # KEGG hsa04662 node -> HGNC symbol(s)
node_of  <- setNames(rep(names(nodes), lengths(nodes)), unlist(nodes))
cand     <- names(node_of)[names(node_of) %in% gene.ref$gene_name]
message("Candidate genes not in the data: ", paste(setdiff(names(node_of), cand), collapse = ", "))

## ---- 2) Instruments: strongest cis-SNPs per gene (F > F_MIN) ---------------------
## For each gene, every SNP within +/- CIS_KB kb of the promoter is tested in
##   log-expression ~ SNP + age + sex + PC1 + PC2 + PC3;
## SNPs with F = t^2 > F_MIN (the usual MR weak-instrument rule) are candidate instruments; they are LD-clumped
## (greedy by |t|, r^2 <= LD_R2 with every SNP already kept) and the TOP_PER_GENE strongest kept SNPs are the
## instruments of the gene; genes without any are removed.
stopifnot(all(donor$sex %in% c("male", "female")))           # SEX is coded male = 1 below
age_bin <- cut(donor$age, breaks = c(-Inf, 30, 40, 50, 60, 70, 80, Inf), right = FALSE, labels = FALSE)
cov_df <- data.frame(age = age_bin, sex = as.numeric(donor$sex == "male"), PC1 = donor$PC1, PC2 = donor$PC2, PC3 = donor$PC3)
top_snp_list <- list(); iv_table <- list(); tval_list <- list()
for (g in cand) {
  gi <- match(g, gene.ref$gene_name)
  snp_idx <- which(countOverlaps(vcf.ref, gene.promoter.ref[gi]) > 0)
  y <- RNA.log[gi, ]
  if (!length(snp_idx) || sd(y) == 0) next
  tval <- vapply(snp_idx, function(s) {
    x <- genotype.mat[s, ]
    if (sd(x) == 0) return(0)
    summary(lm(y ~ x + age + sex + PC1 + PC2 + PC3, data = cov_df))$coefficients[2, 3]
  }, numeric(1))
  ok <- which(tval^2 > F_MIN)
  ## LD clumping among the candidate instruments of the gene: strongest first, keep a SNP only if its squared
  ## correlation with every SNP already kept is <= LD_R2 (one representative of each block of near-identical SNPs)
  ok <- ok[order(-abs(tval[ok]))]
  kept <- integer(0)
  for (s in ok) {
    if (length(kept) >= TOP_PER_GENE) break
    x <- genotype.mat[snp_idx[s], ]
    r2max <- if (length(kept)) max(cor(x, t(genotype.mat[snp_idx[kept], , drop = FALSE]))^2) else 0
    if (r2max <= LD_R2) kept <- c(kept, s)
  }
  iv_table[[g]] <- data.frame(gene = g, node = node_of[[g]], cis_snps = length(snp_idx), strong = length(ok),
                              clumped = length(kept), best_t = max(abs(tval)))
  if (!length(kept)) next
  top_snp_list[[g]] <- vcf.ref$ID[snp_idx[kept]]
  tval_list[[g]] <- setNames(tval[kept], vcf.ref$ID[snp_idx[kept]])
}
iv_table <- do.call(rbind, iv_table)
## A SNP in the overlapping cis windows of two genes can be selected for both: keep it for the gene where |t| is larger
all_sel <- data.frame(gene = rep(names(top_snp_list), lengths(top_snp_list)), snp = unlist(top_snp_list),
                      abs_t = abs(unlist(tval_list)), stringsAsFactors = FALSE, row.names = NULL)
dup <- all_sel$snp[duplicated(all_sel$snp)]
if (length(dup)) {
  message(sprintf("%d SNP(s) selected for more than one gene; each is kept for the gene with the larger |t|: %s", length(unique(dup)), paste(unique(dup), collapse = ", ")))
  all_sel <- all_sel[order(all_sel$snp, -all_sel$abs_t), ]
  all_sel <- all_sel[!duplicated(all_sel$snp), ]
  top_snp_list <- split(all_sel$snp, all_sel$gene)[names(top_snp_list)]
  top_snp_list <- top_snp_list[lengths(top_snp_list) > 0]
}
stopifnot(!anyDuplicated(unlist(top_snp_list)))
## Instruments of different genes in LD: cis windows of neighbouring genes overlap, so an instrument of one gene can
## be correlated with an instrument of another. Pairs with r^2 > LD_R2 are resolved by dropping the SNP with the
## smaller |t| (greedy, strongest first), so that the instruments of different genes are not near-collinear (the
## exclusion restrictions of the known map then concern distinct instruments and Sxx is well conditioned).
all_sel <- all_sel[order(-all_sel$abs_t), ]
G_sel <- t(genotype.mat[match(all_sel$snp, Subsetted_SNP_Names), , drop = FALSE])     # donors x selected SNPs
R2_sel <- cor(G_sel)^2; diag(R2_sel) <- 0
drop <- logical(nrow(all_sel))
for (i in seq_len(nrow(all_sel))) {
  if (drop[i]) next
  later <- which(seq_len(nrow(all_sel)) > i & !drop & all_sel$gene != all_sel$gene[i] & R2_sel[i, ] > LD_R2)
  drop[later] <- TRUE
}
if (any(drop)) {
  message(sprintf("%d instrument(s) dropped because of LD (r^2 > %g) with a stronger instrument of another gene: %s",
                  sum(drop), LD_R2, paste(all_sel$snp[drop], collapse = ", ")))
  all_sel <- all_sel[!drop, ]
  top_snp_list <- split(all_sel$snp, all_sel$gene)[names(top_snp_list)]
  top_snp_list <- top_snp_list[lengths(top_snp_list) > 0]
}
ld_between_genes_max <- if (nrow(all_sel) > 1) { R2 <- cor(t(genotype.mat[match(all_sel$snp, Subsetted_SNP_Names), , drop = FALSE]))^2
  same <- outer(all_sel$gene, all_sel$gene, "=="); max(R2[!same]) } else 0
message(sprintf("Largest r^2 between instruments of different genes after clumping: %.3f", ld_between_genes_max))
symbols  <- names(top_snp_list)                              # genes with >= 1 instrument
labels   <- symbols                                          # figures and tables label genes by symbol
removed_no_iv <- setdiff(cand, symbols)
message(sprintf("Genes kept: %d of %d candidates; removed (no SNP with F > %g): %s",
                length(symbols), length(cand), F_MIN, paste(removed_no_iv, collapse = ", ")))
gene_idx <- match(symbols, gene.ref$gene_name)
snps <- unique(unlist(top_snp_list))
if (!length(snps)) stop("No SNPs selected.")
snp_idx2 <- match(snps, Subsetted_SNP_Names)
if (anyNA(snp_idx2)) stop("Selected SNPs not found: ", paste(head(snps[is.na(snp_idx2)], 20), collapse = ", "))

## ---- 3) Matrices and preprocessing -----------------------------------------------
Y_log <- t(RNA.log[gene_idx, , drop = FALSE]);        colnames(Y_log) <- symbols
X_raw <- t(genotype.mat[snp_idx2, , drop = FALSE]);   colnames(X_raw) <- snps
U_raw <- cbind(SEX = cov_df$sex, AGE = cov_df$age, PC1 = donor$PC1, PC2 = donor$PC2, PC3 = donor$PC3)
if (anyNA(U_raw)) stop("Missing covariate values.")

Y <- standardize_cols(Y_log, "Expression")
X <- standardize_cols(X_raw, "Genotypes")
U <- standardize_cols(U_raw, "Covariates")
if (ncol(Y) != length(symbols)) stop("A gene had zero variance and was dropped; labels would no longer match.")

## ---- 3b) Instrument strength ------------------------------------------------------
## Concentration parameter mu2 = n * R2 / (1 - R2), R2 = partial R^2 of a gene's own SNPs given
## the covariates; the N_GENES genes with the largest mu2 are kept (optional floor MIN_MU2, default 0).
## N_GENES is the largest gene set that the design-matched simulation supports (package_check/check_design_matched.R,
## run on the candidates written by ONEK_STRENGTH_ONLY <- TRUE, see below); default 10.
N_GENES <- if (exists("ONEK_N_GENES")) ONEK_N_GENES else 10   # genes kept (override: ONEK_N_GENES)
MIN_MU2 <- if (exists("ONEK_MIN_MU2")) ONEK_MIN_MU2 else 0    # optional floor (override: ONEK_MIN_MU2)
r2 <- vapply(symbols, function(g) {
  xs <- intersect(top_snp_list[[g]], colnames(X))
  if (!length(xs)) return(0)
  full <- lm(Y[, g] ~ X[, xs, drop = FALSE] + U); base <- lm(Y[, g] ~ U)
  1 - sum(resid(full)^2) / sum(resid(base)^2)
}, numeric(1))
mu2 <- nrow(Y) * r2 / (1 - r2)
## The SNPs are selected in the same data among all cis-SNPs of the gene, so a single SNP with F just above
## F_MIN can be a chance finding: best_p_bonf = Bonferroni-adjusted p-value of the best SNP over the
## cis-SNPs tested (< 0.05: the gene has a real cis-eQTL after accounting for the selection).
cis_n <- iv_table$cis_snps[match(symbols, iv_table$gene)]
best_p <- 2 * pt(-iv_table$best_t[match(symbols, iv_table$gene)], df = nrow(Y) - 7)
F_joint <- vapply(symbols, function(g) {                     # F statistic of the gene's own instruments given the covariates
  m <- length(intersect(top_snp_list[[g]], colnames(X))); if (!m) return(0)
  (r2[[g]] / m) / ((1 - r2[[g]]) / (nrow(Y) - m - ncol(U) - 1))
}, numeric(1))
iv_strength <- data.frame(gene = symbols, node = unname(node_of[symbols]),
                          instruments = lengths(top_snp_list[symbols]), partial_R2 = round(r2, 3), F = round(F_joint, 1), mu2 = round(mu2, 1),
                          cis_snps = cis_n, best_p_bonf = signif(pmin(1, best_p * cis_n), 2),
                          partial_R2_exact = r2)                # unrounded, the design of check_design_matched.R
print(iv_strength[order(-iv_strength$mu2), ], row.names = FALSE)
## A gene is kept only if its best cis-SNP is significant after Bonferroni correction over its cis-SNPs
## (otherwise its instrument may be a chance finding and its outgoing edges would not be identified).
BONF_ALPHA <- if (exists("ONEK_BONF_ALPHA")) ONEK_BONF_ALPHA else 0.05   # override: ONEK_BONF_ALPHA (1 = off)
chance <- if (BONF_ALPHA >= 1) character(0) else symbols[iv_strength$best_p_bonf >= BONF_ALPHA]
if (length(chance)) message(sprintf("Removed (best cis-SNP not significant after Bonferroni, p >= %g): %s", BONF_ALPHA, paste(chance, collapse = ", ")))
strongest <- setdiff(symbols, chance)[order(-mu2[setdiff(symbols, chance)])][seq_len(min(N_GENES, length(setdiff(symbols, chance))))]
weak <- symbols[!(symbols %in% strongest) | mu2 < MIN_MU2]          # pathway order kept for the rest
if (length(weak)) {
  message(sprintf("Removed (not among the %d strongest genes, or mu2 < %g): %s", N_GENES, MIN_MU2, paste(weak, collapse = ", ")))
  keep_g <- !symbols %in% weak
  symbols <- symbols[keep_g]; labels <- labels[keep_g]
  top_snp_list <- top_snp_list[symbols]
  Y <- Y[, symbols, drop = FALSE]
  X <- X[, intersect(colnames(X), unlist(top_snp_list)), drop = FALSE]
}
removed_weak_iv <- weak
if (length(symbols) < 2) stop("Fewer than 2 genes left after the instrument-strength step.")
## Cross-gene associations of the selected SNPs: t-statistic of each instrument on every analysed gene given the
## covariates (marginal regressions; own-gene entries are the selection statistics). These are reported with the
## instruments. A large cross-gene |t| is not by itself evidence of pleiotropy, because an instrument's association
## with another gene can be mediated by its own gene; MR.RGM+ is the analysis that allows direct cross-gene effects.
cross_t <- sapply(colnames(X), function(s) vapply(symbols, function(g) summary(lm(Y[, g] ~ X[, s] + U))$coefficients[2, 3], numeric(1)))
cross_t <- matrix(cross_t, nrow = length(symbols), dimnames = list(symbols, colnames(X)))
own <- vapply(colnames(X), function(s) names(which(vapply(top_snp_list[symbols], function(v) s %in% v, logical(1))))[1], character(1))
off_t <- sapply(colnames(X), function(s) max(abs(cross_t[setdiff(symbols, own[[s]]), s])))
message(sprintf("Cross-gene check: %d of %d instruments have |t| > %.1f on a gene other than their own (largest %.1f)", sum(off_t^2 > F_MIN), ncol(X), sqrt(F_MIN), max(off_t)))
S <- summary_stats(Y, X, U)
message(sprintf("OneK1K: n = %d donors, p = %d genes, k = %d SNPs, l = %d covariates", S$n, ncol(Y), ncol(X), ncol(U)))
## ONEK_STRENGTH_ONLY <- TRUE: no model is fitted; the instrument strengths of every candidate gene that passes the
## Bonferroni check (strongest first) are saved to <OUT_DIR>/output/OneK1K_candidates.rds, the design of
## package_check/check_design_matched.R, which gives the largest supported gene set (ONEK_N_GENES).
if (exists("ONEK_STRENGTH_ONLY") && isTRUE(ONEK_STRENGTH_ONLY)) {
  dir.create(file.path(OUT_DIR, "output"), showWarnings = FALSE, recursive = TRUE)
  cand_ok <- if (BONF_ALPHA >= 1) iv_strength else iv_strength[iv_strength$best_p_bonf < BONF_ALPHA, ]
  cand_ok <- cand_ok[order(-cand_ok$mu2), ]
  cf <- file.path(OUT_DIR, "output", "OneK1K_candidates.rds")
  saveRDS(list(n = nrow(Y), instrument_strength = cand_ok, symbols = cand_ok$gene, all_strength = iv_strength,
               settings = list(top_per_gene = TOP_PER_GENE, f_min = F_MIN, ld_r2 = LD_R2, bonf_alpha = BONF_ALPHA, cis_window = sprintf("+/- %g kb", CIS_KB)),
               package_version = as.character(packageVersion("MR.RGM"))), cf)
  message(sprintf("ONEK_STRENGTH_ONLY: %d candidate genes pass the Bonferroni check; saved %s (no model fitted).", nrow(cand_ok), cf))
  ## scale of the analysed columns before standardization (GrIVET, whose tuning is not scale free, is fitted to the
  ## unstandardized data in OneK1K_Competitors.R)
  sc_file <- file.path(OUT_DIR, "output", "OneK1K_scales.rds")
  saveRDS(list(sd_Y = apply(Y_log[, colnames(Y), drop = FALSE], 2, sd), sd_X = apply(X_raw[, colnames(X), drop = FALSE], 2, sd),
               symbols = colnames(Y), snps = colnames(X)), sc_file)
  message("Saved ", sc_file)
} else {

## ---- 4) Fit -----------------------------------------------------------------------
## MODEL "MR.RGM" (main analysis): each instrument may affect only the gene it was selected for
##   (exclusion restriction; cis-eQTLs act on other genes through their own gene).
## MODEL "MR.RGM+" (sensitivity analysis): every instrument may affect every gene (D = 1).
MODEL <- if (exists("ONEK_MODEL")) ONEK_MODEL else "MR.RGM"
stopifnot(MODEL %in% c("MR.RGM", "MR.RGM+"))
MCMC  <- MCMC_BY_MODEL[[MODEL]]
message(sprintf("Model: %s; %d chains of %s iterations (%s burn-in, thinning %d)", MODEL, length(SEEDS),
                format(MCMC$nIter, big.mark = ",", scientific = FALSE), format(MCMC$nBurnin, big.mark = ",", scientific = FALSE), MCMC$Thin))
D_own <- matrix(0, ncol(Y), ncol(X), dimnames = list(colnames(Y), colnames(X)))
for (g in intersect(names(top_snp_list), colnames(Y))) D_own[g, intersect(top_snp_list[[g]], colnames(X))] <- 1
D <- if (MODEL == "MR.RGM") D_own else matrix(1, ncol(Y), ncol(X))
## The package's warnings (weak instruments, identification, spectral radius, chain halves, acceptance, error-variance
## ratio) are recorded with the fit instead of printed; fit$diagnostics holds the corresponding statistics.
run_chain <- function(seed) {
  set.seed(seed)
  w <- character(0)
  f <- withCallingHandlers(
    RGM(Syy = S$Syy, Syx = S$Syx, Sxx = S$Sxx, Syu = S$Syu, Sxu = S$Sxu, Suu = S$Suu,
        D = D, n = S$n, selectB = (MODEL == "MR.RGM+"),
        nIter = MCMC$nIter, nBurnin = MCMC$nBurnin, Thin = MCMC$Thin,
        prior = "Spike and Slab", SigmaStarModel = "SSSL"),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  f$warnings <- w
  f
}
## Four chains with different seeds from the same starting point; pooled when converged (fit_chains() in
## ../real_data_functions.R: rank-normalized split-R-hat <= 1.01 for the log-likelihood and every edge indicator)
fc <- fit_chains(run_chain, SEEDS)
fit <- fc$fit; convergence <- fc$convergence
if (length(fit$warnings)) message("Package warnings (recorded in the fit, any chain):\n  ", paste(fit$warnings, collapse = "\n  "))
rm(fc)

## ---- 5) Save -------------------------------------------------------------------------
dir.create(file.path(OUT_DIR, "output"), showWarnings = FALSE, recursive = TRUE)
OneK1K_fit <- list(
  fit = fit, convergence = convergence, labels = labels, symbols = symbols, pathway_node = unname(node_of[symbols]),
  removed_no_instrument = removed_no_iv, removed_weak_instrument = removed_weak_iv,
  instrument_strength = iv_strength, instrument_table = iv_table, snps = colnames(X), top_snps = top_snp_list,
  cross_gene_t = cross_t, ld_between_genes_max = ld_between_genes_max, kegg_pathway = KEGG_PATHWAY, kegg_accessed = KEGG_ACCESSED,
  n = S$n, p = ncol(Y), k = ncol(X), l = ncol(U),
  settings = list(top_per_gene = TOP_PER_GENE, f_min = F_MIN, ld_r2 = LD_R2, n_genes = N_GENES, min_mu2 = MIN_MU2, bonf_alpha = BONF_ALPHA, mcmc = MCMC, seed = SEED, seeds = SEEDS, chains = length(SEEDS), cis_window = sprintf("+/- %g kb", CIS_KB), genotype_source = GENOTYPES,
                  expression = "library-size normalised, log1p, centered and scaled",
                  genotypes  = "array dosages (0/1/2; > 5% missing dropped, rest mean-imputed; minor-allele carriers > 5%), centered and scaled",
                  covariates = "SEX (male = 1), AGE (decade bins 1-7), PC1-PC3; centered and scaled",
                  model = MODEL, D = if (MODEL == "MR.RGM") "own-gene instrument map (no selection of the instrument effects)" else "all ones (instrument effects selected)",
                  selectB = fit$settings$selectB, prior = "Spike and Slab", SigmaStarModel = "SSSL"),
  package_version = as.character(packageVersion("MR.RGM")),
  session_info = sessionInfo())
out_file <- file.path(OUT_DIR, "output", if (MODEL == "MR.RGM") "OneK1K_fit_MRRGM.rds" else "OneK1K_fit_MRRGMplus.rds")
saveRDS(OneK1K_fit, out_file)
message("Saved ", out_file)
## The analysed data (standardized expression Y, genotypes X, covariates U) and the own-gene instrument map, for the
## competing methods (OneK1K_Competitors.R) on exactly the same genes, SNPs and donors, with the settings (chain
## length, seeds) of the main analysis; written by the MR.RGM fit
if (MODEL == "MR.RGM")
  saveRDS(list(Y = Y, X = X, U = U, D = D_own, labels = labels, symbols = symbols, n = S$n, settings = OneK1K_fit$settings,
               sd_Y = apply(Y_log[, colnames(Y), drop = FALSE], 2, sd), sd_X = apply(X_raw[, colnames(X), drop = FALSE], 2, sd)),
          file.path(OUT_DIR, "output", "OneK1K_analysis_data.rds"))   # sd_Y, sd_X: scale before standardization (GrIVET)
}
