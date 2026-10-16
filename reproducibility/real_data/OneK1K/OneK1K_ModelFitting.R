###############################################################################
## OneK1K_ModelFitting.R  —  OneK1K B cells, MR.RGM+ with covariates
##
## Prerequisite: OneK1K_PreProcessing.R (same R session)
##
## Steps
##   1) B-cell receptor pathway node labels and the gene symbols used in the data;
##      genes without sufficient cis-SNP coverage are removed (fixed list)
##   2) Instruments: for each gene, the 15 cis-SNPs (promoter +/- 200 kb) with the
##      smallest p-value in  log-expression ~ SNP + age + sex + PC1 + PC2 + PC3
##   3) Preprocessing (Version 4): log1p expression (from preprocessing), genotypes and
##      covariates [SEX (male = 1), AGE (decade bins 1-7), PC1-PC3]; every column
##      centered and scaled
##   4) Fit MR.RGM+ (all SNPs allowed for all genes, D = 1) from summary statistics
##   5) Save the fit and everything needed to interpret it to output/OneK1K_fit.rds
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
TOP_PER_GENE <- 15
MCMC <- list(nIter = 50000, nBurnin = 10000, Thin = 10)    # 4,000 posterior samples; same chain as GTEx and the simulations
SEED <- 12345

## ---- 1) Genes ------------------------------------------------------------------
gene_names_pathway <- c(
  "CD19", "BCR", "AG", "BCAP", "P85", "PI3K", "FGR2B", "SHIP", "DOK3", "PTEN",
  "SHP2", "LYN", "SHP1", "CD22", "CD45", "CBP/PAG", "CSK", "PIR-B", "BAM32",
  "PLCY2", "BLNK", "LAB", "GRB2", "SOS", "SYK", "CBL", "EZRIN", "CLATHRIN",
  "VAV", "RAC", "CDC42", "RASGRP", "RASGAP", "PKC", "AKT", "MTOR", "P70S6K",
  "GSK3", "CARMA1", "TAK1", "BCL10", "MALT1", "DOK1", "IKK", "NFKB", "IKB",
  "RAS", "RAP", "RIAM", "PYK2", "HS1", "MEKK", "MEK", "P38", "JNK", "NFAT",
  "CALCINEURIN", "CAM", "CAMK", "C-RAF", "MEK1/2", "ERK1/2", "CD40", "FOXO",
  "ETS1", "OCT2", "BFL1", "BCL-XL", "ELK1", "BCL6", "EGR1", "JUN", "ATF2",
  "CREB", "MEF2C", "PIP3", "RAPL")
gene_names_dataset <- c(
  "CD19", "BCR", "CD79A", "PIK3AP1", "PIK3R1", "PIK3CA", "FCGR2B", "INPP5D",
  "DOK3", "PTEN", "PTPN11", "LYN", "PTPN6", "CD22", "PTPRC", "PAG1", "CSK",
  "LILRB3", "DAPP1", "PLCG2", "BLNK", "LAT2", "GRB2", "SOS1", "SYK", "CBL",
  "EZR", "CLTB", "VAV1", "RAC1", "CDC42", "RASGRP1", "RASA1", "PRKCA", "AKT1",
  "MTOR", "RPS6KB1", "GSK3B", "CARD11", "MAP3K7", "BCL10", "MALT1", "DOK1",
  "IKBKB", "RELA", "NFKB2", "HRAS", "RAP1A", "APBB1IP", "PTK2B", "YWHAB",
  "MAP3K1", "MAP2K1", "MAPK14", "MAPK8", "NFATC1", "PPP3CA", "CALM1", "CAMK1",
  "RAF1", "MAP2K2", "MAPK1", "CD40", "FOXO1", "ETS1", "POU2F2", "BCL2A1",
  "BCL2L1", "ELK1", "BCL6", "EGR1", "JUN", "ATF2", "CREB1", "MEF2C", "PREX1",
  "RASSF3")
stopifnot(length(gene_names_pathway) == length(gene_names_dataset))
genes_to_remove <- c("CD79A", "DOK3", "PTEN", "PTPN11", "PTPN6", "CDC42", "MTOR", "PPP3CA", "FOXO1", "POU2F2", "ELK1")
keep    <- !(gene_names_dataset %in% genes_to_remove)
symbols <- gene_names_dataset[keep]
labels  <- gene_names_pathway[keep]
gene_idx <- match(symbols, gene.ref$gene_name)
if (anyNA(gene_idx)) stop("Genes not found in gene.ref: ", paste(symbols[is.na(gene_idx)], collapse = ", "))

## ---- 2) Instruments: top cis-SNPs per gene -------------------------------------
cov_df <- data.frame(age = donor$age, sex = donor$sex, PC1 = donor$PC1, PC2 = donor$PC2, PC3 = donor$PC3)
top_snp_list <- setNames(vector("list", length(symbols)), symbols)
for (g in seq_along(symbols)) {
  snp_idx <- which(countOverlaps(vcf.ref, gene.promoter.ref[gene_idx[g]]) > 0)
  if (!length(snp_idx)) next
  y <- RNA.log[gene_idx[g], ]
  pv <- vapply(snp_idx, function(s)
    summary(lm(y ~ genotype.mat[s, ] + age + sex + PC1 + PC2 + PC3, data = cov_df))$coefficients[2, 4], numeric(1))
  ord <- order(pv)
  top_snp_list[[g]] <- vcf.ref$ID[snp_idx[ord[seq_len(min(TOP_PER_GENE, length(ord)))]]]
}
snps <- unique(unlist(top_snp_list))
if (!length(snps)) stop("No SNPs selected.")
snp_idx2 <- match(snps, Subsetted_SNP_Names)
if (anyNA(snp_idx2)) stop("Selected SNPs not found: ", paste(head(snps[is.na(snp_idx2)], 20), collapse = ", "))

## ---- 3) Matrices and preprocessing -----------------------------------------------
Y_log <- t(RNA.log[gene_idx, , drop = FALSE]);        colnames(Y_log) <- symbols
X_raw <- t(genotype.mat[snp_idx2, , drop = FALSE]);   colnames(X_raw) <- snps
age_bin <- cut(donor$age, breaks = c(-Inf, 30, 40, 50, 60, 70, 80, Inf), right = FALSE, labels = FALSE)
U_raw <- cbind(SEX = as.numeric(donor$sex == "male"), AGE = age_bin, PC1 = donor$PC1, PC2 = donor$PC2, PC3 = donor$PC3)
if (anyNA(U_raw)) stop("Missing covariate values.")

Y <- standardize_cols(Y_log, "Expression")
X <- standardize_cols(X_raw, "Genotypes")
U <- standardize_cols(U_raw, "Covariates")
if (ncol(Y) != length(symbols)) stop("A gene had zero variance and was dropped; labels would no longer match.")
S <- summary_stats(Y, X, U)
message(sprintf("OneK1K: n = %d donors, p = %d genes, k = %d SNPs, l = %d covariates", S$n, ncol(Y), ncol(X), ncol(U)))

## ---- 4) Fit MR.RGM+ ---------------------------------------------------------------
set.seed(SEED)
fit <- RGM(Syy = S$Syy, Syx = S$Syx, Sxx = S$Sxx, Syu = S$Syu, Sxu = S$Sxu, Suu = S$Suu,
           D = matrix(1, ncol(Y), ncol(X)), n = S$n,
           nIter = MCMC$nIter, nBurnin = MCMC$nBurnin, Thin = MCMC$Thin,
           prior = "Spike and Slab", SigmaStarModel = "SSSL")

## ---- 5) Save -------------------------------------------------------------------------
dir.create("output", showWarnings = FALSE)
OneK1K_fit <- list(
  fit = fit, labels = labels, symbols = symbols, snps = colnames(X), top_snps = top_snp_list,
  n = S$n, p = ncol(Y), k = ncol(X), l = ncol(U),
  settings = list(top_per_gene = TOP_PER_GENE, mcmc = MCMC, seed = SEED, cis_window = "+/- 200 kb",
                  expression = "library-size normalised, log1p, centered and scaled",
                  genotypes  = "dosage (0/1/2; missing set to 0), centered and scaled",
                  covariates = "SEX (male = 1), AGE (decade bins 1-7), PC1-PC3; centered and scaled",
                  D = "all ones (MR.RGM+)", prior = "Spike and Slab", SigmaStarModel = "SSSL"),
  package_version = as.character(packageVersion("MR.RGM")),
  session_info = sessionInfo())
saveRDS(OneK1K_fit, file.path("output", "OneK1K_fit.rds"))
message("Saved output/OneK1K_fit.rds")
