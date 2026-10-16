###############################################################################
## GTEx_ModelFitting.R  —  GTEx v7 skeletal muscle, MR.RGM+ with covariates
##
## Steps
##   1) Load preprocessed genotype + expression matrices (Zenodo GTEx.zip) and the
##      public GTEx v7 eQTL summary files (Zenodo GTEx_Analysis_v7_eQTL.tar.gz)
##   2) Keep donors with both genotype and expression data
##   3) Map the mTOR-pathway node labels to the gene symbols in the data
##   4) Keep genes listed as eGenes with at least one significant variant-gene pair
##   5) Instruments: for each gene the top 25 cis-eQTL variants (smallest pval_nominal)
##   6) Preprocessing (Version 4):
##        expression  -> log2(x + 1), then centered and scaled per gene
##        genotypes   -> centered and scaled per variant (0/1/2 dosage)
##        covariates  -> SEX (1 = male, 0 = female) and AGE (ordinal decade bins 1-6),
##                       centered and scaled
##      Centering removes the need for an intercept; scaling puts every variable on the
##      unit scale assumed by the prior hyperparameters.
##   7) Fit MR.RGM+ (all variants allowed for all genes, D = 1) from summary statistics
##   8) Save the fit and everything needed to interpret it to output/GTEx_fit.rds
##
## Data layout (after extracting both Zenodo archives into one folder):
##   GTEX_FOLDER_PATH/
##     GTEx/                    (GTEx.zip)
##     GTEx_Analysis_v7_eQTL/   (GTEx_Analysis_v7_eQTL.tar.gz)
##
## Usage: setwd() to this folder, set GTEX_FOLDER_PATH (here or before sourcing), then
##   source("GTEx_ModelFitting.R")      # then source("GTEx_Results.R")
###############################################################################

suppressPackageStartupMessages({
  library(R.matlab); library(readr); library(dplyr); library(MR.RGM)
})
source("../real_data_functions.R")

if (!exists("GTEX_FOLDER_PATH")) GTEX_FOLDER_PATH <- "PATH/TO/YOUR/GTEx_PARENT_FOLDER"
if (!dir.exists(file.path(GTEX_FOLDER_PATH, "GTEx")))
  stop("Set GTEX_FOLDER_PATH to the folder containing GTEx/ and GTEx_Analysis_v7_eQTL/.")

## Settings (reported in the manuscript)
TOP_PER_GENE <- 25
MCMC <- list(nIter = 50000, nBurnin = 10000, Thin = 10)
SEED <- 123

## ---- 1) Data -------------------------------------------------------------------
SNP_data   <- readMat(file.path(GTEX_FOLDER_PATH, "GTEx", "geno450_snp012x_maf_GT_015.mat"))
Expr_data  <- readMat(file.path(GTEX_FOLDER_PATH, "GTEx", "GTEx_v7_4_tissues", "GTEx_v7_Muscle_Skeletal_564_491.mat"))
signif_var <- read_delim(file.path(GTEX_FOLDER_PATH, "GTEx_Analysis_v7_eQTL", "Muscle_Skeletal.v7.signif_variant_gene_pairs.txt.gz"),
                         delim = "\t", show_col_types = FALSE)
egenes     <- read_delim(file.path(GTEX_FOLDER_PATH, "GTEx_Analysis_v7_eQTL", "Muscle_Skeletal.v7.egenes.txt.gz"),
                         delim = "\t", show_col_types = FALSE)
pheno      <- read.table(file.path(GTEX_FOLDER_PATH, "GTEx", "GTEx_v7_4_tissues", "Covariates_New",
                                   "GTEx_v7_Annotations_SubjectPhenotypesDS.txt"),
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## ---- 2) Donors with genotype and expression ------------------------------------
snp_donors  <- read.table(file.path(GTEX_FOLDER_PATH, "GTEx", "donorid_450.txt"),
                          header = FALSE, stringsAsFactors = FALSE)$V1
expr_donors <- unlist(Expr_data$v.donors)
common_ids  <- intersect(snp_donors, expr_donors)
expr_mat    <- Expr_data$v.expr[, match(common_ids, expr_donors)]          # genes x donors
geno_mat    <- SNP_data$geno012[, match(common_ids, snp_donors)]           # variants x donors
pheno       <- pheno[match(common_ids, pheno$SUBJID), ]

## ---- 3) Pathway node labels and the gene symbols used in the data -------------
gene_names_pathway <- c(
  "PTEN", "PI3K", "SHIP1", "RAS", "NF1", "INSULIN", "INSULIN RECEPTOR", "IRS1", "Tel2",
  "PDK1", "SGK", "PHLPP1/2", "AKT/PKB", "MTOR", "MLST8", "PRR5", "B-RAF", "PKCA", "ERK", "RSK1",
  "TSC2", "TSC1", "LKB1", "PKC-Delta", "Folliculin/BHD", "AMPK", "GSK3B", "FOXO1", "REDD1/2",
  "FKBP38", "FKBP12", "PLD1", "PRAS40", "S6K", "ROC1", "HIF1/2A", "VHL", "EIF4A", "4EBP1", "RHEB",
  "PKA", "WNT", "MSIN1", "PA")
gene_names_dataset <- c(
  "PTEN", "PIK3C2B", "INPP5D", "RASD1", "NF1", "INS", "INSR", "IRS1", "TELO2",
  "PDK1", "SGK1", "PHLPP1", "AKT1", "MTOR", "MLST8", "PRR5", "BRAF", "PRKCA", "EPHB2", "RPS6KA1",
  "TSC2", "TSC1", "STK11", "RIPK4", "FLCN", "PRKAA2", "GSK3B", "FOXO1", "DDIT4",
  "FKBP8", "FKBP1A", "PLD1", "AKT1S1", "POLDIP3", "RBX1", "HIF1A", "VHL", "EIF4A1", "EIF4EBP1", "RHEB",
  "PRKACA", "WNT1", "MAPKAP1", "CKM")
stopifnot(length(gene_names_pathway) == length(gene_names_dataset))
all_genes <- unlist(Expr_data$v.genes)

## ---- 4) Genes that are eGenes with >= 1 significant variant-gene pair ----------
##      (order follows the eGenes file, as in the original analysis)
eg <- egenes %>% filter(gene_name %in% intersect(gene_names_dataset, all_genes)) %>% select(gene_id, gene_name)
eg <- eg[eg$gene_id %in% signif_var$gene_id, ]
if (!nrow(eg)) stop("No selected gene has significant variant-gene pairs.")
symbols <- eg$gene_name                                             # final genes (data symbols)
labels  <- gene_names_pathway[match(symbols, gene_names_dataset)]   # final node labels

## ---- 5) Instruments: top cis-eQTL variants per gene ----------------------------
geno_ids <- unlist(SNP_data$vrsid)
top_var <- signif_var %>% filter(gene_id %in% eg$gene_id) %>%
  group_by(gene_id) %>% arrange(pval_nominal, .by_group = TRUE) %>%
  slice_head(n = TOP_PER_GENE) %>% ungroup() %>% filter(variant_id %in% geno_ids)
variant_ids <- unique(top_var$variant_id)
if (!length(variant_ids)) stop("No selected variant found in the genotype matrix.")

## ---- 6) Matrices and preprocessing ----------------------------------------------
Y_raw <- t(expr_mat[match(symbols, all_genes), , drop = FALSE]); colnames(Y_raw) <- symbols
X_raw <- t(geno_mat[match(variant_ids, geno_ids), , drop = FALSE]); colnames(X_raw) <- variant_ids
age_bin <- c("20-29" = 1, "30-39" = 2, "40-49" = 3, "50-59" = 4, "60-69" = 5, "70-79" = 6)
U_raw <- cbind(SEX = as.numeric(pheno$SEX == 1), AGE = unname(age_bin[pheno$AGE]))
if (anyNA(U_raw)) stop("Unrecognised SEX/AGE codes in the phenotype file.")
if (any(Y_raw < 0, na.rm = TRUE)) stop("Negative expression values: log2(x + 1) is not appropriate.")

Y <- standardize_cols(log2(Y_raw + 1), "Expression")
X <- standardize_cols(X_raw, "Genotypes")
U <- standardize_cols(U_raw, "Covariates")
if (ncol(Y) != length(symbols)) stop("A gene had zero variance and was dropped; labels would no longer match.")
S <- summary_stats(Y, X, U)
message(sprintf("GTEx: n = %d donors, p = %d genes, k = %d variants, l = %d covariates",
                S$n, ncol(Y), ncol(X), ncol(U)))

## ---- 7) Fit MR.RGM+ -------------------------------------------------------------
set.seed(SEED)
fit <- RGM(Syy = S$Syy, Syx = S$Syx, Sxx = S$Sxx, Syu = S$Syu, Sxu = S$Sxu, Suu = S$Suu,
           D = matrix(1, ncol(Y), ncol(X)), n = S$n,
           nIter = MCMC$nIter, nBurnin = MCMC$nBurnin, Thin = MCMC$Thin,
           prior = "Spike and Slab", SigmaStarModel = "SSSL")

## ---- 8) Save ---------------------------------------------------------------------
dir.create("output", showWarnings = FALSE)
GTEx_fit <- list(
  fit = fit, labels = labels, symbols = symbols, variants = colnames(X),
  n = S$n, p = ncol(Y), k = ncol(X), l = ncol(U),
  settings = list(top_per_gene = TOP_PER_GENE, mcmc = MCMC, seed = SEED,
                  expression = "log2(x + 1), centered and scaled",
                  genotypes  = "dosage (0/1/2), centered and scaled",
                  covariates = "SEX (male = 1), AGE (decade bins 1-6); centered and scaled",
                  D = "all ones (MR.RGM+)", prior = "Spike and Slab", SigmaStarModel = "SSSL"),
  package_version = as.character(packageVersion("MR.RGM")),
  session_info = sessionInfo())
saveRDS(GTEx_fit, file.path("output", "GTEx_fit.rds"))
message("Saved output/GTEx_fit.rds")
