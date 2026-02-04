###############################################################################
## GTEx v7 (Muscle Skeletal) — Real-data analysis for MR.RGM+
##
## This script:
##  1) Loads preprocessed GTEx genotype + expression matrices (from Zenodo GTEx.zip)
##  2) Loads public GTEx v7 eQTL summary files (from Zenodo GTEx_Analysis_v7_eQTL.tar.gz)
##  3) Matches common donor IDs across genotype/expression/covariates
##  4) Selects a pathway gene set (mapped to GTEx gene symbols available in this dataset)
##  5) For each selected gene, chooses top cis-eQTL variants (by pval_nominal), and builds X
##  6) Builds covariate matrix U = [SEX (binary), AGE (ordinal)]
##  7) Runs MR.RGM+ with confounders via (Syy, Syx, Sxx, Syu, Sxu, Suu) and full D = 1
##
## Assumed folder layout after extracting Zenodo archives into:
##   reproducibility/real_data/GTEx/
##     ├── GTEx/                         (from GTEx.zip)
##     └── GTEx_Analysis_v7_eQTL/        (from GTEx_Analysis_v7_eQTL.tar.gz)
##
###############################################################################

## ---------------------------- ##
## Libraries
## ---------------------------- ##
library(R.matlab)
library(readr)
library(dplyr)
library(tidyr)
library(Matrix)   # large matrices
library(GIGrvg)
library(MR.RGM)

## ---------------------------- ##
## User: set this path
## ---------------------------- ##
## Replace with the folder that CONTAINS both "GTEx/" and "GTEx_Analysis_v7_eQTL/"
GTEX_FOLDER_PATH <- "PATH/TO/YOUR/GTEx_PARENT_FOLDER"  # e.g., "reproducibility/real_data/GTEx"


# SNP Data
SNP_data = readMat("GTEX_FOLDER_PATH/GTEx/geno450_snp012x_maf_GT_015.mat")

# Gene Data
# Read the muscle skeletal (Cell) data
Muscle_Skeletal_Data = readMat("GTEX_FOLDER_PATH/GTEx/GTEx_v7_4_tissues/GTEx_v7_Muscle_Skeletal_564_491.mat")


###################################################################################################

# Read significant variant-gene pairs file (Useful for getting SNPs for genes)
signif_variant_data <- read_delim("GTEX_FOLDER_PATH/GTEx_Analysis_v7_eQTL/Muscle_Skeletal.v7.signif_variant_gene_pairs.txt.gz", delim = "\t")

# Read egenes
egenes_Muscle_Skeletal <-  read_delim("GTEX_FOLDER_PATH/GTEx_Analysis_v7_eQTL/Muscle_Skeletal.v7.egenes.txt.gz", delim = "\t")


###################################################################################################
## Ids
# Get the donor id for SNPs
SNP_donor_ids = as.vector(read.table("GTEX_FOLDER_PATH/GTEx/donorid_450.txt",
                                     header = FALSE, stringsAsFactors = FALSE))

# Get the donor ids for Genes
Gene_donor_ids = unlist(Muscle_Skeletal_Data$v.donors)


# Find the common IDs
common_ids <- intersect(SNP_donor_ids$V1, Gene_donor_ids)


# Get the element-wise indices in SNP_donor_ids for each ID in common_ids
SNP_indices <- match(common_ids, SNP_donor_ids$V1)

# Get the element-wise indices in Gene_donor_ids for each ID in common_ids
Gene_indices <- match(common_ids, Gene_donor_ids)


# Subset the Gene data and SNP data
Muscle_Skeletal_Data_Subset = Muscle_Skeletal_Data$v.expr[, Gene_indices]
SNP_data_Subset = SNP_data$geno012[, SNP_indices]


###################################################################################################

# Read Covariates
Covariates <- read.table("GTEX_FOLDER_PATH/GTEx/GTEx_v7_4_tissues/Covariates_New/GTEx_v7_Annotations_SubjectPhenotypesDS.txt",
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Take subset
Covariates_Subset = Covariates[match(common_ids, Covariates$SUBJID),]

U <- Covariates_Subset %>%
  mutate(
    ## GTEx v7 coding typically: 1 = male, 2 = female
    SEX = ifelse(SEX == 1, 1, 0),
    AGE = case_when(
      AGE == "20-29" ~ 1,
      AGE == "30-39" ~ 2,
      AGE == "40-49" ~ 3,
      AGE == "50-59" ~ 4,
      AGE == "60-69" ~ 5,
      AGE == "70-79" ~ 6,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::select(SEX, AGE)

if (anyNA(U$AGE)) warning("Some AGE entries were not recognized and are NA.")
U <- as.matrix(U)

## ---------------------------- ##
## Gene set (pathway genes) and dataset gene symbols mapping
## ---------------------------- ##
gene_names_pathway <- c(
  "PTEN", "PI3K", "SHIP1", "RAS", "NF1", "INSULIN", "INSULIN RECEPTOR", "IRS1", "Tel2",
  "PDK1", "SGK", "PHLPP1/2", "AKT/PKB", "MTOR", "MLST8", "PRR5", "B-RAF", "PKCA", "ERK", "RSK1",
  "TSC2", "TSC1","LKB1", "PKC-Delta", "Folliculin/BHD", "AMPK", "GSK3B", "FOXO1", "REDD1/2",
  "FKBP38", "FKBP12", "PLD1", "PRAS40", "S6K", "ROC1", "HIF1/2A", "VHL", "EIF4A", "4EBP1", "RHEB",
  "PKA", "WNT", "MSIN1", "PA"
)
gene_names_pathway <- unique(gene_names_pathway)

## Dataset-specific gene symbols available in this GTEx-derived matrix
gene_names_dataset <- c(
  "PTEN", "PIK3C2B", "INPP5D", "RASD1", "NF1", "INS", "INSR", "IRS1", "TELO2",
  "PDK1", "SGK1", "PHLPP1", "AKT1", "MTOR", "MLST8", "PRR5", "BRAF", "PRKCA", "EPHB2", "RPS6KA1",
  "TSC2", "TSC1", "STK11", "RIPK4", "FLCN", "PRKAA2", "GSK3B", "FOXO1", "DDIT4",
  "FKBP8", "FKBP1A", "PLD1", "AKT1S1", "POLDIP3", "RBX1", "HIF1A", "VHL", "EIF4A1", "EIF4EBP1", "RHEB",
  "PRKACA", "WNT1", "MAPKAP1", "CKM"
)

## ---------------------------- ##
## Build Y from expression matrix
##   - Match dataset gene symbols to expression MAT gene list
## ---------------------------- ##
all_genes_in_mat <- unlist(Muscle_Skeletal_Data$v.genes)

RowNo_Genes <- match(gene_names_dataset, all_genes_in_mat)
keep_gene_rows <- which(!is.na(RowNo_Genes))
if (length(keep_gene_rows) == 0) stop("None of the requested genes were found in expression matrix.")

RowNo_Genes <- RowNo_Genes[keep_gene_rows]
Gene_Names  <- unlist(Muscle_Skeletal_Data$v.genes[RowNo_Genes])

## Y: (n samples) x (p genes)
Y <- t(Muscle_Skeletal_Data_Subset[RowNo_Genes, , drop = FALSE])
colnames(Y) <- Gene_Names

## Map to pathway labels for plotting (same order as gene_names_dataset)
gene_names_final <- gene_names_pathway[match(colnames(Y), gene_names_dataset)]

## ---------------------------- ##
## Trim to eGenes table (ensures gene_id mapping exists)
## ---------------------------- ##
gene_ids <- egenes_Muscle_Skeletal %>%
  filter(gene_name %in% colnames(Y)) %>%
  dplyr::select(gene_id, gene_name)

if (nrow(gene_ids) == 0) stop("No overlap between selected genes and eGenes file.")

gene_indices <- match(gene_ids$gene_name, colnames(Y))
gene_indices <- gene_indices[!is.na(gene_indices)]

Y <- Y[, gene_indices, drop = FALSE]
colnames(Y) <- gene_ids$gene_name

gene_names_final <- gene_names_final[match(colnames(Y), Gene_Names)]
gene_names_final <- gene_names_final[!is.na(gene_names_final)]

## ---------------------------- ##
## Keep only genes that have at least one significant variant in signif_variant_data
## ---------------------------- ##
has_sig <- gene_ids$gene_id %in% signif_variant_data$gene_id
gene_ids_final <- gene_ids[has_sig, , drop = FALSE]
if (nrow(gene_ids_final) == 0) stop("None of the selected genes have significant variant-gene pairs in the eQTL file.")

Y <- Y[, has_sig, drop = FALSE]
gene_names_final <- gene_names_final[has_sig]

## ---------------------------- ##
## Select variants: top N variants per gene by pval_nominal, restricted to SNPs present in geno MAT
## ---------------------------- ##
TOP_PER_GENE <- 25

filtered_variants <- signif_variant_data %>%
  filter(gene_id %in% gene_ids_final$gene_id)

selected_SNPs <- unlist(SNP_data$vrsid)

top_variants <- filtered_variants %>%
  group_by(gene_id) %>%
  arrange(pval_nominal, .by_group = TRUE) %>%
  slice_head(n = TOP_PER_GENE) %>%
  ungroup() %>%
  filter(variant_id %in% selected_SNPs)

variant_ids <- unique(top_variants$variant_id)
if (length(variant_ids) == 0) stop("No selected variants found in genotype MAT (SNP_data$vrsid).")

## ---------------------------- ##
## Build X from genotype matrix
##   geno_subset is (SNPs x n donors) -> transpose to (n donors x k variants)
## ---------------------------- ##
X_row_idx <- match(variant_ids, selected_SNPs)
X_row_idx <- X_row_idx[!is.na(X_row_idx)]
if (length(X_row_idx) == 0) stop("Variant IDs could not be matched to genotype SNP IDs.")

X <- t(SNP_data_Subset[X_row_idx, , drop = FALSE])

## ---------------------------- ##
## Compute summary matrices
## ---------------------------- ##
p <- ncol(Y)
k <- ncol(X)
l <- ncol(U)

message("Final dimensions: n=", n, ", p=", p, ", k=", k, ", l=", l)

Syy <- t(Y) %*% Y / n
Syx <- t(Y) %*% X / n
Sxx <- t(X) %*% X / n

Syu <- t(Y) %*% U / n
Sxu <- t(X) %*% U / n
Suu <- t(U) %*% U / n

## ---------------------------- ##
## Run MR.RGM+ with confounders (full D = 1)
## ---------------------------- ##
set.seed(123)

D_full <- matrix(1, nrow = nrow(Syy), ncol = ncol(Sxx))

Output_GTEx <- RGM(
  Syy = Syy, Syx = Syx, Sxx = Sxx,
  Syu = Syu, Sxu = Sxu, Suu = Suu,
  D = D_full, n = n,
  nIter = 50000, nBurnin = 10000, Thin = 10,
  prior = "Spike and Slab", SigmaStarModel = "SSSL"
)

## Output_GTEx contains MR.RGM+ estimates
message("Done. Output object stored in: Output_GTEx")
###############################################################################
