## =========================================================
## OneK1K_ModelFitting.R
##
## Purpose:
##   Run MR.RGM on the OneK1K B-cell dataset.
##
## Prerequisite:
##   Run `OneK1K_PreProcessing.R` first (it should create / load):
##     - donor, gene.ref, vcf.ref, gene.promoter.ref
##     - RNA.count.adj, genotype.mat, Subsetted_SNP_Names
##
## Output:
##   Output_OneK1K (results from RGM)
## =========================================================

library(MR.RGM)
library(GenomicRanges)
library(dplyr)

## ---- 0) Ensure preprocessed objects exist ----
needed <- c("donor", "gene.ref", "vcf.ref", "gene.promoter.ref",
            "RNA.count.adj", "genotype.mat", "Subsetted_SNP_Names")
missing <- needed[!vapply(needed, exists, logical(1))]
if (length(missing) > 0) {
  stop("Missing objects. Run `OneK1K_PreProcessing.R` first.\nMissing: ",
       paste(missing, collapse = ", "))
}

## ---- 1) Gene lists: pathway names + dataset gene symbols ----
gene_names_pathway <- c(
  "CD19", "BCR", "AG", "BCAP", "P85", "PI3K", "FGR2B", "SHIP", "DOK3", "PTEN",
  "SHP2", "LYN", "SHP1", "CD22", "CD45", "CBP/PAG", "CSK", "PIR-B", "BAM32",
  "PLCY2", "BLNK", "LAB", "GRB2", "SOS", "SYK", "CBL", "EZRIN", "CLATHRIN",
  "VAV", "RAC", "CDC42", "RASGRP", "RASGAP", "PKC", "AKT", "MTOR", "P70S6K",
  "GSK3", "CARMA1", "TAK1", "BCL10", "MALT1", "DOK1", "IKK", "NFKB", "IKB",
  "RAS", "RAP", "RIAM", "PYK2", "HS1", "MEKK", "MEK", "P38", "JNK", "NFAT",
  "CALCINEURIN", "CAM", "CAMK", "C-RAF", "MEK1/2", "ERK1/2", "CD40", "FOXO",
  "ETS1", "OCT2", "BFL1", "BCL-XL", "ELK1", "BCL6", "EGR1", "JUN", "ATF2",
  "CREB", "MEF2C", "PIP3", "RAPL"
)

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
  "RASSF3"
)

stopifnot(length(gene_names_pathway) == length(gene_names_dataset))

## ---- 2) Remove genes with insufficient SNP coverage (pre-determined) ----
genes_to_remove <- c(
  "CD79A", "DOK3", "PTEN", "PTPN11", "PTPN6",
  "CDC42", "MTOR", "PPP3CA", "FOXO1", "POU2F2", "ELK1"
)

idx_remove <- which(gene_names_dataset %in% genes_to_remove)

gene_names_dataset_use <- gene_names_dataset[-idx_remove]
gene_names_pathway_use <- gene_names_pathway[-idx_remove]

## ---- 3) Map genes to indices in gene.ref ----
gene_idx <- match(gene_names_dataset_use, gene.ref$gene_name)
if (anyNA(gene_idx)) {
  missing_genes <- gene_names_dataset_use[is.na(gene_idx)]
  stop("These genes were not found in gene.ref$gene_name: ",
       paste(missing_genes, collapse = ", "))
}

## ---- 4) For each gene, select top cis SNPs (up to 15) by association p-value ----
# NOTE: This can be slow because it fits many linear models.
top_snp_list <- vector("list", length(gene_names_dataset_use))
names(top_snp_list) <- gene_names_dataset_use

for (g in gene_names_dataset_use) {

  g_idx <- which(gene.ref$gene_name == g)

  # SNPs overlapping promoter window
  snp_idx <- which(countOverlaps(vcf.ref, gene.promoter.ref[g_idx]) > 0)

  if (length(snp_idx) == 0) {
    top_snp_list[[g]] <- character(0)
    next
  }

  # Compute p-values for SNP ~ expression (adjusting for covariates)
  pvals <- rep(NA_real_, length(snp_idx))

  for (j in seq_along(snp_idx)) {
    fit <- lm(
      RNA.count.adj[g_idx, ] ~ genotype.mat[snp_idx[j], ] +
        donor$age + donor$sex + donor$PC1 + donor$PC2 + donor$PC3
    )
    # p-value for SNP effect (2nd coefficient)
    pvals[j] <- summary(fit)$coefficients[2, 4]
  }

  # Take top 15 by smallest p-value
  ord <- order(pvals)
  k <- min(15, length(ord))
  top_snp_list[[g]] <- vcf.ref$ID[snp_idx[ord[1:k]]]
}

unique_snps <- unique(unlist(top_snp_list))
if (length(unique_snps) == 0) stop("No SNPs selected. Check promoter overlaps / inputs.")

## ---- 5) Build matrices Y (expression) and X (genotypes) ----
# Y: n x p (donors x genes)
Y <- t(RNA.count.adj[gene_idx, , drop = FALSE])

# X: n x q (donors x SNPs)
snp_idx2 <- match(unique_snps, Subsetted_SNP_Names)
if (anyNA(snp_idx2)) {
  missing_snps <- unique_snps[is.na(snp_idx2)]
  stop("These selected SNP IDs were not found in Subsetted_SNP_Names: ",
       paste(head(missing_snps, 30), collapse = ", "),
       if (length(missing_snps) > 30) " ...")
}
X <- t(genotype.mat[snp_idx2, , drop = FALSE])

## ---- 6) Build covariate matrix U ----
# Sex: binary indicator (MALE = 1, FEMALE = 0)
# Age: categorical bins
# PCs: top genotype principal components
U <- donor %>%
  mutate(
    MALE = ifelse(sex == "male", 1, 0),
    AGE = case_when(
      age < 30 ~ 1,
      age >= 30 & age < 40 ~ 2,
      age >= 40 & age < 50 ~ 3,
      age >= 50 & age < 60 ~ 4,
      age >= 60 & age < 70 ~ 5,
      age >= 70 & age < 80 ~ 6,
      age >= 80 ~ 7
    )
  ) %>%
  dplyr::select(MALE, AGE, PC1, PC2, PC3) %>%
  as.matrix()

## ---- 7) Compute summary-statistic matrices ----
n <- nrow(Y)

Syy <- crossprod(Y) / n
Syx <- crossprod(Y, X) / n
Sxx <- crossprod(X) / n
Syu <- crossprod(Y, U) / n
Sxu <- crossprod(X, U) / n
Suu <- crossprod(U) / n

## ---- 8) Run RGM ----
set.seed(12345)

D_mat <- matrix(1, nrow = nrow(Syy), ncol = ncol(Sxx))

Output_OneK1K <- RGM(
  Syy = Syy, Syx = Syx, Sxx = Sxx,
  Syu = Syu, Sxu = Sxu, Suu = Suu,
  D   = D_mat, n = n,
  nIter = 50000, nBurnin = 10000, Thin = 10,
  prior = "Spike and Slab",
  SigmaStarModel = "SSSL"
)

