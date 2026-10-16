###############################################################################
## OneK1K_PreProcessing.R  —  B-cell data from Zenodo (donor_b_cell.rda)
##
## Steps
##   1) Load donor covariates, genotypes, SNP and gene annotation, RNA counts
##   2) cis windows: promoter +/- 200 kb around each gene's TSS
##   3) Donor filtering on relative library size (> 0.25 and <= median + 3 MAD)
##   4) Expression: library-size normalisation, then log1p (Version 4)
##   5) Genotypes: VCF strings -> dosage 0/1/2 (missing './.' set to 0, as in the
##      original analysis); keep SNPs with non-zero dosage in > 5% of donors
##
## Creates: donor, gene.ref, vcf.ref, gene.promoter.ref, RNA.log, genotype.mat,
##          Subsetted_SNP_Names                       -> then run OneK1K_ModelFitting.R
##
## Usage: set data_dir to the folder containing donor_b_cell.rda (Zenodo record
## 10.5281/zenodo.18489687), then source("OneK1K_PreProcessing.R") from this folder.
###############################################################################

suppressPackageStartupMessages(library(GenomicRanges))
if (!exists("data_dir")) data_dir <- "PATH_TO_ZENODO_DATA_FOLDER"
rda <- file.path(data_dir, "donor_b_cell.rda")
if (!file.exists(rda)) stop("Set data_dir to the folder containing donor_b_cell.rda.")
load(rda)

required_objs <- c("donor", "genotype", "vcf.ref", "RNA.count", "gene.ref")
missing <- required_objs[!vapply(required_objs, exists, logical(1))]
if (length(missing)) stop("Missing objects in donor_b_cell.rda: ", paste(missing, collapse = ", "))
message(sprintf("Loaded: %d donors, %d SNPs, %d genes", nrow(donor), nrow(genotype), nrow(RNA.count)))

## ---- 2) cis windows -----------------------------------------------------------
gene.promoter.ref <- promoters(gene.ref, upstream = 200000, downstream = 200000)

## ---- 3) Donor filtering on library size --------------------------------------
lib.size   <- colSums(RNA.count) / median(colSums(RNA.count))
donor.keep <- (lib.size <= median(lib.size) + 3 * mad(lib.size)) & (lib.size > 0.25)
donor      <- donor[donor.keep, ]
genotype   <- genotype[, donor.keep]
RNA.count  <- RNA.count[, donor.keep]
message(sprintf("Donors kept after library-size filtering: %d", sum(donor.keep)))

## ---- 4) Expression: library-size normalisation + log1p --------------------------
lib.size <- colSums(RNA.count) / median(colSums(RNA.count))
RNA.log  <- log1p(sweep(RNA.count, 2, lib.size, "/"))           # genes x donors

## ---- 5) Genotypes -----------------------------------------------------------------
genotype.mat <- matrix(NA_real_, nrow(genotype), ncol(genotype))
genotype.mat[genotype == "./."] <- 0
genotype.mat[genotype == "0/0"] <- 0
genotype.mat[genotype == "0/1"] <- 1
genotype.mat[genotype == "1/1"] <- 2
unexpected <- unique(genotype[is.na(genotype.mat)])
if (length(unexpected)) warning("Unmapped genotype strings: ", paste(head(unexpected, 10), collapse = ", "))
message(sprintf("Missing genotype calls ('./.', set to 0): %.3f%%", 100 * mean(genotype == "./.")))

SNP.keep <- apply(genotype.mat, 1, function(x) mean(x > 0, na.rm = TRUE)) > 0.05
Subsetted_SNP_Names <- vcf.ref$ID[SNP.keep]
genotype.mat <- genotype.mat[SNP.keep, ]
vcf.ref      <- vcf.ref[SNP.keep, ]
message(sprintf("SNPs kept (non-zero dosage in > 5%% of donors): %d", sum(SNP.keep)))
rm(genotype, RNA.count, lib.size, donor.keep, SNP.keep)
message("Preprocessing finished. Next: source('OneK1K_ModelFitting.R')")
