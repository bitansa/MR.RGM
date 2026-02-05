## =========================================================
## OneK1K real data analysis: setup + preprocessing
## - Downloads data from Zenodo (zip), extracts, sets wd
## - Loads donor/genotype/RNA objects
## - Library-size normalization + donor filtering
## - Genotype conversion to numeric + MAF-like filtering
## =========================================================

# Replace with the local path where the Zenodo files were saved
data_dir <- "PATH_TO_ZENODO_DATA_FOLDER"
setwd(data_dir)

library(GenomicRanges)

## ---- Load data ----
load("donor_b_cell_bitan.rda")

# donor      : donor-level covariates (including genotype PCs)
# genotype   : genotype matrix (VCF-style encoding)
# vcf.ref    : SNP information (IDs and genomic coordinates)
# RNA.count  : gene expression counts aggregated by donor
# gene.ref   : gene annotations (symbols and coordinates)

# Expecting these objects to exist after load:
# donor, genotype, vcf.ref, RNA.count, gene.ref
required_objs <- c("donor", "genotype", "vcf.ref", "RNA.count", "gene.ref")
missing <- required_objs[!vapply(required_objs, exists, logical(1))]
if (length(missing) > 0) stop("Missing objects in .rda: ", paste(missing, collapse = ", "))

## Quick checks
message("Dimensions:")
message("  donor:     ", paste(dim(donor), collapse = " x "))
message("  genotype:  ", paste(dim(genotype), collapse = " x "))
message("  RNA.count: ", paste(dim(RNA.count), collapse = " x "))

## ---- 2) Define gene promoter regions (cis window) ----
# Promoters +/- 200kb around TSS
gene.promoter.ref <- promoters(gene.ref, upstream = 200000, downstream = 200000)

## ---- 3) Library size normalization + donor filtering ----
# Library size per donor relative to the median
lib.size <- colSums(RNA.count) / median(colSums(RNA.count))
hist(lib.size, main = "Library size (relative)", xlab = "Relative library size")

# Filter outliers:
# - Keep donors within median + 3*MAD
# - Also remove extremely small library size donors (<= 0.25)
donor.keep <- (lib.size <= median(lib.size) + 3 * mad(lib.size)) & (lib.size > 0.25)

# Apply filter consistently across all donor-indexed objects
donor     <- donor[donor.keep, ]
genotype  <- genotype[, donor.keep]
RNA.count <- RNA.count[, donor.keep]
rm(donor.keep)

# Recompute library size after filtering
lib.size <- colSums(RNA.count) / median(colSums(RNA.count))
hist(lib.size, main = "Library size after filtering", xlab = "Relative library size")

# Normalize RNA counts by library size (simple scaling)
RNA.count.adj <- sweep(RNA.count, 2, lib.size, "/")

## ---- 4) Convert genotype strings to numeric dosage ----
# genotype is assumed to be a character matrix with values like:
# './.' , '0/0', '0/1', '1/1'
genotype.mat <- matrix(NA_real_, nrow = nrow(genotype), ncol = ncol(genotype))

# Map common genotype encodings
genotype.mat[genotype == "./."] <- 0   # missing treated as 0 (keep your original behavior)
genotype.mat[genotype == "0/0"] <- 0
genotype.mat[genotype == "0/1"] <- 1
genotype.mat[genotype == "1/1"] <- 2

# Optional sanity check: any unexpected encodings?
unexpected <- unique(genotype[is.na(genotype.mat)])
unexpected <- unexpected[!is.na(unexpected)]
if (length(unexpected) > 0) {
  warning("Unexpected genotype strings encountered (not mapped): ",
          paste(head(unexpected, 10), collapse = ", "),
          if (length(unexpected) > 10) " ...")
}

## ---- 5) SNP filtering (keep SNPs with >5% non-zero dosage) ----
# Your original logic: keep SNPs where mean(x>0) > 0.05
SNP.keep <- apply(genotype.mat, 1, function(x) mean(x > 0, na.rm = TRUE)) > 0.05
message("SNPs kept: ", sum(SNP.keep), " / ", length(SNP.keep))

Subsetted_SNP_Names <- vcf.ref$ID[SNP.keep]

# Subset all SNP-indexed objects consistently
genotype     <- genotype[SNP.keep, ]
genotype.mat <- genotype.mat[SNP.keep, ]
vcf.ref      <- vcf.ref[SNP.keep, ]
rm(SNP.keep)

## Done
message("Finished preprocessing.")
