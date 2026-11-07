###############################################################################
## OneK1K_PreProcessing.R  —  B-cell data from Zenodo (donor_b_cell.rda)
##
## Steps
##   1) Load donor covariates, genotypes, SNP and gene annotation, RNA counts
##   2) cis windows: +/- CIS_KB kb (default 1000, the standard eQTL window) around each gene's transcription start
##      site (gene annotation: Ensembl 75, GRCh37; the array positions must be on the same build)
##   3) Donor filtering on relative library size (> 0.25 and <= median + 3 MAD)
##   4) Expression: library-size normalisation, then log1p
##   5) Genotypes (subset option): VCF strings -> dosage 0/1/2 (missing './.' set to 0); keep SNPs with
##      minor allele carried by > 5% of donors
##      Default (ONEK_GENOTYPES = "chip"): every genotyped SNP of the OneK1K array (PLINK files in
##      data_dir/science_data/raw_genotype; about 760,000 SNPs, of which donor_b_cell.rda holds a
##      96,663-SNP subset) in the cis windows of the pathway genes; SNPs with > 5% missing calls are dropped,
##      the remaining missing calls are mean-imputed, and SNPs whose minor allele is carried by <= 5% of the
##      donors (genotype 1 or 2 copies of the minor allele) are dropped. Needs snpStats.
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
CIS_KB <- if (exists("ONEK_CIS_KB")) ONEK_CIS_KB else 1000   # cis window half-width in kb (override: ONEK_CIS_KB)
## windows reaching beyond a chromosome end are trimmed to the chromosome
gene.promoter.ref <- trim(suppressWarnings(promoters(gene.ref, upstream = CIS_KB * 1000, downstream = CIS_KB * 1000)))

## ---- 2b) Genotypes: all array SNPs (default) or the 96,663-SNP subset (ONEK_GENOTYPES <- "subset") ----
GENOTYPES <- if (exists("ONEK_GENOTYPES")) ONEK_GENOTYPES else "chip"     # "chip" (all array SNPs, final analysis) or "subset"
stopifnot(GENOTYPES %in% c("subset", "chip"))
if (GENOTYPES == "chip") {
  suppressPackageStartupMessages(library(snpStats))
  gdir <- file.path(data_dir, "science_data", "raw_genotype")
  if (!file.exists(file.path(gdir, "OneK1K_AllChr.bed"))) stop("Array genotypes (OneK1K_AllChr.bed/.bim/.fam) not found in ", gdir)
  bim <- read.table(file.path(gdir, "OneK1K_AllChr.bim"), stringsAsFactors = FALSE,
                    col.names = c("chr", "id", "cm", "pos", "a1", "a2"), colClasses = c("character", "character", "numeric", "integer", "character", "character"))
  bim$chr[bim$chr == "23"] <- "X"
  bim <- bim[bim$chr %in% c(as.character(1:22), "X") & !duplicated(bim$id) & !(bim$id %in% bim$id[duplicated(bim$id)]), ]
  source("OneK1K_nodes.R")                                  # pathway genes: read only their cis windows
  cand_win <- gene.promoter.ref[gene.ref$gene_name %in% unlist(nodes)]
  bim <- bim[overlapsAny(GRanges(bim$chr, IRanges(bim$pos, width = 1)), cand_win), ]
  if (!nrow(bim)) stop("No array SNP falls in the cis windows of the pathway genes: the genome build of the array (.bim positions) ",
                       "and of the gene annotation (GRCh37, Ensembl 75) must agree, and the chromosome names must match.")
  chip <- read.plink(file.path(gdir, "OneK1K_AllChr.bed"), file.path(gdir, "OneK1K_AllChr.bim"),
                     file.path(gdir, "OneK1K_AllChr.fam"), select.snps = bim$id)
  G <- as(chip$genotypes, "numeric")                         # donors x SNPs (copies of allele 2; NA = missing)
  ids <- paste(chip$fam$pedigree, chip$fam$member, sep = "_")
  idx <- match(colnames(genotype), ids)
  if (anyNA(idx)) stop(sum(is.na(idx)), " donors of donor_b_cell.rda are not in the array genotype file.")
  G <- G[idx, , drop = FALSE]
  call_ok <- colMeans(is.na(G)) <= 0.05                    # drop SNPs with > 5% missing calls (incl. all-missing)
  message(sprintf("Array genotypes: %d SNPs in the cis windows of %d pathway genes; %d dropped (> 5%% missing); missing calls in the rest %.3f%% (mean-imputed)",
                  ncol(G), length(cand_win), sum(!call_ok), 100 * mean(is.na(G[, call_ok]))))
  G <- G[, call_ok, drop = FALSE]
  miss <- which(is.na(G), arr.ind = TRUE)
  if (nrow(miss)) G[miss] <- colMeans(G, na.rm = TRUE)[miss[, 2]]
  map <- bim[match(colnames(G), bim$id), ]
  vcf.ref <- GRanges(map$chr, IRanges(map$pos, width = 1), ID = map$id)
  genotype.chip <- t(G); colnames(genotype.chip) <- colnames(genotype)   # SNPs x donors
  rm(G, chip, bim, map, miss, idx, ids, cand_win, call_ok)
}

## ---- 3) Donor filtering on library size --------------------------------------
lib.size   <- colSums(RNA.count) / median(colSums(RNA.count))
donor.keep <- (lib.size <= median(lib.size) + 3 * mad(lib.size)) & (lib.size > 0.25)
donor      <- donor[donor.keep, ]
genotype   <- genotype[, donor.keep]
if (GENOTYPES == "chip") genotype.chip <- genotype.chip[, donor.keep]
RNA.count  <- RNA.count[, donor.keep]
message(sprintf("Donors kept after library-size filtering: %d", sum(donor.keep)))

## ---- 4) Expression: library-size normalisation + log1p --------------------------
lib.size <- colSums(RNA.count) / median(colSums(RNA.count))
RNA.log  <- log1p(sweep(RNA.count, 2, lib.size, "/"))           # genes x donors

## ---- 5) Genotypes -----------------------------------------------------------------
if (GENOTYPES == "subset") {
  genotype.mat <- matrix(NA_real_, nrow(genotype), ncol(genotype))
  genotype.mat[genotype == "./."] <- 0
  genotype.mat[genotype == "0/0"] <- 0
  genotype.mat[genotype == "0/1"] <- 1
  genotype.mat[genotype == "1/1"] <- 2
  unexpected <- unique(genotype[is.na(genotype.mat)])
  if (length(unexpected)) warning("Unmapped genotype strings: ", paste(head(unexpected, 10), collapse = ", "))
  message(sprintf("Missing genotype calls ('./.', set to 0): %.3f%%", 100 * mean(genotype == "./.")))
  ## minor-allele carriers (the same rule as for the array genotypes below)
  SNP.keep <- apply(genotype.mat, 1, function(x) { r <- round(x); f <- mean(r, na.rm = TRUE) / 2; if (f <= 0.5) mean(r >= 1, na.rm = TRUE) else mean(r <= 1, na.rm = TRUE) }) > 0.05
} else {
  genotype.mat <- genotype.chip; rm(genotype.chip)
  ## minor-allele carriers: donors with at least one copy of the less frequent allele (allele 2 if its frequency <= 0.5)
  SNP.keep <- apply(genotype.mat, 1, function(x) { r <- round(x); f <- mean(r) / 2; if (f <= 0.5) mean(r >= 1) else mean(r <= 1) }) > 0.05
  SNP.keep <- SNP.keep %in% TRUE
}
Subsetted_SNP_Names <- vcf.ref$ID[SNP.keep]
genotype.mat <- genotype.mat[SNP.keep, ]
vcf.ref      <- vcf.ref[SNP.keep, ]
message(sprintf("SNPs kept (%s genotypes; minor allele carried by > 5%% of donors): %d", GENOTYPES, sum(SNP.keep)))
rm(genotype, RNA.count, lib.size, donor.keep, SNP.keep)
message("Preprocessing finished. Next: source('OneK1K_ModelFitting.R')")
