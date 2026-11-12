# OneK1K B cells (KEGG B cell receptor signaling pathway, hsa04662)

Reproduces the OneK1K analysis of the manuscript: MR.RGM with covariates (MR.RGM+ as sensitivity analysis) on the
B-cell receptor pathway genes that have at least one strong cis-SNP instrument.

| Script | Purpose |
|---|---|
| `OneK1K_PreProcessing.R` | Load data, cis windows, donor filtering, log expression, genotype dosages |
| `OneK1K_nodes.R` | Pathway genes: KEGG hsa04662 (82 named genes, grouped by KEGG orthology) |
| `OneK1K_ModelFitting.R` | Instrument selection, preprocessing, MR.RGM (or MR.RGM+) fit, four chains pooled → `output/OneK1K_fit_MRRGM*.rds` |
| `OneK1K_Competitors.R` | Every competing method of the simulation study on the same genes, SNPs and donors → `output/OneK1K_competitors.rds` |
| `OneK1K_PriorSensitivity.R` | MR.RGM under a uniform Beta(1, 1) prior on the causal inclusion probabilities → `output/OneK1K_fit_MRRGM_Beta11.rds`, `tables/onek_prior_sensitivity.{csv,tex}` |
| `OneK1K_Bootstrap.R` | Parametric bootstrap: the fitted MR.RGM network taken as the truth, expression generated for the same donors from their genotypes and covariates (100 replicates), MR.RGM refitted → `output/OneK1K_bootstrap.rds`, `tables/onek_bootstrap{,_edges}.{csv,tex}` |
| `OneK1K_Results.R` | All tables (`tables/`) and figures (`figures/`) of the manuscript, generated from the fits |
| `OneK1K_RunAll.R` | Step 3 below in one call, with the model of every step set explicitly |
| `../real_data_functions.R` | Shared helpers (preprocessing, tables, network figures, modules, motifs, diagnostics) |

## Data

Zenodo, DOI **10.5281/zenodo.18489687** (it resolves to the latest version of the record)
<https://doi.org/10.5281/zenodo.18489687>

* `donor_b_cell.rda`: donor covariates, B-cell RNA counts, gene annotation (Ensembl 75, GRCh37) and a 96,663-SNP
  genotype subset. Download it into a folder of your choice (`data_dir`) and check it:
  `tools::md5sum("donor_b_cell.rda")` must give `63a105565f84a45cf654e7860db8b5e9`.
* Array genotypes (default, `ONEK_GENOTYPES = "chip"`, the setting of the manuscript): the PLINK files
  `OneK1K_AllChr.bed/.bim/.fam` of the OneK1K genotyping array (about 760,000 SNPs, GRCh37) in
  `data_dir/science_data/raw_genotype/`, from the public release of the OneK1K Consortium (Yazar et al., Science 2022; <https://onek1k.org>); they are not part of the Zenodo record. Needs the Bioconductor
  package `snpStats`. With `ONEK_GENOTYPES <- "subset"` the analysis runs from `donor_b_cell.rda` alone (fewer cis-SNPs), which does not
  reproduce the paper.

## Run (one R session; a few hours: SNP-wise regressions, four chains per model, the competing methods)

```r
setwd("path/to/MR.RGM/reproducibility/real_data/OneK1K")
data_dir <- "path/to/folder/with/donor_b_cell.rda"
source("OneK1K_PreProcessing.R")
## 1. candidate genes and their instrument strengths (no fit) -> output/OneK1K_candidates.rds
ONEK_STRENGTH_ONLY <- TRUE; source("OneK1K_ModelFitting.R"); rm(ONEK_STRENGTH_ONLY)
## 2. largest supported gene set (design-matched simulation; reads output/OneK1K_candidates.rds); prints
##    "largest supported gene set = 10"
setwd("../../package_check"); source("check_design_matched.R"); setwd("../real_data/OneK1K")
## 3. MR.RGM and MR.RGM+ fits (four chains each), every competing method, prior sensitivity, parametric bootstrap,
##    tables and figures
ONEK_N_GENES <- 10
source("OneK1K_RunAll.R")
```
`OneK1K_RunAll.R` runs `OneK1K_ModelFitting.R` (MR.RGM, then MR.RGM+: `output/OneK1K_fit_MRRGM.rds`,
`OneK1K_fit_MRRGMplus.rds`, `OneK1K_analysis_data.rds`), `OneK1K_Competitors.R` (`output/OneK1K_competitors.rds`),
`OneK1K_PriorSensitivity.R` (`output/OneK1K_fit_MRRGM_Beta11.rds`), `OneK1K_Bootstrap.R` (`output/OneK1K_bootstrap.rds`) and `OneK1K_Results.R` for both models (`tables/`, `figures/`), setting the model of each step itself.

`ONEK_OUT_DIR` (default: this folder) sets where `output/`, `tables/` and `figures/` are written.

Checkpoints in the console output: 891 of 961 donors kept; 35 genes with an instrument with F > 10; 14 candidate genes
pass the Bonferroni check; the design-matched simulation supports 10 genes; 48 instruments in total.

## Analysis settings

* **Donors:** relative library size > 0.25 and ≤ median + 3 MAD.
* **Expression:** library-size normalized, log1p.
* **Genotypes:** array dosages 0/1/2; SNPs with > 5% missing calls dropped, remaining missing calls mean-imputed;
  SNPs whose minor allele is carried by > 5% of donors; only SNPs in the cis windows of the pathway genes.
* **Genes:** the 82 named genes of KEGG pathway hsa04662, B cell receptor signaling (`OneK1K_nodes.R`, retrieved
  from the KEGG REST API on 26 October 2026; genes sharing a KEGG orthology group form one node); genes without an
  instrument are removed. Figures and tables label genes by symbol; `tables/onek_label_map.tex` gives the node.
* **Instruments:** for each gene, the SNPs within ±1 Mb of the transcription start site (`ONEK_CIS_KB`, default
  1000) with F > 10 in log-expression ~ SNP + sex + age (decade bins) + PC1 + PC2 + PC3 (the covariates of the
  model), LD-clumped in R (strongest first; a SNP is kept only if its r² with every SNP already kept is ≤ 0.5,
  `ONEK_LD_R2`), at most 10 per gene (`ONEK_TOP_PER_GENE`). A SNP selected for two genes is kept for the gene with
  the larger |t|, and instruments of different genes with r² > 0.5 are resolved the same way, so that the
  instruments of different genes are not near-collinear (the largest remaining between-gene r² is saved). Each SNP
  may affect only its own gene and its effect is not subject to selection (MR.RGM, main analysis);
  `ONEK_MODEL <- "MR.RGM+"` lets every SNP affect every gene with spike-and-slab selection (sensitivity analysis).
  The t-statistics of every instrument on every other analyzed gene are saved and reported (`cross_gene_t`); a
  large cross-gene association can be mediated by the instrument's own gene, so it is not by itself evidence
  against the exclusion restriction.
* **Gene selection:** for every gene the partial R² of its own instruments given the covariates, the joint
  first-stage F statistic and μ² = n R² / (1 − R²) are computed. A gene is kept only if its best cis-SNP is
  significant after Bonferroni correction over its cis-SNPs (`ONEK_BONF_ALPHA`, default 0.05; the instruments are
  selected in the same data), and the m genes with the largest μ² are analyzed (`ONEK_N_GENES`, default 10; optional floor
  `ONEK_MIN_MU2`). The number of genes is chosen beforehand by the design-matched simulation: `ONEK_STRENGTH_ONLY <- TRUE`
  writes the instrument strengths of every candidate gene (`output/OneK1K_candidates.rds`, no fit), and
  `package_check/check_design_matched.R` gives the largest gene set that MR.RGM recovers at these strengths and n
  (mean AUC ≥ 0.80, mean chain-half PIP correlation ≥ 0.90, error-variance warning in at most 10% of the fits).
* **Preprocessing:** covariates SEX (male = 1), AGE (decade bins 1–7), PC1–PC3; every column of expression,
  genotypes and covariates centered and scaled.
* **Model:** spike-and-slab prior, SSSL error covariance; 100,000 iterations, 20,000 burn-in, thinning 10
  (8,000 posterior samples per chain) for MR.RGM and MR.RGM_NoConf; MR.RGM+, which estimates every SNP-gene effect,
  500,000 iterations, 100,000 burn-in, thinning 50 (8,000 per chain), so that a lack of convergence cannot be
  attributed to the chain length. Four chains (seeds 12345, 54321, 1001, 2002, the same starting point) are run
  for every model (MR.RGM, MR.RGM+ and MR.RGM_NoConf); convergence: rank-normalized split-R̂ ≤ 1.01 for the
  log-likelihood and every causal edge indicator, bulk effective sample size and the PIPs of every pair of chains
  (`convergence` in the saved fit; package `posterior`). Converged chains are pooled (32,000 draws); for a model whose
  chains do not converge the network, modules and motifs are not drawn and the PIPs of every chain are reported.
  The package's warnings (weak instruments, identification, spectral radius, chain halves, acceptance,
  error-variance ratio) and diagnostics are recorded with the fit and written to
  `tables/onek*_mcmc_diagnostics.csv` and `onek*_package_warnings.txt`. Requires MR.RGM ≥ 0.1.4.
* **Tables for the paper** (`tables/onek_candidates.tex`, `onek_instruments_list.tex`, `onek_all_edges.tex`,
  `realdata_diagnostics.tex`, `onek_model_comparison.tex`, `onek_methods_comparison.tex`, `onek_edge_support.tex`)
  are written by `OneK1K_Results.R` from the fitted objects.
* **Competing methods** (`OneK1K_Competitors.R`): every method of the simulation study that applies to a network of
  genes (MR.RGM_NoConf, MVMR, MR2G, MrDAG one gene at a time as outcome, GrIVET on the original scale of the data, IVW,
  simple and weighted median, OneSampleMR, mrbayes; MR2 needs an exposure/outcome split) on the same genes, SNPs and donors, with the wrappers and package
  settings of the simulation study; covariates are removed from expression and genotypes before the summary-level
  and pairwise methods (they have no covariate argument); MCMC methods with the chain length of the MR.RGM fit.
* **Reporting:** causal edges with PIP > 0.5 (median probability model, as in the simulations) and confounding
  links with PIP > 0.5 (`ONEK_THR_DIR`, `ONEK_THR_CONF`), the same threshold in every table. Prior inclusion
  probabilities are 0.75 (causal) and 0.5 (confounding); `tables/onek_threshold_sweep.csv` gives edge counts and
  Bayes factors at several thresholds. Module panels are the Louvain communities (≥ 3 genes) of the displayed network.

Original OneK1K data were generated by the OneK1K Consortium; the Zenodo record redistributes derived data for
methodological reproducibility only.
