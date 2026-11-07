# Real-data analysis

MR.RGM (spike-and-slab graph prior, SSSL error covariance, measured covariates; each cis-eQTL instrument may affect
only its own gene and its effect is estimated without selection) applied to B-cell expression data from OneK1K,
with MR.RGM+ (every instrument may affect every gene, instrument effects selected) as a sensitivity analysis. With
cis-eQTL instruments, which act mainly on their own gene and explain little expression variance, MR.RGM+ is only
weakly identified (its MCMC diagnostics show this); see the methods paper.

```text
real_data/
├── real_data_functions.R   shared helpers: preprocessing, tables, network figures, modules, motifs, MCMC diagnostics
└── OneK1K/                 OneK1K B cells, KEGG B cell receptor signaling pathway (hsa04662)
```

## Requirements

* MR.RGM >= 0.1.3 (`remotes::install_github("bitansa/MR.RGM")`, or from CRAN); the scripts stop otherwise.
* R packages: `dplyr`, `tidyr`, `ggplot2`, `ggraph`, `igraph`, `tidygraph`, `scales`; `GenomicRanges` and `snpStats` (Bioconductor);
  `posterior` (convergence diagnostics of the chains);
  for `OneK1K_Competitors.R` and `../package_check/check_design_matched.R` also the packages of the simulation study
  (`pROC`, `MendelianRandomization`, `OneSampleMR`, `mrbayes`, `rjags` with JAGS, `MR2G`, `MrDAG`, `grivet`; see
  `../simulation/README.md`).
* Data: Zenodo record 10.5281/zenodo.18489687, file `donor_b_cell.rda` (md5 `63a105565f84a45cf654e7860db8b5e9`), and the
  array genotypes of the public OneK1K release (see `OneK1K/README.md`). The only path to set is the data folder at the
  top of the scripts (`data_dir`).

## Pipeline

1. **Data and instruments.** Donors with genotype and expression data; the genes of KEGG pathway hsa04662; for each
   gene the cis-SNPs with F > 10 given the covariates, LD-clumped (r² ≤ 0.5, within and between genes), at most 10
   per gene, each SNP assigned to one gene; genes whose best cis-SNP is not significant after Bonferroni correction
   over the cis-SNPs are removed; the m genes with the strongest instruments are kept, m = the largest gene set supported by the design-matched simulation (concentration parameter
   μ² = n R²/(1 − R²) and the joint first-stage F, reported for every gene).
2. **Preprocessing.** Expression library-size normalized and log-transformed; covariates SEX, AGE (decade bins) and
   PC1–PC3; every column of expression, genotypes and covariates centered and scaled (the model has no intercept,
   and the prior hyperparameters assume unit scale).
3. **Fit** (`OneK1K_ModelFitting.R`). MR.RGM (default) or MR.RGM+ (`ONEK_MODEL <- "MR.RGM+"`) from the sufficient
   statistics; 100,000 iterations, 20,000 burn-in, thinning 10 (8,000 posterior draws per chain; MR.RGM+ five times
   longer, 500,000 / 100,000 / 50), four chains
   with different seeds from the same starting point, checked for convergence (rank-normalized split-R̂ ≤ 1.01 for the
   log-likelihood and every edge indicator, bulk effective sample size, between-chain PIPs; package `posterior`) and
   pooled (32,000 draws); the results of a model whose chains do not converge are not interpreted (per-chain PIPs are
   reported instead); the package's warnings and diagnostics are recorded; saves
   `output/OneK1K_fit_MRRGM.rds` / `output/OneK1K_fit_MRRGMplus.rds` with the settings and package version.
   The number of genes is the largest gene set supported by the design-matched simulation
   (`../package_check/check_design_matched.R`, run from `package_check/` on the candidate genes written with
   `ONEK_STRENGTH_ONLY <- TRUE`).
4. **Every competing method on the same data** (`OneK1K_Competitors.R`): MR.RGM_NoConf, MVMR, MR2G, MrDAG (one gene
   at a time as outcome), GrIVET (on the original scale of the data, as in the simulations) and the pairwise methods
   (IVW, simple and weighted median, OneSampleMR, mrbayes), with the wrappers and package settings of the simulation
   study, on covariate-adjusted data; saves `output/OneK1K_competitors.rds`.
   **Prior sensitivity** (`OneK1K_PriorSensitivity.R`): MR.RGM refitted with a uniform Beta(1, 1) prior on the causal
   inclusion probabilities (same data, chain length and seeds); saves `output/OneK1K_fit_MRRGM_Beta11.rds` and
   `tables/onek_prior_sensitivity.{csv,tex}`.
   **Parametric bootstrap** (`OneK1K_Bootstrap.R`): the fitted MR.RGM network taken as the truth, expression generated
   100 times for the same donors from their genotypes and covariates, and MR.RGM refitted (simulation-study chain
   length); saves `output/OneK1K_bootstrap.rds` and `tables/onek_bootstrap{,_edges}.{csv,tex}` (network AUC, false
   discovery proportions, recovery of every edge and 2-cycle, warnings). About 7 minutes on 9 cores.
5. **Results** (`OneK1K_Results.R`, needs only the saved fits). Written to `tables/` and `figures/`:
   * convergence and identification diagnostics: causal PIPs from the two halves of the chain and from the second
     chain, acceptance rate, log-likelihood trace, first-stage F, error-variance ratios, spectral radius, and the
     package warnings that fired (`realdata_diagnostics.tex` when both fits exist)
   * all causal and confounding PIPs, edge counts with Bayes factors at several thresholds
   * reported causal edges and confounding links, label map, candidate-gene and instrument tables
   * network figure (plus one panel per data-driven module) and the most probable motif of each type
   * MR.RGM vs MR.RGM+ comparison (`onek_model_comparison.tex`) and every method on the same data
     (`onek_methods_comparison.tex`: edges, overlap with MR.RGM, sign agreement; `onek_edge_support.tex`: which
     methods report each edge of MR.RGM / MR.RGM+)

Reporting thresholds: causal PIP > 0.50 (median probability model, as in the simulations), confounding > 0.50;
prior inclusion probabilities are 0.75 (causal edge) and 0.5 (confounding link).

See `OneK1K/README.md` for the data layout, every setting and the exact commands.
