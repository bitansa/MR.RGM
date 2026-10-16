# Real-data analyses

MR.RGM+ (spike-and-slab graph prior, SSSL error covariance, measured covariates, every instrument allowed
for every gene) applied to two expression data sets.

```text
real_data/
├── real_data_functions.R   shared helpers: preprocessing, tables, network figures, modules, motifs, MCMC diagnostics
├── GTEx/                   GTEx v7 skeletal muscle, mTOR pathway (18 genes)       -> Main text
└── OneK1K/                 OneK1K B cells, B-cell receptor pathway (66 genes)     -> Supplement
```

## Requirements

* MR.RGM >= 0.1.2 (`remotes::install_github("bitansa/MR.RGM")` until 0.1.2 is on CRAN); the scripts stop otherwise.
* R packages: `dplyr`, `tidyr`, `ggplot2`, `ggraph`, `igraph`, `tidygraph`, `scales`;
  GTEx also `R.matlab`, `readr`; OneK1K also `GenomicRanges` (Bioconductor).
* Data: Zenodo records listed in each subfolder README (GTEx: 10.5281/zenodo.18476622; OneK1K: 10.5281/zenodo.18489687).
  The only path to set is the data folder at the top of each script (`GTEX_FOLDER_PATH`, `data_dir`).

## Pipeline (identical for both data sets)

1. **Data and instruments.** Donors with genotype and expression data; pathway genes; top cis-eQTL variants per gene.
2. **Preprocessing.** Expression log-transformed; covariates SEX (male = 1), AGE (ordinal decade bins)
   (+ PC1–PC3 for OneK1K); every column of expression, genotypes and covariates centered and scaled
   (the model has no intercept, and the prior hyperparameters assume unit scale).
3. **Fit** (`*_ModelFitting.R`). MR.RGM+ from summary statistics; 50,000 iterations, 10,000 burn-in,
   thinning 10 (4,000 posterior draws); saves `output/*_fit.rds` with the settings and package version.
4. **Results** (`*_Results.R`, needs only the saved fit). Written to `tables/` and `figures/`:
   * MCMC diagnostics: causal PIPs from the two halves of the chain, acceptance rates, log-likelihood trace
   * all causal and confounding PIPs, and edge counts with Bayes factors at several thresholds
   * reported causal edges (explicit direction) and confounding links, label map (node label -> gene symbol)
   * network figure (OneK1K also one panel per data-driven module)
   * every feedback loop, feed-forward loop and cascade among the reported edges, scored with `NetworkMotif()`,
     and a figure of the most probable motif of each type

Reporting thresholds: GTEx causal PIP >= 0.85, confounding >= 0.50; OneK1K causal >= 0.50, confounding >= 0.50.
Prior inclusion probabilities are 0.75 (causal edge) and 0.5 (confounding link).

See `GTEx/README.md` and `OneK1K/README.md` for the data layout and the exact commands.
