# Simulation study

This folder reproduces every simulation figure and table in the manuscript and its Supplement.

| File | Purpose |
|---|---|
| `sim_functions.R` | Settings, data generation, competing methods, MR.RGM fits, metrics |
| `run_simulations.R` | Runs all settings (20 replicates each, in parallel); one result file per setting in `results/` |
| `make_figures_tables.R` | Builds every figure (`figures/`) and Supplement table (`tables/*.tex`) from `results/` |
| `Runtime Analysis/RunTimeAnalysis.R` | Runtime benchmark (scale-free, n = 30,000, p ∈ {2, 5, 10, 20}, 20 runs per method) |

## Requirements

* MR.RGM >= 0.1.2 (`remotes::install_github("bitansa/MR.RGM")` until 0.1.2 is on CRAN); the scripts stop otherwise.
* R packages: `MASS`, `igraph`, `pROC`, `MendelianRandomization`, `OneSampleMR`, `mrbayes`, `rjags`
  (+ `ggplot2`, `scales` for the figures; `microbenchmark`, `ggplot2` for the runtime analysis).
* `rjags` needs the JAGS program itself (https://mcmc-jags.sourceforge.io), installed before `install.packages("rjags")`.
* Parallel runs use `parallel::mclapply` (forking): on Windows set `N_CORES <- 1`.

## How to run

```r
setwd("path/to/MR.RGM/reproducibility/simulation")

SMOKE_TEST <- TRUE;  source("run_simulations.R")    # a few minutes: every method and metric once
SMOKE_TEST <- FALSE; source("run_simulations.R")    # full study (restartable: finished settings are skipped)
source("make_figures_tables.R")                     # figures/ and tables/
source("Runtime Analysis/RunTimeAnalysis.R")        # runtime table and figure
```

Set `N_CORES` before sourcing to choose the number of cores (default: all but one). A subset can be run with
`SCENARIOS`, `P_VALUES`, `N_VALUES` (e.g. `SCENARIOS <- "scalefree"`).

## Design

* **Scenarios.** Case 1: scale-free network (Barabási–Albert tree, every edge reciprocal). Case 2: small-world
  network (Watts–Strogatz, neighbourhood 2, rewiring probability 0.1; each edge oriented at random). Case 3: as
  Case 2 plus horizontal pleiotropy (one shared instrument per consecutive trait pair).
* **Sizes.** p ∈ {5, 10}; n ∈ {500, 1000, 10000, 30000}; 20 replicates per setting.
* **Truth.** Non-zero entries of **A** are ±0.1 (equal probability) on the edges of the graph; each trait has 3
  unique instruments with effect 1; 5 latent confounders with loadings ±1, half of which are set to 0.
  In Case 3 the instrument map given to every method assigns each shared instrument to the first of its two
  traits (traits are arbitrarily ordered); **MR.RGM+** instead allows every instrument for every trait.
* **Data.** X, U ~ N(0, 1); E ~ N(0, 9 I); **Y** = (I − **A**)⁻¹(**B X** + **C U** + **E**).
* **Replicates.** In each replicate a new graph, **A** and **C** are drawn (seed = replicate number, so a
  replicate has the same true network at every n), and the data are drawn with a separate seed unique to the
  setting. All methods are fitted to the same data set.
* **Methods.** MendelianRandomization (simple median, weighted median, IVW; 10,000 bootstrap draws),
  OneSampleMR (two-stage residual inclusion), mrbayes (Bayesian IVW; 2 chains, 10,000 burn-in, 50,000 iterations),
  MR.RGM (spike-and-slab, SSSL error covariance), MR.RGM_NoConf (spike-and-slab, diagonal error covariance) and,
  in Case 3, MR.RGM+ (all instruments allowed for all traits). MR.RGM family: 50,000 iterations, 10,000 burn-in,
  thinning 10. Competitors are pairwise (exposure *i* → outcome *j* using the instruments of *i*).
* **Metrics.** Graph recovery: AUC (scores 1 − p-value or posterior inclusion probability), TPR, FDR, MCC
  (edge if p < 0.05, or PIP > 0.5). Effect estimation: maximum, mean absolute and mean squared deviation over the
  off-diagonal entries of **A**. Confounding: truth = |Σ\*| (off-diagonal), divided by its maximum and thresholded
  at its mean; AUC and TPR/FDR/MCC (link if PIP > 0.5). Instrument selection (MR.RGM+): AUC of the inclusion
  probabilities against **B** ≠ 0. Tables report mean (SD) over replicates.

`run_simulations.R` ends with a check that every method has 20 valid replicates in every setting and lists any
pairwise fit that errored.
