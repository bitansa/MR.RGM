# Reproducibility

This directory contains all code and instructions required to reproduce the results reported in the methods paper.

**Method names.** In the paper, the methods are called `IV-RGM` (the R package `MR.RGM` with a known
instrument-response map, `selectB = FALSE`), `IV-RGM+` (all instruments candidates for every response,
`selectB = TRUE`) and `IV-RGM_NoConf` (a diagonal error covariance). The scripts and saved result files use the
internal keys `MR.RGM`, `MR.RGM+` and `MR.RGM_NoConf`; tables and figures display the names of the paper
(function `display_text()`).

```text
reproducibility/
├── package_check/   checks of the installed MR.RGM on simulated data, and the design-matched simulation
├── simulation/      simulation study: all experiments, competitors, figures and tables
└── real_data/       OneK1K B-cell analysis: preprocessing, model fitting, figures and tables
```

**Package version.** The results in the paper were produced with MR.RGM 0.1.3
(`remotes::install_github("bitansa/MR.RGM")`, or from CRAN). Check with `packageVersion("MR.RGM")`; every script
stops with an older version (a 0.1.3 build is accepted only if `RGM()` has the argument `gibbsB`), and every script
saves the package version with its results.

**Software versions.** R 4.4.3; MendelianRandomization 0.10.0, OneSampleMR 0.1.8, mrbayes 0.5.3 (rjags 4.17),
MR2G 0.1.0, MR2 0.2.0, MrDAG 0.1.3 and grivet 1.0 (GitHub commit `1777dbf`). `simulation/cluster/install_packages.R`
installs MendelianRandomization, mrbayes, OneSampleMR and grivet at these versions; MR2G (GitHub) and MR2 and MrDAG
(r-universe) are installed at their current versions, so check `packageVersion()` against the list above.

**Data.** The real-data scripts read the preprocessed OneK1K B-cell data archived on Zenodo
(<https://doi.org/10.5281/zenodo.18489687>, file `donor_b_cell.rda`, md5 `63a105565f84a45cf654e7860db8b5e9`; check
with `tools::md5sum()`), and the array genotypes of the public OneK1K release (details in `real_data/OneK1K/README.md`).
Set the data folder at the top of the real-data scripts (`data_dir`); on a cluster, the paths and partitions are in
`simulation/cluster/config.sh`.

**Results.** The result files (per-replicate results, stored estimates, error logs, runtime benchmark, design-matched
simulation and OneK1K fits) are archived on Zenodo, <https://doi.org/10.5281/zenodo.24043120>, and are not in the
repository. To regenerate every table and figure without refitting any model, unzip

* `simulation_results.zip` into `simulation/` (it contains the `results_<experiment>/` folders),
* `runtime_results.zip` into `simulation/Runtime Analysis/`,
* `design_matched_results.zip` into `package_check/`,
* `onek1k_results.zip` into `real_data/OneK1K/output/`,

and run the figure and table scripts below (with `STRENGTH_PLOT_ONLY <- TRUE` and `RUNTIME_PLOT_ONLY <- TRUE` for the
instrument-strength and runtime figures).

## Run order

```r
## 1. package checks (about 10 minutes)
setwd("reproducibility/package_check")
source("check_recovery.R")            # all 24 sampler variants, MR.RGM and MR.RGM+, on one known network
## 2. simulations (simulation/README.md). On one machine this takes several weeks; the runtime benchmark alone takes
##    about two weeks, since its runs are sequential so that the timings are undisturbed. Recommended on a SLURM
##    cluster (simulation/cluster/README.md): one task per setting, the same result files.
##    On one machine, from a Terminal (not RStudio):  cd reproducibility/simulation; Rscript run_all.R > run_all.log 2>&1 &
##    run_all.R runs every experiment, the instrument-strength study, the figures and tables and the runtime benchmark,
##    without GrIVET; GrIVET is then added in a second pass with a limit of 48 hours per fit (see "GrIVET" below).
##    Before the run, check the competitor wrappers once:
setwd("../simulation")
source("pilot_competitors.R")
## 3. real data (real_data/OneK1K/README.md)
setwd("../real_data/OneK1K"); data_dir <- "..."
source("OneK1K_PreProcessing.R")
ONEK_STRENGTH_ONLY <- TRUE; source("OneK1K_ModelFitting.R")    # candidate genes and instrument strengths, no fit
## 4. design-matched simulation: the largest gene set that MR.RGM recovers at the observed instrument strengths
setwd("../../package_check"); source("check_design_matched.R")  # prints "largest supported gene set = 10"
## 5. real-data fits with that gene set, every competing method on the same data, tables and figures
setwd("../real_data/OneK1K"); rm(ONEK_STRENGTH_ONLY); ONEK_N_GENES <- 10
source("OneK1K_RunAll.R")                                       # MR.RGM and MR.RGM+ (four chains each), competitors,
                                                                #   prior sensitivity, parametric bootstrap, results
```

## package_check/

* `check_recovery.R` fits all 24 sampler variants (two graph priors × three error models × with/without
  instruments × with/without covariates) to one simulated data set with a known cyclic network and reports AUC and
  mean absolute error of the causal effects, for MR.RGM and (with instruments) MR.RGM+ (about 10 minutes).
* `check_design_matched.R`: simulations with the instrument strengths and n of the OneK1K candidate genes (the 14
  genes whose strongest cis-SNP passes the Bonferroni check), for gene sets of increasing size in decreasing order of
  the concentration parameter, m in {5, 8, 10, 12, 14}, with 100 replicates and the MCMC settings of the simulation
  study. The largest gene set that meets the pre-specified criteria for itself and every smaller set (mean AUC ≥ 0.80,
  mean chain-half PIP correlation ≥ 0.90, error-variance warning in at most 10% of the fits) is the gene set of the
  real-data analysis. Run it from `package_check/`; it writes `design_matched_OneK1K.tex` there.

## simulation/

Seven experiments: the six of `run_simulations.R` (main cases 1–3; error variance; strength of confounding; large n;
hyperparameter sensitivity; exposure → outcome design with known roles) and the instrument-strength study of
`run_instrument_strength.R`, together with the runtime benchmark. They compare MR.RGM, MR.RGM_NoConf and MR.RGM+ with
MR2G, MVMR, MR2, MrDAG, GrIVET and the pairwise estimators (IVW, simple and weighted median, OneSampleMR, mrbayes),
every competitor with its package's default settings except the √n input scaling and chain length of MR2 and MrDAG,
the tuning grid and unstandardized input of GrIVET, and a correction of MR2G's ScreenAug step for igraph ≥ 2.0 (all
explained in `simulation/README.md`).

## real_data/

The OneK1K B-cell analysis of the KEGG B cell receptor signaling pathway: data download, preprocessing, instrument
selection (LD clumping), the number of genes from the design-matched simulation, MR.RGM and MR.RGM+ fits, every
competing method on the same data, posterior inclusion probabilities, network figures, motifs and diagnostics.
See `real_data/README.md` and `real_data/OneK1K/README.md`.

## Notes

- All scripts are written in R. Execution order matters and is documented in the README of each folder.
- Simulations run the replicates in forked R processes (`run_parallel()` in `simulation/sim_functions.R`, macOS and
  Linux); on Windows, or with `N_CORES <- 1`, they run sequentially (slow).
- Every run is restartable and saves its settings and session information with the results.

## Where each result of the paper comes from

| Paper | Output file | Script |
|---|---|---|
| Figure 2 (main study, p = 10); main-study figures and tables of Sections OA.2–OA.3 | `main_summary_p=10.pdf`, `feedback_*.tex`, `overview_*.tex`, `failures.tex`, `main_diagnostics.tex` | `simulation/make_figures_tables.R` |
| Further experiments (exposure-outcome, error variance, confounding, large samples, prior sensitivity) | `roles*`, `roles_auc.pdf`, `noise*`, `confscale*`, `largen*`, `largen_seconds.tex`, `sensitivity*` | `simulation/make_figures_tables.R` |
| Mean time per fit in the main study | `timing.tex` | `simulation/make_figures_tables.R` |
| Figure 3 and the instrument-strength study | `strength_*.pdf`, `strength_*.tex` | `simulation/run_instrument_strength.R` |
| Computation time against the number of responses | `runtime_vs_p_n30000.pdf` | timed by `simulation/cluster/run_runtime_task.R` and collected by `collect_runtime.R`, drawn by `simulation/Runtime Analysis/RunTimeAnalysis.R` |
| Choice of the number of genes (design-matched simulation) | `design_matched_OneK1K.tex` (in `package_check/`) | `package_check/check_design_matched.R` |
| Figure 4 and the OneK1K results of Section OA.6 | `onek1k_network.pdf`, `onek1k_cluster_*.pdf`, `motif_onek_all.pdf`, `onek_loglik_trace.pdf`, `onek_plus_loglik_trace.pdf`, `onek_chains_PIPs.pdf`, `onek_plus_chains_PIPs.pdf`, `onek_all_edges.tex`, `onek_candidates.tex`, `onek_instruments*.tex`, `onek_methods_comparison.tex`, `onek_edge_support.tex`, `realdata_diagnostics.tex` | `real_data/OneK1K/OneK1K_RunAll.R` (fits, competitors, results) |
| Parametric bootstrap of the OneK1K network | `onek_bootstrap.tex`, `onek_bootstrap_edges.tex` | `real_data/OneK1K/OneK1K_Bootstrap.R` |
| Prior sensitivity of the OneK1K fit | `onek_prior_sensitivity.tex` | `real_data/OneK1K/OneK1K_PriorSensitivity.R` |

The figures in the paper are the figures written by these scripts, redrawn from the saved results at their printed
size (base font size 9 pt); no model is refitted for that step.

**GrIVET.** Every fit of GrIVET ends in one of four ways: a result; a failure of its peeling step; a numerical error;
or no result within 48 hours. GrIVET is fitted in a second pass that adds it to the existing result files with a limit
of 48 hours per fit (the first pass skips it through `SKIP_METHODS`, since a replicate that exceeds its time limit is
stopped with every method in it):

* on a cluster: `REP_TIMEOUT_MIN=2880 REP_RETRIES=0 bash cluster/submit_add.sh` (from `simulation/`); the 48-hour job
  limit (`SIM_TIME` in `cluster/config.sh`) also stops the fits at n = 10^6;
* on one machine: `REP_TIMEOUT_MIN <- 720; REP_RETRIES <- 0; source("cluster/add_methods.R")` for n = 10^6
  (`add_methods.R` allows four times `REP_TIMEOUT_MIN` at n ≥ 10^6), and `REP_TIMEOUT_MIN <- 2880` for the other
  result files (`ADD_FILES`).

A fit without a result at the limit is recorded as not finished. The tables report the number of fits of each kind
(`failures.tex`, `strength_failures.tex`); the causes are read from the error logs `results_*/errors.log`.
