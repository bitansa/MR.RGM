# Reproducibility

This directory contains all code and instructions required to reproduce the results reported in this study.

It is organized into two main components:

```text
reproducibility/
├── package_check/   quick check that the installed MR.RGM recovers a known network (run first)
├── simulation/
└── real_data/
```

**Package version.** The results in the manuscript (Version 4) were produced with MR.RGM 0.1.2
(`remotes::install_github("bitansa/MR.RGM")` until 0.1.2 is on CRAN). Check with `packageVersion("MR.RGM")`.
Every script saves the package version with its results.

**Data.** The real-data scripts read the data from the Zenodo records given in each subfolder README;
set the data folder at the top of each script (`GTEX_FOLDER_PATH`, `data_dir`). Nothing else is path-specific.

---

## package_check/

`check_recovery.R` fits all 24 sampler variants (two graph priors x three error models x with/without
instruments x with/without covariates) to one simulated data set with a known cyclic network and reports
AUC and mean absolute error of the causal effects (about 10 minutes).

---

## simulation/

This folder contains scripts for **simulation studies** used to evaluate the performance of the proposed methods under controlled settings.  
It includes code to generate synthetic data, run **MR.RGM / MR.RGM+**, and compare their performance against existing methods such as **Inverse-Variance Weighting (IVW)**, **Simple Median**, **Weighted Median**, **One-Sample MR**, and **MR-Bayes** across a range of simulation scenarios.

The simulation framework is designed to assess accuracy, robustness to confounding and pleiotropy, and overall estimation performance of **MR.RGM / MR.RGM+** relative to these baseline approaches.

Please refer to the README inside the `simulation/` folder for details on how to run each experiment and reproduce the corresponding figures and tables.

---

## real_data/

This folder contains **real data analyses**.  
It includes two dataset-specific subfolders (e.g., GTEx and OneK1K), each providing:

- data download instructions  
- preprocessing scripts  
- **MR.RGM / MR.RGM+** model fitting  
- extraction of posterior inclusion probabilities (PIPs) for causal and confounding effects  
- network visualization and motif analysis

Each subfolder is self-contained and includes its own README with step-by-step instructions.

---

## Notes

- All scripts are written in R unless otherwise noted.
- Execution order matters and is documented in the README files within each subfolder.
- Some analyses may be computationally intensive.

This structure ensures full reproducibility of both the simulation results and the real-data analyses presented in the paper.
