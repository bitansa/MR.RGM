# Reproducibility

This directory contains all code and instructions required to reproduce the results reported in this study.

It is organized into two main components:

```text
reproducibility/
├── simulation/
└── real_data/
```

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
