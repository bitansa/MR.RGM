# Confounding Recovery

This folder contains code to evaluate **confounding recovery performance** for **MR.RGM** and **MR.RGM+** under multiple network-generating scenarios.

---

## Main script

### `ConfoundingRecoveryParallel.R`
Runs confounding-recovery simulations and produces:
1. **Performance tables** (mean ± SD) across replicates  
2. **AUC boxplots** (by sample size)

---

## What is being evaluated?

For each replicate, the script compares:

- **Estimated confounding structure**: `ZEst` (from MR.RGM / MR.RGM+)  
- **Ground-truth confounding structure**: `Sigma_True`

Only the **upper-triangular off-diagonal entries** are used for evaluation.

### Metrics
- **AUC**: computed using *continuous* scores (`ZEst_upper`) against a binary ground truth derived from `Sigma_True`
- **TPR / FDR / MCC**: computed using *thresholded* predictions (`ZEst_upper > 0.5`)

The script uses:
- `metrics_binary_graph()` for **TPR/FDR/MCC** 

---

## Simulation scenarios

### 1) Scale-free network (MR.RGM only)
- Network sizes: `p ∈ {5, 10}`
- Sample sizes: `n ∈ {500, 1000, 10000, 30000}`
- Output:
  - Summary table: mean & SD of (AUC, TPR, FDR, MCC) by `(p, n)`
  - Two AUC boxplots: one for `p=5`, one for `p=10`

Key function:
- `ScaleFree_Conf(run, Sample_Size, Network_Size)`

---

### 2) Small-world network (MR.RGM only)
- Network sizes: `p ∈ {5, 10}`
- Sample sizes: `n ∈ {500, 1000, 10000, 30000}`
- Output:
  - Summary table: mean & SD of (AUC, TPR, FDR, MCC) by `(p, n)`
  - Two AUC boxplots: one for `p=5`, one for `p=10`

Key function:
- `SmallWorld_Conf(run, Sample_Size, Network_Size)`

---

### 3) Small-world network with horizontal pleiotropy (MR.RGM vs MR.RGM+)
In this setting, both methods are run and compared:
- **MR.RGM**: uses structured `D`
- **MR.RGM+**: uses **full D** (`D = 1` matrix)

- Network sizes: `p ∈ {5, 10}`
- Sample sizes: `n ∈ {500, 1000, 10000, 30000}`
- Output:
  - Summary table: mean & SD of (AUC, TPR, FDR, MCC) by `(Method, p, n)`
  - Two **comparison** plots (MR.RGM vs MR.RGM+):
    - One plot for `p=5`
    - One plot for `p=10`
    - Each plot is faceted by sample size

Key function:
- `SmallWorld_Pleio_Conf(run, Sample_Size, Network_Size)`

---

