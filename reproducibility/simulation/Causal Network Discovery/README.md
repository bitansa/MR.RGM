# Causal Network Discovery — Simulation Reproducibility

This folder contains scripts to reproduce the **causal network recovery** simulation results
reported in the MR.RGM paper. The simulations evaluate graph recovery performance under
different network structures and data-generating assumptions.

Performance is assessed using the following metrics:
- Area Under the ROC Curve (**AUC**)
- True Positive Rate (**TPR**)
- False Discovery Rate (**FDR**)
- Matthews Correlation Coefficient (**MCC**)

---

## Contents of this folder

### 1. `GraphRecoverMetricsCalculation.R`

**Purpose:**  
Defines helper functions to compute graph recovery performance metrics.

- Computes: **AUC, TPR, FDR, MCC**
- Inputs:
  - True causal graph (binary vector or adjacency)
  - Estimated causal graph (binary or probabilistic)
- This file **must be sourced first** before running any simulation script.

All simulation scripts in this folder depend on the functions defined here.

---

### 2. Simulation scripts (performance computation)

These scripts generate simulated data, apply different causal discovery methods,
and compute graph recovery performance.

Each script allows the user to vary:
- **Network size:** `p = 5` or `p = 10`
- **Sample size:** `n ∈ {500, 1000, 10000, 30000}`

Parallel computation is used where applicable.

#### a) `ScaleFree_CausalNetwork_Parallel.R`

- Scenario: **Scale-free networks**
- Methods evaluated include:
  - MR.RGM
  - MR.RGM (no confounder)
  - Classical MR baselines (e.g., IVW, weighted median, simple median, OneSampleMR, MrBayes)
- Outputs performance metrics for causal network recovery.

---

#### b) `SmallWorld_CausalNetwork_Parallel.R`

- Scenario: **Small-world networks without horizontal pleiotropy**
- Same set of methods as the scale-free case
---

#### c) `SmallWorld_CausalNetwork_Pleiotropy_Parallel.R`

- Scenario: **Small-world networks with horizontal pleiotropy**
- Includes an additional method:
  - **MR.RGM+** (FullD variant)

---

### 3. `Network_Dicovery_Performnace_comparison_Plot.R`

**Purpose:**  
Generates comparison plots for causal network recovery performance.

- Produces boxplots (primarily AUC) across:
  - Methods
  - Sample sizes
  - Network sizes
- Handles:
  - Scale-free networks
  - Small-world networks without pleiotropy
  - Small-world networks with horizontal pleiotropy (including MR.RGM+)

This script assumes that the simulation scripts above have already been run and
that the corresponding result objects/files are available.

---

## Recommended execution order

1. **Set working directory** to the repository root.
2. **Source the metrics file**:
   - `GraphRecoverMetricsCalculation.R`
3. **Run one simulation script** depending on the desired scenario:
   - Scale-free
   - Small-world (no pleiotropy)
   - Small-world (with pleiotropy)
4. **Run the plotting script** to generate performance comparison figures.

---

## Notes

- All scripts are intended for **GitHub-based reproducibility**.
- Users may need to adjust parallel settings depending on available cores.
- Any locally hard-coded paths should be replaced with relative paths
  under the `reproducibility/` directory when running elsewhere.
