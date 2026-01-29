# Causal Effect Estimation — Simulation Reproducibility

This folder contains scripts to reproduce the **causal effect estimation accuracy** simulation results reported in the MR.RGM paper.  
The simulations evaluate how accurately different Mendelian Randomization (MR) methods estimate *pairwise causal effects* under different network structures and data-generating assumptions.

Performance is assessed using the following **causal-effect error metrics** (computed off-diagonal only):
- **Maximum Absolute Deviation** (MaxAD)
- **Mean Absolute Deviation** (MAD)
- **Mean Squared Error** (MSE)

---

## Contents of this folder

### 1. Simulation scripts (metric computation)

These scripts generate simulated data, apply different causal effect estimation methods, and compute **MaxAD / MAD / MSE** of the estimated causal effect matrix against the true causal effect matrix.

Each script allows the user to vary:
- **Network size:** `p = 5` or `p = 10`
- **Sample size:** `n ∈ {500, 1000, 10000, 30000}`

Parallel computation is used where applicable.

#### a) `ScaleFree_CausalEffect_Parallel.R`

- Scenario: **Scale-free networks**
- Purpose:
  - Simulate a scale-free causal graph and generate genotype/exposure/outcome data
  - Estimate causal effects using multiple MR methods
  - Compute causal-effect accuracy metrics: **MaxAD, MAD, MSE**
  - Print summary **means and SDs** across simulation replicates
- Methods evaluated typically include:
  - MR.RGM
  - MR.RGM (No Conf)
  - MR-SimpleMedian
  - MR-WeightedMedian
  - MR-IVW
  - OneSampleMR
  - MrBayes

> Output objects:
> Each method produces a matrix `T_*` where columns are:
> `col 1 = Maximum Absolute Deviation`, `col 2 = Mean Absolute Deviation`, `col 3 = Mean Squared Error`.

---

#### b) `SmallWorld_CausalEffect_Parallel.R`

- Scenario: **Small-world networks without horizontal pleiotropy**
- Same purpose and outputs as the scale-free case, but under a small-world causal graph structure.
- Methods evaluated are the same as above (no MR.RGM+ in this non-pleiotropic setting).

---

#### c) `SmallWorld_CausalEffect_Pleiotropy_Parallel.R`

- Scenario: **Small-world networks with horizontal pleiotropy**
- In addition to the baseline methods, this scenario includes:
  - **MR.RGM+ (Full-D variant)** to handle horizontal pleiotropy
- Outputs the same causal-effect accuracy metrics (**MaxAD, MAD, MSE**) and summary statistics.

> Output objects:
> In this pleiotropic setting you may additionally see:
> - `T_MR.RGM_FullD` (treated as **MR.RGM+** in plotting)

---

### 2. `CausalEffectEstimation_Performance_Comparison_Plot.R`

**Purpose:**  
Generates comparison plots for **causal effect estimation accuracy**, focusing primarily on **MAD (Mean Absolute Deviation)**.

- Produces boxplots comparing methods across:
  - Sample sizes
  - Network sizes
  - Scenarios:
    - Scale-free
    - Small-world (no pleiotropy)
    - Small-world (with horizontal pleiotropy; includes MR.RGM+)
- This script assumes that you have already run one of the simulation scripts above and that the corresponding `T_*` result objects have been loaded (or read from saved files).

**Key detail:**  
The plotting script uses **column 2** of each `T_*` matrix:
- `T_*[, 2]` = **Mean Absolute Deviation (MAD)**

It also includes example code to:
- store simulation outputs (e.g., `saveRDS(...)`)
- reload stored results (e.g., `readRDS(...)`)
- append results across sample sizes / methods
- generate the final figures

---

## Recommended execution order

1. **Set working directory** to the repository root.
2. **Run one simulation script** depending on the desired scenario:
   - `ScaleFree_CausalEffect_Parallel.R`
   - `SmallWorld_CausalEffect_Parallel.R`
   - `SmallWorld_CausalEffect_Pleiotropy_Parallel.R`
3. **Save or load the simulation outputs** (`T_*` objects) as needed.
4. **Run the plotting script**:
   - `CausalEffectEstimation_Performance_Comparison_Plot.R`
   to generate MAD comparison boxplots.

---

## Notes

- Each simulation script computes error metrics comparing the **true causal effect matrix** to the **estimated causal effect matrix**, excluding diagonal entries.
- Parallel settings (number of cores, cluster backend) may need to be adjusted depending on your machine/HPC environment.
- Any locally hard-coded paths should be replaced with **relative paths** under the `reproducibility/` directory to ensure portability.
- The plotting script is intentionally modular:
  - run simulations → save results → reload/append → plot
  which makes it easy to regenerate figures without rerunning computationally expensive simulations.
