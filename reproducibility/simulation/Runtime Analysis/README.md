# Runtime Analysis — `RunTimeAnalysis.R`

## Overview
`RunTimeAnalysis.R` benchmarks the **computational runtime** of six Mendelian randomization methods as the number of response variables (**p**) increases, while keeping the **sample size (n)** fixed.

The script produces:
1. A **runtime table** (median runtime across repetitions)
2. A **plot** showing median runtime versus `p`

---

## Methods Benchmarked
The following **six methods** are evaluated:

- **MR.RGM**
- **MR.RGM+** (Full-D variant)
- **MR.RGM_NoConf**
- **MendelianRandomization**
- **OneSampleMR**
- **mrbayes**

Each method is executed through a **wrapper function** to ensure consistent benchmarking.

---

## Simulation Setup
For each value of `p`, the script:
- Simulates a **scale-free causal network**
- Generates instrument variables (3 per response)
- Adds latent confounders
- Simulates genotype (`X`) and response (`Y`) data
- Computes summary statistics:
  - `Syy`, `Syx`, `Sxx`

The same simulated dataset is used across all methods for fair runtime comparison.

---

## Runtime Measurement
- Runtime is measured using `microbenchmark()`
- Each method is executed **20 times** per setting (`RUNTIME_TIMES`)
- The **median runtime** (seconds) is reported
- Data generator and method settings are those of the simulation study (`../sim_functions.R`):
  n = 30,000; MR.RGM family 50,000 iterations, 10,000 burn-in, thinning 10; mrbayes 2 chains,
  10,000 burn-in, 50,000 iterations; MendelianRandomization 10,000 bootstrap draws

## Usage
From the `simulation` folder: `source("Runtime Analysis/RunTimeAnalysis.R")`.
Outputs `Runtime Analysis/runtime_vs_p_n30000.csv` and `figures/runtime_vs_p_n30000.pdf`.

---




