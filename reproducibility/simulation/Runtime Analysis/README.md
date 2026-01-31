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
- Each method is executed **multiple times** per setting
- The **median runtime** is reported
- Runtimes are converted to **seconds** for plotting

---




