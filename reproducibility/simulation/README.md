# Simulation Studies — Reproducibility

This directory contains all **simulation studies** used to reproduce the numerical results and figures in the paper.  
Each subfolder corresponds to a **distinct evaluation task**, focusing on different aspects of model performance.

All simulations are implemented in **R**, use **synthetic data generated under controlled network structures**, and are designed to be **fully reproducible**.

---

## Folder Structure

The `simulation/` directory contains the following five subfolders:

 - `Causal Effect Estimation`
 - `Causal Network Discovery`
 - `Confounding Recovery`
 - `IV-Response Recovery Performance`
 - `Runtime Analysis`


Each folder includes:
- One or more **R scripts** for running simulations
- A **folder-specific README** describing the scripts in detail
- Code to **store results** and **generate manuscript figures**

This top-level README provides a **conceptual overview only**.

---

## 1. Causal Network Discovery

**Purpose:**  
Evaluate how well different methods recover the **causal network structure** among response variables.

**Focus:**  
- Graph recovery performance under different network types  
- Comparison across multiple MR-based and Bayesian methods

**Typical outputs:**  
- AUC, TPR, FDR, MCC  
- Boxplots and comparison figures across sample sizes and network sizes

---

## 2. Causal Effect Estimation

**Purpose:**  
Assess the accuracy of **causal effect estimation** between response variables.

**Focus:**  
- Estimation error under different data-generating mechanisms  

**Typical outputs:**  
- Maximum Absolute Deviation  
- Mean Absolute Deviation  
- Mean Squared Error  
- Summary tables and manuscript-ready plots

---

## 3. Confounding Recovery

**Purpose:**  
Evaluate how well models recover **latent confounding structure**.

**Focus:**  
- Confounding recovery using posterior dependence estimates  
- Comparison between:
  - **MR.RGM**
  - **MR.RGM+**

**Typical outputs:**  
- AUC, TPR, FDR, MCC for confounding recovery  
- Plots comparing performance across sample sizes and network sizes

---

## 4. IV–Response Recovery Performance

**Purpose:**  
Assess recovery of **instrument–response (IV–trait) links**.

**Focus:**  
- Instrument selection accuracy under **horizontal pleiotropy**
- Evaluation of **MR.RGM+** using posterior IV weights

**Typical outputs:**  
- AUC for IV recovery  
- Boxplots of AUC vs sample size  
- Separate figures for different network sizes

---

## 5. Runtime Analysis

**Purpose:**  
Compare **computational runtime** across different MR methods.

**Focus:**  
- Scalability with respect to the number of response variables (**p**)  
- Fixed sample size benchmarking

**Methods compared:**  
- MR.RGM  
- MR.RGM+  
- MR.RGM_NoConf  
- MendelianRandomization  
- OneSampleMR  
- mrbayes  

**Typical outputs:**  
- Runtime tables (median runtime)  
- A manuscript-ready runtime comparison plot

---

## General Notes

- All simulations use **synthetic data** with known ground truth.
- Parallel computation is used where appropriate.
- Random seeds are controlled at the replicate level.
- Folder-level READMEs provide **exact execution order and parameter details**.

For reproducing specific figures or tables, please refer to the **README inside the corresponding subfolder**.

---

## Intended Use

This directory is designed to support:
- **Full reproducibility** of simulation results
- **Transparent comparison** of methods
- Direct generation of **manuscript figures**

No additional preprocessing or external datasets are required.
