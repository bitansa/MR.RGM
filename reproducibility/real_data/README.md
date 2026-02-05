# Reproducibility: Real Data Analyses

This directory contains reproducible workflows for applying **MR.RGM / MR.RGM+** to real biological datasets.  
It is organized into two dataset-specific subfolders, each of which includes data-download instructions, preprocessing, model fitting, and post-analysis visualization.

---

## Directory structure

```text
realdata/
├── GTEx/
└── OneK1K/
```


Each folder is **self-contained** and can be run independently.

---

## Overview

Both **GTEx** and **OneK1K** workflows follow the same high-level pipeline:

1. **Data acquisition & preprocessing**
   - Download raw data from a public repository
   - Perform quality control and normalization
   - Construct inputs required by MR.RGM / MR.RGM+

2. **Model fitting**
   - Fit the MR.RGM+ model
   - Estimate causal effects, confounding structure, and covariance components
   - Store posterior samples and posterior inclusion probabilities (PIPs)

3. **Post-analysis**
   - Query PIPs for:
     - directed causal edges
     - undirected confounding links
   - Visualize inferred gene networks
   - Plot biologically motivated subnetworks

4. **Network motif analysis**
   - Apply the `NetworkMotif()` function in **MR.RGM**
   - Compute posterior probabilities for specific motifs
     (e.g., feedback loops, feedforward loops, cascades)

---

## GTEx folder

The `GTEx/` directory contains scripts and instructions to reproduce analyses on GTEx gene expression data.

Typical contents include:
- data download instructions for GTEx
- preprocessing scripts
- MR.RGM / MR.RGM+ model fitting
- extraction of causal and confounding PIPs
- network plots and motif probability calculations

Please see `GTEx/README.md` for dataset-specific details and step-by-step instructions.

---

## OneK1K folder

The `OneK1K/` directory contains scripts and instructions to reproduce analyses on the OneK1K B-cell dataset.

This workflow includes:
- downloading the OneK1K data from Zenodo
- preprocessing donor, genotype, and expression data
- fitting MR.RGM+ to infer causal and confounding gene networks
- querying posterior inclusion probabilities (PIPs)
- plotting subnetworks used in the manuscript
- computing posterior probabilities of predefined network motifs using `NetworkMotif()`

Please see `OneK1K/README.md` for full instructions, script order, and expected outputs.

---

## Key outputs across both datasets

Running the full pipelines produces:
- fitted MR.RGM+ model objects
- posterior inclusion probabilities (PIPs) for:
  - causal edges
  - confounding links
- network and subnetwork visualizations
- posterior probabilities for biologically relevant network motifs

---

## Notes

- Each dataset folder contains its **own README** with precise instructions.
- Script execution order is important and documented within each folder.
- Some steps (e.g., cis-SNP selection) may be computationally intensive.

---

This structure is designed to make all real-data results in the paper fully reproducible, from raw data download to final figures and motif probability estimates.
