# OneK1K (B cells) — Reproducibility scripts (MR.RGM+)

This folder contains the complete workflow to reproduce the **OneK1K B-cell** real-data analysis results reported in the manuscript, including:
- preprocessing the Zenodo dataset,
- fitting **MR.RGM+**,
- querying/plotting posterior inclusion probabilities (PIPs) for **causal** and **confounding** links, and
- computing posterior probabilities of selected **network motifs**.

---

## Contents

- `OneK1K_PreProcessing.R`  
  Loads the OneK1K B-cell dataset (`donor_b_cell_bitan.rda`) and performs:
  - promoter-region construction (±200 kb),
  - donor filtering using library-size outlier rules,
  - library-size normalization of expression,
  - conversion of genotype strings to numeric dosages,
  - SNP filtering using a non-zero dosage frequency threshold.

- `OneK1K_ModelFitting.R`  
  Runs **MR.RGM+** on the processed dataset:
  - selects cis-SNPs (top up to 15 per gene by association p-value),
  - constructs matrices **Y** (expression), **X** (genotypes), **U** (covariates),
  - computes summary-statistic matrices (Syy, Syx, Sxx, Syu, Sxu, Suu),
  - fits the model via `MR.RGM::RGM()` and creates `Output_OneK1K`.

- `OneK1K_Plot.R`  
  Uses the fitted output to:
  - extract/query **PIPs** for any **causal edge** (`GammaEst[to, from]`) and any **confounding link** (`ZEst[g1, g2]`),
  - plot the **9 cluster subnetworks** shown in the manuscript figure using default thresholds  
    **causal PIP ≥ 0.50** and **confounding PIP ≥ 0.48** (editable in the script).

- `OneK1K_NetworkMotif.R`  
  Plots three canonical motifs (saved as PDFs) and computes the posterior probability of each motif using  
  `MR.RGM::NetworkMotif()` with posterior samples `Output_OneK1K$GammaPst`:
  - a 3-cycle feedback loop,
  - a feedforward loop,
  - a 3-node cascade.

---

## Data download (Zenodo)

The OneK1K B-cell dataset used by these scripts is available on Zenodo:

- DOI: **10.5281/zenodo.18489688**  
- Download link (DOI landing page): **https://doi.org/10.5281/zenodo.18489688**

Download `donor_b_cell_bitan.rda` from the Zenodo record and save it into a local folder on your machine.

---

## How to run

### 0) Set your local data path (required)
In `OneK1K_PreProcessing.R`, set:

```r
# Replace with the local path where the Zenodo files were saved
data_dir <- "PATH_TO_ZENODO_DATA_FOLDER"
setwd(data_dir)
```

That folder **must contain**:

- `donor_b_cell_bitan.rda`

### 1) Run preprocessing

```r
source("OneK1K_PreProcessing.R")
```

This creates key objects used downstream, including:

- `RNA.count.adj` 
- `genotype.mat`  
- `Subsetted_SNP_Names`  
- `gene.promoter.ref`  
- `updated donor`, `vcf.ref`, etc.

---

### 2) Fit MR.RGM+  

```r
source("OneK1K_ModelFitting.R")
```

This fits the model and creates:

- `Output_OneK1K`  
  (contains `AEst`, `GammaEst`, `ZEst`, `SigmaEst`, and posterior samples such as `GammaPst`)
- `gene_names_pathway_use`  
  (the final gene list used for plotting/motifs)

---

### 3) Query PIPs + plot 9 subnetworks 

```r
source("OneK1K_Plot.R")
```

This script:

- shows examples for extracting PIPs of a causal edge and a confounding pair,  
- generates the 9 cluster plots (thresholds can be changed inside the script).

**Interpretation reminder**

- Causal edge PIP: `GammaEst[to, from]` corresponds to `from → to`  
- Confounding link PIP: `ZEst[g1, g2]` (symmetric / undirected)

---

### 4) Plot motifs + compute motif posterior probabilities 

```r
source("OneK1K_NetworkMotif.R")
```

This writes three PDF files in the working directory:

- `motif_onek_feedback.pdf`  
- `motif_onek_feedforward.pdf` 
- `motif_onek_cascade.pdf`  

and prints the posterior probabilities of these motifs (computed from `GammaPst`).

---

**Notes**

- **Order matters:**  
  `OneK1K_PreProcessing.R` and `OneK1K_ModelFitting.R` must be run before  
  `OneK1K_Plot.R` and `OneK1K_NetworkMotif.R`.

- **Packages:**  
  These scripts require R packages including  
  `MR.RGM`, `GenomicRanges`, and (for plotting)  
  `ggraph`, `tidygraph`, `scales`, `ggplot2`, `dplyr`, `igraph`.

- **Runtime:**  
  The cis-SNP selection step in `OneK1K_ModelFitting.R` can be slow because it fits many linear models (up to 15 per gene across genes).

---

**Expected outputs**

After running all scripts, you should have:

- a fitted model object: `Output_OneK1K`
- PIP querying + cluster network plots  
  (displayed / optionally saved if you set `out_file`)
- motif PDFs:
  - `motif_onek_feedback.pdf`  
  - `motif_onek_feedforward.pdf`  
  - `motif_onek_cascade.pdf`  
- printed motif posterior probabilities from `NetworkMotif()`
