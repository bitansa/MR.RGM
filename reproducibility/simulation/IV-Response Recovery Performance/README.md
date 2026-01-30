# IV Recovery (Instrument–trait selection) — Small-world + Horizontal Pleiotropy (MR.RGM+)

This code evaluates **instrument–trait (IV–exposure) recovery** performance under a **small-world network** with **horizontal pleiotropy**, using **MR.RGM+** (i.e., MR.RGM with **Full D = all ones**).

## What the simulation does

For each combination of:

- **Network size**: `p ∈ {5, 10}`
- **Sample size**: `n ∈ {500, 1000, 10000, 30000}`
- **Replicates**: `n_runs` (e.g., 20)

the script runs a synthetic data-generating process and computes an **AUC** score for recovering the true **instrument–trait links**.

### Data generation (per replicate)

For a given `(p, n)` and replicate seed `run`:

1. **Generate a small-world graph** on `p` nodes and convert it to a **directed** graph.
2. Construct a sparse **structural matrix** `A` aligned with the directed adjacency (edge weights ±0.1).
3. Construct the **true instrument–trait matrix `B`**:
   - Start with **3 instruments per trait** (block-diagonal structure).
   - Add **pleiotropic instruments** that affect *two traits* at once.
4. Generate **confounders** via a sparse matrix `C` and latent factors `U`.
5. Simulate:
   - Instruments `X`
   - Traits `Y` from the linear SEM using `(A, B, C, U)` and Gaussian noise.

### Method fit (MR.RGM+)

Using summary statistics:

- `S_YY = YᵀY/n`, `S_YX = YᵀX/n`, `S_XX = XᵀX/n`

the code fits **MR.RGM+** using:

- `D = 1` (a matrix of ones of dimension `p × k`),
- Spike-and-slab prior with `SigmaStarModel = "SSSL"`.

The fitted output includes an estimated matrix:

- `PsiEst` (probabilities for instrument–trait links).

## Metric: AUC for IV recovery

For each replicate, the script computes:

- **Truth**: `as.vector((B != 0) * 1)`  
- **Score**: `as.vector(Output$PsiEst)`

and reports:

- **AUC = AUC(truth, score)** using `pROC`.

So the AUC measures how well `PsiEst` ranks true instrument–trait edges above non-edges.

## Parallelization

For each `(p, n)` scenario, replicates are run in parallel via:

- `parallel::mclapply(..., mc.cores = n_cores)`

## Outputs

### 1) Raw AUC table (wide format)

The script stores AUC values in a wide table `AUC_results` with columns:

- `SampleSize`, `NetworkSize`, and `AUC1 ... AUC{n_runs}`

Each row corresponds to one `(p, n)` setting.

### 2) Manuscript boxplots (AUC vs sample size)

The wide results are converted to a long format (`df_long`) and used to generate **two manuscript plots**:

- **p = 5**: AUC boxplots across sample sizes  
- **p = 10**: AUC boxplots across sample sizes

Each plot shows:

- **x-axis**: sample size (`n`)
- **y-axis**: replicate AUC distribution
