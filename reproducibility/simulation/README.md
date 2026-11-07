# Simulation study

This folder reproduces every simulation figure and table of the manuscript and its online appendix.

**Method names.** In the manuscript, the methods are called `IV-RGM` (the R package `MR.RGM` with a known
instrument-response map, `selectB = FALSE`), `IV-RGM+` (all instruments candidates for every response,
`selectB = TRUE`) and `IV-RGM_NoConf` (a diagonal error covariance). The scripts and saved result files use the
internal keys `MR.RGM`, `MR.RGM+` and `MR.RGM_NoConf`; tables and figures display the manuscript names
(function `display_text()`).

| File | Purpose |
|---|---|
| `sim_functions.R` | Settings (single source of truth), data generation and standardization, the MR.RGM family, every competing method, metrics |
| `run_all.R` | The complete study in one run (from a Terminal; see below) |
| `pilot_competitors.R` | One-off check of the competitor wrappers on one replicate (prints raw outputs next to the truth; run once before the study) |
| `run_simulations.R` | Runs one experiment (100 replicates per setting, in parallel); one result file per setting in `results_<experiment>/` |
| `make_figures_tables.R` | Builds every figure (`figures/`) and table (`tables/*.tex`) from the result files; `tables/simulation_summary.csv` holds every mean, SD and replicate count |
| `run_instrument_strength.R` | Every method under strong, moderate and weak instruments (runs, figures `strength_*`, tables incl. diagnostics) |
| `Runtime Analysis/RunTimeAnalysis.R` | Runtime benchmark (scale-free, n = 30,000, p ∈ {2, 5, 10, 20}, 100 runs per method) |
| `cluster/` | The same study on a SLURM cluster: one task per setting, runtime benchmark in single-core tasks (`cluster/README.md`) |

## Requirements

* MR.RGM >= 0.1.3 (`remotes::install_github("bitansa/MR.RGM")`, or from CRAN); the scripts stop otherwise. Its MCMC draws the instrument effects `B` exactly from their normal full conditional (Gibbs)
  and updates the causal effects `A` by random-walk Metropolis with proposal variances tuned during the burn-in
  (adaptive, then frozen); the chain starts at the least-squares instrument effects. 0.1.3 distinguishes **MR.RGM** (instrument map `D` given; the allowed instrument effects are estimated without
  selection, `selectB = FALSE`) from **MR.RGM+** (`D = 1`, every instrument effect subject to spike-and-slab
  selection, `selectB = TRUE`).
* R packages: `igraph`, `pROC`, `rjags` (+ `ggplot2`, `scales` for the figures; `microbenchmark` for the runtime
  analysis), and the competing methods at the versions of the study (list in `../README.md`); a method whose package
  is missing is skipped and reported:
  ```r
  remotes::install_version("MendelianRandomization", "0.10.0")                    # IVW, medians and MVMR
  remotes::install_version("mrbayes", "0.5.3")
  remotes::install_version("OneSampleMR", "0.1.8")                                # archived on CRAN
  remotes::install_github("ZhaotongL/MR2G")                                       # MR2G 0.1.0
  install.packages(c("MR2", "MrDAG"), repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
  install.packages(c("cvTools", "mnormt")); remotes::install_github("chunlinli/grivet@1777dbf", subdir = "grivet")   # GrIVET
  ```
* `rjags` needs the JAGS program itself (https://mcmc-jags.sourceforge.io), installed before `install.packages("rjags")`.
* Parallel runs use `run_parallel()` (sim_functions.R): one forked R process per replicate (`parallel::mcparallel`;
  macOS and Linux). On Windows (no forking), or with `N_CORES <- 1`, the replicates run sequentially (slow).
* The LaTeX tables use the packages `adjustbox` and `multirow`.

## How to run

**Recommended: the whole study from a Terminal** (not inside RStudio: the replicates run in forked R processes, which
R's documentation discourages inside a GUI, and which may not return there):

```bash
cd path/to/MR.RGM/reproducibility/simulation
caffeinate -i Rscript run_all.R > run_all.log 2>&1 &      # macOS; on Linux: nohup Rscript run_all.R > run_all.log 2>&1 &
tail -f run_all.log                                       # progress ("done in ... min"); Ctrl+C stops only the display
```

`run_all.R` runs the experiments and the instrument-strength study (main, roles, strength, noise, confounding,
sensitivity, large n), then the figures and tables, then the runtime benchmark; it is restartable (finished settings
are skipped), and a step that fails is reported and the next step still runs. Result files carry a sampler tag: a
result file without the current tag is not reused (the scripts stop and ask for it to be moved out of `results_*/`).
A replicate that runs longer than `REP_TIMEOUT_MIN` (default 240 minutes, twice that at n = 10^6) is terminated and
restarted once, so that no single replicate can block a setting. Since such a replicate is stopped with every method
in it, `run_all.R` skips GrIVET (`SKIP_METHODS`), which is then added in a second pass with a limit of 48 hours per fit
(`cluster/add_methods.R`; commands in `../README.md`).

Step by step instead:

```r
setwd("path/to/MR.RGM/reproducibility/simulation")

source("pilot_competitors.R")                                   # once: confirms the orientation of every competitor's output
for (e in c("main", "noise", "confounding", "largen", "sensitivity", "roles")) {
  SMOKE_TEST <- TRUE; EXPERIMENT <- e; source("run_simulations.R")   # minutes: every method and metric once (results_smoke_*/)
}
SMOKE_TEST <- FALSE
for (e in c("main", "noise", "confounding", "largen", "sensitivity", "roles")) {
  EXPERIMENT <- e; source("run_simulations.R")                  # full experiment (restartable: finished settings are skipped)
}
source("make_figures_tables.R")                                 # figures/ and tables/ from results_*/
source("run_instrument_strength.R")                             # instrument-strength study (several days on 9 cores; restartable)
source("Runtime Analysis/RunTimeAnalysis.R")                    # runtime table and figure
```

**On a SLURM cluster** (hours instead of days): see `cluster/README.md`.

Set `N_CORES` before sourcing to choose the number of cores (default: all but one). A subset of an experiment can
be run with `SCENARIOS`, `P_VALUES`, `N_VALUES` and `SETTINGS` (indices).

## Design

* **Scenarios.** Case 1: scale-free network (Barabási–Albert tree, every edge reciprocal). Case 2: small-world
  network (Watts–Strogatz, neighbourhood 2, rewiring probability 0.1; each edge oriented at random; for p = 5 the
  ring lattice with neighbourhood 2 is already complete, so the graph is a random orientation of all 10 pairs). Case 3:
  as Case 2 (same graph, **A** and **C** in each replicate) plus invalid instruments: one shared instrument per
  pair of responses (1,2), (3,4), …, and for odd p an extra one for (1, p), acting on both responses of its pair with
  the effect of an own instrument. Experiment `roles`: 4 exposures and
  4 outcomes with known roles (edges exposure → outcome, a random DAG among the outcomes, effects ±0.3), the design
  built for MR2 and MrDAG (which are given the true split into exposures and outcomes there).
* **Sizes.** p ∈ {5, 10} (8 in `roles`); n ∈ {500, 1000, 10000, 30000} in `main` and in the instrument-strength study,
  n ∈ {1000, 10000} in `noise`, `confounding`, `sensitivity` and `roles`, n ∈ {10⁵, 10⁶} in `largen`; 100 replicates per
  setting everywhere.
* **Truth.** Non-zero entries of **A** are ±0.1 (equal probability) on the edges of the graph; each response has 3
  unique instruments with effect 1 on the generating scale; 5 latent confounders with loadings ±1 (times the
  confounder scale, 1 in the main study), half of them (rounded down) set to 0. In Case 3 the instrument map given to every
  method assigns each shared instrument to the first of its two responses; MR.RGM+ instead allows every instrument for
  every response.
* **Data.** X, U ~ N(0, 1); E ~ N(0, σ² I) with σ² = 9 in the main study; **Y** = (I − **A**)⁻¹(**B X** + **C U** + **E**).
  Every column of X and Y is then centered and scaled (as recommended for the SSSL prior and done for the real
  data); the truth is expressed on the same scale (a*ᵢⱼ = aᵢⱼ sⱼ/sᵢ, b*ᵢₕ = bᵢₕ tₕ/sᵢ, Σ* rescaled). All methods see the
  standardized data. The instrument strength of the simulated responses is reported as the population R² of each
  response on its own instruments (median ≈ 0.2 in the main study), with the implied μ² and F.
* **Replicates and seeds.** In each replicate a new graph, **A** and **C** are drawn (seed = replicate number, so a
  replicate has the same true network at every n); the data are drawn with a seed unique to (scenario, p, n,
  replicate) and shared by the settings of an experiment (paired comparison; the "base" setting of every
  experiment reproduces the data of the main study). Every method is fitted to the same data set after its own
  `set.seed(data seed + method offset)`, so no method's random stream depends on which other methods were run.
* **Methods.** MR.RGM (spike-and-slab graph prior, SSSL error covariance, known instrument map, no selection of
  the instrument effects), MR.RGM_NoConf (as MR.RGM with a diagonal error covariance) and MR.RGM+ (D = 1, instrument
  effects selected); MR2G (Lin, Pan & Xue; package defaults: `IV_screen = "ScreenAug"`, 100 perturbations; inputs =
  per-SNP marginal effects, SEs, the phenotypic correlation matrix, the candidate instruments of each pair according
  to D and the LD matrix); MVMR (`MendelianRandomization::mr_mvivw`, package defaults; one outcome at a time with all
  other responses as exposures); MR2 (`MR2::MR2`, package defaults including `EVgamma = c(2, 2)`, which with 4 exposures
  is a Beta(1, 1) prior on the inclusion probability and is proper only with at least 4 exposures) and MrDAG
  (`MrDAG::MrDAG`, package defaults `pp = 0.05`, `tempMax = 10`), both of which need a split of the responses into
  exposures and outcomes: in the `roles` experiment they are given the true split; in the designs without roles MrDAG
  is applied one outcome at a time with all other responses as exposures, exactly as MVMR (`by_outcome()`; row j of
  the estimated network comes from the fit in which response j is the outcome), while MR2 cannot be (its sampler
  conditions each response on the others and needs at least two responses), so MR2 is evaluated in the `roles`
  experiment only; GrIVET (Chen, Li, Shen & Pan 2024; package `grivet`): causal discovery for a directed acyclic
  graph with hidden confounding, its data-driven pipeline (`cv.intdag.pmle.diff.aic()`, `topological_order()`,
  `cv.intdag.coe.penalized()`) with every instrument a candidate for every response; the package has no default tuning
  parameters, and the grid of its own scripts (tau 0.2–0.3, gamma 0.1–1), tuned to effects of order 1, selects nothing
  here, so the cross-validation grid is extended to the effect sizes of this study (`SIM_SETTINGS$grivet`); its
  penalties are not scale free, so it is fitted to the unstandardized data (on the standardized data it selects
  almost nothing), and its estimated effects are converted to the standardized scale as the truth is; edge score =
  |effect|, edge = its own selection. Every GrIVET fit ends in one of four ways: a result; a peeling failure (its peeling step
  cannot order every response, which a directed cycle can prevent); a numerical error; or no result within 48 hours.
  Under R >= 4.3 the package's internal selection of the tuning parameters stops with a coercion error, which
  `grivet_compat()` corrects by running the package's own fallback; and the pairwise baselines
  (MendelianRandomization simple median, weighted median, IVW with 10,000 bootstrap draws; OneSampleMR two-stage
  residual inclusion; mrbayes Bayesian IVW, 3 chains as in the package default) for every ordered pair of responses.
  **Every method is run in every experiment and in the instrument-strength study** (MR2: `roles` only, see above), except `sensitivity`, which varies
  the hyperparameters of the MR.RGM family only (MR.RGM, MR.RGM_NoConf, MR.RGM+; in the ν₂ settings MR.RGM+ only,
  since the others do not select instrument effects): the competitors do not use these hyperparameters, and the
  sensitivity data are those of the main study (same seeds), so their results are the main-study results. **Every MCMC method uses the same chain length:** 50,000 iterations in total,
  of which 10,000 are burn-in, thinning 10 (4,000 draws) — the MR.RGM family, MR2 and MrDAG (whose packages have no
  default chain length), and mrbayes (10,000 burn-in + 40,000 iterations per chain; no thinning argument). Every other
  tuning parameter and prior hyperparameter of every competitor is its package default. MR.RGM family: package
  defaults for every hyperparameter (`aRho = 3, bRho = 1, nu1 = 1e-3, nu2 = 1e-4, tauC = 2,
  omega1 = 5, omega2 = 0.005, varpi = 0.5, lambda = 1`), varied only in the `sensitivity` experiment.
  Two input details, both documented in `sim_functions.R`: (i) MR2G 0.1.0 fails with igraph >= 2.0 in its ScreenAug
  step (the diagonal of an adjacency matrix is NaN); `mr2g_screenaug()` is `MR2G()` with `IV_screen = "ScreenAug"`
  line by line, with only those NA entries set to FALSE. (ii) The default priors of MR2 and MrDAG assume summary
  statistics of unit order; the marginal effects of the standardized data are therefore multiplied by sqrt(n)
  (approximately z-statistics), which leaves their direct and causal effects unchanged.
  A fit that errors (any method) is stored (competitors: as "no network found"; `seconds = NA`, `failed_fits = 1`)
  and its message written to `results_<experiment>/errors.log`. Every table and figure summarises a method over its
  successful fits; a method without failures is summarised over all 100 replicates. Where a method has fewer than
  100, the table entries are marked † and the caption gives the number of successful fits per setting; with fewer
  than 10 successful fits the method is not summarised in that setting (table entry "failed (k/100)", left out of
  the figure, listed in the figure notes); `tables/figure_notes.txt` gives the same
  for every figure (for the captions in the paper), `tables/failures.tex` lists all failed fits, and
  `tables/simulation_all_replicates.csv` flags them (column `failed_fit`). Typical case: MR2G needs a genome-wide
  significant instrument (p < 5e-8, fixed in the package) for every response, which a small n or weak instruments may
  not provide.
  Every competitor is run through a wrapper in `sim_functions.R` that documents its inputs, settings and the
  orientation of its output (confirmed on one replicate by `pilot_competitors.R`).
* **Metrics.** Graph recovery: AUC (scores 1 − p-value or posterior inclusion probability), TPR, FDR, MCC
  (edge if p < 0.05, or PIP > 0.5). Effect estimation: maximum, mean absolute and mean squared deviation over the
  off-diagonal entries of **A** (entries a method does not estimate are ignored). Confounding: truth = |Σ\*| on the generating
  scale (off-diagonal), divided by its maximum and thresholded at its mean; AUC and TPR/FDR/MCC (link if PIP > 0.5); in
  the `roles` experiment restricted to the outcomes, and MR2's residual correlation is scored the same way.
  Instrument selection (MR.RGM+): AUC of the inclusion probabilities against **B** ≠ 0. In the `roles` experiment
  the graph and effect metrics are also reported for the exposure → outcome block and the outcome → outcome block
  separately. Diagnostics of every MR.RGM fit: error-variance ratio, chain-half PIP correlation and largest difference,
  spectral radius, acceptance rate of A, smallest first-stage F, and which of the package's warnings fired
  (`warn_*`; the warnings are recorded, not printed); `fits_*.rds` also stores the tuned proposal variances
  (`PropVarA`), the package's `diagnostics` and the warning messages;
  competitor fits that error are stored as "no network found", flagged and counted (see above). Wall-clock seconds per fit are stored.
  Tables report mean (SD) over replicates.

`run_simulations.R` ends with a check that every method has 100 valid replicates in every setting and lists any
competitor fit that errored.

## Instrument-strength study

The main study uses strong instruments (population R² ≈ 0.2 per response). `run_instrument_strength.R` repeats Case 2
(small-world, valid instruments) and Case 3 (small-world, invalid instruments) with p ∈ {5, 10} and n ∈ {500, 1,000, 10,000, 30,000},
for instrument effects 1 (strong), 0.5 (moderate) and 0.3 (weak; population R² ≈ 0.02–0.03 per response; the OneK1K
genes have sample partial R² between 0.03 and 0.19), each with 3 instruments per response, and the weak setting with
10 instruments per response (weak, 10 instruments: the same total strength spread over ten instruments), fitting every method of the main study (weak instruments are the regime of the real
data; one replicate is `run_replicate()` of the main study, `strength_replicate()`). Truth and data seeds are those of the main study and
U, E are drawn before X, so only the instrument effects differ between strength levels. It reports the main-study
metrics (network AUC, MAD, confounding AUC, instrument-selection AUC) plus the diagnostics (population R², μ² and F
of the simulated responses, largest error-variance ratio and the share of fits that trigger the package warning,
chain-half correlation of the causal PIPs, spectral radius) in `tables/strength_*.tex` (`strength_network_auc_p=*`,
`strength_effect_mad_p=*`: every method; `strength_failures.tex`: failed fits).
