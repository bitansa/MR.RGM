# MR.RGM 0.1.4

## Model change

* **Point-mass spike for the instrument effects under selection (MR.RGM+).** With instrument selection
  (`selectB = TRUE`, the default when `D` has no zeros) and the spike-and-slab prior, an entry of `B` is now exactly
  zero when its inclusion indicator is zero: `b ~ phi N(0, eta) + (1 - phi) delta_0`. In 0.1.3 the spike was the
  narrow normal distribution `N(0, nu2 * eta)`. The identification and posterior consistency results of the
  accompanying paper for MR.RGM+ are stated for the point-mass prior. Fits with a known instrument map
  (`selectB = FALSE`, MR.RGM), fits without instruments and fits with the threshold prior are unchanged.

## Sampler

* With selection under the spike-and-slab prior, each inclusion indicator and its entry of `B` are drawn jointly
  from their exact full conditional: the indicator with the entry integrated out under the slab, then the entry
  given the indicator. The chain starts with all indicators of `B` at zero. `gibbsB = FALSE` is not available in
  this case and gives a warning.

## Arguments

* `nu2` is no longer used. The argument is kept so that earlier code runs unchanged.

# MR.RGM 0.1.3

## Defaults

* **Default chain length** `nIter = 50000`, `nBurnin = 10000`, `Thin = 10` (was 10,000 / 2,000 / 1): the settings
  validated in the simulation study of the accompanying paper, so that a call with the defaults (for example in a
  comparison with other methods) runs a chain of the recommended length. The help page of `RGM()` recommends several
  chains with different seeds and a convergence check (split-R-hat) for a reported analysis.

## Model change

* **Known instrument map (MR.RGM) versus instrument selection (MR.RGM+).** The instrument effects in `B` are
  now selected by the spike-and-slab (or threshold) prior **only when the map `D` has no zeros** (`D = 1`,
  MR.RGM+). When `D` contains zeros (MR.RGM), the entries with `D = 0` are fixed at zero as before and the
  allowed entries receive the slab prior `N(0, eta)` with their inclusion fixed at 1: every allowed instrument
  stays in the model, as the exclusion restrictions of a known map imply. In 0.1.2 and earlier the allowed
  entries were also subject to selection, so a weak own instrument could be shrunk into the spike. The new
  argument `selectB` (default `NULL` = `all(D == 1)`) overrides the default; the value used is returned in
  `settings$selectB`. All 12 samplers with instruments (both priors, three error models, with and without
  covariates) implement the switch; with the threshold prior and no selection the threshold `tB` is fixed at 0.
  The model of MR.RGM+ (`D = 1`) is unchanged (its results change only through the new sampler below).

## Sampler

* **Exact Gibbs update of the instrument effects.** Given `A`, the error covariance (and `C`), the likelihood is
  quadratic in `B` and its prior is normal given the inclusion indicators, so every row of `B` is now drawn from its
  multivariate normal full conditional instead of by entry-wise random-walk Metropolis (spike-and-slab prior; with
  the threshold prior when the instrument effects are not selected). No tuning, no rejections. The new argument
  `gibbsB = FALSE` restores Metropolis updates of `B` (now with adaptive proposals). A row with more than 100
  allowed instruments is drawn entry by entry (each entry from its exact scalar full conditional), so the cost stays
  linear in the number of instruments.
* **Adaptive random-walk proposals for `A` (and for `B` when it is updated by Metropolis).** Every entry has its own
  proposal variance, started at 2.38^2 times an approximation of its posterior variance (`PropVarA = NULL`,
  `PropVarB = NULL`, the new defaults; a number gives a fixed starting variance) and tuned during the burn-in, in
  batches of 50 iterations, towards an acceptance rate of 0.44 (Roberts and Rosenthal, 2009); the variances are
  frozen after the burn-in, so the kept draws come from a fixed Markov kernel. The fixed default variance 0.01
  of 0.1.2 gave acceptance rates near 1% for large n. The tuned variances are
  returned as `PropVarA` (and `PropVarB`).
* **Starting values.** The chain starts at `A = 0`, at the least-squares effects of each response on its allowed
  instruments (known map; `B = 0` with selection) and at the corresponding residual covariance (its diagonal under
  the diagonal error model), instead of a random inverse-Wishart draw.
* **Acceptance rates** (`AccptA`, `AccptB`, `AccpttA`, `AccpttB`) are now reported over the iterations after the
  burn-in; `AccptB` is omitted when `B` is drawn by Gibbs.
* **Warm-up of the causal-edge indicators (spike-and-slab prior, MR.RGM).** When the instrument map is known
  (`D` with zeros, `selectB = FALSE`), the inclusion indicators of `A` are held at 1 during the first half of the burn-in, so that
  every causal effect first moves under its slab prior. Previously the chain started at `A = 0` and the first Gibbs
  update of each indicator, made while its effect was still 0, almost always set it to 0; the narrow spike could
  then hold a well-supported effect near 0 for the rest of the run. With a known map the likelihood identifies `A` while all indicators are 1, so the warm-up is safe; it is
  not used for MR.RGM+ (`selectB = TRUE`), where `A` is identified only through the sparsity of `B` (there a warm-up lets `A` drift to an observationally equivalent solution with large effects), nor without instruments. The
  warm-up changes only where the chain starts; the kept draws target the same posterior. The threshold prior is
  unchanged.

## New arguments

* The hyperparameters that were hard-coded in the samplers are now arguments of `RGM()`, with the same
  defaults: `tauC = 2` (column variance of the matrix-normal prior on the covariate effects `C`),
  `omega1 = 5` and `omega2 = 0.005` (slab and spike variances of the SSSL prior on the error covariances),
  `varpi = 0.5` (prior inclusion probability of a confounding link) and `lambda = 1` (the error variances
  have an exponential prior with rate `lambda / 2`, the parametrization of Wang, 2015).

## Interface and checks

* **Diagnostics and warnings.** `RGM()` returns `diagnostics` and warns when: the first-stage F statistic of a
  response on its allowed instruments is below 10 (known map); the spectral radius of `AEst` is 1 or larger; the
  causal inclusion probabilities differ between the two halves of the kept draws by more than 0.2 and more than four
  batch-means Monte Carlo standard errors (with at least 200 kept draws); the acceptance
  rate of `A` after burn-in is below 10%; an estimated error variance exceeds five times the second moment of its
  response (as before). The returned values are `instrumentF`, `spectralRadiusA`, `halvesMaxDiff` with its Monte Carlo
  standard error `halvesMCSE`, `acceptanceA` and `errorVarianceRatio`.
* `X` without `Y` is an error (it was ignored); a map `D` without any 1 is an error; with `Beta` input, `Sxx` is
  symmetrized and `D` checked once.
* The examples and the README no longer use `MASS` (removed from Suggests); the first README example standardizes its data.
* `src/Makevars`: unused OpenMP flags removed; `RcppDist` is no longer needed (LinkingTo).
* Integer overflow in the acceptance-rate denominators for large p and long chains fixed.
* `RGM()` returns a `settings` component with every setting used (prior, error model, `selectB`,
  dimensions, MCMC settings, hyperparameters, package version).
* `SigmaEst` is a vector of length `p` when `SigmaStarModel = "diagonal"` (it was a `1 x p` matrix).
* `PsiEst`, `tBEst` and `AccpttB` are returned only when the instrument effects are selected.
* Covariate summary matrices must be complete: with instruments, `Syu`, `Sxu` and `Suu` are all required
  (previously a missing `Sxu` silently dropped the covariates).
* **Marginal-regression input (`Sxx`, `Beta`, `SigmaHat`).** This input does not contain the covariances between
  the responses, which are reconstructed under independent errors; it now requires `SigmaStarModel = "diagonal"`
  and no covariates (both were silently accepted), and the documentation states that `SigmaHat` holds the residual
  sums of squares divided by `n`, that the input is exact only for uncorrelated instruments, and that it needs an
  instrument that affects only one response for every response (so it is not available for MR.RGM+). The
  reconstruction of `A` from the marginal coefficients now works for two responses (it failed for `p = 2`), and
  inconsistent `Beta`/`SigmaHat` (non-positive implied error variances) are an error.
* Mixed inputs are errors instead of being silently ignored: individual-level data together with summary matrices,
  `Syy` or `Syx` together with `Beta`, `Beta` without `SigmaHat`, and `Syx` without `Sxx`. Fewer than two responses
  are an error. Without instruments, a supplied `D` or `selectB` gives a warning that it is not used.
* Argument checks: every hyperparameter must be a single strictly positive number (0 was accepted), and `varpi`
  must lie in (0, 1); `omega1 <= omega2` is an error under the
  SSSL model (and is no longer checked for the other error models); `nu1` or `nu2` of 1 or more gives a warning;
  `nIter`, `nBurnin` and `Thin` must be finite integers within R's integer range; a `PropVarB` that is not used
  (Gibbs update of `B`) gives a warning; a logical map `D` is accepted.
* Outputs carry the names of the responses, instruments and covariates when the input has them (column names of
  `Y`, `X`, `U` or the dimnames of the summary matrices), and the nodes of `Graph` are labeled by the response names.
  `LLPst` is a numeric vector (it was a one-column matrix). `settings$input` records the input type.
* `NetworkMotif()`: a logical motif matrix is accepted; a motif with self-loops, non-binary draws and an empty set
  of draws are errors; `Gamma` must be binary, and a motif without edges has probability 1. The documentation
  states the orientation (`Gamma[i, j] = 1` for the edge from response `j` to response `i`).
* The samplers can be interrupted by the user (`Rcpp::checkUserInterrupt()` every 100 iterations).
* The inverse-Wishart draws use the package's own implementation of the Bartlett decomposition (same algorithm and
  random-number stream as before, so the results are unchanged).
* Examples rewritten: fast (well under a second each), standardized data, and covering the known map, MR.RGM+,
  second-moment and marginal-regression input, covariates, the confounding output `ZEst`, the diagnostics and
  `NetworkMotif()`. New tests for every input type, the argument checks, the output format and `NetworkMotif()`.
* DESCRIPTION: shorter title, rewritten description with references.
* `Syy`, `Sxx` and `Suu` are symmetrized before use, with a warning when the input is clearly asymmetric
  or not positive semi-definite (the incremental likelihood updates of the samplers assume symmetry).
* Individual-level input that is not centered triggers a warning (the model has no intercept).
* The warning that a response has no instrument of its own is issued only when the instrument effects
  are not selected (it does not apply to MR.RGM+), and it names the responses concerned.
* The documentation of `RGM()` was rewritten (model, the two versions, inputs, every argument and output);
  the README plots of the true graph now use `t(A != 0)` so that the arrows follow the convention
  `A[i, j]` = effect of `Y_j` on `Y_i`; a new example contrasts MR.RGM and MR.RGM+.
* Internally, `RGM()` dispatches to the 24 samplers through one code path (the 24 near-identical call and
  return blocks were replaced); the outputs are unchanged apart from the points above.

# MR.RGM 0.1.2

Released on GitHub only; not submitted to CRAN.

## Bug fixes

* **Update of A (all 24 samplers).** Each proposal for an entry of `A` uses the log-determinant
  `log|det(I - A)|`, which is updated with the matrix determinant lemma and the inverse
  `(I - A)^{-1}`, which is updated with the Sherman-Morrison formula. The Sherman-Morrison update
  was computed after the current value had been overwritten by the accepted value, so the change it
  applied was zero and `(I - A)^{-1}` was not updated within a sweep. After the first accepted move
  in a sweep, the Jacobian term `n * log|det(I - A)|` of every later proposal used an out-of-date
  inverse. The error grows with `n`. Both quantities are now updated correctly; they are still
  recomputed exactly at the start of every sweep.
* **Samplers with covariates** (`U` / `Syu`, `Sxu`, `Suu` supplied): the full conditional of the
  error covariance now includes the determinant contribution of the covariate-effect prior
  `C ~ MN(0, Sigma, TAU * I_l)`. Previously only its trace term was included.
  - `SigmaStarModel = "SSSL"`: GIG shape `1 - (n + l)/2` (was `1 - n/2`)
  - `SigmaStarModel = "IW"`: inverse-Wishart degrees of freedom `p + 1 + n + l` (was `p + 1 + n`)
  - `SigmaStarModel = "diagonal"`: inverse-gamma shape `(n + l)/2 + aSigma` (was `n/2 + aSigma`)

  Here `l` is the number of covariates.
* **Thinning.** Posterior draws are now kept at iterations `nBurnin + Thin, nBurnin + 2 * Thin, ...`.
  Previously draws were kept when `i %% Thin == 0`, which stored one draw more than allocated (and
  stopped with an index error) when `Thin` did not divide `nBurnin`. Results are unchanged when
  `Thin` divides `nBurnin`.
* **Threshold prior.** After a threshold (`tA`, `tB`) is accepted, the edge indicators (`Gamma`,
  `Phi`) are refreshed, so `GammaEst`, `zAEst`, `GammaPst` (and `PhiEst`, `zBEst`) agree with `AEst` and `BEst`.
* **Inclusion probabilities** of the spike-and-slab indicators (`Gamma`, `Phi`) and of the SSSL
  indicators (`Z`) are computed on the log scale. Previously a `0/0` (both densities underflowing)
  was set to 0 or 0.5 depending on the sampler. An entry of `A` or `B` that never leaves its starting
  value 0 lets its variance parameter drift until it underflows to 0 or overflows; the indicator
  probability is now defined in that limit, so no `NaN` can enter `GammaEst`, `TauEst`, `RhoEst`
  (or `PhiEst`, `EtaEst`, `PsiEst`).
* **Diagonal error model:** the error variances start at the residual variances of the starting values instead of a draw from the vague `IG(aSigma, bSigma)` prior, which could be
  numerically 0 or infinite.
* `LLPst`: the Gaussian constant is now `-(n p / 2) log(2 pi)` (was `-(n / 2) log(2 pi)`). Only the
  reported values change.

## Interface and documentation

* `RGM()` warns when the estimated error variances exceed the second moments of the responses by more than
  a factor of 5, a sign of weak identification (typically `D = 1`, i.e. MR.RGM+, with instruments that explain
  little variance, or a response without any instrument). The help page explains what to check and suggests
  a more restrictive instrument map `D`.
* `prior` and `SigmaStarModel` are matched with `match.arg()`; the defaults are unchanged
  (`"Spike and Slab"`, `"SSSL"`).
* `RGM()` warns when `SigmaStarModel = "SSSL"` and the responses are far from unit scale (second
  moments below 0.01 or above 100). The SSSL prior uses fixed spike and slab variances (0.005 and 5),
  so on such scales the inclusion probabilities of the error covariances stay at their prior value 0.5.
* The documentation states that the model has no intercept (center the data, or add a column of ones
  to `U`) and that the `(Sxx, Beta, SigmaHat)` input treats marginal instrument effects as joint
  effects, which is exact for uncorrelated instruments.
* Each C++ sampler file starts with a header that describes its model, priors and update order.
