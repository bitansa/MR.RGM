# MR.RGM 0.1.2

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
  was set to 0 or 0.5 depending on the sampler.
* **Diagonal error model:** the error precisions start at 1 instead of a draw from the vague
  `IG(aSigma, bSigma)` prior, which could be numerically 0 or infinite.
* `LLPst`: the Gaussian constant is now `-(n p / 2) log(2 pi)` (was `-(n / 2) log(2 pi)`). Only the
  reported values change.

## Interface and documentation

* `prior` and `SigmaStarModel` are matched with `match.arg()`; the defaults are unchanged
  (`"Spike and Slab"`, `"SSSL"`).
* `RGM()` warns when `SigmaStarModel = "SSSL"` and the responses are far from unit scale (second
  moments below 0.01 or above 100). The SSSL prior uses fixed spike and slab variances (0.005 and 5),
  so on such scales the inclusion probabilities of the error covariances stay at their prior value 0.5.
* The documentation states that the model has no intercept (center the data, or add a column of ones
  to `U`) and that the `(Sxx, Beta, SigmaHat)` input treats marginal instrument effects as joint
  effects, which is exact for uncorrelated instruments.
* Each C++ sampler file starts with a header that describes its model, priors and update order.
