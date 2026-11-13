## Fast checks that every one of the 24 samplers runs and returns consistent output.
## (2 priors x 3 error models x {instruments, no instruments} x {covariates, no covariates})

sim_small <- function(n = 400, p = 3, k = 6, l = 2, seed = 1) {
  set.seed(seed)
  A <- matrix(0, p, p); A[2, 1] <- 0.3; A[3, 2] <- -0.3
  D <- matrix(0, p, k); for (i in 1:p) D[i, (2 * i - 1):(2 * i)] <- 1
  X <- matrix(rnorm(n * k), n, k); U <- matrix(rnorm(n * l), n, l)
  E <- matrix(rnorm(n * p), n, p)
  Y <- (X %*% t(D) + U %*% matrix(0.5, l, p) + E) %*% t(solve(diag(p) - A))
  list(X = scale(X), Y = scale(Y), U = scale(U), D = D, A = A, n = n, p = p, k = k, l = l)
}

d <- sim_small()
mcmc <- list(nIter = 300, nBurnin = 101, Thin = 2)   # Thin deliberately does not divide nBurnin
n_keep <- (mcmc$nIter - mcmc$nBurnin) %/% mcmc$Thin

grid <- expand.grid(prior = c("Spike and Slab", "Threshold"), Sigma = c("SSSL", "IW", "diagonal"),
                    iv = c(TRUE, FALSE), cov = c(TRUE, FALSE), stringsAsFactors = FALSE)

for (r in seq_len(nrow(grid))) {
  g <- grid[r, ]
  test_that(sprintf("RGM runs: %s / %s / instruments = %s / covariates = %s", g$prior, g$Sigma, g$iv, g$cov), {
    args <- list(Y = d$Y, prior = g$prior, SigmaStarModel = g$Sigma,
                 nIter = mcmc$nIter, nBurnin = mcmc$nBurnin, Thin = mcmc$Thin)
    if (g$iv)  { args$X <- d$X; args$D <- d$D }
    if (g$cov) args$U <- d$U
    set.seed(10 + r)
    fit <- do.call(RGM, args)

    expect_equal(dim(fit$AEst), c(d$p, d$p))
    expect_true(all(is.finite(fit$AEst)))
    ## no NaN/Inf anywhere in the returned estimates
    for (nm in intersect(c("GammaEst", "TauEst", "RhoEst", "PhiEst", "EtaEst", "PsiEst", "BEst", "CEst", "SigmaEst", "ZEst", "LLPst"), names(fit)))
      expect_true(all(is.finite(fit[[nm]])), info = nm)
    expect_true(all(is.finite(fit$GammaPst)))
    expect_equal(dim(fit$GammaPst), c(d$p, d$p, n_keep))
    expect_true(all(fit$GammaEst >= 0 & fit$GammaEst <= 1))
    expect_true(all(diag(fit$GammaEst) == 0))
    expect_equal(fit$zAEst, (fit$GammaEst > 0.5) * 1)
    expect_equal(length(fit$LLPst), n_keep)
    expect_s3_class(fit$Graph, "igraph")
    expect_equal(fit$settings$prior, g$prior)
    expect_equal(fit$settings$SigmaStarModel, g$Sigma)
    if (g$cov) expect_equal(dim(fit$CEst), c(d$p, d$l)) else expect_null(fit$CEst)
    if (g$iv) {
      expect_equal(dim(fit$BEst), c(d$p, d$k))
      expect_true(all(fit$BEst[d$D == 0] == 0))
      ## D has zeros: MR.RGM, no selection of the instrument effects
      expect_false(fit$settings$selectB)
      expect_equal(fit$PhiEst, d$D)
      expect_equal(fit$zBEst, d$D)
      expect_null(fit$PsiEst)
      expect_null(fit$tBEst)
    } else {
      expect_null(fit$BEst)
      expect_true(is.na(fit$settings$selectB))
    }
    if (g$Sigma == "SSSL") {
      expect_true(all(fit$ZEst[upper.tri(fit$ZEst)] >= 0 & fit$ZEst[upper.tri(fit$ZEst)] <= 1))
    } else {
      expect_null(fit$ZEst)
    }
    if (g$Sigma == "diagonal") {
      expect_false(is.matrix(fit$SigmaEst))
      expect_equal(length(fit$SigmaEst), d$p)
      expect_true(all(fit$SigmaEst > 0))
    } else {
      expect_equal(dim(fit$SigmaEst), c(d$p, d$p))
      expect_equal(fit$SigmaEst, t(fit$SigmaEst))
    }
    if (g$prior == "Threshold") {
      ## the stored edge indicators are binary
      expect_true(all(fit$GammaPst[, , n_keep] %in% c(0, 1)))
      expect_true(is.numeric(fit$tAEst))
    } else {
      expect_null(fit$tAEst)
    }
  })
}

test_that("MR.RGM+ (D = 1) selects the instrument effects", {
  D1 <- matrix(1, d$p, d$k)
  for (pr in c("Spike and Slab", "Threshold")) {
    set.seed(5)
    fit <- suppressWarnings(RGM(X = d$X, Y = d$Y, D = D1, prior = pr, nIter = mcmc$nIter, nBurnin = mcmc$nBurnin, Thin = mcmc$Thin))
    expect_true(fit$settings$selectB)
    expect_true(all(fit$PhiEst >= 0 & fit$PhiEst <= 1))
    expect_equal(fit$zBEst, (fit$PhiEst > 0.5) * 1)
    if (pr == "Spike and Slab") expect_equal(dim(fit$PsiEst), c(d$p, d$k)) else expect_true(is.numeric(fit$tBEst))
    ## point-mass spike: an instrument effect that is never selected has posterior mean exactly zero
    if (pr == "Spike and Slab") {
      expect_true(all(fit$BEst[fit$PhiEst == 0] == 0))
      expect_true(fit$settings$gibbsB)
    }
  }
})

test_that("with instrument selection under the spike-and-slab prior, gibbsB = FALSE gives a warning", {
  set.seed(6)
  expect_warning(RGM(X = d$X, Y = d$Y, D = matrix(1, d$p, d$k), gibbsB = FALSE, nIter = mcmc$nIter, nBurnin = mcmc$nBurnin,
                     Thin = mcmc$Thin), "gibbsB = FALSE is not used")
})

test_that("the point-mass update of B runs under every error model, with and without covariates", {
  for (sg in c("SSSL", "IW", "diagonal")) for (cv in c(FALSE, TRUE)) {
    args <- list(X = d$X, Y = d$Y, D = d$D, selectB = TRUE, SigmaStarModel = sg,
                 nIter = mcmc$nIter, nBurnin = mcmc$nBurnin, Thin = mcmc$Thin)
    if (cv) args$U <- d$U
    set.seed(40)
    fit <- suppressWarnings(do.call(RGM, args))
    info <- paste(sg, if (cv) "with covariates" else "without covariates")
    expect_true(fit$settings$selectB, info = info)
    expect_true(all(is.finite(fit$BEst)), info = info)
    expect_true(all(fit$PhiEst >= 0 & fit$PhiEst <= 1), info = info)
    ## entries outside the map are never selected and stay at zero
    expect_true(all(fit$PhiEst[d$D == 0] == 0), info = info)
    expect_true(all(fit$BEst[d$D == 0] == 0), info = info)
    ## the own instruments (true effect 1) are selected
    expect_true(all(fit$PhiEst[d$D == 1] > 0.5), info = info)
  }
})

test_that("selectB overrides the default given by D", {
  set.seed(6)
  ## selection restricted to the allowed entries
  fit <- RGM(X = d$X, Y = d$Y, D = d$D, selectB = TRUE, nIter = mcmc$nIter, nBurnin = mcmc$nBurnin, Thin = mcmc$Thin)
  expect_true(fit$settings$selectB)
  expect_true(all(fit$PhiEst[d$D == 0] == 0))
  expect_equal(dim(fit$PsiEst), c(d$p, d$k))
  ## full map without selection: every entry has the slab prior
  fit2 <- suppressWarnings(RGM(X = d$X, Y = d$Y, D = matrix(1, d$p, d$k), selectB = FALSE, nIter = mcmc$nIter, nBurnin = mcmc$nBurnin, Thin = mcmc$Thin))
  expect_false(fit2$settings$selectB)
  expect_equal(fit2$PhiEst, matrix(1, d$p, d$k))
  expect_null(fit2$PsiEst)
  expect_error(RGM(X = d$X, Y = d$Y, D = d$D, selectB = "yes", nIter = 20, nBurnin = 10), "selectB")
})

test_that("identification warning only without selection", {
  D_bad <- d$D; D_bad[1, ] <- D_bad[2, ]   # response 1 has no instrument of its own
  ## response 1 also has no strong instrument of its own: the weak-instrument warning fires as well
  expect_warning(expect_warning(RGM(X = d$X, Y = d$Y, D = D_bad, nIter = 20, nBurnin = 10, Thin = 1), "Weak instruments"),
                 "affects only response")
  expect_no_warning(RGM(X = d$X, Y = d$Y, D = matrix(1, d$p, d$k), nIter = 20, nBurnin = 10, Thin = 1), message = "affects only response")
})

test_that("input checks", {
  S <- list(Syy = crossprod(d$Y) / d$n, Syx = crossprod(d$Y, d$X) / d$n, Sxx = crossprod(d$X) / d$n,
            Syu = crossprod(d$Y, d$U) / d$n, Sxu = crossprod(d$X, d$U) / d$n, Suu = crossprod(d$U) / d$n)
  ## incomplete covariate summaries are an error, not silently dropped
  expect_error(RGM(Syy = S$Syy, Syx = S$Syx, Sxx = S$Sxx, Syu = S$Syu, Suu = S$Suu, D = d$D, n = d$n, nIter = 20, nBurnin = 10), "Sxu")
  ## asymmetric input is symmetrized with a warning
  Sxx_bad <- S$Sxx; Sxx_bad[1, 2] <- Sxx_bad[1, 2] + 0.05
  expect_warning(RGM(Syy = S$Syy, Syx = S$Syx, Sxx = Sxx_bad, D = d$D, n = d$n, nIter = 20, nBurnin = 10), "not symmetric")
  ## uncentered individual-level data
  expect_warning(RGM(X = d$X, Y = d$Y + 3, D = d$D, nIter = 20, nBurnin = 10), "not centered")
  ## SSSL scale
  expect_warning(RGM(X = d$X, Y = d$Y * 100, D = d$D, nIter = 20, nBurnin = 10, Thin = 1), "unit scale")
  ## hyperparameters
  expect_error(RGM(X = d$X, Y = d$Y, D = d$D, nIter = 20, nBurnin = 10, varpi = 1.5), "varpi")
  expect_error(RGM(X = d$X, Y = d$Y, D = d$D, nIter = 20, nBurnin = 10, aRho = 0), "aRho")
  expect_error(RGM(X = d$X, Y = d$Y, D = d$D, nIter = 20, nBurnin = 10, tauC = -1), "tauC")
  expect_error(RGM(X = d$X, Y = d$Y, D = d$D[, 1:3], nIter = 20, nBurnin = 10), "D must")
  expect_error(RGM(Syy = S$Syy, Syx = S$Syx, Sxx = S$Sxx, D = d$D, nIter = 20, nBurnin = 10), "n \\(number of observations\\)")
  ## incomplete instrument summaries are an error (not a responses-only fit)
  expect_error(RGM(Syy = S$Syy, Syx = S$Syx, D = d$D, n = d$n, nIter = 20, nBurnin = 10), "Syx and Sxx")
  ## without instruments a supplied selectB also gives the warning that it is not used
  expect_error(suppressWarnings(RGM(Y = d$Y, selectB = "yes", nIter = 20, nBurnin = 10)), "selectB")
})

test_that("hyperparameters are passed through and recorded", {
  set.seed(7)
  fit <- RGM(X = d$X, Y = d$Y, U = d$U, D = d$D, nIter = mcmc$nIter, nBurnin = mcmc$nBurnin, Thin = mcmc$Thin,
             omega1 = 2, omega2 = 0.01, varpi = 0.3, lambda = 2, tauC = 1, aRho = 1, bRho = 1)
  expect_equal(fit$settings[c("omega1", "omega2", "varpi", "lambda", "tauC", "aRho")], list(omega1 = 2, omega2 = 0.01, varpi = 0.3, lambda = 2, tauC = 1, aRho = 1))
  expect_equal(fit$settings$package_version, as.character(utils::packageVersion("MR.RGM")))
})

test_that("NetworkMotif returns the probability that all motif edges are present", {
  set.seed(3)
  fit <- RGM(X = d$X, Y = d$Y, D = d$D, nIter = 200, nBurnin = 100, Thin = 1, SigmaStarModel = "diagonal")
  G <- matrix(0, d$p, d$p); G[2, 1] <- 1
  v <- NetworkMotif(G, fit$GammaPst)
  expect_true(v >= 0 && v <= 1)
  expect_equal(v, mean(fit$GammaPst[2, 1, ] == 1))
  G2 <- G; G2[3, 2] <- 1
  expect_equal(NetworkMotif(G2, fit$GammaPst), mean(fit$GammaPst[2, 1, ] == 1 & fit$GammaPst[3, 2, ] == 1))
  expect_equal(NetworkMotif(matrix(0, d$p, d$p), fit$GammaPst), 1)
  expect_error(NetworkMotif(G * 2, fit$GammaPst), "0 and 1")
})

test_that("known truth is recovered (MR.RGM and MR.RGM+)", {
  skip_on_cran()
  big <- sim_small(n = 2000, seed = 11)
  set.seed(12)
  fit <- RGM(X = big$X, Y = big$Y, D = big$D, nIter = 3000, nBurnin = 1000, Thin = 2)
  expect_true(fit$GammaEst[2, 1] > 0.5 && fit$GammaEst[3, 2] > 0.5)
  expect_true(fit$GammaEst[1, 2] < 0.5 && fit$GammaEst[1, 3] < 0.5 && fit$GammaEst[3, 1] < 0.5)
  ## the data are scaled, so the true effect of Y_j on Y_i on that scale is A[i, j] * sd(Y_j) / sd(Y_i)
  sdY <- attr(big$Y, "scaled:scale")
  truth <- c(big$A[2, 1] * sdY[1] / sdY[2], big$A[3, 2] * sdY[2] / sdY[3])
  expect_true(all(abs(c(fit$AEst[2, 1], fit$AEst[3, 2]) - truth) < 0.15))
  set.seed(13)
  fit_plus <- suppressWarnings(RGM(X = big$X, Y = big$Y, D = matrix(1, big$p, big$k), nIter = 3000, nBurnin = 1000, Thin = 2))
  expect_true(all(fit_plus$PhiEst[big$D == 1] > 0.5))
  expect_true(fit_plus$GammaEst[2, 1] > 0.5 && fit_plus$GammaEst[3, 2] > 0.5)
})

test_that("Gibbs update of B, adaptive proposals and acceptance output", {
  set.seed(21)
  fit <- RGM(X = d$X, Y = d$Y, D = d$D, nIter = mcmc$nIter, nBurnin = mcmc$nBurnin, Thin = mcmc$Thin)
  ## B by Gibbs: no acceptance rate or proposal variance for B
  expect_true(fit$settings$gibbsB)
  expect_null(fit$AccptB); expect_null(fit$PropVarB)
  ## tuned proposal variances of A: p x p, zero diagonal, positive off-diagonal
  expect_equal(dim(fit$PropVarA), c(d$p, d$p))
  expect_true(all(diag(fit$PropVarA) == 0))
  expect_true(all(fit$PropVarA[row(fit$PropVarA) != col(fit$PropVarA)] > 0))
  expect_true(fit$AccptA >= 0 && fit$AccptA <= 100)
  expect_equal(fit$settings$PropVarA, "auto")
  ## Metropolis for B (gibbsB = FALSE): acceptance rate and proposal variances (zero where D = 0)
  set.seed(21)
  fit_mh <- RGM(X = d$X, Y = d$Y, D = d$D, gibbsB = FALSE, PropVarA = 0.01, PropVarB = 0.01,
                nIter = mcmc$nIter, nBurnin = mcmc$nBurnin, Thin = mcmc$Thin)
  expect_false(fit_mh$settings$gibbsB)
  expect_true(fit_mh$AccptB >= 0 && fit_mh$AccptB <= 100)
  expect_equal(dim(fit_mh$PropVarB), c(d$p, d$k))
  expect_true(all(fit_mh$PropVarB[d$D == 0] == 0) && all(fit_mh$PropVarB[d$D == 1] > 0))
  expect_equal(fit_mh$settings$PropVarA, 0.01)
  ## threshold prior with selection: B by Metropolis even with gibbsB = TRUE
  set.seed(22)
  fit_t <- suppressWarnings(RGM(X = d$X, Y = d$Y, D = matrix(1, d$p, d$k), prior = "Threshold",
                                nIter = mcmc$nIter, nBurnin = mcmc$nBurnin, Thin = mcmc$Thin))
  expect_false(fit_t$settings$gibbsB)
  expect_true(is.numeric(fit_t$AccptB))
  ## selection restricted to the allowed entries: B stays 0 where D = 0
  set.seed(23)
  fit_s <- RGM(X = d$X, Y = d$Y, D = d$D, selectB = TRUE, nIter = mcmc$nIter, nBurnin = mcmc$nBurnin, Thin = mcmc$Thin)
  expect_true(all(fit_s$BEst[d$D == 0] == 0))
  ## diagnostics are returned
  expect_true(all(c("instrumentF", "spectralRadiusA", "acceptanceA", "errorVarianceRatio") %in% names(fit$diagnostics)))   # halvesMaxDiff needs >= 200 kept draws
  expect_equal(length(fit$diagnostics$instrumentF), d$p)
  expect_equal(length(fit$diagnostics$errorVarianceRatio), d$p)
  expect_error(RGM(X = d$X, Y = d$Y, D = d$D, gibbsB = NA, nIter = 20, nBurnin = 10), "gibbsB")
  expect_error(RGM(X = d$X, Y = d$Y, D = d$D, PropVarA = -1, nIter = 20, nBurnin = 10), "PropVarA")
  expect_error(RGM(X = d$X, Y = d$Y, D = d$D, PropVarA = c(0.1, 0.2), nIter = 20, nBurnin = 10), "PropVarA")
})

test_that("same seed, same chain", {
  set.seed(31); f1 <- RGM(X = d$X, Y = d$Y, D = d$D, nIter = 200, nBurnin = 100, Thin = 1)
  set.seed(31); f2 <- RGM(X = d$X, Y = d$Y, D = d$D, nIter = 200, nBurnin = 100, Thin = 1)
  expect_identical(f1$AEst, f2$AEst); expect_identical(f1$BEst, f2$BEst); expect_identical(f1$GammaPst, f2$GammaPst)
})

test_that("further input checks", {
  expect_error(RGM(X = d$X, D = d$D, nIter = 20, nBurnin = 10), "X was given without Y")
  expect_error(RGM(X = d$X, Y = d$Y, D = d$D * 0, nIter = 20, nBurnin = 10), "no entry equal to 1")
})

test_that("Gibbs and Metropolis updates of B target the same posterior", {
  skip_on_cran()
  big <- sim_small(n = 2000, seed = 41)
  set.seed(42); fg <- RGM(X = big$X, Y = big$Y, D = big$D, nIter = 6000, nBurnin = 2000, Thin = 2)
  set.seed(43); fm <- RGM(X = big$X, Y = big$Y, D = big$D, gibbsB = FALSE, nIter = 6000, nBurnin = 2000, Thin = 2)
  expect_true(max(abs(fg$BEst - fm$BEst)) < 0.05)
  expect_true(max(abs(fg$AEst - fm$AEst)) < 0.05)
  expect_true(max(abs(fg$GammaEst - fm$GammaEst)) < 0.2)
})

test_that("rows of B with more than 100 allowed instruments (entry-wise Gibbs) work", {
  set.seed(51)
  n <- 500; p <- 3; k <- 105
  X <- scale(matrix(rnorm(n * k), n, k)); Y <- scale(X[, 1:3] %*% diag(3) + matrix(rnorm(n * p), n, p))
  D <- matrix(0, p, k); D[1, ] <- 1; D[2, c(2, 50)] <- 1; D[3, c(3, 60)] <- 1   # row 1: 105 allowed instruments
  fit <- suppressWarnings(RGM(X = X, Y = Y, D = D, nIter = 200, nBurnin = 100, Thin = 1))
  expect_true(all(is.finite(fit$BEst)))
  expect_true(all(fit$BEst[D == 0] == 0))
  expect_true(abs(fit$BEst[1, 1]) > 0.3)   # the true instrument of response 1 is found among the 105
})
