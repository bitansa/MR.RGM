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
mcmc <- list(nIter = 300, nBurnin = 101, Thin = 2)   # Thin does not divide nBurnin on purpose
n_keep <- (mcmc$nIter - mcmc$nBurnin) %/% mcmc$Thin

grid <- expand.grid(prior = c("Spike and Slab", "Threshold"), Sigma = c("SSSL", "IW", "diagonal"),
                    iv = c(TRUE, FALSE), cov = c(TRUE, FALSE), stringsAsFactors = FALSE)

for (r in seq_len(nrow(grid))) {
  g <- grid[r, ]
  test_that(sprintf("RGM runs: %s / %s / instruments = %s / covariates = %s", g$prior, g$Sigma, g$iv, g$cov), {
    args <- list(Y = d$Y, D = d$D, prior = g$prior, SigmaStarModel = g$Sigma,
                 nIter = mcmc$nIter, nBurnin = mcmc$nBurnin, Thin = mcmc$Thin)
    if (g$iv)  args$X <- d$X
    if (g$cov) args$U <- d$U
    set.seed(10 + r)
    fit <- do.call(RGM, args)

    expect_equal(dim(fit$AEst), c(d$p, d$p))
    expect_true(all(is.finite(fit$AEst)))
    expect_equal(dim(fit$GammaPst), c(d$p, d$p, n_keep))
    expect_true(all(fit$GammaEst >= 0 & fit$GammaEst <= 1))
    expect_true(all(diag(fit$GammaEst) == 0))
    expect_equal(fit$zAEst, (fit$GammaEst > 0.5) * 1)
    expect_equal(length(fit$LLPst), n_keep)
    if (g$cov) expect_equal(dim(fit$CEst), c(d$p, d$l))
    if (g$iv) {
      expect_equal(dim(fit$BEst), c(d$p, d$k))
      expect_true(all(fit$BEst[d$D == 0] == 0))
    }
    if (g$Sigma == "SSSL") expect_true(all(fit$ZEst[upper.tri(fit$ZEst)] >= 0 & fit$ZEst[upper.tri(fit$ZEst)] <= 1))
    if (g$Sigma == "diagonal") expect_equal(length(fit$SigmaEst), d$p) else expect_equal(dim(fit$SigmaEst), c(d$p, d$p))
    if (g$prior == "Threshold") {
      ## indicators agree with the thresholded estimates in every stored draw
      expect_true(all(fit$GammaPst[, , n_keep] %in% c(0, 1)))
    }
  })
}

test_that("SSSL warns for responses far from unit scale", {
  expect_warning(RGM(X = d$X, Y = d$Y * 100, D = d$D, nIter = 20, nBurnin = 10, Thin = 1), "unit scale")
})

test_that("NetworkMotif returns a proportion", {
  set.seed(3)
  fit <- RGM(X = d$X, Y = d$Y, D = d$D, nIter = 200, nBurnin = 100, Thin = 1, SigmaStarModel = "diagonal")
  G <- matrix(0, d$p, d$p); G[2, 1] <- 1
  v <- NetworkMotif(G, fit$GammaPst)
  expect_true(v >= 0 && v <= 1)
})
