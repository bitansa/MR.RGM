## Input modes, argument checks, output format and NetworkMotif checks.

sim_in <- function(n = 400, p = 3, k = 6, l = 2, seed = 1) {
  set.seed(seed)
  A <- matrix(0, p, p); A[2, 1] <- 0.3; A[3, 2] <- -0.3
  D <- matrix(0, p, k); for (i in 1:p) D[i, (2 * i - 1):(2 * i)] <- 1
  X <- matrix(rnorm(n * k), n, k); U <- matrix(rnorm(n * l), n, l)
  Y <- (X %*% t(D) + U %*% matrix(0.5, l, p) + matrix(rnorm(n * p), n, p)) %*% t(solve(diag(p) - A))
  list(X = scale(X), Y = scale(Y), U = scale(U), D = D, n = n, p = p, k = k, l = l)
}
d <- sim_in()
S <- with(d, list(Syy = crossprod(Y) / n, Syx = crossprod(Y, X) / n, Sxx = crossprod(X) / n,
                  Syu = crossprod(Y, U) / n, Sxu = crossprod(X, U) / n, Suu = crossprod(U) / n))
short <- list(nIter = 200, nBurnin = 100, Thin = 1)

test_that("second-moment input reproduces the individual-level fit", {
  set.seed(61); f1 <- do.call(RGM, c(list(X = d$X, Y = d$Y, U = d$U, D = d$D), short))
  set.seed(61); f2 <- do.call(RGM, c(list(Syy = S$Syy, Syx = S$Syx, Sxx = S$Sxx, Syu = S$Syu, Sxu = S$Sxu,
                                          Suu = S$Suu, D = d$D, n = d$n), short))
  expect_equal(f1$GammaPst, f2$GammaPst)
  expect_equal(f1$AEst, f2$AEst)
  expect_equal(f1$CEst, f2$CEst)
  expect_equal(f1$settings$input, "XY"); expect_equal(f2$settings$input, "S")
  set.seed(62); g1 <- do.call(RGM, c(list(Y = d$Y), short))
  set.seed(62); g2 <- do.call(RGM, c(list(Syy = S$Syy, n = d$n), short))
  expect_equal(g1$AEst, g2$AEst)
})

test_that("marginal-regression input: diagonal errors only, no covariates, p = 2 works", {
  for (pp in 2:3) {
    set.seed(pp)
    n <- 2000
    X <- scale(matrix(rnorm(n * pp), n, pp)); Y <- scale(X + matrix(rnorm(n * pp), n, pp))
    Bt <- SH <- matrix(0, pp, pp)
    for (i in 1:pp) for (j in 1:pp) {
      b <- sum(X[, j] * Y[, i]) / sum(X[, j]^2); Bt[i, j] <- b; SH[i, j] <- sum((Y[, i] - b * X[, j])^2) / n
    }
    fit <- suppressWarnings(RGM(Sxx = crossprod(X) / n, Beta = Bt, SigmaHat = SH, D = diag(pp), n = n,
                                SigmaStarModel = "diagonal", nIter = 200, nBurnin = 100, Thin = 1))
    expect_true(all(is.finite(fit$AEst)))
    expect_equal(fit$settings$input, "Beta")
  }
  args <- list(Sxx = diag(3), Beta = matrix(0.1, 3, 3), SigmaHat = matrix(1, 3, 3), D = diag(3), n = 100, nIter = 20, nBurnin = 10)
  expect_error(do.call(RGM, args), "diagonal")
  expect_error(do.call(RGM, c(args[names(args) != "D"], list(D = matrix(1, 3, 3), SigmaStarModel = "diagonal"))), "only that response")
  expect_error(do.call(RGM, c(args, list(SigmaStarModel = "diagonal", Suu = diag(1), Syu = matrix(0, 3, 1),
                                         Sxu = matrix(0, 3, 1)))), "Covariates cannot")
  expect_error(RGM(Sxx = diag(3), Beta = matrix(0.1, 3, 3), D = diag(3), n = 100, nIter = 20, nBurnin = 10), "together")
})

test_that("mixed inputs and invalid arguments are errors", {
  expect_error(do.call(RGM, c(list(Y = d$Y, Syx = S$Syx, Sxx = S$Sxx, D = d$D, n = d$n), short)), "not both")
  expect_error(do.call(RGM, c(list(X = d$X, Y = d$Y, Syy = S$Syy, D = d$D), short)), "not both")
  expect_error(do.call(RGM, c(list(Y = d$Y[, 1, drop = FALSE]), short)), "at least two responses")
  a <- list(X = d$X, Y = d$Y, D = d$D)
  expect_error(do.call(RGM, c(a, nIter = 20, nBurnin = 20)), "nBurnin")
  expect_error(do.call(RGM, c(a, nIter = 20, nBurnin = 10, Thin = 11)), "Thin")
  expect_error(do.call(RGM, c(a, nIter = 20.5, nBurnin = 10)), "nIter")
  expect_error(do.call(RGM, c(a, nIter = NA_real_, nBurnin = 10)), "nIter")
  expect_error(do.call(RGM, c(a, list(aRho = c(1, 2)), nIter = 20, nBurnin = 10)), "aRho")
  expect_error(do.call(RGM, c(a, omega1 = 0.001, nIter = 20, nBurnin = 10)), "omega1")
  expect_warning(do.call(RGM, c(a, nu1 = 2, short)), "nu1")
  expect_warning(do.call(RGM, c(a, PropVarB = 0.1, short)), "PropVarB")
  expect_warning(do.call(RGM, c(list(Y = d$Y, D = d$D), short)), "not used")
  expect_error(RGM(Syy = S$Syy, Syx = S$Syx, Sxx = S$Sxx, U = d$U, D = d$D, n = d$n, nIter = 20, nBurnin = 10), "cannot be computed")
  expect_error(RGM(X = d$X, Y = d$Y, U = d$U, Suu = S$Suu, D = d$D, nIter = 20, nBurnin = 10), "not both")
  expect_error(RGM(X = d$X, Y = d$Y, U = d$U[-1, ], D = d$D, nIter = 20, nBurnin = 10), "rows of U")
  expect_error(RGM(Syy = S$Syy, Syx = S$Syx[, -1], Sxx = S$Sxx, D = d$D, n = d$n, nIter = 20, nBurnin = 10), "Syx")
})

test_that("a logical map D is accepted", {
  set.seed(5); f1 <- do.call(RGM, c(list(X = d$X, Y = d$Y, D = d$D), short))
  set.seed(5); f2 <- do.call(RGM, c(list(X = d$X, Y = d$Y, D = d$D == 1), short))
  expect_equal(f1$AEst, f2$AEst)
})

test_that("output format: names from the input, vector log-likelihood, graph", {
  Y <- d$Y; X <- d$X; U <- d$U
  colnames(Y) <- c("g1", "g2", "g3"); colnames(X) <- paste0("snp", 1:d$k); colnames(U) <- c("age", "sex")
  set.seed(71); f <- RGM(X = X, Y = Y, U = U, D = d$D, nIter = 600, nBurnin = 200, Thin = 2)
  expect_equal(dimnames(f$AEst), list(colnames(Y), colnames(Y)))
  expect_equal(dimnames(f$BEst), list(colnames(Y), colnames(X)))
  expect_equal(dimnames(f$CEst), list(colnames(Y), colnames(U)))
  expect_equal(dimnames(f$GammaPst)[1:2], list(colnames(Y), colnames(Y)))
  expect_true(is.null(dim(f$LLPst)) && is.numeric(f$LLPst))
  expect_true(is.numeric(f$diagnostics$halvesMaxDiff))
  expect_true(is.numeric(f$diagnostics$halvesMCSE))
  expect_equal(sort(igraph::V(f$Graph)$name), sort(colnames(Y)))
  expect_equal(igraph::ecount(f$Graph), sum(f$zAEst))
  ## without names in the input the outputs have no dimnames
  set.seed(72); g <- do.call(RGM, c(list(X = d$X, Y = d$Y, D = d$D), short))
  expect_null(dimnames(g$AEst))
})

test_that("NetworkMotif input checks", {
  set.seed(3); fit <- do.call(RGM, c(list(X = d$X, Y = d$Y, D = d$D, SigmaStarModel = "diagonal"), short))
  G <- matrix(0, d$p, d$p); G[2, 1] <- 1
  expect_equal(NetworkMotif(G == 1, fit$GammaPst), NetworkMotif(G, fit$GammaPst))
  expect_error(NetworkMotif(diag(d$p), fit$GammaPst), "zero diagonal")
  expect_error(NetworkMotif(G, fit$GammaPst * 2), "only 0 and 1")
  expect_error(NetworkMotif(G, matrix(0, 3, 3)), "three dimensions")
  expect_error(NetworkMotif(matrix(0, 2, 2), fit$GammaPst), "square matrix")
  expect_error(NetworkMotif(G, array(0, c(3, 2, 5))), "rows and columns")
})
