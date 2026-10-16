###############################################################################
## sim_functions.R
##
## Shared code for the simulation study of MR.RGM:
##   - simulation settings (single source of truth, reported in the Supplement)
##   - generation of the true network and of the data
##   - competing methods (MendelianRandomization, OneSampleMR, mrbayes)
##   - MR.RGM, MR.RGM_NoConf and MR.RGM+ fits
##   - performance metrics
##
## Design (one replicate):
##   1. set.seed(truth seed)  -> draw graph, A and C            (gen_truth)
##   2. set.seed(data seed)   -> draw X, U and E, compute Y      (gen_data)
##   3. every method is fitted to this SAME data set
##   The truth seed depends only on the replicate, so replicate r has the same
##   true network at every sample size n; the data seed is unique to
##   (scenario, p, n, r), so no random numbers are shared between data sets.
###############################################################################

suppressPackageStartupMessages({
  library(MASS)
  library(igraph)
  library(pROC)
  library(MendelianRandomization)
  library(OneSampleMR)
  library(mrbayes)
  library(rjags)
  library(MR.RGM)
})

## Version 4 results require MR.RGM >= 0.1.2 (corrected samplers; see NEWS.md)
if (packageVersion("MR.RGM") < "0.1.2")
  stop("MR.RGM ", packageVersion("MR.RGM"), " is installed; these scripts need >= 0.1.2: remotes::install_github(\"bitansa/MR.RGM\")")


## ---------------------------------------------------------------------------
## Settings
## ---------------------------------------------------------------------------
SIM_SETTINGS <- list(
  scenarios     = c("scalefree", "smallworld", "smallworld_pleio"),
  p_values      = c(5, 10),
  n_values      = c(500, 1000, 10000, 30000),
  n_rep         = 20,
  ivs_per_trait = 3,                 # unique instruments per trait, effect 1
  l_star        = 5,                 # latent confounders U (each N(0,1))
  err_var       = 9,                 # error variance: E ~ N(0, 9 I_p)
  a_values      = c(-0.1, 0.1),      # non-zero entries of A
  sw_nei        = 2,                 # small-world: neighbourhood size
  sw_rewire     = 0.1,               # small-world: rewiring probability
  mcmc          = list(nIter = 50000, nBurnin = 10000, Thin = 10),   # MR.RGM family
  mr_boot       = 10000,             # MendelianRandomization: bootstrap draws for median SEs (package default)
  mrbayes       = list(n.chains = 2, n.burn = 10000, n.iter = 50000),
  alpha         = 0.05,              # competitors: edge if p-value < alpha (no multiplicity correction)
  pip_thr       = 0.5                # MR.RGM family: edge / confounding link if PIP > pip_thr
)

## Seeds -----------------------------------------------------------------------
truth_seed <- function(rep) as.integer(rep)
data_seed  <- function(scenario, p, n, rep, set = SIM_SETTINGS) {
  si <- match(scenario, set$scenarios); pi <- match(p, set$p_values); ni <- match(n, set$n_values)
  as.integer(1e6 + 1e5 * si + 1e4 * pi + 1e3 * ni + rep)
}

## ---------------------------------------------------------------------------
## Truth: graph, A (p x p, A[to, from]), B (p x k), D (instrument map), C (p x l)
## ---------------------------------------------------------------------------
gen_truth <- function(scenario, p, set = SIM_SETTINGS) {

  if (scenario == "scalefree") {
    ## Barabasi-Albert tree; every edge made reciprocal (directed 2-cycles)
    g <- sample_pa(n = p, directed = FALSE)
    g <- as.directed(g, mode = "mutual")
  } else {
    ## Watts-Strogatz small-world graph; each edge gets a random direction
    sw <- sample_smallworld(dim = 1, size = p, nei = set$sw_nei, p = set$sw_rewire)
    el <- as_edgelist(sw, names = FALSE)
    el <- t(apply(el, 1, function(e) if (runif(1) > 0.5) e else rev(e)))
    g  <- make_graph(as.vector(t(el)), n = p, directed = TRUE)
  }
  adj <- as_adjacency_matrix(g, sparse = FALSE)          # adj[from, to]
  adj <- as.matrix(adj); dimnames(adj) <- NULL

  A <- matrix(sample(set$a_values, p^2, replace = TRUE), p, p) * t(adj)   # A[to, from]
  diag(A) <- 0

  ## Instruments: ivs_per_trait unique instruments per trait, effect 1
  k0 <- set$ivs_per_trait * p
  B  <- matrix(0, p, k0)
  for (i in 1:p) B[i, ((i - 1) * set$ivs_per_trait + 1):(i * set$ivs_per_trait)] <- 1
  D  <- B

  if (scenario == "smallworld_pleio") {
    ## One shared (pleiotropic) instrument per consecutive trait pair (1,2), (3,4), ...
    ## (if p is odd, an extra one for (1, p)). The shared instrument affects both traits;
    ## the instrument map D given to methods keeps only the first trait of each pair.
    starts <- if (p %% 2 == 0) seq(1, p, by = 2) else seq(1, p - 1, by = 2)
    pleio  <- sapply(starts, function(i) { v <- rep(0, p); v[c(i, i + 1)] <- 1; v })
    if (p %% 2 == 1) { v <- rep(0, p); v[c(1, p)] <- 1; pleio <- cbind(pleio, v) }
    pleio  <- matrix(pleio, nrow = p)
    first_only <- pleio
    for (j in seq_len(ncol(pleio))) {
      r1 <- which(pleio[, j] == 1)
      if (length(r1) > 1) first_only[r1[-1], j] <- 0
    }
    B <- cbind(B, pleio)
    D <- cbind(D, first_only)
  }

  ## Latent confounders: l_star factors, loadings +-1, half of them set to 0
  C <- matrix(sample(c(-1, 1), p * set$l_star, replace = TRUE), p, set$l_star)
  C[sample(which(C != 0), length(which(C != 0)) * 1 / 2)] <- 0

  list(scenario   = scenario, p = p,
       A          = A,
       graph      = (A != 0) * 1,
       B          = B,
       D          = D,
       C          = C,
       Sigma_star = C %*% t(C) + set$err_var * diag(p))   # covariance of C U + E
}

## ---------------------------------------------------------------------------
## Data: Y = (I - A)^{-1} (B X + C U + E), rows = individuals
## ---------------------------------------------------------------------------
gen_data <- function(truth, n, set = SIM_SETTINGS) {
  p <- truth$p; k <- ncol(truth$B); l <- ncol(truth$C)
  X <- matrix(rnorm(n * k), n, k)
  U <- matrix(rnorm(n * l), n, l)
  E <- matrix(rnorm(n * p, sd = sqrt(set$err_var)), n, p)
  M <- solve(diag(p) - truth$A)
  Y <- (X %*% t(truth$B) + U %*% t(truth$C) + E) %*% t(M)
  list(X = X, Y = Y)
}

## ---------------------------------------------------------------------------
## Metrics
## ---------------------------------------------------------------------------
## TPR / FDR / MCC (unchanged from the original GraphRecoveryMetricsCalculation.R)
metrics_binary_graph <- function(true_graph, estimate_graph) {
  stopifnot(length(true_graph) == length(estimate_graph))
  TP <- sum(true_graph == 1 & estimate_graph == 1)
  FP <- sum(true_graph == 0 & estimate_graph == 1)
  FN <- sum(true_graph == 1 & estimate_graph == 0)
  TN <- sum(true_graph == 0 & estimate_graph == 0)
  TPR <- TP / max(1, sum(true_graph))
  FDR <- FP / max(1, (TP + FP))
  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  MCC <- if (denom == 0) 0 else (TP * TN - FP * FN) / denom
  c(TPR = TPR, FDR = FDR, MCC = MCC)
}

## AUC with a fixed direction: larger score = more evidence for an edge
auc_safe <- function(truth, score) {
  if (length(unique(truth)) < 2 || all(is.na(score))) return(NA_real_)
  tryCatch(as.numeric(pROC::auc(pROC::roc(truth, score, levels = c(0, 1),
                                          direction = "<", quiet = TRUE))),
           error = function(e) NA_real_)
}

off_diag <- function(M) M[row(M) != col(M)]

## Graph recovery; all matrices oriented [to, from] like A
graph_metrics <- function(truth_graph, decision, score) {
  tg <- off_diag(truth_graph)
  c(AUC = auc_safe(tg, off_diag(score)), metrics_binary_graph(tg, off_diag(decision) * 1))
}

## Causal effect estimation; A_hat oriented [to, from]
effect_metrics <- function(A, A_hat) {
  d <- off_diag(abs(A - A_hat))
  c(MaxAbsDev = max(d, na.rm = TRUE), MeanAbsDev = mean(d, na.rm = TRUE), MeanSqDev = mean(d^2, na.rm = TRUE))
}

## Confounding structure: truth = |Sigma*| (upper triangle) / max, thresholded at its mean
conf_metrics <- function(Sigma_star, Z_hat, thr = SIM_SETTINGS$pip_thr) {
  ut <- upper.tri(Sigma_star)
  s  <- abs(Sigma_star[ut]); s <- s / max(s)
  tg <- (s > mean(s)) * 1
  z  <- Z_hat[ut]
  c(AUC = auc_safe(tg, z), metrics_binary_graph(tg, (z > thr) * 1))
}

## Instrument-trait selection (MR.RGM+): truth = B != 0, score = PsiEst
iv_metrics <- function(B, Psi_hat) c(AUC = auc_safe(as.vector((B != 0) * 1), as.vector(Psi_hat)))

## ---------------------------------------------------------------------------
## Competing methods. Pairwise: exposure i -> outcome j using the instruments of
## trait i (D[i, ]). Returned matrices are [exposure, outcome]; transpose to
## compare with A[to, from].
## ---------------------------------------------------------------------------
snp_betas <- function(y, G) {
  t(sapply(seq_len(ncol(G)), function(m) {
    co <- summary(lm(y ~ G[, m]))$coefficients
    c(beta = co[2, 1], se = co[2, 2])
  }))
}

fit_mr_package <- function(X, Y, D, set = SIM_SETTINGS) {
  p <- ncol(Y)
  meth <- c(SimpleMedian = "Simple median", WeightedMedian = "Weighted median", IVW = "IVW")
  est <- lapply(meth, function(m) matrix(0, p, p))
  pv  <- lapply(meth, function(m) matrix(1, p, p))
  failed <- 0
  for (i in 1:p) {
    ivs <- which(D[i, ] == 1)
    if (length(ivs) < 3) next
    G  <- X[, ivs, drop = FALSE]
    bx <- snp_betas(Y[, i], G)
    for (j in 1:p) {
      if (i == j) next
      by  <- snp_betas(Y[, j], G)
      obj <- mr_input(bx = bx[, 1], bxse = bx[, 2], by = by[, 1], byse = by[, 2],
                      snps = paste0("rs", ivs))
      res <- tryCatch(mr_allmethods(obj, method = "main", iterations = set$mr_boot),
                      error = function(e) NULL)
      if (is.null(res)) { failed <- failed + 1; next }
      v <- res@Values
      for (m in names(meth)) {
        est[[m]][i, j] <- v[v$Method == meth[[m]], "Estimate"]
        pv[[m]][i, j]  <- v[v$Method == meth[[m]], "P-value"]
      }
    }
  }
  list(est = est, pval = pv, failed = failed)
}

fit_onesamplemr <- function(X, Y, D) {
  p <- ncol(Y)
  est <- matrix(0, p, p); pv <- matrix(1, p, p); failed <- 0
  for (i in 1:p) {
    ivs <- which(D[i, ] == 1)
    if (length(ivs) < 1) next
    G <- X[, ivs, drop = FALSE]
    nm <- paste0("z", ivs); colnames(G) <- nm
    for (j in 1:p) {
      if (i == j) next
      df  <- data.frame(y = Y[, j], x = Y[, i], G)
      fit <- tryCatch(tsri(y ~ x, reformulate(nm), df, link = "identity"), error = function(e) NULL)
      ok <- FALSE
      if (!is.null(fit) && !is.null(fit$fit)) {
        co <- summary(fit$fit)$coefficients
        if ("x" %in% rownames(co)) {
          est[i, j] <- co["x", "Estimate"]
          pv[i, j]  <- co["x", "Pr(>|t|)"]
          ok <- TRUE
        }
      }
      if (!ok) failed <- failed + 1
    }
  }
  list(est = est, pval = pv, failed = failed)
}

fit_mrbayes <- function(X, Y, D, seed, set = SIM_SETTINGS) {
  p <- ncol(Y)
  est <- matrix(0, p, p); pv <- matrix(1, p, p); failed <- 0
  has_seed <- "seed" %in% names(formals(mrbayes::mr_ivw_rjags))
  for (i in 1:p) {
    ivs <- which(D[i, ] == 1)
    if (length(ivs) < 2) next
    G  <- X[, ivs, drop = FALSE]
    bx <- snp_betas(Y[, i], G)
    for (j in 1:p) {
      if (i == j) next
      by <- snp_betas(Y[, j], G)
      dat <- mr_format(rsid = paste0("SNP", ivs), xbeta = bx[, 1], ybeta = by[, 1],
                       xse = bx[, 2], yse = by[, 2])
      args <- list(object = dat, prior = "default", n.chains = set$mrbayes$n.chains,
                   n.burn = set$mrbayes$n.burn, n.iter = set$mrbayes$n.iter)
      if (has_seed) args$seed <- seed + 100 * i + j + seq_len(set$mrbayes$n.chains) - 1
      fit <- tryCatch(do.call(mr_ivw_rjags, args), error = function(e) NULL)
      if (is.null(fit) && has_seed) {          # retry without seed if the seed format was rejected
        args$seed <- NULL
        fit <- tryCatch(do.call(mr_ivw_rjags, args), error = function(e) NULL)
      }
      if (is.null(fit)) { failed <- failed + 1; next }
      s <- do.call(rbind.data.frame, fit$samples)$Estimate
      est[i, j] <- mean(fit$CausalEffect)
      pv[i, j]  <- 2 * min(mean(s > 0), mean(s < 0))
    }
  }
  list(est = est, pval = pv, failed = failed)
}

## ---------------------------------------------------------------------------
## MR.RGM family (summary-level input; simulated data have mean zero)
## ---------------------------------------------------------------------------
fit_rgm <- function(X, Y, D, model = c("SSSL", "diagonal"), set = SIM_SETTINGS) {
  model <- match.arg(model)
  n <- nrow(Y)
  RGM(Syy = crossprod(Y) / n, Syx = crossprod(Y, X) / n, Sxx = crossprod(X) / n,
      D = D, n = n,
      nIter = set$mcmc$nIter, nBurnin = set$mcmc$nBurnin, Thin = set$mcmc$Thin,
      prior = "Spike and Slab", SigmaStarModel = model)
}

## ---------------------------------------------------------------------------
## One replicate: all methods on the same data; returns long data.frame
## ---------------------------------------------------------------------------
run_replicate <- function(rep, scenario, p, n, set = SIM_SETTINGS) {

  set.seed(truth_seed(rep));                      truth <- gen_truth(scenario, p, set)
  dseed <- data_seed(scenario, p, n, rep, set);   set.seed(dseed)
  dat   <- gen_data(truth, n, set)
  X <- dat$X; Y <- dat$Y; D <- truth$D

  out <- list()
  add <- function(method, family, vals, secs) {
    out[[length(out) + 1]] <<- data.frame(scenario = scenario, p = p, n = n, rep = rep,
                                          method = method, family = family,
                                          metric = names(vals), value = unname(vals),
                                          seconds = secs, stringsAsFactors = FALSE)
  }
  timed <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr; list(v = v, s = proc.time()[["elapsed"]] - t0) }
  add_pairwise <- function(name, fit, secs) {
    add(name, "graph",  graph_metrics(truth$graph, t(fit$pval < set$alpha), t(1 - fit$pval)), secs)
    add(name, "effect", effect_metrics(truth$A, t(fit$est)), secs)
    add(name, "diagnostic", c(failed_pairs = fit$failed), secs)   # pairs where the fit errored (kept at estimate 0, p = 1)
  }
  add_rgm <- function(name, fit, secs, conf = FALSE, iv = FALSE) {
    add(name, "graph",  graph_metrics(truth$graph, fit$GammaEst > set$pip_thr, fit$GammaEst), secs)
    add(name, "effect", effect_metrics(truth$A, fit$AEst), secs)
    if (conf) add(name, "confounding", conf_metrics(truth$Sigma_star, fit$ZEst, set$pip_thr), secs)
    if (iv)   add(name, "iv_selection", iv_metrics(truth$B, fit$PsiEst), secs)
  }
  safely <- function(name, expr_fun) {
    r <- tryCatch(expr_fun(), error = function(e) { message(sprintf("[%s p=%d n=%d rep=%d] %s failed: %s",
                                                            scenario, p, n, rep, name, conditionMessage(e))); NULL })
    invisible(r)
  }

  ## Competitors
  safely("MendelianRandomization", function() {
    r <- timed(fit_mr_package(X, Y, D, set))
    for (m in names(r$v$est)) add_pairwise(m, list(est = r$v$est[[m]], pval = r$v$pval[[m]], failed = r$v$failed), r$s)
  })
  safely("OneSampleMR", function() { r <- timed(fit_onesamplemr(X, Y, D));          add_pairwise("OneSampleMR", r$v, r$s) })
  safely("mrbayes",     function() { r <- timed(fit_mrbayes(X, Y, D, dseed, set));  add_pairwise("mrbayes", r$v, r$s) })

  ## MR.RGM family
  safely("MR.RGM",        function() { r <- timed(fit_rgm(X, Y, D, "SSSL", set));     add_rgm("MR.RGM", r$v, r$s, conf = TRUE) })
  safely("MR.RGM_NoConf", function() { r <- timed(fit_rgm(X, Y, D, "diagonal", set)); add_rgm("MR.RGM_NoConf", r$v, r$s) })
  if (scenario == "smallworld_pleio") {
    Dfull <- matrix(1, nrow = p, ncol = ncol(X))
    safely("MR.RGM+", function() { r <- timed(fit_rgm(X, Y, Dfull, "SSSL", set)); add_rgm("MR.RGM+", r$v, r$s, conf = TRUE, iv = TRUE) })
  }

  do.call(rbind, out)
}
