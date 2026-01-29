###############################################################################
## Block 1 — MendelianRandomization baselines (summary-data MR)
## Methods:
##   • Simple median
##   • Weighted median
##   • IVW
##
## What this block does:
##   For each ordered trait pair (i -> j), use the instruments in D[i, ] to
##   estimate the causal effect of Y_i on Y_j via mr_allmethods().
##
## Effect estimate used:
##   res_df[Method, "Estimate"] from mr_allmethods().
##
## Evaluation metric (causal effect accuracy):
##   Compare the estimated causal effect matrix (transposed to match A's
##   direction convention) against the true direct-effect matrix A:
##     error(i, j) = | A[i, j] - t(causal_hat)[i, j] |
##   Report (off-diagonal only):
##     • Maximum Absolute Deviation
##     • Mean Absolute Deviation
##     • Mean Squared Error
###############################################################################

# Mendelian Randomization
library(doParallel)
library(foreach)
library(MendelianRandomization)
library(parallel)
library(MASS)
library(igraph)

# Number of simulations and parallel cores
n_runs <- 20
n_cores <- detectCores() - 1

# Parallel execution
results <- mclapply(1:n_runs, function(run) {

  set.seed(run)

  # Sample size (n), number of response variables (p), number of instrument variables (k) and number of confounders (l_star)
  n = 30000
  p =  10
  k = 3 * p
  l_star = 5

  # Variance
  var = 9

  # Create scale free netwrok
  g <- sample_pa(n = p, directed = FALSE)
  g <- as.directed(g, mode = "mutual")  # Make edges bi-directional -> cycles introduced

  # Assuming 'g' is your directed graph
  adj_matrix <- as_adjacency_matrix(g, sparse = FALSE)

  # Transpose so that [i,j] means edge from j to i
  adj_matrix_reversed <- t(adj_matrix)

  # Generate A matrix
  A = matrix(sample(c(-0.1, 0.1), p^2, replace = TRUE), p, p)

  # Make the diagonals of A to be 0
  A = A * adj_matrix_reversed

  B <- matrix(0, nrow = p, ncol = 3 * p)  # initialize zero matrix

  block_size <- 3

  for (i in 1:p) {
    start_col <- (i - 1) * block_size + 1
    end_col <- start_col + block_size - 1
    B[i, start_col:end_col] <- 1
  }


  # Generate C matrix i.e. influence of confounders on response variables, each row corresponds to a particular response
  C = matrix(sample(c(-1, 1), p * l_star, replace = TRUE), p, l_star)

  # Make the network sparse
  C[sample(which(C!=0), length(which(C!=0)) * 1 / 2)] = 0


  U = matrix(rnorm(n * l_star, 0, 1), nrow = n, ncol = l_star)

  # Calculate variance-covariance matrix
  Sigma = var * diag(p)

  # Calculate (I_p - A)^(-1)
  Mult_Mat = solve(diag(p) - A)

  Variance = Mult_Mat %*% Sigma %*% t(Mult_Mat)

  # True variance-covariance matrix
  Sigma_True = C %*% t(C) + Sigma

  # Set seed
  set.seed(run)

  # Data generation
  X = matrix(rnorm(n * k, 0, 1), nrow = n, ncol = k)

  Y = matrix(0, nrow = n, ncol = p)

  for (i in 1:n) {

    Y[i, ] = mvrnorm(n = 1, Mult_Mat %*% (B %*% X[i, ] + C %*% U[i, ]), Variance)

  }


  library(MendelianRandomization)

  # Assume Y, X, and D are already defined
  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(X)

  D = B


  # Initialize matrices
  init_mat <- function() matrix(0, nrow = p, ncol = p)
  causal_matrix_SimpleMedian <- init_mat()
  causal_matrix_WeightedMedian <- init_mat()
  causal_matrix_IVW <- init_mat()

  # MR function
  run_mr_allmethods <- function(X, Y, D, exposure_index, outcome_index) {
    ivs <- which(D[exposure_index, ] == 1)
    if (length(ivs) < 3) return(NULL)

    G <- X[, ivs, drop = FALSE]
    snp_names <- paste0("rs", ivs)
    exposure <- Y[, exposure_index]
    outcome <- Y[, outcome_index]

    get_betas <- function(y, G) {
      t(sapply(1:ncol(G), function(j) {
        fit <- summary(lm(y ~ G[, j]))
        coef <- fit$coefficients[2, ]
        c(beta = coef["Estimate"], se = coef["Std. Error"])
      }))
    }

    bx_data <- get_betas(exposure, G)
    by_data <- get_betas(outcome, G)

    mr_input_obj <- mr_input(
      bx = bx_data[, 1],
      bxse = bx_data[, 2],
      by = by_data[, 1],
      byse = by_data[, 2],
      snps = snp_names
    )

    mr_allmethods(mr_input_obj, method = "main", iterations = 50000)
  }

  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) next
      result <- run_mr_allmethods(X, Y, D, i, j)
      if (is.null(result)) next
      res_df <- result@Values

      get_val <- function(method, metric) {
        val <- res_df[res_df$Method == method, metric]
        if (length(val) == 0) NA else val
      }

      causal_matrix_SimpleMedian[i, j] <- get_val("Simple median", "Estimate")
      causal_matrix_WeightedMedian[i, j] <- get_val("Weighted median", "Estimate")
      causal_matrix_IVW[i, j] <- get_val("IVW", "Estimate")
    }
  }

  # Metric computation (off-diagonal only)
  causal_list <- list(
    causal_matrix_SimpleMedian,
    causal_matrix_WeightedMedian,
    causal_matrix_IVW
  )

  metric_result <- lapply(causal_list, function(causal_mat) {
    diff_mat <- abs(A - t(causal_mat))
    off_diag <- row(diff_mat) != col(diff_mat)
    max_err <- max(diff_mat[off_diag], na.rm = TRUE)
    mean_abs_err <- mean(diff_mat[off_diag], na.rm = TRUE)
    mean_sq_err <- mean((A - t(causal_mat))[off_diag]^2, na.rm = TRUE)
    c(max_err, mean_abs_err, mean_sq_err)
  })

  do.call(rbind, metric_result)

}, mc.cores = n_cores)

# Combine and compute column means
results_array <- simplify2array(results)
mean_metrics <- apply(results_array, c(1, 2), mean, na.rm = TRUE)
sd_metrics <- apply(results_array, c(1, 2), sd, na.rm = TRUE)


# Display results
rownames(mean_metrics) <- rownames(sd_metrics) <- c(
  "SimpleMedian", "WeightedMedian", "IVW"
)
colnames(mean_metrics) <- colnames(sd_metrics) <-  c("Maximum Absolute Deviation", "Mean Absolute Deviation", "Mean Squared Error")
print(mean_metrics)
print(sd_metrics)

T_SimpleMedian = t(results_array[1, , ])
T_WeightedMedian = t(results_array[2, , ])
T_IVW = t(results_array[3, , ])



##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
###############################################################################
## Block 2 — OneSampleMR + MrBayes baselines (individual-level MR)
## Methods:
##   • OneSampleMR (TSRI via tsri)
##   • MrBayes / Bayesian IVW via JAGS (mr_ivw_rjags)
##
## What this block does:
##   For each ordered trait pair (i -> j), use instruments in D[i, ] and
##   estimate the causal effect of Y_i on Y_j using:
##     (a) OneSampleMR (tsri)
##     (b) Bayesian IVW (mr_ivw_rjags)
##
## Effect estimate used:
##   • OneSampleMR: coefficient of x in the TSRI regression (Estimate)
##   • MrBayes: posterior mean causal effect (mean(bayes_model$CausalEffect))
##
## Evaluation metric (causal effect accuracy):
##   Compare estimated causal effect matrices (transposed to match A) to A
##   and report (off-diagonal only):
##     • Maximum Absolute Deviation
##     • Mean Absolute Deviation
##     • Mean Squared Error
##
## Note:
##   This block is computationally heavy because it fits per-pair models
##   and the Bayesian method runs MCMC for each pair.
###############################################################################

# OneSampleMR + MrBAYES
library(parallel)
library(OneSampleMR)
library(ivreg)
library(gmm)
library(mrbayes)
library(rjags)

# Function to perform a single simulation run
run_simulation_OSMR_MRBAYES <- function(run) {

  set.seed(run)

  # Sample size (n), number of response variables (p), number of instrument variables (k) and number of confounders (l_star)
  n = 30000
  p =  10
  k = 3 * p
  l_star = 5

  # Variance
  var = 9

  # Create scale free netwrok
  g <- sample_pa(n = p, directed = FALSE)
  g <- as.directed(g, mode = "mutual")  # Make edges bi-directional -> cycles introduced

  # Assuming 'g' is your directed graph
  adj_matrix <- as_adjacency_matrix(g, sparse = FALSE)

  # Transpose so that [i,j] means edge from j to i
  adj_matrix_reversed <- t(adj_matrix)

  # Generate A matrix
  A = matrix(sample(c(-0.1, 0.1), p^2, replace = TRUE), p, p)

  # Make the diagonals of A to be 0
  A = A * adj_matrix_reversed

  B <- matrix(0, nrow = p, ncol = 3 * p)  # initialize zero matrix

  block_size <- 3

  for (i in 1:p) {
    start_col <- (i - 1) * block_size + 1
    end_col <- start_col + block_size - 1
    B[i, start_col:end_col] <- 1
  }


  # Generate C matrix i.e. influence of confounders on response variables, each row corresponds to a particular response
  C = matrix(sample(c(-1, 1), p * l_star, replace = TRUE), p, l_star)

  # Make the network sparse
  C[sample(which(C!=0), length(which(C!=0)) * 1 / 2)] = 0


  U = matrix(rnorm(n * l_star, 0, 1), nrow = n, ncol = l_star)

  # Calculate variance-covariance matrix
  Sigma = var * diag(p)

  # Calculate (I_p - A)^(-1)
  Mult_Mat = solve(diag(p) - A)

  Variance = Mult_Mat %*% Sigma %*% t(Mult_Mat)

  # True variance-covariance matrix
  Sigma_True = C %*% t(C) + Sigma

  # Set seed
  set.seed(run)

  # Data generation
  X = matrix(rnorm(n * k, 0, 1), nrow = n, ncol = k)

  Y = matrix(0, nrow = n, ncol = p)

  for (i in 1:n) {

    Y[i, ] = mvrnorm(n = 1, Mult_Mat %*% (B %*% X[i, ] + C %*% U[i, ]), Variance)

  }


  # Assume Y, X, and D are already defined
  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(X)

  D = B


  ### OneSampleMR ###
  pval_threshold <- 0.05
  adjacency_matrix <- matrix(0, p, p)
  effect_matrix <- matrix(0, p, p)

  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) next
      iv_idx <- which(D[i, ] == 1)
      if (length(iv_idx) < 1) next
      G <- X[, iv_idx, drop = FALSE]
      colnames(G) <- paste0("z", iv_idx)
      exposure <- Y[, i]
      outcome <- Y[, j]
      df <- data.frame(y = outcome, x = exposure, G)
      instruments_formula <- reformulate(paste0("z", iv_idx))
      formula <- y ~ x

      fit <- tryCatch({
        tsri(formula = formula,
             instruments = instruments_formula,
             data = df,
             link = "identity")
      }, error = function(e) NULL)

      if (!is.null(fit) && !is.null(fit$fit)) {
        coef_table <- summary(fit$fit)$coefficients
        if ("x" %in% rownames(coef_table)) {
          effect_matrix[i, j] <- coef_table["x", "Estimate"]
        }
      }
    }
  }

  # Compute OneSampleMR metrics (off-diagonal only)
  diff1 <- abs(A - t(effect_matrix))
  off_diag <- row(diff1) != col(diff1)
  T1 <- c(
    max(diff1[off_diag], na.rm = TRUE),
    mean(diff1[off_diag], na.rm = TRUE),
    mean(diff1[off_diag]^2, na.rm = TRUE)
  )

  ### MrBAYES ###
  causal_matrix <- matrix(0, nrow = p, ncol = p)

  get_snp_effects <- function(expr, snp_matrix) {
    betas <- se <- numeric(ncol(snp_matrix))
    for (m in seq_len(ncol(snp_matrix))) {
      model <- summary(lm(expr ~ snp_matrix[, m]))
      betas[m] <- coef(model)[2, 1]
      se[m] <- coef(model)[2, 2]
    }
    list(beta = betas, se = se)
  }

  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) next
      snp_indices <- which(D[i, ] == 1)
      if (length(snp_indices) < 2) next

      X_sub <- X[, snp_indices, drop = FALSE]
      exp_effects <- get_snp_effects(Y[, i], X_sub)
      out_effects <- get_snp_effects(Y[, j], X_sub)

      stan_data <- mr_format(
        rsid = paste0("SNP", snp_indices),
        xbeta = exp_effects$beta,
        ybeta = out_effects$beta,
        xse = exp_effects$se,
        yse = out_effects$se
      )

      bayes_model <- tryCatch({
        mr_ivw_rjags(
          object = stan_data,
          prior = "default",
          n.chains = 2,
          n.burn = 10000,
          n.iter = 50000
        )
      }, error = function(e) NULL)

      if (!is.null(bayes_model)) {
        causal_matrix[i, j] <- mean(bayes_model$CausalEffect)
      }
    }
  }

  # Compute MrBAYES metrics (off-diagonal only)
  diff2 <- abs(A - t(causal_matrix))
  off_diag <- row(diff2) != col(diff2)
  T2 <- c(
    max(diff2[off_diag], na.rm = TRUE),
    mean(diff2[off_diag], na.rm = TRUE),
    mean(diff2[off_diag]^2, na.rm = TRUE)
  )

  list(T1 = T1, T2 = T2)
}

# Run simulations in parallel
n_runs <- 20
n_cores <- detectCores() - 1
all_results <- mclapply(1:n_runs, run_simulation_OSMR_MRBAYES, mc.cores = n_cores)

# Extract and combine results
T1 <- do.call(rbind, lapply(all_results, function(x) x$T1))
T2 <- do.call(rbind, lapply(all_results, function(x) x$T2))

# Final output: average errors across runs
colnames(T1) <- colnames(T2) <- c("Maximum Absolute Deviation", "Mean Absolute Deviation", "Mean Squared Error")
cat("OneSampleMR:\n")
print(colMeans(T1))
apply(T1, MARGIN = 2, FUN = sd) / sqrt(nrow(T1)) * sqrt(nrow(T1) - 1)

cat("MrBAYES:\n")
print(colMeans(T2))
apply(T2, MARGIN = 2, FUN = sd) / sqrt(nrow(T2)) * sqrt(nrow(T2) - 1)

T_OneSampleMR = T1
T_MrBayes = T2



##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
###############################################################################
## Block 3 — MR.RGM with correlated residuals (allows unmeasured confounding)
## Model choice:
##   • prior = "Spike and Slab"
##   • SigmaStarModel = "SSSL"
##     (permits correlated residuals among traits, capturing latent/confounding
##      dependence not explained by instruments)
##
## What this block does:
##   Runs MR.RGM on summary matrices (Syy, Syx, Sxx) with instrument design D,
##   producing posterior estimates of the direct-effect matrix among traits.
##
## Effect estimate used:
##   • AEst (estimated direct-effect matrix among traits)
##
## Evaluation metric (causal effect accuracy):
##   Compare AEst to the true direct-effect matrix A and report (off-diagonal only):
##     • Maximum Absolute Deviation
##     • Mean Absolute Deviation
##     • Mean Squared Error
###############################################################################

# Load necessary libraries
library(MR.RGM)

# Function to perform a single simulation run
run_simulation_RGMCONF <- function(run) {

  set.seed(run)

  # Sample size (n), number of response variables (p), number of instrument variables (k) and number of confounders (l_star)
  n = 30000
  p =  10
  k = 3 * p
  l_star = 5

  # Variance
  var = 9

  # Create scale free netwrok
  g <- sample_pa(n = p, directed = FALSE)
  g <- as.directed(g, mode = "mutual")  # Make edges bi-directional -> cycles introduced

  # Assuming 'g' is your directed graph
  adj_matrix <- as_adjacency_matrix(g, sparse = FALSE)

  # Transpose so that [i,j] means edge from j to i
  adj_matrix_reversed <- t(adj_matrix)

  # Generate A matrix
  A = matrix(sample(c(-0.1, 0.1), p^2, replace = TRUE), p, p)

  # Make the diagonals of A to be 0
  A = A * adj_matrix_reversed

  B <- matrix(0, nrow = p, ncol = 3 * p)  # initialize zero matrix

  block_size <- 3

  for (i in 1:p) {
    start_col <- (i - 1) * block_size + 1
    end_col <- start_col + block_size - 1
    B[i, start_col:end_col] <- 1
  }


  # Generate C matrix i.e. influence of confounders on response variables, each row corresponds to a particular response
  C = matrix(sample(c(-1, 1), p * l_star, replace = TRUE), p, l_star)

  # Make the network sparse
  C[sample(which(C!=0), length(which(C!=0)) * 1 / 2)] = 0


  U = matrix(rnorm(n * l_star, 0, 1), nrow = n, ncol = l_star)

  # Calculate variance-covariance matrix
  Sigma = var * diag(p)

  # Calculate (I_p - A)^(-1)
  Mult_Mat = solve(diag(p) - A)

  Variance = Mult_Mat %*% Sigma %*% t(Mult_Mat)

  # True variance-covariance matrix
  Sigma_True = C %*% t(C) + Sigma

  # Set seed
  set.seed(run)

  # Data generation
  X = matrix(rnorm(n * k, 0, 1), nrow = n, ncol = k)

  Y = matrix(0, nrow = n, ncol = p)

  for (i in 1:n) {

    Y[i, ] = mvrnorm(n = 1, Mult_Mat %*% (B %*% X[i, ] + C %*% U[i, ]), Variance)

  }


  # Assume Y, X, and D are already defined
  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(X)

  D = B

  ### MR.RGM with confounding ###
  ## Calculate S_YY, S_YX, S_XX
  S_YY = t(Y) %*% Y / n
  S_YX = t(Y) %*% X / n
  S_XX = t(X) %*% X / n

  # Run MR.RGM with confounders
  Output = RGM(Syy = S_YY, Syx = S_YX, Sxx = S_XX,
               D = D, n = n, nIter = 50000, nBurnin = 10000, Thin = 10,
               prior = "Spike and Slab", SigmaStarModel = "SSSL")


  # Compute MR.RGM metrics (off-diagonal only)
  diff <- abs(A - Output$AEst)
  off_diag <- row(diff) != col(diff)
  T <- c(
    max(diff[off_diag], na.rm = TRUE),
    mean(diff[off_diag], na.rm = TRUE),
    mean(diff[off_diag]^2, na.rm = TRUE)
  )

  list(T = T)
}

# Run simulations in parallel
n_runs <- 20
n_cores <- detectCores() - 1
all_results <- mclapply(1:n_runs, run_simulation_RGMCONF, mc.cores = n_cores)

# Extract and combine results
T <- do.call(rbind, lapply(all_results, function(x) x$T))

# Final output: average errors across runs
colnames(T) <- c("Maximum Absolute Deviation", "Mean Absolute Deviation", "Mean Squared Error")
cat("MR.RGM:\n")
print(colMeans(T))
apply(T, MARGIN = 2, FUN = sd) / sqrt(nrow(T)) * sqrt(nrow(T) - 1)


T_MR.RGM = T



##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
###############################################################################
## Block 4 — MR.RGM with independent residuals (no residual correlation)
## Model choice:
##   • prior = "Spike and Slab"
##   • SigmaStarModel = "diagonal"
##     (assumes residual errors across traits are independent)
##
## Interpretation:
##   Baseline MR.RGM variant without modeling residual correlation among traits,
##   i.e., no unmeasured confounding represented through correlated residuals.
##
## Effect estimate used:
##   • AEst (estimated direct-effect matrix among traits)
##
## Evaluation metric (causal effect accuracy):
##   Compare AEst to the true direct-effect matrix A and report (off-diagonal only):
##     • Maximum Absolute Deviation
##     • Mean Absolute Deviation
##     • Mean Squared Error
###############################################################################

# Load necessary libraries
library(MR.RGM)

# Function to perform a single simulation run
run_simulation_RGMNOCONF <- function(run) {

  set.seed(run)

  # Sample size (n), number of response variables (p), number of instrument variables (k) and number of confounders (l_star)
  n = 30000
  p =  10
  k = 3 * p
  l_star = 5

  # Variance
  var = 9

  # Create scale free netwrok
  g <- sample_pa(n = p, directed = FALSE)
  g <- as.directed(g, mode = "mutual")  # Make edges bi-directional -> cycles introduced

  # Assuming 'g' is your directed graph
  adj_matrix <- as_adjacency_matrix(g, sparse = FALSE)

  # Transpose so that [i,j] means edge from j to i
  adj_matrix_reversed <- t(adj_matrix)

  # Generate A matrix
  A = matrix(sample(c(-0.1, 0.1), p^2, replace = TRUE), p, p)

  # Make the diagonals of A to be 0
  A = A * adj_matrix_reversed

  B <- matrix(0, nrow = p, ncol = 3 * p)  # initialize zero matrix

  block_size <- 3

  for (i in 1:p) {
    start_col <- (i - 1) * block_size + 1
    end_col <- start_col + block_size - 1
    B[i, start_col:end_col] <- 1
  }


  # Generate C matrix i.e. influence of confounders on response variables, each row corresponds to a particular response
  C = matrix(sample(c(-1, 1), p * l_star, replace = TRUE), p, l_star)

  # Make the network sparse
  C[sample(which(C!=0), length(which(C!=0)) * 1 / 2)] = 0


  U = matrix(rnorm(n * l_star, 0, 1), nrow = n, ncol = l_star)

  # Calculate variance-covariance matrix
  Sigma = var * diag(p)

  # Calculate (I_p - A)^(-1)
  Mult_Mat = solve(diag(p) - A)

  Variance = Mult_Mat %*% Sigma %*% t(Mult_Mat)

  # True variance-covariance matrix
  Sigma_True = C %*% t(C) + Sigma

  # Set seed
  set.seed(run)

  # Data generation
  X = matrix(rnorm(n * k, 0, 1), nrow = n, ncol = k)

  Y = matrix(0, nrow = n, ncol = p)

  for (i in 1:n) {

    Y[i, ] = mvrnorm(n = 1, Mult_Mat %*% (B %*% X[i, ] + C %*% U[i, ]), Variance)

  }


  # Assume Y, X, and D are already defined
  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(X)

  D = B


  ### MR.RGM without confounding ###
  ## Calculate S_YY, S_YX, S_XX
  S_YY = t(Y) %*% Y / n
  S_YX = t(Y) %*% X / n
  S_XX = t(X) %*% X / n

  # Run MR.RGM without confounder
  Output = RGM(Syy = S_YY, Syx = S_YX, Sxx = S_XX,
                D = D, n = n, nIter = 50000, nBurnin = 10000, Thin = 10,
                prior = "Spike and Slab", SigmaStarModel = "diagonal")


  # Compute MR.RGM metrics (off-diagonal only)
  diff <- abs(A - Output$AEst)
  off_diag <- row(diff) != col(diff)
  T <- c(
    max(diff[off_diag], na.rm = TRUE),
    mean(diff[off_diag], na.rm = TRUE),
    mean(diff[off_diag]^2, na.rm = TRUE)
  )

  list(T = T)
}

# Run simulations in parallel
n_runs <- 20
n_cores <- detectCores() - 1
all_results <- mclapply(1:n_runs, run_simulation_RGMNOCONF, mc.cores = n_cores)

# Extract and combine results
T <- do.call(rbind, lapply(all_results, function(x) x$T))

# Final output: average errors across runs
colnames(T) <- c("Maximum Absolute Deviation", "Mean Absolute Deviation", "Mean Squared Error")
cat("MR.RGM(No Conf.):\n")
print(colMeans(T))
apply(T, MARGIN = 2, FUN = sd) / sqrt(nrow(T)) * sqrt(nrow(T) - 1)

T_MR.RGM_NoConf = T



