###############################################################################
## Block 1 — MendelianRandomization package
## Methods:
##   • Simple median
##   • Weighted median
##   • IVW
##
## Implementation detail:
##   For each ordered pair (i -> j), regress Y_i ~ G and Y_j ~ G across SNPs in D[i,]
##   and call mr_allmethods(). Then threshold p-values to form an adjacency matrix.
###############################################################################

# Mendelian Randomization
library(doParallel)
library(foreach)
library(MendelianRandomization)
library(pROC)
library(igraph)
library(MASS)

# Setup
num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Set number of simulations
num_sim <- 20

# Run in parallel
results <- foreach(run = 1:num_sim, .packages = c("MendelianRandomization", "pROC", "MASS", "igraph")) %dopar% {

  set.seed(run)

  # Sample size (n), number of response variables (p), number of instrument variables (k) and number of confounders (l_star)
  n = 30000
  p =  10
  k = 3 * p
  l_star = 5

  # Variance
  var = 9


  # Generate small-world undirected graph
  sw_undirected <- sample_smallworld(dim = 1, size = p, nei = 2, p = 0.1)

  # Convert to directed graph
  edge_list <- as_edgelist(sw_undirected)

  directed_edges <- t(apply(edge_list, 1, function(e) {
    if (runif(1) > 0.5) e else rev(e)
  }))

  # Create directed graph
  sw_directed <- graph_from_edgelist(directed_edges, directed = TRUE)


  # Assuming 'g' is your directed graph
  adj_matrix <- as_adjacency_matrix(sw_directed, sparse = FALSE)

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


  # Create D
  D <- B

  # Again Calculate k
  k = ncol(B)



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

  # Get the true grpah
  True_Graph = ((A) != 0) * 1


  library(MendelianRandomization)

  # Assume Y, X, and D are already defined
  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(X)

  causal_matrix_SimpleMedian <- matrix(0, p, p)
  pval_matrix_SimpleMedian <- matrix(1, p, p)
  causal_matrix_WeightedMedian <- matrix(0, p, p)
  pval_matrix_WeightedMedian <- matrix(1, p, p)
  causal_matrix_IVW <- matrix(0, p, p)
  pval_matrix_IVW <- matrix(1, p, p)

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

    mr_input_obj <- mr_input(bx = bx_data[, 1], bxse = bx_data[, 2],
                             by = by_data[, 1], byse = by_data[, 2],
                             snps = snp_names)
    tryCatch(mr_allmethods(mr_input_obj, method = "main", iterations = 50000),
             error = function(e) NULL)
  }

  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) next
      result <- run_mr_allmethods(X, Y, D, i, j)
      if (is.null(result)) next

      res_df <- result@Values

      causal_matrix_SimpleMedian[i, j] <- res_df[res_df$Method == "Simple median", "Estimate"]
      causal_matrix_WeightedMedian[i, j] <- res_df[res_df$Method == "Weighted median", "Estimate"]
      causal_matrix_IVW[i, j] <- res_df[res_df$Method == "IVW", "Estimate"]

      pval_matrix_SimpleMedian[i, j] <- res_df[res_df$Method == "Simple median", "P-value"]
      pval_matrix_WeightedMedian[i, j] <- res_df[res_df$Method == "Weighted median", "P-value"]
      pval_matrix_IVW[i, j] <- res_df[res_df$Method == "IVW", "P-value"]
    }
  }

  off_diag <- which(row(True_Graph) != col(True_Graph))
  TG <- as.vector((True_Graph))[off_diag]

  # Define metric output for each method
  make_metrics <- function(pval_matrix) {
    predicted <- as.vector((t(pval_matrix) < 0.05))[off_diag] * 1
    scores <- as.vector((t(1 - pval_matrix)))[off_diag]
    met <- metrics_binary_graph(TG, predicted)
    AUC <- tryCatch(pROC::auc(pROC::roc(TG, scores)), error = function(e) NA)
    met[1] <- AUC
    met
  }

  list(
    T1 = make_metrics(pval_matrix_SimpleMedian),
    T2 = make_metrics(pval_matrix_WeightedMedian),
    T3 = make_metrics(pval_matrix_IVW)
  )
}

# Convert results to matrices
T_SimpleMedian <- do.call(rbind, lapply(results, `[[`, "T1"))
T_WeightedMedian <- do.call(rbind, lapply(results, `[[`, "T2"))
T_IVW <- do.call(rbind, lapply(results, `[[`, "T3"))

# Stop cluster
stopCluster(cl)

# Final averages
colMeans(T_SimpleMedian)
colMeans(T_WeightedMedian)
colMeans(T_IVW)


# Get SD
apply(T_SimpleMedian, MARGIN = 2, FUN = sd) / sqrt(nrow(T_SimpleMedian)) * sqrt(nrow(T_SimpleMedian) - 1)
apply(T_WeightedMedian, MARGIN = 2, FUN = sd) / sqrt(nrow(T_WeightedMedian)) * sqrt(nrow(T_WeightedMedian) - 1)
apply(T_IVW, MARGIN = 2, FUN = sd) / sqrt(nrow(T_IVW)) * sqrt(nrow(T_IVW) - 1)






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
## Block 2 — OneSampleMR + MrBayes baselines
## Methods:
##   • OneSampleMR (TSRI / tsri)
##   • MrBayes / Bayesian IVW via JAGS (mr_ivw_rjags)
##
## Note:
##   This block is computationally heavy because it runs per-pair models (i -> j),
##   and the Bayesian method runs MCMC for each pair.
###############################################################################

# Load necessary libraries
library(parallel)
library(OneSampleMR)
library(pROC)
library(mrbayes)
library(rjags)

# Number of cores to use
num_cores <- detectCores() - 1  # Use one fewer than total to avoid freezing your system
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Set number of simulations
num_sim <- 20


# Define the simulation function
run_simulation_OSMR_MRBAYES <- function(run) {

  set.seed(run)

  # Sample size (n), number of response variables (p), number of instrument variables (k) and number of confounders (l_star)
  n = 30000
  p =  10
  k = 3 * p
  l_star = 5

  # Variance
  var = 9


  # Generate small-world undirected graph
  sw_undirected <- sample_smallworld(dim = 1, size = p, nei = 2, p = 0.1)

  # Convert to directed graph
  edge_list <- as_edgelist(sw_undirected)

  directed_edges <- t(apply(edge_list, 1, function(e) {
    if (runif(1) > 0.5) e else rev(e)
  }))

  # Create directed graph
  sw_directed <- graph_from_edgelist(directed_edges, directed = TRUE)


  # Assuming 'g' is your directed graph
  adj_matrix <- as_adjacency_matrix(sw_directed, sparse = FALSE)

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



  # Create D
  D <- B

  # Again Calculate k
  k = ncol(B)



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


  # Get the true grpah
  True_Graph = ((A) != 0) * 1

  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(X)

  # Get off-diagonal indices
  off_diag_indices <- which(row(True_Graph) != col(True_Graph))

  ## =============== METHOD 1: OneSampleMR ================
  adjacency_matrix <- matrix(0, p, p)
  pval_matrix1 <- matrix(1, p, p)

  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) next
      iv_idx <- which(D[i, ] == 1)
      if (length(iv_idx) < 1) next

      G <- X[, iv_idx, drop = FALSE]
      instrument_names <- paste0("z", iv_idx)
      colnames(G) <- instrument_names

      df <- data.frame(y = Y[, j], x = Y[, i], G)
      instruments_formula <- reformulate(termlabels = instrument_names)
      formula <- y ~ x

      fit <- tryCatch({
        tsri(formula, instruments_formula, df, link = "identity")
      }, error = function(e) NULL)

      if (!is.null(fit) && !is.null(fit$fit)) {
        coef_table <- summary(fit$fit)$coefficients
        if ("x" %in% rownames(coef_table)) {
          pval <- coef_table["x", "Pr(>|t|)"]
          pval_matrix1[i, j] <- pval
          if (!is.na(pval) && pval < 0.05) {
            adjacency_matrix[i, j] <- 1
          }
        }
      }
    }
  }

  metrics1_res1 <- metrics_binary_graph(
    as.vector((True_Graph))[off_diag_indices],
    as.vector((t(pval_matrix1) < 0.05))[off_diag_indices] * 1
  )
  auc1 <- pROC::auc(pROC::roc(
    as.vector((True_Graph)[off_diag_indices]),
    as.vector((t(1 - pval_matrix1))[off_diag_indices])
  ))

  ## =============== METHOD 2: Bayesian MR ================
  pval_matrix2 <- matrix(1, p, p)

  get_snp_effects <- function(expr_vector, snp_matrix) {
    betas <- se <- numeric(ncol(snp_matrix))
    for (m in seq_len(ncol(snp_matrix))) {
      model <- summary(lm(expr_vector ~ snp_matrix[, m]))
      betas[m] <- coef(model)[2, 1]
      se[m] <- coef(model)[2, 2]
    }
    return(list(beta = betas, se = se))
  }

  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) next
      snp_indices <- which(D[i, ] == 1)
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

      bayes_model <- mr_ivw_rjags(
        object = stan_data,
        prior = "default",
        n.chains = 2,
        n.burn = 10000,
        n.iter = 50000
      )

      samples_df <- do.call(rbind.data.frame, bayes_model$samples)
      p_value <- 2 * min(mean(samples_df$Estimate > 0), mean(samples_df$Estimate < 0))
      pval_matrix2[i, j] <- p_value
    }
  }

  metrics1_res2 <- metrics_binary_graph(
    as.vector((True_Graph))[off_diag_indices],
    as.vector((t(pval_matrix2) < 0.05))[off_diag_indices] * 1
  )
  auc2 <- pROC::auc(pROC::roc(
    as.vector((True_Graph)[off_diag_indices]),
    as.vector((t(1 - pval_matrix2))[off_diag_indices])
  ))

  return(list(
    T1 = c(auc1, metrics1_res1[-1]),
    T2 = c(auc2, metrics1_res2[-1])
  ))
}

# Run in parallel
results <- mclapply(1:num_sim, run_simulation_OSMR_MRBAYES, mc.cores = num_cores)


# Stop cluster
stopCluster(cl)

# Combine results
T_OneSampleMR <- do.call(rbind, lapply(results, function(x) x$T1))
T_MrBayes <- do.call(rbind, lapply(results, function(x) x$T2))

# Summary
colMeans(T_OneSampleMR)
colMeans(T_MrBayes)

# Get SD
apply(T_OneSampleMR, MARGIN = 2, FUN = sd) / sqrt(nrow(T_OneSampleMR)) * sqrt(nrow(T_OneSampleMR) - 1)
apply(T_MrBayes, MARGIN = 2, FUN = sd) / sqrt(nrow(T_MrBayes)) * sqrt(nrow(T_MrBayes) - 1)










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
## Block 3 — MR.RGM with correlated errors (unmeasured confounding allowed)
## Model choice:
##   SigmaStarModel = "SSSL"  (allows correlated residuals among responses)
##
## Graph estimate used:
##   GammaEst (posterior edge probabilities); threshold at 0.5 for adjacency.
###############################################################################

# Load necessary libraries
library(MR.RGM)

# Number of cores to use
num_cores <- detectCores() - 1  # Use one fewer than total to avoid freezing your system
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Set number of simulations
num_sim <- 20


# Define the simulation function
run_simulation_RGMCONF <- function(run) {

  set.seed(run)

  # Sample size (n), number of response variables (p), number of instrument variables (k) and number of confounders (l_star)
  n = 30000
  p =  10
  k = 3 * p
  l_star = 5

  # Variance
  var = 9


  # Generate small-world undirected graph
  sw_undirected <- sample_smallworld(dim = 1, size = p, nei = 2, p = 0.1)

  # Convert to directed graph
  edge_list <- as_edgelist(sw_undirected)

  directed_edges <- t(apply(edge_list, 1, function(e) {
    if (runif(1) > 0.5) e else rev(e)
  }))

  # Create directed graph
  sw_directed <- graph_from_edgelist(directed_edges, directed = TRUE)


  # Assuming 'g' is your directed graph
  adj_matrix <- as_adjacency_matrix(sw_directed, sparse = FALSE)

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



  # Create D
  D <- B

  # Again Calculate k
  k = ncol(B)



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

  # Get the true grpah
  True_Graph = ((A) != 0) * 1



  # Get off-diagonal indices
  off_diag_indices <- which(row(True_Graph) != col(True_Graph))


  ## =============== METHOD: MR.RGM with confounders ================
  ## Calculate S_YY, S_YX, S_XX
  S_YY = t(Y) %*% Y / n
  S_YX = t(Y) %*% X / n
  S_XX = t(X) %*% X / n

  # Run MR.RGM with confounders
  Output = RGM(Syy = S_YY, Syx = S_YX, Sxx = S_XX,
               D = D, n = n, nIter = 50000, nBurnin = 10000, Thin = 10,
               prior = "Spike and Slab", SigmaStarModel = "SSSL")

  metrics1_res <- metrics_binary_graph(
    as.vector((True_Graph))[off_diag_indices],
    as.vector((Output$GammaEst > 0.5))[off_diag_indices] * 1
  )
  auc <- pROC::auc(pROC::roc(
    as.vector((True_Graph)[off_diag_indices]),
    as.vector((Output$GammaEst)[off_diag_indices])
  ))


  return(list(
    T = c(auc, metrics1_res[-1])
  ))
}

# Run in parallel
results <- mclapply(1:num_sim, run_simulation_RGMCONF, mc.cores = num_cores)


# Stop cluster
stopCluster(cl)

# Combine results
T_MR.RGM <- do.call(rbind, lapply(results, function(x) x$T))

# Summary
colMeans(T_MR.RGM)


# Get SD
apply(T_MR.RGM, MARGIN = 2, FUN = sd) / sqrt(nrow(T_MR.RGM)) * sqrt(nrow(T_MR.RGM) - 1)








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
## Block 4 — MR.RGM without correlated errors (independent residuals)
## Model choice:
##   SigmaStarModel = "diagonal"
##
## Interpretation:
##   This is the “no unmeasured confounding through residual correlation” baseline.
###############################################################################

# Load necessary libraries
library(MR.RGM)


# Number of cores to use
num_cores <- detectCores() - 1  # Use one fewer than total to avoid freezing your system
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Set number of simulations
num_sim <- 20


# Define the simulation function
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

  # Get the true grpah
  True_Graph = ((A) != 0) * 1

  D = B

  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(X)

  # Get off-diagonal indices
  off_diag_indices <- which(row(True_Graph) != col(True_Graph))



  ## =============== METHOD 3: MR.RGM without confounders ================
  ## Calculate S_YY, S_YX, S_XX
  S_YY = t(Y) %*% Y / n
  S_YX = t(Y) %*% X / n
  S_XX = t(X) %*% X / n

  # Run MR.RGM without confounder
  Output =  RGM(Syy = S_YY, Syx = S_YX, Sxx = S_XX, D = D,
                n = n, nIter = 50000, nBurnin = 10000, Thin = 10,
                prior = "Spike and Slab", SigmaStarModel = "diagonal")

  metrics1_res <- metrics_binary_graph(
    as.vector((True_Graph))[off_diag_indices],
    as.vector((Output$GammaEst > 0.5))[off_diag_indices] * 1
  )
  auc <- pROC::auc(pROC::roc(
    as.vector((True_Graph)[off_diag_indices]),
    as.vector((Output$GammaEst)[off_diag_indices])
  ))

  return(list(
    T = c(auc, metrics1_res[-1])
  ))
}

# Run in parallel
results <- mclapply(1:num_sim, run_simulation_RGMNOCONF, mc.cores = num_cores)


# Stop cluster
stopCluster(cl)

# Extract T entries and bind them
T_MR.RGM_NoConf <- do.call(rbind, lapply(results, function(x) x$T))

# Compute summary statistics
colMeans(T_MR.RGM_NoConf)


# Get SD
apply(T_MR.RGM_NoConf, MARGIN = 2, FUN = sd) / sqrt(nrow(T_MR.RGM_NoConf)) * sqrt(nrow(T_MR.RGM_NoConf) - 1)




