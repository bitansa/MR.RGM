# Load necessary library
library(pROC)

# ------------------------------------------------------------------------------
# Evaluate performance metrics for binary graph recovery
#
# Inputs:
#   true_graph     : binary vector (0/1), ground-truth edges
#   estimate_graph : binary vector (0/1), estimated edges
#
# Returns:
#   Named numeric vector with:
#     - AUC : Area Under the ROC Curve
#     - TPR : True Positive Rate (Recall)
#     - FDR : False Discovery Rate
#     - MCC : Matthews Correlation Coefficient
#
# ------------------------------------------------------------------------------
metrics_binary_graph <- function(true_graph, estimate_graph) {

  ## ---- Sanity checks ----
  if (length(true_graph) != length(estimate_graph)) {
    stop("Input vectors must be of the same length")
  }

  if (!all(true_graph %in% c(0, 1)) || !all(estimate_graph %in% c(0, 1))) {
    stop("Inputs must be binary vectors (0/1)")
  }

  ## ---- Confusion matrix ----
  TP <- sum(true_graph == 1 & estimate_graph == 1)
  FP <- sum(true_graph == 0 & estimate_graph == 1)
  FN <- sum(true_graph == 1 & estimate_graph == 0)
  TN <- sum(true_graph == 0 & estimate_graph == 0)

  ## ---- Metrics ----
  TPR <- TP / max(1, sum(true_graph))          # True Positive Rate
  FDR <- FP / max(1, (TP + FP))                # False Discovery Rate

  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  MCC <- if (denom == 0) 0 else (TP * TN - FP * FN) / denom

  ## ---- AUC ----
  AUC <- tryCatch({
    roc_obj <- pROC::roc(true_graph, estimate_graph, quiet = TRUE)
    as.numeric(pROC::auc(roc_obj))
  }, error = function(e) {
    NA_real_
  })

  c(AUC = AUC, TPR = TPR, FDR = FDR, MCC = MCC)
}



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
## CONF. RECOVERY — SCALE-FREE NETWORKS (MR.RGM only)
##   - Evaluate confounding recovery using MR.RGM's ZEst vs true Sigma_True
##   - Metrics: AUC (continuous ZEst), TPR/FDR/MCC (thresholded ZEst > 0.5)
##   - Run for p ∈ {5,10} and n ∈ {500,1000,10000,30000}
##   - Output:
##       (1) Table of mean & SD for AUC/TPR/FDR/MCC by (p,n)
##       (2) AUC boxplots across sample sizes, separately for each p
###############################################################################

## ----------------------------
## Libraries
## ----------------------------
library(MR.RGM)
library(igraph)
library(MASS)
library(pROC)
library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)  # for comma()

## ----------------------------
## One replicate: scale-free confounding recovery
## ----------------------------
ScaleFree_Conf = function(run, Sample_Size, Network_Size) {

  # Set seed
  set.seed(run)


  # Sample size (n), number of response variables (p), number of instrument variables (k) and number of confounders (l_star)
  n = Sample_Size
  p =  Network_Size
  k = 3 * p
  l_star = 5

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

  # Variance
  var = 9


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

  ## Run MR.RGM
  ## Calculate S_YY, S_YX, S_XX
  S_YY = t(Y) %*% Y / n
  S_YX = t(Y) %*% X / n
  S_XX = t(X) %*% X / n

  # Run MR.RGM with confounders
  Output = RGM(Syy = S_YY, Syx = S_YX, Sxx = S_XX,
               D = D, n = n, nIter = 50000, nBurnin = 10000, Thin = 10,
               prior = "Spike and Slab", SigmaStarModel = "SSSL")


  #############################################
  # Calculation with Z
  Z_Est = Output$ZEst

  Z_Est_UpperTriangular = Z_Est[upper.tri(Z_Est, diag = FALSE)]
  Sigma_True_UpperTriangular = Sigma_True[upper.tri(Sigma_True, diag = FALSE)]


  # Take absolute values and normalize it
  Sigma_True_Modified = abs(Sigma_True_UpperTriangular) / max(abs(Sigma_True_UpperTriangular))

  Conf = metrics_binary_graph(as.vector((Sigma_True_Modified > mean(Sigma_True_Modified)) * 1), as.vector((Z_Est_UpperTriangular > 0.5) * 1))

  Conf[1] = pROC::auc(pROC::roc((Sigma_True_Modified > mean(Sigma_True_Modified)) * 1, Z_Est_UpperTriangular))

  return(list(Conf = Conf))

}

## ----------------------------
## Run full grid (p x n) and collect results
## ----------------------------
sample_sizes  <- c(500, 1000, 10000, 30000)
network_sizes <- c(5, 10)
n_runs  <- 20
n_cores <- max(1, parallel::detectCores() - 1)

grid <- expand.grid(
  NetworkSize = network_sizes,
  SampleSize  = sample_sizes,
  stringsAsFactors = FALSE
)

results_all <- lapply(seq_len(nrow(grid)), function(idx) {

  p <- grid$NetworkSize[idx]
  n <- grid$SampleSize[idx]

  # Run n_runs replicates in parallel; each replicate returns list(Conf = c(AUC,TPR,FDR,MCC))
  res_list <- parallel::mclapply(
    X = 1:n_runs,
    FUN = function(run) ScaleFree_Conf(run = run, Sample_Size = n, Network_Size = p),
    mc.cores = n_cores
  )

  # Extract numeric vectors and row-bind into a matrix
  met_mat <- do.call(rbind, lapply(res_list, function(x) x$Conf))

  # Ensure column names exist (critical for summarise() later)
  met_mat <- as.matrix(met_mat)
  colnames(met_mat) <- c("AUC", "TPR", "FDR", "MCC")

  df <- as.data.frame(met_mat)
  df$Run <- seq_len(n_runs)
  df$NetworkSize <- p
  df$SampleSize  <- n

  df
})

Conf_df <- dplyr::bind_rows(results_all)

## ----------------------------
## (A) Table: mean & SD by (p,n)
## ----------------------------
Conf_summary <- Conf_df %>%
  dplyr::group_by(NetworkSize, SampleSize) %>%
  dplyr::summarise(
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_sd   = sd(AUC, na.rm = TRUE),
    TPR_mean = mean(TPR, na.rm = TRUE),
    TPR_sd   = sd(TPR, na.rm = TRUE),
    FDR_mean = mean(FDR, na.rm = TRUE),
    FDR_sd   = sd(FDR, na.rm = TRUE),
    MCC_mean = mean(MCC, na.rm = TRUE),
    MCC_sd   = sd(MCC, na.rm = TRUE),
    .groups = "drop"
  )

print(Conf_summary)

## ----------------------------
## (B) Plot: AUC boxplots across sample sizes (MR.RGM only), facet by p
## ----------------------------
## =========================================
## Scale-free (p = 5): Confounding recovery AUC — MR.RGM
## =========================================

library(dplyr)
library(ggplot2)
library(scales)  # for comma()

AUC_mr_sf_p5 <- Conf_df %>%
  filter(NetworkSize == 5) %>%
  mutate(SampleSize = as.integer(SampleSize))

x_lvls   <- sort(unique(AUC_mr_sf_p5$SampleSize))
x_labels <- setNames(paste0("n = ", comma(x_lvls)), x_lvls)

p_sf_conf_p5 <- ggplot(
  AUC_mr_sf_p5 %>% mutate(SampleSize = factor(SampleSize, levels = x_lvls)),
  aes(x = SampleSize, y = AUC)
) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             linewidth = 0.4, color = "grey60") +
  geom_boxplot(
    width = 0.6,
    fill = "#0072B2",
    color = "grey15",
    outlier.shape = 16,
    outlier.size = 1.8,
    alpha = 0.88
  ) +
  geom_jitter(width = 0.08, height = 0,
              size = 1.4, alpha = 0.5, color = "grey20") +
  scale_x_discrete(labels = x_labels) +
  coord_cartesian(ylim = c(0.4, 1.0)) +
  labs(
    title = "Confounding recovery AUC — MR.RGM across sample sizes (scale-free, p = 5)",
    x     = "Sample Size",
    y     = "Area Under Curve (AUC)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title.position = "panel",
    plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 8)),
    axis.title.y = element_text(size = 24, face = "bold", margin = margin(r = 8)),
    axis.text.x  = element_text(size = 24),
    axis.text.y  = element_text(size = 24)
  )

p_sf_conf_p5


## =========================================
## Scale-free (p = 10): Confounding recovery AUC — MR.RGM
## =========================================

AUC_mr_sf_p10 <- Conf_df %>%
  filter(NetworkSize == 10) %>%
  mutate(SampleSize = as.integer(SampleSize))

x_lvls   <- sort(unique(AUC_mr_sf_p10$SampleSize))
x_labels <- setNames(paste0("n = ", comma(x_lvls)), x_lvls)

p_sf_conf_p10 <- ggplot(
  AUC_mr_sf_p10 %>% mutate(SampleSize = factor(SampleSize, levels = x_lvls)),
  aes(x = SampleSize, y = AUC)
) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             linewidth = 0.4, color = "grey60") +
  geom_boxplot(
    width = 0.6,
    fill = "#0072B2",
    color = "grey15",
    outlier.shape = 16,
    outlier.size = 1.8,
    alpha = 0.88
  ) +
  geom_jitter(width = 0.08, height = 0,
              size = 1.4, alpha = 0.5, color = "grey20") +
  scale_x_discrete(labels = x_labels) +
  coord_cartesian(ylim = c(0.4, 1.0)) +
  labs(
    title = "Confounding recovery AUC — MR.RGM across sample sizes (scale-free, p = 10)",
    x     = "Sample Size",
    y     = "Area Under Curve (AUC)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title.position = "panel",
    plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 8)),
    axis.title.y = element_text(size = 24, face = "bold", margin = margin(r = 8)),
    axis.text.x  = element_text(size = 24),
    axis.text.y  = element_text(size = 24)
  )

p_sf_conf_p10



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
## CONF. RECOVERY — SMALL-WORLD NETWORKS (MR.RGM only)
##   - Evaluate confounding recovery using MR.RGM's ZEst vs true Sigma_True
##   - Metrics: AUC (continuous ZEst), TPR/FDR/MCC (thresholded ZEst > 0.5)
##   - Run for p ∈ {5,10} and n ∈ {500,1000,10000,30000}
##   - Output:
##       (1) Table of mean & SD for AUC/TPR/FDR/MCC by (p,n)
##       (2) AUC boxplots across sample sizes, separately for each p
##
###############################################################################

## ----------------------------
## Libraries
## ----------------------------
library(MR.RGM)
library(igraph)
library(MASS)
library(pROC)
library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)  # comma()

## ----------------------------
## One replicate: small-world confounding recovery
## ----------------------------
SmallWorld_Conf = function(run, Sample_Size, Network_Size) {

  set.seed(run)

  # Sample size (n), number of response variables (p), number of instrument variables (k) and number of confounders (l_star)
  n = Sample_Size
  p =  Network_Size
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

  # Create D and B for horizontal pleiotropy
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

  # Assume Y, X, and D are already defined
  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(X)

  ## Run MR.RGM
  ## Calculate S_YY, S_YX, S_XX
  S_YY = t(Y) %*% Y / n
  S_YX = t(Y) %*% X / n
  S_XX = t(X) %*% X / n

  # Run MR.RGM with confounders
  Output = RGM(Syy = S_YY, Syx = S_YX, Sxx = S_XX,
               D = D, n = n, nIter = 50000, nBurnin = 10000, Thin = 10,
               prior = "Spike and Slab", SigmaStarModel = "SSSL")

  #############################################
  # Calculation with Z
  Z_Est = Output$ZEst

  Z_Est_UpperTriangular = Z_Est[upper.tri(Z_Est, diag = FALSE)]
  Sigma_True_UpperTriangular = Sigma_True[upper.tri(Sigma_True, diag = FALSE)]

  # Take absolute values and normalize it
  Sigma_True_Modified = abs(Sigma_True_UpperTriangular) / max(abs(Sigma_True_UpperTriangular))

  Conf = metrics_binary_graph(as.vector((Sigma_True_Modified > mean(Sigma_True_Modified)) * 1),
                              as.vector((Z_Est_UpperTriangular > 0.5) * 1))

  Conf[1] = pROC::auc(pROC::roc((Sigma_True_Modified > mean(Sigma_True_Modified)) * 1, Z_Est_UpperTriangular))

  return(list(Conf = Conf))
}

## ----------------------------
## Run full grid (p x n) and collect results
## ----------------------------
sample_sizes  <- c(500, 1000, 10000, 30000)
network_sizes <- c(5, 10)

n_runs  <- 20
n_cores <- max(1, parallel::detectCores() - 1)

grid <- expand.grid(
  NetworkSize = network_sizes,
  SampleSize  = sample_sizes,
  stringsAsFactors = FALSE
)

results_all <- lapply(seq_len(nrow(grid)), function(idx) {

  p <- grid$NetworkSize[idx]
  n <- grid$SampleSize[idx]

  # Run n_runs replicates in parallel; keep only Conf1 (MR.RGM)
  res_list <- parallel::mclapply(
    X = 1:n_runs,
    FUN = function(run) SmallWorld_Conf(run = run, Sample_Size = n, Network_Size = p),
    mc.cores = n_cores
  )

  # Extract Conf1 (each is a numeric vector: AUC,TPR,FDR,MCC) and bind to matrix
  met_mat <- do.call(rbind, lapply(res_list, function(x) x$Conf))

  # Ensure column names exist (critical for summarise() later)
  met_mat <- as.matrix(met_mat)
  colnames(met_mat) <- c("AUC", "TPR", "FDR", "MCC")

  df <- as.data.frame(met_mat)
  df$Run <- seq_len(n_runs)
  df$NetworkSize <- p
  df$SampleSize  <- n

  df
})

Conf_df <- dplyr::bind_rows(results_all)

## ----------------------------
## (A) Table: mean & SD by (p,n)
## ----------------------------
Conf_summary <- Conf_df %>%
  dplyr::group_by(NetworkSize, SampleSize) %>%
  dplyr::summarise(
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_sd   = sd(AUC, na.rm = TRUE),
    TPR_mean = mean(TPR, na.rm = TRUE),
    TPR_sd   = sd(TPR, na.rm = TRUE),
    FDR_mean = mean(FDR, na.rm = TRUE),
    FDR_sd   = sd(FDR, na.rm = TRUE),
    MCC_mean = mean(MCC, na.rm = TRUE),
    MCC_sd   = sd(MCC, na.rm = TRUE),
    .groups = "drop"
  )

print(Conf_summary)

## ----------------------------
## (B) Plot: AUC boxplots across sample sizes (MR.RGM only)
## ----------------------------

## =========================================
## Small-world (p = 5): Confounding recovery AUC — MR.RGM
## =========================================
AUC_mr_sw_p5 <- Conf_df %>%
  filter(NetworkSize == 5) %>%
  mutate(SampleSize = as.integer(SampleSize))

x_lvls   <- sort(unique(AUC_mr_sw_p5$SampleSize))
x_labels <- setNames(paste0("n = ", comma(x_lvls)), x_lvls)

p_sw_conf_p5 <- ggplot(
  AUC_mr_sw_p5 %>% mutate(SampleSize = factor(SampleSize, levels = x_lvls)),
  aes(x = SampleSize, y = AUC)
) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             linewidth = 0.4, color = "grey60") +
  geom_boxplot(
    width = 0.6,
    fill = "#0072B2",
    color = "grey15",
    outlier.shape = 16,
    outlier.size = 1.8,
    alpha = 0.88
  ) +
  geom_jitter(width = 0.08, height = 0,
              size = 1.4, alpha = 0.5, color = "grey20") +
  scale_x_discrete(labels = x_labels) +
  coord_cartesian(ylim = c(0.4, 1.0)) +
  labs(
    title = "Confounding recovery AUC — MR.RGM across sample sizes (small-world, p = 5)",
    x     = "Sample Size",
    y     = "Area Under Curve (AUC)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title.position = "panel",
    plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 8)),
    axis.title.y = element_text(size = 24, face = "bold", margin = margin(r = 8)),
    axis.text.x  = element_text(size = 24),
    axis.text.y  = element_text(size = 24)
  )

p_sw_conf_p5


## =========================================
## Small-world (p = 10): Confounding recovery AUC — MR.RGM
## =========================================
AUC_mr_sw_p10 <- Conf_df %>%
  filter(NetworkSize == 10) %>%
  mutate(SampleSize = as.integer(SampleSize))

x_lvls   <- sort(unique(AUC_mr_sw_p10$SampleSize))
x_labels <- setNames(paste0("n = ", comma(x_lvls)), x_lvls)

p_sw_conf_p10 <- ggplot(
  AUC_mr_sw_p10 %>% mutate(SampleSize = factor(SampleSize, levels = x_lvls)),
  aes(x = SampleSize, y = AUC)
) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             linewidth = 0.4, color = "grey60") +
  geom_boxplot(
    width = 0.6,
    fill = "#0072B2",
    color = "grey15",
    outlier.shape = 16,
    outlier.size = 1.8,
    alpha = 0.88
  ) +
  geom_jitter(width = 0.08, height = 0,
              size = 1.4, alpha = 0.5, color = "grey20") +
  scale_x_discrete(labels = x_labels) +
  coord_cartesian(ylim = c(0.4, 1.0)) +
  labs(
    title = "Confounding recovery AUC — MR.RGM across sample sizes (small-world, p = 10)",
    x     = "Sample Size",
    y     = "Area Under Curve (AUC)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title.position = "panel",
    plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 8)),
    axis.title.y = element_text(size = 24, face = "bold", margin = margin(r = 8)),
    axis.text.x  = element_text(size = 24),
    axis.text.y  = element_text(size = 24)
  )

p_sw_conf_p10




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
## CONF. RECOVERY — SMALL-WORLD NETWORKS WITH HORIZONTAL PLEIOTROPY
##   - Compare MR.RGM vs MR.RGM+ (Full D) for confounding recovery
##   - Metrics: AUC (continuous ZEst), TPR/FDR/MCC (thresholded ZEst > 0.5)
##   - Run for p ∈ {5,10} and n ∈ {500,1000,10000,30000}
##   - Output:
##       (1) Table of mean & SD for AUC/TPR/FDR/MCC by (p,n,Method)
##       (2) AUC comparison boxplots across sample sizes, separately for p=5 and p=10
##
## Notes:
##   - SmallWorld_Pleio_Conf returns Conf1 (MR.RGM with correct D) and Conf2 (MR.RGM+ with Full D).
###############################################################################

## ----------------------------
## Libraries
## ----------------------------
library(MR.RGM)
library(igraph)
library(MASS)
library(pROC)
library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)  # comma()

## ----------------------------
## One replicate: small-world confounding recovery with horizontal pleiotropy
## ----------------------------
SmallWorld_Pleio_Conf = function(run, Sample_Size, Network_Size) {

  set.seed(run)

  # Sample size (n), number of response variables (p), number of instrument variables (k) and number of confounders (l_star)
  n = Sample_Size
  p =  Network_Size
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

  # Adjacency
  adj_matrix <- as_adjacency_matrix(sw_directed, sparse = FALSE)
  adj_matrix_reversed <- t(adj_matrix)

  # Generate A matrix
  A = matrix(sample(c(-0.1, 0.1), p^2, replace = TRUE), p, p)
  A = A * adj_matrix_reversed

  # Base B: 3 IVs per exposure
  B <- matrix(0, nrow = p, ncol = 3 * p)
  block_size <- 3
  for (i in 1:p) {
    start_col <- (i - 1) * block_size + 1
    end_col <- start_col + block_size - 1
    B[i, start_col:end_col] <- 1
  }

  # Add pleiotropic IVs
  pleio_cols <- list()

  if (p %% 2 == 0) {
    for (i in seq(1, p, by = 2)) {
      new_col <- rep(0, p)
      new_col[i] <- 1
      new_col[i + 1] <- 1
      pleio_cols[[length(pleio_cols) + 1]] <- new_col
    }
  } else {
    for (i in seq(1, p - 1, by = 2)) {
      new_col <- rep(0, p)
      new_col[i] <- 1
      new_col[i + 1] <- 1
      pleio_cols[[length(pleio_cols) + 1]] <- new_col
    }
    new_col <- rep(0, p)
    new_col[1] <- 1
    new_col[p] <- 1
    pleio_cols[[length(pleio_cols) + 1]] <- new_col
  }

  pleio_matrix <- do.call(cbind, pleio_cols)

  # Subset matrix: retain only one '1' per column (prioritize top-to-bottom)
  subset_matrix <- pleio_matrix
  for (j in 1:ncol(pleio_matrix)) {
    rows_with_1 <- which(pleio_matrix[, j] == 1)
    if (length(rows_with_1) > 1) {
      subset_matrix[rows_with_1[-1], j] <- 0
    }
  }

  # D for MR.RGM (subset), B for data generation (full pleio)
  D <- cbind(B, subset_matrix)
  B <- cbind(B, pleio_matrix)
  k = ncol(B)

  # Confounders
  C = matrix(sample(c(-1, 1), p * l_star, replace = TRUE), p, l_star)
  C[sample(which(C!=0), length(which(C!=0)) * 1 / 2)] = 0
  U = matrix(rnorm(n * l_star, 0, 1), nrow = n, ncol = l_star)

  # Variances
  Sigma = var * diag(p)
  Mult_Mat = solve(diag(p) - A)
  Variance = Mult_Mat %*% Sigma %*% t(Mult_Mat)
  Sigma_True = C %*% t(C) + Sigma

  set.seed(run)

  # Data generation
  X = matrix(rnorm(n * k, 0, 1), nrow = n, ncol = k)
  Y = matrix(0, nrow = n, ncol = p)
  for (i in 1:n) {
    Y[i, ] = mvrnorm(n = 1, Mult_Mat %*% (B %*% X[i, ] + C %*% U[i, ]), Variance)
  }

  n <- nrow(Y); p <- ncol(Y); k <- ncol(X)

  # Summary stats
  S_YY = t(Y) %*% Y / n
  S_YX = t(Y) %*% X / n
  S_XX = t(X) %*% X / n

  # MR.RGM (D = subset)
  Output1 = RGM(Syy = S_YY, Syx = S_YX, Sxx = S_XX,
                D = D, n = n, nIter = 50000, nBurnin = 10000, Thin = 10,
                prior = "Spike and Slab", SigmaStarModel = "SSSL")

  # MR.RGM+ (Full D)
  Output2 = RGM(Syy = S_YY, Syx = S_YX, Sxx = S_XX,
                D = matrix(1, nrow = p, ncol = k), n = n, nIter = 50000, nBurnin = 10000, Thin = 10,
                prior = "Spike and Slab", SigmaStarModel = "SSSL")

  # Z1
  Z1_Est = Output1$ZEst
  Z1_Est_UpperTriangular = Z1_Est[upper.tri(Z1_Est, diag = FALSE)]
  Sigma_True_UpperTriangular = Sigma_True[upper.tri(Sigma_True, diag = FALSE)]
  Sigma_True_Modified = abs(Sigma_True_UpperTriangular) / max(abs(Sigma_True_UpperTriangular))

  Conf1 = metrics_binary_graph(as.vector((Sigma_True_Modified > mean(Sigma_True_Modified)) * 1),
                               as.vector((Z1_Est_UpperTriangular > 0.5) * 1))
  Conf1[1] = pROC::auc(pROC::roc((Sigma_True_Modified > mean(Sigma_True_Modified)) * 1, Z1_Est_UpperTriangular))

  # Z2
  Z2_Est = Output2$ZEst
  Z2_Est_UpperTriangular = Z2_Est[upper.tri(Z2_Est, diag = FALSE)]

  Conf2 = metrics_binary_graph(as.vector((Sigma_True_Modified > mean(Sigma_True_Modified)) * 1),
                               as.vector((Z2_Est_UpperTriangular > 0.5) * 1))
  Conf2[1] = pROC::auc(pROC::roc((Sigma_True_Modified > mean(Sigma_True_Modified)) * 1, Z2_Est_UpperTriangular))

  return(list(Conf1 = Conf1, Conf2 = Conf2))
}

## ----------------------------
## Run full grid (p x n) and collect results for BOTH methods
## ----------------------------
sample_sizes  <- c(500, 1000, 10000, 30000)
network_sizes <- c(5, 10)

n_runs  <- 20
n_cores <- max(1, parallel::detectCores() - 1)

grid <- expand.grid(
  NetworkSize = network_sizes,
  SampleSize  = sample_sizes,
  stringsAsFactors = FALSE
)

results_all <- lapply(seq_len(nrow(grid)), function(idx) {

  p <- grid$NetworkSize[idx]
  n <- grid$SampleSize[idx]

  res_list <- parallel::mclapply(
    X = 1:n_runs,
    FUN = function(run) SmallWorld_Pleio_Conf(run = run, Sample_Size = n, Network_Size = p),
    mc.cores = n_cores
  )

  # Build two matrices (MR.RGM, MR.RGM+)
  met1 <- do.call(rbind, lapply(res_list, function(x) x$Conf1))
  met2 <- do.call(rbind, lapply(res_list, function(x) x$Conf2))

  met1 <- as.matrix(met1); colnames(met1) <- c("AUC", "TPR", "FDR", "MCC")
  met2 <- as.matrix(met2); colnames(met2) <- c("AUC", "TPR", "FDR", "MCC")

  df1 <- as.data.frame(met1) %>%
    mutate(Run = seq_len(n_runs),
           NetworkSize = p,
           SampleSize  = n,
           Method = "MR.RGM")

  df2 <- as.data.frame(met2) %>%
    mutate(Run = seq_len(n_runs),
           NetworkSize = p,
           SampleSize  = n,
           Method = "MR.RGM+")

  bind_rows(df1, df2)
})

Conf_df <- bind_rows(results_all)

## ----------------------------
## (A) Table: mean & SD by (p,n,Method)
## ----------------------------
Conf_summary <- Conf_df %>%
  group_by(Method, NetworkSize, SampleSize) %>%
  summarise(
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_sd   = sd(AUC, na.rm = TRUE),
    TPR_mean = mean(TPR, na.rm = TRUE),
    TPR_sd   = sd(TPR, na.rm = TRUE),
    FDR_mean = mean(FDR, na.rm = TRUE),
    FDR_sd   = sd(FDR, na.rm = TRUE),
    MCC_mean = mean(MCC, na.rm = TRUE),
    MCC_sd   = sd(MCC, na.rm = TRUE),
    .groups = "drop"
  )

print(Conf_summary)

## ----------------------------
## (B) Plot: AUC comparison (MR.RGM vs MR.RGM+) across sample sizes
##   - One plot for p = 5
##   - One plot for p = 10
## ----------------------------

## Helper: pretty facet labels
sample_labels <- c(
  "500"   = "Sample Size: 500",
  "1000"  = "Sample Size: 1000",
  "10000" = "Sample Size: 10,000",
  "30000" = "Sample Size: 30,000"
)

## =========================================
## Small-world + pleiotropy (p = 5): MR.RGM vs MR.RGM+
## =========================================
AUC_df_p5 <- Conf_df %>%
  filter(NetworkSize == 5) %>%
  transmute(
    SampleSize = as.integer(SampleSize),
    Method = Method,
    AUC = AUC
  )

AUC_df_p5$Method <- factor(AUC_df_p5$Method, levels = c("MR.RGM", "MR.RGM+"))
AUC_df_p5$SampleSize <- factor(AUC_df_p5$SampleSize, levels = sample_sizes)

p_sw_pleio_p5 <- ggplot(AUC_df_p5, aes(x = Method, y = AUC, fill = Method)) +
  geom_boxplot(width = 0.6, outlier.shape = 16, outlier.size = 2, alpha = 0.85) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 1.5, color = "grey20") +
  facet_wrap(~ SampleSize, scales = "free_x",
             labeller = labeller(SampleSize = sample_labels)) +
  guides(fill = "none") +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x  = element_text(size = 24, angle = 25, hjust = 1, face = "bold"),
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_text(size = 26, face = "bold"),
    axis.title.y = element_text(size = 26, face = "bold"),
    strip.text   = element_text(size = 22, face = "bold"),
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Confounding recovery AUC: MR.RGM vs MR.RGM+ across sample sizes\n(small-world network with horizontal pleiotropy, network size = 5)",
    x = "Methods",
    y = "Area Under Curve (AUC)"
  )

p_sw_pleio_p5


## =========================================
## Small-world + pleiotropy (p = 10): MR.RGM vs MR.RGM+
## =========================================
AUC_df_p10 <- Conf_df %>%
  filter(NetworkSize == 10) %>%
  transmute(
    SampleSize = as.integer(SampleSize),
    Method = Method,
    AUC = AUC
  )

AUC_df_p10$Method <- factor(AUC_df_p10$Method, levels = c("MR.RGM", "MR.RGM+"))
AUC_df_p10$SampleSize <- factor(AUC_df_p10$SampleSize, levels = sample_sizes)

p_sw_pleio_p10 <- ggplot(AUC_df_p10, aes(x = Method, y = AUC, fill = Method)) +
  geom_boxplot(width = 0.6, outlier.shape = 16, outlier.size = 2, alpha = 0.85) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 1.5, color = "grey20") +
  facet_wrap(~ SampleSize, scales = "free_x",
             labeller = labeller(SampleSize = sample_labels)) +
  guides(fill = "none") +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x  = element_text(size = 24, angle = 25, hjust = 1, face = "bold"),
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_text(size = 26, face = "bold"),
    axis.title.y = element_text(size = 26, face = "bold"),
    strip.text   = element_text(size = 22, face = "bold"),
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Confounding recovery AUC: MR.RGM vs MR.RGM+ across sample sizes\n(small-world network with horizontal pleiotropy, network size = 10)",
    x = "Methods",
    y = "Area Under Curve (AUC)"
  )

p_sw_pleio_p10

