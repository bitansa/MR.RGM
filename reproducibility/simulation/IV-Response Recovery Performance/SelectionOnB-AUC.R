###############################################################################
## IV RECOVERY (Instrument–trait selection) — Small-world + horizontal pleiotropy
##   - Method: MR.RGM+ (Full D = 1 matrix)
##   - Metric: AUC for recovering true B (instrument–trait links)
##            * truth  : (B != 0)
##            * score  : Output$Psi_Est
##   - Two manuscript plots:
##       (1) p = 5  (AUC vs sample size)
##       (2) p = 10 (AUC vs sample size)
##   - Parallel across replicates using mclapply()
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
## Settings (edit as needed)
## ----------------------------
sample_sizes  <- c(500, 1000, 10000, 30000)
network_sizes <- c(5, 10)

n_runs  <- 20
n_cores <- max(1, parallel::detectCores() - 1)

## -------------------------------------------------------------------------
## One replicate:
##   Returns a single AUC for IV recovery under MR.RGM+ (Full D)
## -------------------------------------------------------------------------
SmallWorld_Pleio_IV_AUC <- function(run, Sample_Size, Network_Size) {

  set.seed(run)

  n <- Sample_Size
  p <- Network_Size
  k <- 3 * p
  l_star <- 5

  ## variance
  var <- 9

  nei_sw  <- 2        # small-world parameter
  rewire_p <- 0.1     # small-world parameter

  ## --- small-world graph (directed) ---
  sw_undirected <- igraph::sample_smallworld(dim = 1, size = p, nei = nei_sw, p = rewire_p)
  edge_list <- igraph::as_edgelist(sw_undirected)

  directed_edges <- t(apply(edge_list, 1, function(e) {
    if (runif(1) > 0.5) e else rev(e)
  }))

  sw_directed <- igraph::graph_from_edgelist(directed_edges, directed = TRUE)

  adj_matrix <- igraph::as_adjacency_matrix(sw_directed, sparse = FALSE)
  adj_matrix_reversed <- t(adj_matrix)

  ## --- A matrix (structural graph) ---
  A <- matrix(sample(c(-0.1, 0.1), p^2, replace = TRUE), p, p)
  A <- A * adj_matrix_reversed

  ## --- Base B: 3 IVs per exposure ---
  B <- matrix(0, nrow = p, ncol = 3 * p)
  block_size <- 3
  for (i in 1:p) {
    start_col <- (i - 1) * block_size + 1
    end_col   <- start_col + block_size - 1
    B[i, start_col:end_col] <- 1
  }

  ## --- Add pleiotropic IVs ---
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

  ## NOTE: For MR.RGM+ we use full D anyway, so we do not need subset_matrix here.
  ##       But we keep D construction conceptually consistent:
  ##       - Data generation uses full B (includes pleio_matrix)
  B <- cbind(B, pleio_matrix)
  k <- ncol(B)

  ## --- Confounders ---
  C <- matrix(sample(c(-1, 1), p * l_star, replace = TRUE), p, l_star)
  C[sample(which(C != 0), length(which(C != 0)) * 1 / 2)] <- 0
  U <- matrix(rnorm(n * l_star, 0, 1), nrow = n, ncol = l_star)

  ## --- Noise covariance ---
  Sigma <- var * diag(p)

  Mult_Mat <- solve(diag(p) - A)
  Variance <- Mult_Mat %*% Sigma %*% t(Mult_Mat)

  ## --- Data generation ---
  set.seed(run)
  X <- matrix(rnorm(n * k, 0, 1), nrow = n, ncol = k)
  Y <- matrix(0, nrow = n, ncol = p)

  for (i in 1:n) {
    Y[i, ] <- MASS::mvrnorm(
      n = 1,
      mu = Mult_Mat %*% (B %*% X[i, ] + C %*% U[i, ]),
      Sigma = Variance
    )
  }

  ## --- Summary stats ---
  S_YY <- t(Y) %*% Y / n
  S_YX <- t(Y) %*% X / n
  S_XX <- t(X) %*% X / n

  ## --- MR.RGM+ : Full D ---
  Output <- RGM(
    Syy = S_YY, Syx = S_YX, Sxx = S_XX,
    D = matrix(1, nrow = p, ncol = k),
    n = n, nIter = 50000, nBurnin = 10000, Thin = 10,
    prior = "Spike and Slab", SigmaStarModel = "SSSL"
  )

  ## --- IV recovery AUC ---
  ## truth: (B != 0), score: Psi_Est
  auc_val <- as.numeric(
    pROC::auc(pROC::roc(as.vector((B != 0) * 1), as.vector(Output$PsiEst), quiet = TRUE))
  )

  return(auc_val)
}

## -------------------------------------------------------------------------
## Run grid (p x n) and store AUCs (wide format for your existing workflow)
## -------------------------------------------------------------------------
AUC_results <- data.frame()

append_auc <- function(df, sample_size, network_size, auc_vec) {
  auc_list <- as.list(auc_vec)
  names(auc_list) <- paste0("AUC", seq_along(auc_vec))

  new_row <- data.frame(
    SampleSize  = sample_size,
    NetworkSize = network_size,
    auc_list,
    stringsAsFactors = FALSE
  )

  rbind(df, new_row)
}

for (p in network_sizes) {
  for (n in sample_sizes) {

    auc_vec <- unlist(parallel::mclapply(
      X = 1:n_runs,
      FUN = function(run) SmallWorld_Pleio_IV_AUC(run = run, Sample_Size = n, Network_Size = p),
      mc.cores = n_cores
    ))

    AUC_results <- append_auc(AUC_results, sample_size = n, network_size = p, auc_vec = auc_vec)
    message("Done: p = ", p, ", n = ", n)
  }
}


## -------------------------------------------------------------------------
## Tidy + manuscript plots (ONLY p = 5 and p = 10)
## -------------------------------------------------------------------------
tidy_auc_plus <- function(df_wide) {
  stopifnot(all(c("SampleSize", "NetworkSize") %in% names(df_wide)))
  auc_cols <- grep("^AUC", names(df_wide), value = TRUE)

  df_wide %>%
    pivot_longer(cols = all_of(auc_cols), names_to = "RunID", values_to = "AUC") %>%
    mutate(
      SampleSize  = as.integer(SampleSize),
      NetworkSize = as.integer(NetworkSize),
      Method      = "MR.RGM_Plus"
    )
}

df_long <- tidy_auc_plus(AUC_results)

plot_snp_auc_plus <- function(df_long, p_value, title_text) {

  d <- df_long %>% filter(NetworkSize == p_value)

  x_levels <- sort(unique(d$SampleSize))
  x_labels <- setNames(paste0("n = ", scales::comma(x_levels)), x_levels)

  ggplot(
    d %>% mutate(SampleSize = factor(SampleSize, levels = x_levels)),
    aes(x = SampleSize, y = AUC)
  ) +
    geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.4, color = "grey60") +
    geom_boxplot(
      width = 0.6, fill = "#0072B2", color = "grey15",
      outlier.shape = 16, outlier.size = 1.8, alpha = 0.88
    ) +
    geom_jitter(width = 0.08, height = 0, size = 1.4, alpha = 0.5, color = "grey20") +
    scale_x_discrete(labels = x_labels) +
    coord_cartesian(ylim = c(0.85, 1.0)) +
    labs(
      title = title_text,
      x     = "Sample size",
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
      axis.text.x  = element_text(size = 18),
      axis.text.y  = element_text(size = 18)
    )
}

## --- Manuscript plots ---
p5 <- plot_snp_auc_plus(
  df_long, p_value = 5,
  title_text = "Instrument–trait selection AUC (MR.RGM+, SW + horizontal pleiotropy, p = 5)"
)

p10 <- plot_snp_auc_plus(
  df_long, p_value = 10,
  title_text = "Instrument–trait selection AUC (MR.RGM+, SW + horizontal pleiotropy, p = 10)"
)

p5
p10
