## ============================================================
## Runtime benchmark: 6 methods vs number of responses p
## Fixed n = 30000, p in {2,5,10,20}
## ============================================================

## ----------------------------
## Libraries
## ----------------------------
library(igraph)
library(MASS)
library(microbenchmark)
library(reshape2)
library(dplyr)
library(ggplot2)
library(scales)
library(MR.RGM)
library(MendelianRandomization)
library(mrbayes)
library(OneSampleMR)
library(rjags)  # needed by mr_ivw_rjags

## ----------------------------
## Wrapper functions (RETURN OUTPUT)
## ----------------------------
Fn_MR_RGM <- function(Syy, Syx, Sxx, D, n) {
  out <- RGM(
    Syy = Syy, Syx = Syx, Sxx = Sxx,
    D = D, n = n,
    nIter = 50000, nBurnin = 10000, Thin = 2,
    prior = "Spike and Slab", SigmaStarModel = "SSSL"
  )
  return(out)
}

Fn_MR_RGM_Plus <- function(Syy, Syx, Sxx, D_full, n) {
  out <- RGM(
    Syy = Syy, Syx = Syx, Sxx = Sxx,
    D = D_full, n = n,
    nIter = 50000, nBurnin = 10000, Thin = 2,
    prior = "Spike and Slab", SigmaStarModel = "SSSL"
  )
  return(out)
}

Fn_MR_RGM_NoConf <- function(Syy, Syx, Sxx, D, n) {
  out <- RGM(
    Syy = Syy, Syx = Syx, Sxx = Sxx,
    D = D, n = n,
    nIter = 50000, nBurnin = 10000, Thin = 2,
    prior = "Spike and Slab", SigmaStarModel = "diagonal"
  )
  return(out)
}

## ---------- MendelianRandomization wrapper ----------
## Uses mr_allmethods() for each (i,j). .
# Function to compute MRinput object
run_mr_allmethods <- function(X, Y, D, exposure_index, outcome_index) {
  ivs <- which(D[exposure_index, ] == 1)
  if (length(ivs) < 3) {
    stop("Need at least 3 instruments for stable MR estimates")
  }

  G <- X[, ivs, drop = FALSE]
  snp_names <- paste0("rs", ivs)
  exposure <- Y[, exposure_index]
  outcome <- Y[, outcome_index]

  # Estimate SNP-exposure and SNP-outcome associations
  get_betas <- function(y, G) {
    t(sapply(1:ncol(G), function(j) {
      fit <- summary(lm(y ~ G[, j]))
      coef <- fit$coefficients[2, ]  # slope
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

  # Run all MR methods
  res <- mr_allmethods(mr_input_obj, method = "main", iterations = 50000)
  return(res)
}

Fn_MendelianRandomization <- function(X, Y, D) {
  p <- ncol(Y)
  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) next
      tmp <- run_mr_allmethods(X, Y, D = D, exposure_index = i, outcome_index = j)
    }
  }
  return(TRUE)
}

## ---------- mrbayes wrapper ----------
get_snp_effects <- function(expr_vector, snp_matrix) {
  betas <- se <- numeric(ncol(snp_matrix))
  for (m in seq_len(ncol(snp_matrix))) {
    model <- summary(lm(expr_vector ~ snp_matrix[, m]))
    betas[m] <- coef(model)[2, 1]
    se[m]    <- coef(model)[2, 2]
  }
  list(beta = betas, se = se)
}

Fn_mrbayes <- function(X, Y, D) {
  p <- ncol(Y)
  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) next
      snp_indices <- which(D[i, ] == 1)
      if (length(snp_indices) < 1) next

      X_sub <- X[, snp_indices, drop = FALSE]
      exp_effects <- get_snp_effects(Y[, i], X_sub)
      out_effects <- get_snp_effects(Y[, j], X_sub)

      df <- data.frame(
        bx   = exp_effects$beta,
        bxse = exp_effects$se,
        by   = out_effects$beta,
        byse = out_effects$se
      )

      stan_data <- mr_format(
        rsid  = paste0("SNP", snp_indices),
        xbeta = df$bx, ybeta = df$by,
        xse   = df$bxse, yse = df$byse
      )

      bayes_model <- mr_ivw_rjags(
        object   = stan_data,
        prior    = "default",
        n.chains = 2,
        n.burn   = 50000,
        n.iter   = 10000
      )
    }
  }
  return(TRUE)
}

## ---------- OneSampleMR wrapper ----------
Fn_OneSampleMR <- function(X, Y, D) {
  p <- ncol(Y)

  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) next

      iv_idx <- which(D[i, ] == 1)
      if (length(iv_idx) < 1) next

      G <- X[, iv_idx, drop = FALSE]
      instrument_names <- paste0("z", iv_idx)
      colnames(G) <- instrument_names

      exposure <- Y[, i]
      outcome  <- Y[, j]

      df <- data.frame(y = outcome, x = exposure, G)
      instruments_formula <- reformulate(termlabels = instrument_names)
      formula <- y ~ x

      fit <- tryCatch(
        tsri(
          formula     = formula,
          instruments = instruments_formula,
          data        = df,
          link        = "identity"
        ),
        error = function(e) NULL
      )
    }
  }
  return(TRUE)
}

## ----------------------------
## Runtime table vs p
## ----------------------------
p_vals <- c(2, 5, 10, 20)
n <- 30000

results_df <- data.frame(
  Method = c("MR.RGM", "MR.RGM_Plus", "MR.RGM_NoConf",
             "MendelianRandomization", "OneSampleMR", "MrBayes"),
  stringsAsFactors = FALSE
)

for (p in p_vals) {

  k <- 3 * p
  l_star <- 5

  ## --- scale-free graph ---
  g <- sample_pa(n = p, directed = FALSE)
  g <- as.directed(g, mode = "mutual")

  adj_matrix <- as_adjacency_matrix(g, sparse = FALSE)
  adj_matrix_reversed <- t(adj_matrix)

  ## --- A ---
  A <- matrix(sample(c(-0.1, 0.1), p^2, replace = TRUE), p, p)
  A <- A * adj_matrix_reversed

  ## --- B and D ---
  B <- matrix(0, nrow = p, ncol = 3 * p)
  block_size <- 3
  for (i in 1:p) {
    start_col <- (i - 1) * block_size + 1
    end_col <- start_col + block_size - 1
    B[i, start_col:end_col] <- 1
  }
  D <- B

  ## --- confounders ---
  var <- 9
  C <- matrix(sample(c(-1, 1), p * l_star, replace = TRUE), p, l_star)
  C[sample(which(C != 0), length(which(C != 0)) / 2)] <- 0
  U <- matrix(rnorm(n * l_star, 0, 1), nrow = n, ncol = l_star)

  ## --- covariance pieces ---
  Sigma <- var * diag(p)
  Mult_Mat <- solve(diag(p) - A)
  Variance <- Mult_Mat %*% Sigma %*% t(Mult_Mat)

  ## --- simulate X, Y ---
  X <- matrix(rnorm(n * k, 0, 1), nrow = n, ncol = k)
  Y <- matrix(0, nrow = n, ncol = p)
  for (ii in 1:n) {
    Y[ii, ] <- mvrnorm(n = 1, Mult_Mat %*% (B %*% X[ii, ] + C %*% U[ii, ]), Variance)
  }

  ## --- summary stats ---
  Syy <- t(Y) %*% Y / n
  Syx <- t(Y) %*% X / n
  Sxx <- t(X) %*% X / n

  ## full-D for MR.RGM+
  D_full <- matrix(1, nrow = nrow(Syy), ncol = nrow(Sxx))

  ## --- microbenchmark ---
  mb <- microbenchmark(
    MR.RGM = Fn_MR_RGM(Syy = Syy, Syx = Syx, Sxx = Sxx, D = D, n = n),
    MR.RGM_Plus = Fn_MR_RGM_Plus(Syy = Syy, Syx = Syx, Sxx = Sxx, D_full = D_full, n = n),
    MR.RGM_NoConf = Fn_MR_RGM_NoConf(Syy = Syy, Syx = Syx, Sxx = Sxx, D = D, n = n),
    MendelianRandomization = Fn_MendelianRandomization(X = X, Y = Y, D = D),
    OneSampleMR = Fn_OneSampleMR(X = X, Y = Y, D = D),
    MrBayes = Fn_mrbayes(X = X, Y = Y, D = D),
    times = 20
  )

  medians <- tapply(mb$time, mb$expr, median) / 1e6
  results_df[[paste0("p=", p)]] <- as.numeric(medians[
    c("MR.RGM", "MR.RGM_Plus", "MR.RGM_NoConf",
      "MendelianRandomization", "OneSampleMR", "MrBayes")
  ])
}

cat(paste0("# Benchmark Runtime Table (n = ", n, ")\n"))
print(results_df, row.names = FALSE)

## ----------------------------
## Plot
## ----------------------------
results_df$Method[which(results_df$Method == "MR.RGM_Plus")] <- "MR.RGM+"

results_long <- melt(results_df, id.vars = "Method",
                     variable.name = "PSize", value.name = "Runtime_ms")
results_long$PSize <- as.numeric(gsub("p=", "", results_long$PSize))
results_long <- results_long %>% mutate(Runtime_sec = Runtime_ms / 1000)

results_long$Method[results_long$Method == "MrBayes"] <- "mrbayes"

method_order <- c("MR.RGM", "MR.RGM+", "MR.RGM_NoConf",
                  "MendelianRandomization", "OneSampleMR", "mrbayes")
results_long <- results_long %>% mutate(Method = factor(Method, levels = method_order))

okabe_ito <- c("#0072B2", "#E69F00", "#009E73",
               "#56B4E9", "#000000", "#CC79A7")

end_pts <- results_long %>%
  group_by(Method) %>%
  filter(PSize == max(PSize, na.rm = TRUE)) %>%
  ungroup()

lab_fun <- scales::label_number(scale_cut = scales::cut_short_scale())

method_label_size <- 5.3
title_size        <- 24
subtitle_size     <- 18
axis_title_size   <- 24
axis_text_size    <- 24

p_plot <- ggplot(results_long, aes(x = PSize, y = Runtime_sec, color = Method)) +
  geom_line(aes(group = Method), linewidth = 1.2) +
  geom_point(size = 2.6, stroke = 0.5) +
  geom_text(
    data = end_pts,
    aes(label = Method),
    hjust = -0.05, vjust = 0.5,
    size = method_label_size, fontface = "bold",
    show.legend = FALSE
  ) +
  scale_color_manual(values = okabe_ito, guide = "none") +
  scale_x_continuous(
    breaks = sort(unique(results_long$PSize)),
    expand = expansion(mult = c(0.02, 0.22))
  ) +
  scale_y_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  labs(
    title    = "Median runtime versus number of responses (p)",
    subtitle = paste0("Fixed sample size n = ", format(n, big.mark=","), " · medians across replicates"),
    x = "Number of Responses (p)",
    y = "Median Runtime (seconds)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title.position = "plot",
    plot.title   = element_text(face = "bold", size = title_size),
    plot.subtitle= element_text(margin = margin(b = 6), size = subtitle_size),
    axis.title   = element_text(face = "bold", size = axis_title_size),
    axis.text    = element_text(size = axis_text_size),
    axis.ticks.length = unit(3, "pt"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    plot.margin = margin(8, 42, 8, 8)
  ) +
  coord_cartesian(clip = "off")

print(p_plot)
