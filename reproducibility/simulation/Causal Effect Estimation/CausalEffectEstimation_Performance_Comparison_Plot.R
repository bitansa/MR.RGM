###############################################################################
## PLOTTING CODE — SCALE-FREE NETWORKS and SMALL-WORLD NETWORKS
## (NO HORIZONTAL PLEIOTROPY)  [Causal effect accuracy metrics]
##
## This block is used for:
##   • Scale-free networks
##   • Small-world networks WITHOUT horizontal pleiotropy
##
## Methods compared:
##   MR.RGM
##   MR.RGM (No Conf)
##   MR-SimpleMedian
##   MR-WeightedMedian
##   MR-IVW
##   OneSampleMR
##   MrBayes
##
## NOTE:
##   There is NO MR.RGM+ variant in this scenario because pleiotropy is absent.
##
## Expected inputs per sample size (loaded inside the loop):
##   T_MR.RGM, T_MR.RGM_NoConf,
##   T_SimpleMedian, T_WeightedMedian, T_IVW,
##   T_OneSampleMR, T_MrBayes
##
## IMPORTANT (matches your causal-effect scripts):
##   Each T_* is a matrix with columns:
##     col 1 = Maximum Absolute Deviation
##     col 2 = Mean Absolute Deviation      <-- we plot this (MAD)
##     col 3 = Mean Squared Error
###############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

## ----------------------------
## User controls
## ----------------------------
net_size     <- 10
sample_sizes <- c(500, 1000, 10000, 30000)

method_levels <- c(
  "MR.RGM", "MR.RGM_NoConf",
  "SimpleMedian", "WeightedMedian", "IVW",
  "OneSampleMR", "MrBayes"
)

method_labels <- c(
  "MR.RGM"         = "MR.RGM",
  "MR.RGM_NoConf"  = "MR.RGM (No Conf)",
  "SimpleMedian"   = "MR-SimpleMedian",
  "WeightedMedian" = "MR-WeightedMedian",
  "IVW"            = "MR-IVW",
  "OneSampleMR"    = "OneSampleMR",
  "MrBayes"        = "MrBayes"
)

sample_labels <- c(
  "500"   = "Sample Size: 500",
  "1000"  = "Sample Size: 1000",
  "10000" = "Sample Size: 10000",
  "30000" = "Sample Size: 30000"
)

## ----------------------------
## Helper: append one method row
## ----------------------------
append_mad <- function(df, sample_size, network_size, method_name, mad_vec, max_len) {
  if (length(mad_vec) < max_len) {
    mad_vec <- c(mad_vec, rep(NA_real_, max_len - length(mad_vec)))
  }
  mad_list <- as.list(mad_vec)
  names(mad_list) <- paste0("MADrep", seq_along(mad_vec))

  new_row <- data.frame(
    SampleSize  = sample_size,
    NetworkSize = network_size,
    Method      = method_name,
    mad_list,
    stringsAsFactors = FALSE
  )
  rbind(df, new_row)
}

## ----------------------------
## Build MAD_results
## ----------------------------
MAD_results <- data.frame()

for (sample_size in sample_sizes) {

  ## Load per-sample-size outputs here (EDIT PATHS):
  ## load(file.path("results", "CausalEffect_NoPleio",
  ##                paste0("net", net_size, "_n", sample_size, ".RData")))

  required_objs <- c("T_SimpleMedian","T_WeightedMedian","T_IVW",
                     "T_OneSampleMR","T_MrBayes","T_MR.RGM","T_MR.RGM_NoConf")
  missing <- required_objs[!vapply(required_objs, exists, logical(1))]
  if (length(missing) > 0) stop("Missing objects in workspace: ", paste(missing, collapse = ", "))

  ## Column 2 = Mean Absolute Deviation (MAD)
  max_len <- max(
    nrow(T_MR.RGM),
    nrow(T_MR.RGM_NoConf),
    nrow(T_SimpleMedian),
    nrow(T_WeightedMedian),
    nrow(T_IVW),
    nrow(T_OneSampleMR),
    nrow(T_MrBayes)
  )

  MAD_results <- append_mad(MAD_results, sample_size, net_size, "MR.RGM",        T_MR.RGM[, 2],        max_len)
  MAD_results <- append_mad(MAD_results, sample_size, net_size, "MR.RGM_NoConf", T_MR.RGM_NoConf[, 2], max_len)

  MAD_results <- append_mad(MAD_results, sample_size, net_size, "SimpleMedian",   T_SimpleMedian[, 2],   max_len)
  MAD_results <- append_mad(MAD_results, sample_size, net_size, "WeightedMedian", T_WeightedMedian[, 2], max_len)
  MAD_results <- append_mad(MAD_results, sample_size, net_size, "IVW",            T_IVW[, 2],            max_len)
  MAD_results <- append_mad(MAD_results, sample_size, net_size, "OneSampleMR",    T_OneSampleMR[, 2],    max_len)
  MAD_results <- append_mad(MAD_results, sample_size, net_size, "MrBayes",        T_MrBayes[, 2],        max_len)
}

## ----------------------------
## Plot
## ----------------------------
MAD_long <- MAD_results %>%
  pivot_longer(cols = starts_with("MADrep"), names_to = "Replicate", values_to = "MAD") %>%
  filter(!is.na(MAD)) %>%
  mutate(
    SampleSize  = factor(SampleSize),
    NetworkSize = factor(NetworkSize),
    Method      = factor(Method, levels = method_levels),
    MethodLabel = recode(as.character(Method), !!!method_labels)
  )

MAD_long$MethodLabel <- factor(MAD_long$MethodLabel, levels = unname(method_labels[method_levels]))

p <- ggplot(MAD_long, aes(x = MethodLabel, y = MAD, fill = MethodLabel)) +
  geom_boxplot() +
  facet_wrap(~ SampleSize, scales = "free_x", labeller = labeller(SampleSize = sample_labels)) +
  guides(fill = "none") +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x  = element_text(size = 24, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_text(size = 26, face = "bold"),
    axis.title.y = element_text(size = 26, face = "bold"),
    strip.text   = element_text(size = 22, face = "bold"),
    plot.title   = element_text(size = 30, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = paste0("Mean Absolute Deviation (MAD) comparison (no pleiotropy, network size = ", net_size, ")"),
    y = "Mean Absolute Deviation (MAD)",
    x = "Methods"
  )

print(p)



###############################################################################
## PLOTTING CODE — SMALL-WORLD NETWORK WITH HORIZONTAL PLEIOTROPY
## (PLEIOTROPY PRESENT)  [Causal effect accuracy metrics]
##
## This block is used ONLY for:
##   • Small-world network simulation WITH horizontal pleiotropy
##
## Additional method included:
##   MR.RGM+  (Full-D)
##
## Methods compared:
##   MR.RGM
##   MR.RGM+ (Full-D)
##   MR.RGM (No Conf)
##   MR-SimpleMedian
##   MR-WeightedMedian
##   MR-IVW
##   OneSampleMR
##   MrBayes
##
## Expected inputs per sample size (loaded inside the loop):
##   T_MR.RGM, T_MR.RGM_FullD, T_MR.RGM_NoConf,
##   T_SimpleMedian, T_WeightedMedian, T_IVW,
##   T_OneSampleMR, T_MrBayes
##
## IMPORTANT (matches your causal-effect scripts):
##   Each T_* is a matrix with columns:
##     col 1 = Maximum Absolute Deviation
##     col 2 = Mean Absolute Deviation      <-- we plot this (MAD)
##     col 3 = Mean Squared Error
###############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

net_size     <- 10
sample_sizes <- c(500, 1000, 10000, 30000)

method_levels <- c(
  "MR.RGM", "MR.RGM+",
  "MR.RGM_NoConf",
  "SimpleMedian", "WeightedMedian", "IVW",
  "OneSampleMR", "MrBayes"
)

method_labels <- c(
  "MR.RGM"         = "MR.RGM",
  "MR.RGM+"        = "MR.RGM+",
  "MR.RGM_NoConf"  = "MR.RGM (No Conf)",
  "SimpleMedian"   = "MR-SimpleMedian",
  "WeightedMedian" = "MR-WeightedMedian",
  "IVW"            = "MR-IVW",
  "OneSampleMR"    = "OneSampleMR",
  "MrBayes"        = "MrBayes"
)

sample_labels <- c(
  "500"   = "Sample Size: 500",
  "1000"  = "Sample Size: 1000",
  "10000" = "Sample Size: 10000",
  "30000" = "Sample Size: 30000"
)

append_mad <- function(df, sample_size, network_size, method_name, mad_vec, max_len) {
  if (length(mad_vec) < max_len) {
    mad_vec <- c(mad_vec, rep(NA_real_, max_len - length(mad_vec)))
  }
  mad_list <- as.list(mad_vec)
  names(mad_list) <- paste0("MADrep", seq_along(mad_vec))

  new_row <- data.frame(
    SampleSize  = sample_size,
    NetworkSize = network_size,
    Method      = method_name,
    mad_list,
    stringsAsFactors = FALSE
  )
  rbind(df, new_row)
}

MAD_results <- data.frame()

for (sample_size in sample_sizes) {

  ## Load per-sample-size pleiotropy outputs here (EDIT PATHS):
  ## load(file.path("results", "CausalEffect_SmallWorldPleio",
  ##                paste0("net", net_size, "_n", sample_size, ".RData")))

  required_objs <- c("T_SimpleMedian","T_WeightedMedian","T_IVW",
                     "T_OneSampleMR","T_MrBayes",
                     "T_MR.RGM","T_MR.RGM_FullD","T_MR.RGM_NoConf")
  missing <- required_objs[!vapply(required_objs, exists, logical(1))]
  if (length(missing) > 0) stop("Missing objects in workspace: ", paste(missing, collapse = ", "))

  max_len <- max(
    nrow(T_MR.RGM),
    nrow(T_MR.RGM_FullD),
    nrow(T_MR.RGM_NoConf),
    nrow(T_SimpleMedian),
    nrow(T_WeightedMedian),
    nrow(T_IVW),
    nrow(T_OneSampleMR),
    nrow(T_MrBayes)
  )

  MAD_results <- append_mad(MAD_results, sample_size, net_size, "MR.RGM",        T_MR.RGM[, 2],       max_len)
  MAD_results <- append_mad(MAD_results, sample_size, net_size, "MR.RGM+",       T_MR.RGM_FullD[, 2], max_len)
  MAD_results <- append_mad(MAD_results, sample_size, net_size, "MR.RGM_NoConf", T_MR.RGM_NoConf[, 2],max_len)

  MAD_results <- append_mad(MAD_results, sample_size, net_size, "SimpleMedian",   T_SimpleMedian[, 2],   max_len)
  MAD_results <- append_mad(MAD_results, sample_size, net_size, "WeightedMedian", T_WeightedMedian[, 2], max_len)
  MAD_results <- append_mad(MAD_results, sample_size, net_size, "IVW",            T_IVW[, 2],            max_len)
  MAD_results <- append_mad(MAD_results, sample_size, net_size, "OneSampleMR",    T_OneSampleMR[, 2],    max_len)
  MAD_results <- append_mad(MAD_results, sample_size, net_size, "MrBayes",        T_MrBayes[, 2],        max_len)
}

MAD_long <- MAD_results %>%
  pivot_longer(cols = starts_with("MADrep"), names_to = "Replicate", values_to = "MAD") %>%
  filter(!is.na(MAD)) %>%
  mutate(
    SampleSize  = factor(SampleSize),
    NetworkSize = factor(NetworkSize),
    Method      = factor(Method, levels = method_levels),
    MethodLabel = recode(as.character(Method), !!!method_labels)
  )

MAD_long$MethodLabel <- factor(MAD_long$MethodLabel, levels = unname(method_labels[method_levels]))

p <- ggplot(MAD_long, aes(x = MethodLabel, y = MAD, fill = MethodLabel)) +
  geom_boxplot() +
  facet_wrap(~ SampleSize, scales = "free_x", labeller = labeller(SampleSize = sample_labels)) +
  guides(fill = "none") +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x  = element_text(size = 22, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 18),
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    strip.text   = element_text(size = 18, face = "bold"),
    plot.title   = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = paste0("Mean Absolute Deviation (MAD) comparison (small-world + pleiotropy, network size = ", net_size, ")"),
    y = "Mean Absolute Deviation (MAD)",
    x = "Methods"
  )

print(p)
