###############################################################################
## PLOTTING CODE — SCALE-FREE NETWORKS and SMALL-WORLD NETWORKS
## (NO HORIZONTAL PLEIOTROPY)
##
## This block is used for:
##   • Scale-free networks
##   • Small-world networks WITHOUT horizontal pleiotropy
##
## In these settings, the following methods are compared:
##   MR.RGM
##   MR.RGM (No Conf)
##   MR-SimpleMedian
##   MR-WeightedMedian
##   MR-IVW
##   OneSampleMR
##   MrBayes
##
## NOTE:
##   There is NO MR.RGM+ variant in this scenario because
##   horizontal pleiotropy is not present.
##
## Only run THIS block for non-pleiotropic simulations.
###############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

## ----------------------------
## User controls
## ----------------------------
net_size <- 10                 # <- change to 5 or 10
sample_sizes <- c(500, 1000, 10000, 30000)

## Methods to include (order matters)
method_levels <- c(
  "MR.RGM", "MR.RGM_NoConf",
  "SimpleMedian", "WeightedMedian", "IVW",
  "OneSampleMR", "MrBayes"
)

## Nicer labels (optional)
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
append_auc <- function(df, sample_size, network_size, method_name, auc_vec, max_len) {
  stopifnot(is.numeric(auc_vec))
  if (length(auc_vec) < max_len) {
    auc_vec <- c(auc_vec, rep(NA_real_, max_len - length(auc_vec)))
  }
  auc_list <- as.list(auc_vec)
  names(auc_list) <- paste0("AUC", seq_along(auc_vec))

  new_row <- data.frame(
    SampleSize  = sample_size,
    NetworkSize = network_size,
    Method      = method_name,
    auc_list,
    stringsAsFactors = FALSE
  )

  rbind(df, new_row)
}

## ----------------------------
## Build AUC_results
## ----------------------------
AUC_results <- data.frame()

## IMPORTANT:
## This loop assumes that the method-specific AUC objects (T_*) correspond to the
## current sample_size. In practice you should:
##   (a) load per-sample-size results inside this loop (recommended), OR
##   (b) store results in a list keyed by sample_size and index it here.
##
## If you do NOT do (a) or (b), you will append the SAME AUC vectors for every sample_size.

for (sample_size in sample_sizes) {

  ## TODO (recommended): load the correct T_* objects here for this sample_size.
  ## Example pattern:
  ## T_SimpleMedian <- readRDS(paste0("results/net", net_size, "/n", sample_size, "_SimpleMedian.rds"))
  ## ... etc.

  ## Safety: check required objects exist in the workspace
  required_objs <- c("T_SimpleMedian","T_WeightedMedian","T_IVW","T_OneSampleMR","T_MrBayes","T_MR.RGM","T_MR.RGM_NoConf")
  missing <- required_objs[!vapply(required_objs, exists, logical(1))]
  if (length(missing) > 0) stop("Missing objects in workspace: ", paste(missing, collapse = ", "))

  ## Determine padding length
  max_len <- max(
    length(T_SimpleMedian[, 1]),
    length(T_WeightedMedian[, 1]),
    length(T_IVW[, 1]),
    length(T_OneSampleMR[, 1]),
    length(T_MrBayes[, 1]),
    length(T_MR.RGM[, 1]),
    length(T_MR.RGM_NoConf[, 1])
  )

  ## Append in desired order
  AUC_results <- append_auc(AUC_results, sample_size, net_size, "MR.RGM",        T_MR.RGM[, 1],        max_len)
  AUC_results <- append_auc(AUC_results, sample_size, net_size, "MR.RGM_NoConf", T_MR.RGM_NoConf[, 1], max_len)

  AUC_results <- append_auc(AUC_results, sample_size, net_size, "SimpleMedian",   T_SimpleMedian[, 1],   max_len)
  AUC_results <- append_auc(AUC_results, sample_size, net_size, "WeightedMedian", T_WeightedMedian[, 1], max_len)
  AUC_results <- append_auc(AUC_results, sample_size, net_size, "IVW",            T_IVW[, 1],            max_len)
  AUC_results <- append_auc(AUC_results, sample_size, net_size, "OneSampleMR",    T_OneSampleMR[, 1],    max_len)
  AUC_results <- append_auc(AUC_results, sample_size, net_size, "MrBayes",        T_MrBayes[, 1],        max_len)
}

## ----------------------------
## Plot
## ----------------------------
AUC_long <- AUC_results %>%
  pivot_longer(cols = starts_with("AUC"), names_to = "Replicate", values_to = "AUC") %>%
  filter(!is.na(AUC)) %>%  # <- important: drop padding NAs
  mutate(
    SampleSize  = factor(SampleSize),
    NetworkSize = factor(NetworkSize),
    Method      = factor(Method, levels = method_levels),
    MethodLabel = recode(as.character(Method), !!!method_labels)
  )

## sanity: ensure no missing labels
if (any(is.na(AUC_long$MethodLabel))) {
  bad <- unique(as.character(AUC_long$Method[is.na(AUC_long$MethodLabel)]))
  stop("Some methods were not mapped in method_labels: ", paste(bad, collapse = ", "))
}
AUC_long$MethodLabel <- factor(AUC_long$MethodLabel, levels = unname(method_labels[method_levels]))

p <- ggplot(AUC_long, aes(x = MethodLabel, y = AUC, fill = MethodLabel)) +
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
    title = paste0("AUC comparison (network size = ", net_size, ")"),
    y = "Area Under Curve (AUC)",
    x = "Methods"
  )

print(p)







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
## PLOTTING CODE — SMALL-WORLD NETWORK WITH HORIZONTAL PLEIOTROPY
##
## This block is used ONLY for the small-world network
## simulation WITH horizontal pleiotropy.
##
## In this scenario, an additional method is included:
##   MR.RGM+  (Full-D)
##
## Methods compared in this setting:
##   MR.RGM
##   MR.RGM+
##   MR.RGM (No Conf)
##   MR-SimpleMedian
##   MR-WeightedMedian
##   MR-IVW
##   OneSampleMR
##   MrBayes
##
## NOTE:
##   MR.RGM+ is important in the presence of
##   horizontal pleiotropy
##
## Only run THIS block for pleiotropic simulations.
###############################################################################



library(dplyr)
library(tidyr)
library(ggplot2)

## ----------------------------
## User controls
## ----------------------------
net_size     <- 10                  # <- set to 5 or 10
sample_sizes <- c(500, 1000, 10000, 30000)

## Methods to keep (order matters)
method_levels <- c(
  "MR.RGM", "MR.RGM+",
  "MR.RGM_NoConf",
  "SimpleMedian", "WeightedMedian", "IVW",
  "OneSampleMR", "MrBayes"
)

## Pretty labels for plotting
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

## ----------------------------
## Helper: append one method row
## ----------------------------
append_auc <- function(df, sample_size, network_size, method_name, auc_vec, max_len) {
  if (length(auc_vec) < max_len) {
    auc_vec <- c(auc_vec, rep(NA_real_, max_len - length(auc_vec)))
  }
  auc_list <- as.list(auc_vec)
  names(auc_list) <- paste0("AUC", seq_along(auc_vec))

  new_row <- data.frame(
    SampleSize  = sample_size,
    NetworkSize = network_size,
    Method      = method_name,
    auc_list,
    stringsAsFactors = FALSE
  )

  rbind(df, new_row)
}

## ======================================================================
## OPTION A:
##   Build AUC_results from the vectors produced by your simulation script.
##
##   IMPORTANT:
##   For this loop to work, the objects T_* must correspond to each sample_size.
##   Common approaches:
##     (1) you load each sample_size result into T_* inside the loop, or
##     (2) you run simulation for one sample_size at a time and then append.
## ======================================================================

AUC_results <- data.frame()

for (sample_size in sample_sizes) {

  ## ---- You must ensure T_* corresponds to this sample_size ----
  ## Example pattern (pseudo):
  ## load(paste0("results_SWpleio_n", sample_size, "_net", net_size, ".RData"))
  ##
  ## After loading, you should have:
  ## T_SimpleMedian, T_WeightedMedian, T_IVW, T_OneSampleMR, T_MrBayes,
  ## T_MR.RGM, T_MR.RGM_FullD, T_MR.RGM_NoConf

  max_len <- max(
    length(T_SimpleMedian[, 1]),
    length(T_WeightedMedian[, 1]),
    length(T_IVW[, 1]),
    length(T_OneSampleMR[, 1]),
    length(T_MrBayes[, 1]),
    length(T_MR.RGM[, 1]),
    length(T_MR.RGM_FullD[, 1]),     # <- MR.RGM+
    length(T_MR.RGM_NoConf[, 1])
  )

  ## Add methods in the desired order
  AUC_results <- append_auc(AUC_results, sample_size, net_size, "MR.RGM",
                            T_MR.RGM[, 1], max_len)

  AUC_results <- append_auc(AUC_results, sample_size, net_size, "MR.RGM+",
                            T_MR.RGM_FullD[, 1], max_len)

  AUC_results <- append_auc(AUC_results, sample_size, net_size, "MR.RGM_NoConf",
                            T_MR.RGM_NoConf[, 1], max_len)

  AUC_results <- append_auc(AUC_results, sample_size, net_size, "SimpleMedian",
                            T_SimpleMedian[, 1], max_len)

  AUC_results <- append_auc(AUC_results, sample_size, net_size, "WeightedMedian",
                            T_WeightedMedian[, 1], max_len)

  AUC_results <- append_auc(AUC_results, sample_size, net_size, "IVW",
                            T_IVW[, 1], max_len)

  AUC_results <- append_auc(AUC_results, sample_size, net_size, "OneSampleMR",
                            T_OneSampleMR[, 1], max_len)

  AUC_results <- append_auc(AUC_results, sample_size, net_size, "MrBayes",
                            T_MrBayes[, 1], max_len)
}

## ======================================================================
## OPTION B (also fine):
##   Load a pre-built AUC_results and just plot.
##   Keep GitHub-friendly relative paths (no /Users/...).
## ======================================================================
# AUC_results <- readRDS(file.path("reproducibility", "simulation", "results",
#                                 paste0("AUC_results-SmallWorldPleio_p_", net_size, ".rds")))

## ----------------------------
## Plot
## ----------------------------
AUC_long <- AUC_results %>%
  pivot_longer(cols = starts_with("AUC"),
               names_to = "Replicate", values_to = "AUC") %>%
  mutate(
    SampleSize  = factor(SampleSize),
    NetworkSize = factor(NetworkSize),
    Method      = factor(Method, levels = method_levels),
    MethodLabel = factor(recode(as.character(Method), !!!method_labels),
                         levels = unname(method_labels[method_levels]))
  )

p <- ggplot(AUC_long, aes(x = MethodLabel, y = AUC, fill = MethodLabel)) +
  geom_boxplot() +
  facet_wrap(~ SampleSize, scales = "free_x",
             labeller = labeller(SampleSize = sample_labels)) +
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
    title = paste0("AUC comparison (small-world + horizontal pleiotropy, network size = ", net_size, ")"),
    y = "Area Under Curve (AUC)",
    x = "Methods"
  )

print(p)

