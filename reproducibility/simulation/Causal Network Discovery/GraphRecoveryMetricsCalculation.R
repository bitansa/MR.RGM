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
# Note:
#   AUC is computed using binary predictions; it should be interpreted
#   accordingly (ties may occur).
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
