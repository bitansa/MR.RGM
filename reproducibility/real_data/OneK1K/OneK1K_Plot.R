################################################################################
## OneK1K Network: extract edge PIPs + plot 9 cluster subnetworks
##
## Purpose:
##   (1) Query posterior inclusion probabilities (PIPs) for:
##        - causal edges via GammaEst  (directed)
##        - confounding links via ZEst (undirected / symmetric)
##   (2) Plot the 9 cluster subnetworks used in the manuscript figure
##
## Prerequisites (run first):
##   - OneK1K_PreProcessing.R   : builds/loads processed inputs
##   - OneK1K_ModelFitting.R    : fits the model and creates:
##        Output_OneK1K           (list with AEst, GammaEst, ZEst, SigmaEst)
##        gene_names_pathway_use  (final gene list in the figure order)
##
## Notes on interpretation:
##   - GammaEst: T2[to, from] is PIP for causal edge  from -> to
##   - ZEst:     T3[g1, g2]   is PIP for confounding link between g1 and g2
################################################################################

## ---- Packages ----
library(dplyr)
library(ggplot2)
library(ggraph)
library(tidygraph)
library(scales)  # rescale()
library(grid)    # unit()

## ---- 1) Pull estimated matrices from model output ----
T1 <- Output_OneK1K$AEst
T2 <- Output_OneK1K$GammaEst
T3 <- Output_OneK1K$ZEst
T4 <- Output_OneK1K$SigmaEst

## Ensure dimnames are consistent
gene_names <- gene_names_pathway_use

rownames(T1) <- rownames(T2) <- rownames(T3) <- rownames(T4) <- gene_names
colnames(T1) <- colnames(T2) <- colnames(T3) <- colnames(T4) <- gene_names

## ---- 2) Helper functions: query PIPs for edges ----
# Causal edge PIP: from -> to  (GammaEst is indexed as [to, from])
get_pip_causal <- function(T2, from, to) {
  stopifnot(from %in% colnames(T2), to %in% rownames(T2))
  as.numeric(T2[to, from])
}

# Confounding link PIP: undirected pair (ZEst is symmetric)
get_pip_confound <- function(T3, g1, g2) {
  stopifnot(g1 %in% colnames(T3), g2 %in% rownames(T3))
  as.numeric(T3[g1, g2])
}

## Example queries (edit genes as needed)
# Causal: NFKB -> CD40
get_pip_causal(T2, from = "NFKB", to = "CD40")

# Confounding: SHIP -- FGR2B
get_pip_confound(T3, g1 = "SHIP", g2 = "FGR2B")


################################################################################
## 3) Plotting function: draw directed causal + undirected confounding edges
##
## Defaults:
##   - causal threshold  thr_dir  = 0.50
##   - conf threshold    thr_conf = 0.48
################################################################################
plot_onexk_network <- function(
    T2, T3, genes,
    thr_dir  = 0.50,
    thr_conf = 0.48,
    # colors
    col_causal = "#0072B2",
    col_conf   = "#D55E00",
    # scaling
    w_range_causal = c(0.55, 1.55),
    w_range_conf   = c(0.35, 1.05),
    a_range_causal = c(0.65, 1.00),
    a_range_conf   = c(0.45, 0.90),
    # geometry
    conf_curvature = 0.45,
    arrow_len_mm   = 3.2,
    # labels
    label_size      = 4.8,
    label_color     = "grey10",
    label_fill      = "white",
    label_pad_lines = 0.18,
    label_border    = 0.25,
    # layout
    layout_alg = "kk",      # "kk" or "stress" (stress needs {graphlayouts})
    genes_keep = NULL,
    seed = 1,
    # output (optional)
    out_file = NULL,
    width = 8, height = 6
) {
  set.seed(seed)

  # Ensure dimnames exist
  if (is.null(rownames(T2)) || is.null(colnames(T2))) rownames(T2) <- colnames(T2) <- genes
  if (is.null(rownames(T3)) || is.null(colnames(T3))) rownames(T3) <- colnames(T3) <- genes

  # Optional: subset to a cluster
  if (!is.null(genes_keep)) {
    genes_keep <- intersect(genes_keep, genes)
    T2 <- T2[genes_keep, genes_keep, drop = FALSE]
    T3 <- T3[genes_keep, genes_keep, drop = FALSE]
    genes <- genes_keep
  }

  ## ---- Causal edges (directed): include i != j ----
  dir_idx <- which(T2 >= thr_dir & row(T2) != col(T2), arr.ind = TRUE)

  edges_causal <- if (nrow(dir_idx) > 0) {
    tibble(
      from = colnames(T2)[dir_idx[, 2]],
      to   = rownames(T2)[dir_idx[, 1]],
      prob = as.numeric(T2[dir_idx]),
      type = "causal"
    ) %>%
      mutate(
        edge_w = rescale(prob, to = w_range_causal),
        edge_a = rescale(prob, to = a_range_causal)
      )
  } else {
    tibble(from=character(), to=character(), prob=numeric(),
           type=character(), edge_w=numeric(), edge_a=numeric())
  }

  ## ---- Confounding edges (undirected): draw each pair ONCE (upper triangle) ----
  conf_idx <- which(T3 >= thr_conf & upper.tri(T3), arr.ind = TRUE)

  edges_conf <- if (nrow(conf_idx) > 0) {
    tibble(
      from = colnames(T3)[conf_idx[, 2]],
      to   = rownames(T3)[conf_idx[, 1]],
      prob = as.numeric(T3[conf_idx]),
      type = "conf"
    ) %>%
      mutate(
        edge_w = rescale(prob, to = w_range_conf),
        edge_a = rescale(prob, to = a_range_conf)
      )
  } else {
    tibble(from=character(), to=character(), prob=numeric(),
           type=character(), edge_w=numeric(), edge_a=numeric())
  }

  nodes <- tibble(name = genes)
  all_edges <- bind_rows(edges_causal, edges_conf)

  g_tbl <- tbl_graph(nodes = nodes, edges = all_edges, directed = TRUE)

  chosen_layout <- layout_alg
  if (layout_alg == "stress" && !requireNamespace("graphlayouts", quietly = TRUE)) {
    message("`graphlayouts` not installed; falling back to `kk` layout.")
    chosen_layout <- "kk"
  }
  lay <- create_layout(g_tbl, layout = chosen_layout)

  cap_mm <- 4.5  # helps keep arrowheads away from node label boxes

  p <- ggraph(lay) +
    geom_node_label(
      aes(label = name),
      size = label_size,
      label.size = label_border,
      label.padding = unit(label_pad_lines, "lines"),
      fill = label_fill,
      colour = label_color
    ) +
    # Confounding: curved red arcs (once per pair)
    geom_edge_arc(
      aes(edge_width = edge_w, edge_alpha = edge_a, filter = type == "conf"),
      curvature = conf_curvature,
      colour = col_conf,
      show.legend = FALSE,
      start_cap = circle(cap_mm, "mm"),
      end_cap   = circle(cap_mm, "mm"),
      lineend = "round",
      linejoin = "round"
    ) +
    # Causal: directed blue edges
    geom_edge_link(
      aes(edge_width = edge_w, edge_alpha = edge_a, filter = type == "causal"),
      colour = col_causal,
      show.legend = FALSE,
      arrow = arrow(length = unit(arrow_len_mm, "mm"), type = "closed"),
      start_cap = circle(cap_mm, "mm"),
      end_cap   = circle(cap_mm, "mm"),
      lineend = "round",
      linejoin = "mitre"
    ) +
    scale_edge_width(range = c(0.25, 1.6)) +
    scale_edge_alpha(range = c(0.8, 1.0)) +
    theme_void(base_size = 12) +
    theme(
      plot.margin  = margin(3, 3, 3, 3),
      plot.title   = element_text(face = "bold", size = 16),
      plot.subtitle= element_text(size = 12),
      panel.border = element_rect(colour = "grey85", fill = NA, linewidth = 0.4)
    ) +
    labs(
      title   = paste0("OneK1K B-cell network  (Causal ≥ ", thr_dir,
                       ", Confounding ≥ ", thr_conf, ")"),
      caption = "Blue: causal (double-headed if bidirectional). Red: confounding (curved)."
    )

  if (!is.null(out_file)) {
    if (!dir.exists(dirname(out_file))) dir.create(dirname(out_file), recursive = TRUE)
    ggsave(out_file, plot = p, width = width, height = height, units = "in",
           device = grDevices::cairo_pdf)
  }

  p
}


################################################################################
## 4) Plot the 9 clusters used in the figure
##    Default thresholds: causal ≥ 0.50, confounding ≥ 0.48
##    You can edit thr_dir_vec / thr_conf_vec if you want panel-specific cutoffs.
################################################################################

clusters <- list(
  c("CD19","PI3K","PIP3","AKT","PLCY2"),
  c("SYK","PYK2","CBL","DOK1","PLCY2"),
  c("MEKK","MEK1/2","ERK1/2","JNK","JUN"),
  c("CD40","NFKB","IKB","PKC","CAMK"),
  c("RIAM","RAP","CREB","JNK"),
  c("VAV","RAC","EZRIN","HS1","PYK2"),
  c("PLCY2","VAV","PI3K","PIP3"),
  c("SHIP","FGR2B","PI3K","PIP3","AKT"),
  c("CREB","MEF2C","CAMK","JNK")
)

titles <- c(
  "Proximal PI3K axis",
  "SYK & negative regulation",
  "MAPK core",
  "NF-κB & Ca2+ cross-talk",
  "Integrin → TF bridge",
  "Cytoskeleton/adhesion",
  "PI3K–PLCY2–VAV bridge",
  "Inhibitory receptor module",
  "TFs & kinases (compact)"
)

# Default: same thresholds for all panels (edit if needed)
thr_dir_vec  <- rep(0.50, length(clusters))
thr_conf_vec <- rep(0.48, length(clusters))

for (i in seq_along(clusters)) {
  thrD <- thr_dir_vec[i]
  thrC <- thr_conf_vec[i]

  p <- plot_onexk_network(
    T2, T3, genes = gene_names,
    thr_dir  = thrD,
    thr_conf = thrC,
    genes_keep = clusters[[i]],
    seed = 1
  ) +
    labs(
      title    = paste("OneK1K B cells —", titles[i]),
      subtitle = sprintf("Causal PIP ≥ %.2f   |   Confounding PIP ≥ %.2f", thrD, thrC)
    )

  print(p)
}
