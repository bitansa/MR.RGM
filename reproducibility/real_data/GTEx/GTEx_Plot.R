###############################################################################
## GTEx (MR.RGM+) — Extract posterior inclusion probabilities (PIPs)
## and plot the inferred network (causal + confounding).
##
## Assumes you already ran RGM() and have:
##   Output_GTEx  : RGM output object
##   gene_names_final : length-p vector of gene labels (used for plotting/lookup)
##
## What this script does:
##   1) Extracts posterior estimates from Output_GTEx:
##        - AEst      : (optional) causal effect estimate (kept for completeness)
##        - GammaEst  : causal edge inclusion probabilities (directed)
##        - ZEst      : confounding inclusion probabilities (undirected, symmetric)
##        - SigmaEst  : (optional) covariance estimate (kept for completeness)
##   2) Shows how to query PIPs for:
##        - a directed causal edge Gene_From -> Gene_To  (from GammaEst)
##        - an undirected confounding link Gene1 -- Gene2 (from ZEst)
##   3) Plots a network:
##        - Blue arrows: causal edges with PIP >= thr_dir
##        - Red curved arcs: confounding links with PIP >= thr_conf
###############################################################################

## ---------------------------- ##
## 1) Extract posterior matrices
## ---------------------------- ##
T1 <- Output_GTEx$AEst
T2 <- Output_GTEx$GammaEst   # causal PIPs (directed)
T3 <- Output_GTEx$ZEst       # confounding PIPs (undirected; treat as symmetric)
T4 <- Output_GTEx$SigmaEst

## Attach gene labels for easier indexing / plotting
rownames(T1) <- colnames(T1) <- gene_names_final
rownames(T2) <- colnames(T2) <- gene_names_final
rownames(T3) <- colnames(T3) <- gene_names_final
rownames(T4) <- colnames(T4) <- gene_names_final

## ---------------------------- ##
## 2) Examples: query PIPs
## ---------------------------- ##

## ---- (A) Causal edge inclusion probability: Gene_From -> Gene_To ----
##   - rows correspond to "to"
##   - cols correspond to "from"
Gene_From <- "MTOR"
Gene_To   <- "S6K"

pip_causal <- T2[Gene_To, Gene_From]
pip_causal

## Change Gene_From / Gene_To to query any other directed pair.

## ---- (B) Confounding inclusion probability: Gene1 -- Gene2 ----
## Confounding is undirected; using either [Gene1, Gene2] or [Gene2, Gene1]
Gene1 <- "MTOR"
Gene2 <- "TSC2"

pip_conf <- T3[Gene1, Gene2]
pip_conf

## Change Gene1 / Gene2 to query any other undirected pair.

## ---------------------------- ##
## 3) Plot network (causal + confounding)
## ---------------------------- ##
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggraph)
library(igraph)
library(tidygraph)
library(scales)
library(grid)

plot_gtex_network <- function(
    T2, T3, genes,
    thr_dir  = 0.85,
    thr_conf = 0.50,

    ## colors
    col_causal = "#0072B2",  # blue
    col_conf   = "#D55E00",  # red

    ## edge aesthetics scaling
    w_range_causal = c(0.55, 1.55),
    w_range_conf   = c(0.35, 1.05),
    a_range_causal = c(0.65, 1.00),
    a_range_conf   = c(0.45, 0.90),

    ## curvature for confounding arcs
    conf_curvature = 0.25,

    ## arrows & labels
    arrow_len_mm    = 3.2,
    label_size      = 4.8,
    label_color     = "grey10",
    label_fill      = "white",
    label_pad_lines = 0.18,   # padding inside label box (in "lines")
    label_border    = 0.25,   # outline width of label box

    ## text sizing
    title_size   = 24,
    caption_size = 18,

    ## layout & subsetting
    layout_alg = "kk",
    genes_keep = NULL,
    seed = 1,

    ## output
    out_file = NULL,
    width = 8,
    height = 6
) {
  set.seed(seed)

  ## Ensure row/colnames exist
  if (is.null(rownames(T2)) || is.null(colnames(T2))) rownames(T2) <- colnames(T2) <- genes
  if (is.null(rownames(T3)) || is.null(colnames(T3))) rownames(T3) <- colnames(T3) <- genes

  ## Optional: subset to a smaller set of genes for a cleaner figure
  if (!is.null(genes_keep)) {
    genes_keep <- intersect(genes_keep, genes)
    T2 <- T2[genes_keep, genes_keep, drop = FALSE]
    T3 <- T3[genes_keep, genes_keep, drop = FALSE]
    genes <- genes_keep
  }

  ## ---------------------------- ##
  ## Build edge list: causal (directed)
  ## ---------------------------- ##
  dir_idx <- which(T2 >= thr_dir & row(T2) != col(T2), arr.ind = TRUE)

  if (nrow(dir_idx) > 0) {
    edges_causal_all <- tibble(
      from = colnames(T2)[dir_idx[, 2]],
      to   = rownames(T2)[dir_idx[, 1]],
      prob = as.numeric(T2[dir_idx])
    ) %>%
      ## mark bidirectional pairs (A->B and B->A both present)
      mutate(key_undir = paste(pmin(from, to), pmax(from, to), sep = "|"))

    pair_counts <- count(edges_causal_all, key_undir, name = "n")

    edges_causal_all <- edges_causal_all %>%
      left_join(pair_counts, by = "key_undir") %>%
      mutate(type = if_else(n == 2, "causal_bi", "causal_uni")) %>%
      select(-n, -key_undir) %>%
      mutate(
        edge_w = rescale(prob, to = w_range_causal),
        edge_a = rescale(prob, to = a_range_causal)
      )
  } else {
    edges_causal_all <- tibble(
      from = character(), to = character(), prob = numeric(),
      type = character(), edge_w = numeric(), edge_a = numeric()
    )
  }

  ## ---------------------------- ##
  ## Build edge list: confounding (undirected)
  ## We only take upper triangle to avoid duplicates, then mirror it.
  ## ---------------------------- ##
  conf_idx <- which(T3 >= thr_conf & upper.tri(T3), arr.ind = TRUE)

  if (nrow(conf_idx) > 0) {
    conf_ud <- tibble(
      from = colnames(T3)[conf_idx[, 2]],
      to   = rownames(T3)[conf_idx[, 1]],
      prob = as.numeric(T3[conf_idx])
    )

    edges_conf <- bind_rows(
      conf_ud,
      transmute(conf_ud, from = to, to = from, prob = prob)
    ) %>%
      mutate(
        type   = "conf",
        edge_w = rescale(prob, to = w_range_conf),
        edge_a = rescale(prob, to = a_range_conf)
      )
  } else {
    edges_conf <- tibble(
      from = character(), to = character(), prob = numeric(),
      type = character(), edge_w = numeric(), edge_a = numeric()
    )
  }

  ## Nodes + combined edges
  nodes <- tibble(name = genes)
  all_edges <- bind_rows(edges_causal_all, edges_conf)
  g_tbl <- tbl_graph(nodes = nodes, edges = all_edges, directed = TRUE)

  ## Layout choice
  chosen_layout <- layout_alg
  if (layout_alg == "stress" && !requireNamespace("graphlayouts", quietly = TRUE)) {
    message("`graphlayouts` not installed; falling back to `kk` layout.")
    chosen_layout <- "kk"
  }
  lay <- create_layout(g_tbl, layout = chosen_layout)

  ## Stand-off so arrows don't poke into label boxes (global conservative setting)
  cap_mm <- 4.0

  ## ---------------------------- ##
  ## Plot
  ## ---------------------------- ##
  p <- ggraph(lay) +
    ## Node labels (auto-sized boxes)
    geom_node_label(
      aes(label = name),
      size = label_size,
      label.size = label_border,
      label.padding = unit(label_pad_lines, "lines"),
      fill = label_fill,
      colour = label_color
    ) +
    ## Confounding edges: curved red arcs
    geom_edge_arc(
      aes(edge_width = edge_w, edge_alpha = edge_a, filter = type == "conf"),
      curvature = conf_curvature,
      colour = col_conf, show.legend = FALSE,
      start_cap = circle(cap_mm, "mm"), end_cap = circle(cap_mm, "mm"),
      lineend = "round", linejoin = "round"
    ) +
    ## Causal edges: blue arrows (bidirectional edges appear as two arrows)
    geom_edge_link(
      aes(edge_width = edge_w, edge_alpha = edge_a,
          filter = type %in% c("causal_uni", "causal_bi")),
      colour = col_causal, show.legend = FALSE,
      arrow = arrow(length = unit(arrow_len_mm, "mm"), type = "closed"),
      start_cap = circle(cap_mm, "mm"), end_cap = circle(cap_mm, "mm"),
      lineend = "round", linejoin = "mitre"
    ) +
    scale_edge_width(range = c(0.25, 1.6)) +
    scale_edge_alpha(range = c(0.8, 1.0)) +
    theme_void(base_size = 12) +
    theme(
      plot.margin  = margin(3, 3, 3, 3),
      plot.title   = element_text(face = "bold", size = title_size),
      plot.caption = element_text(size = caption_size),
      panel.border = element_rect(colour = "grey85", fill = NA, linewidth = 0.4)
    ) +
    labs(
      title = paste0(
        "GTEx skeletal muscle network (causal inc. prob. ≥ ", thr_dir,
        ", confounding inc. prob. ≥ ", thr_conf, ")"
      ),
      caption = "Blue: causal (double-headed if bidirectional); red: confounding (curved)."
    )

  ## Save (optional)
  if (!is.null(out_file)) {
    if (!dir.exists(dirname(out_file))) dir.create(dirname(out_file), recursive = TRUE)
    ggsave(out_file, plot = p, width = width, height = height, units = "in",
           device = grDevices::cairo_pdf)
  }

  return(p)
}

## ---------------------------- ##
## Example: full network plot
## ---------------------------- ##
p_gtex_full <- plot_gtex_network(
  T2, T3, genes = gene_names_final,
  thr_dir = 0.85,
  thr_conf = 0.50,
  conf_curvature = 0.45,
  label_size = 5,
  title_size = 12,
  caption_size = 12
)

p_gtex_full
