###############################################################################
## GTEx (MR.RGM+) — Motif figures + motif posterior probabilities
##
## This script is used AFTER fitting the GTEx model (see: GTEx_ModelFitting.R).
##
## Required objects (created in GTEx_ModelFitting.R):
##   - Output_GTEx       : output from MR.RGM+ fit
##   - gene_names_final  : character vector of final gene labels used in the fit
##
## What this script does:
##   (A) Creates three clean PDF figures (used in the paper) for specific motifs:
##       1) Feedback loop (3-cycle)
##       2) Feedforward loop
##       3) Cascade
##
##   (B) Computes the posterior probability of each motif using:
##       NetworkMotif(Gamma0, GammaPst)
##       where:
##         - Gamma0   : binary adjacency (p x p) encoding the motif edges
##         - GammaPst : posterior samples from the GTEx fit (Output_GTEx$GammaPst)
##
## Notes:
##   - In MR.RGM output, Gamma matrices use the convention:
##       rows = "to", columns = "from"
##     i.e., Gamma[to, from] corresponds to edge from -> to.
###############################################################################

## ---------------------------- ##
## Packages
## ---------------------------- ##
suppressPackageStartupMessages({
  library(igraph)
})

## ---------------------------- ##
## Helper: build directed motif graph (for plotting)
## ---------------------------- ##
make_motif_graph <- function(nodes, edges_from, edges_to) {
  stopifnot(length(edges_from) == length(edges_to))
  graph_from_data_frame(
    data.frame(from = edges_from, to = edges_to, stringsAsFactors = FALSE),
    directed = TRUE,
    vertices = data.frame(name = nodes, stringsAsFactors = FALSE)
  )
}

## ---------------------------- ##
## Helper: save a motif PDF (no clipping)
## ---------------------------- ##
save_motif_pdf <- function(
    g, layout_mat, file,
    width = 4.2, height = 4.2,
    pad = 0.35,
    vertex_size = 30,
    label_cex = 0.95,
    edge_width = 2.0,
    arrow_size = 0.6,
    edge_curved = 0.0
) {
  ## Expand plot limits so nothing is clipped
  xlim <- range(layout_mat[, 1]) + c(-pad, pad)
  ylim <- range(layout_mat[, 2]) + c(-pad, pad)

  pdf(file, width = width, height = height, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)

  par(mar = c(1.6, 1.6, 1.6, 1.6), xpd = NA)

  plot(
    g,
    layout = layout_mat,
    rescale = FALSE,
    xlim = xlim,
    ylim = ylim,
    asp = 0,  ## avoid clipping due to aspect enforcement

    vertex.size = vertex_size,
    vertex.color = "white",
    vertex.frame.color = "black",
    vertex.label = V(g)$name,
    vertex.label.cex = label_cex,
    vertex.label.color = "black",
    vertex.label.family = "serif",

    edge.width = edge_width,
    edge.color = "black",
    edge.arrow.size = arrow_size,
    edge.curved = edge_curved
  )

  box()
}

## ---------------------------- ##
## Canonical layouts (named by node)
## ---------------------------- ##

## Triangle (equilateral-ish): good for feedback loops
layout_triangle <- function(nodes_order) {
  coords <- matrix(
    c(0,   1,
      -0.87, -0.5,
      0.87, -0.5),
    ncol = 2, byrow = TRUE
  )
  rownames(coords) <- nodes_order
  coords
}

## Feedforward triangle: A(left), B(top), C(right)
layout_feedforward <- function(nodes_order) {
  coords <- matrix(
    c(-1, 0,
      0, 1,
      1, 0),
    ncol = 2, byrow = TRUE
  )
  rownames(coords) <- nodes_order
  coords
}

## Cascade (L-shape): reduces overlap in small figure boxes
layout_cascade_L <- function(nodes_order) {
  coords <- matrix(
    c(0,   1,
      0,  -0.2,
      1,  -0.2),
    ncol = 2, byrow = TRUE
  )
  rownames(coords) <- nodes_order
  coords
}

## Convert named coords -> igraph vertex order
layout_for_graph <- function(g, coords_named) {
  coords_named[V(g)$name, , drop = FALSE]
}

## ---------------------------- ##
## Helper: build Gamma0 motif adjacency (p x p, binary)
## using the MR.RGM convention: rows = to, cols = from
## ---------------------------- ##
make_gamma0 <- function(gene_names, edges_from, edges_to) {
  p <- length(gene_names)
  Gamma0 <- matrix(0, nrow = p, ncol = p)
  rownames(Gamma0) <- colnames(Gamma0) <- gene_names

  stopifnot(length(edges_from) == length(edges_to))
  for (e in seq_along(edges_from)) {
    from <- edges_from[e]
    to   <- edges_to[e]
    Gamma0[to, from] <- 1
  }
  Gamma0
}

## ---------------------------- ##
## Output filenames (paper figures)
## ---------------------------- ##
out_gtex_feedback    <- "motif_gtex_feedback.pdf"
out_gtex_feedforward <- "motif_gtex_feedforward.pdf"
out_gtex_cascade     <- "motif_gtex_cascade.pdf"

## ---------------------------- ##
## Sanity checks: required objects exist
## ---------------------------- ##
stopifnot(exists("Output_GTEx"))
stopifnot(exists("gene_names_final"))
stopifnot(!is.null(Output_GTEx$GammaPst))

## ============================================================
## 1) Feedback loop (3-cycle): PDK1 -> VHL -> PKCA -> PDK1
## ============================================================
nodes_gtex_fb <- c("PDK1", "VHL", "PKCA")

## (A) Plot motif (for paper)
g_gtex_fb <- make_motif_graph(
  nodes      = nodes_gtex_fb,
  edges_from = c("PDK1", "VHL", "PKCA"),
  edges_to   = c("VHL",  "PKCA", "PDK1")
)
coords_gtex_fb <- layout_triangle(nodes_gtex_fb)

save_motif_pdf(
  g = g_gtex_fb,
  layout_mat = layout_for_graph(g_gtex_fb, coords_gtex_fb),
  file = out_gtex_feedback,
  vertex_size = 30,
  pad = 0.38
)

## (B) Compute motif posterior probability
Gamma0_fb <- make_gamma0(
  gene_names = gene_names_final,
  edges_from = c("PDK1", "VHL", "PKCA"),
  edges_to   = c("VHL",  "PKCA", "PDK1")
)
prob_fb <- NetworkMotif(Gamma0_fb, Output_GTEx$GammaPst)
prob_fb

## ============================================================
## 2) Feedforward loop: PDK1 -> VHL, PDK1 -> S6K, VHL -> S6K
## (A=PDK1 left, B=VHL top, C=S6K right)
## ============================================================
nodes_gtex_ff <- c("PDK1", "VHL", "S6K")

## (A) Plot motif (for paper)
g_gtex_ff <- make_motif_graph(
  nodes      = nodes_gtex_ff,
  edges_from = c("PDK1", "PDK1", "VHL"),
  edges_to   = c("VHL",  "S6K",  "S6K")
)
coords_gtex_ff <- layout_feedforward(nodes_gtex_ff)

save_motif_pdf(
  g = g_gtex_ff,
  layout_mat = layout_for_graph(g_gtex_ff, coords_gtex_ff),
  file = out_gtex_feedforward,
  vertex_size = 30,
  pad = 0.38
)

## (B) Compute motif posterior probability
Gamma0_ff <- make_gamma0(
  gene_names = gene_names_final,
  edges_from = c("PDK1", "PDK1", "VHL"),
  edges_to   = c("VHL",  "S6K",  "S6K")
)
prob_ff <- NetworkMotif(Gamma0_ff, Output_GTEx$GammaPst)
prob_ff

## ============================================================
## 3) Cascade: VHL -> PKCA -> MTOR
## (use L-shape layout to reduce overlap)
## ============================================================
nodes_gtex_cas <- c("VHL", "PKCA", "MTOR")

## (A) Plot motif (for paper)
g_gtex_cas <- make_motif_graph(
  nodes      = nodes_gtex_cas,
  edges_from = c("VHL", "PKCA"),
  edges_to   = c("PKCA", "MTOR")
)
coords_gtex_cas <- layout_cascade_L(nodes_gtex_cas)

save_motif_pdf(
  g = g_gtex_cas,
  layout_mat = layout_for_graph(g_gtex_cas, coords_gtex_cas),
  file = out_gtex_cascade,
  vertex_size = 26,
  label_cex = 0.9,
  pad = 0.45
)

## (B) Compute motif posterior probability
Gamma0_cas <- make_gamma0(
  gene_names = gene_names_final,
  edges_from = c("VHL", "PKCA"),
  edges_to   = c("PKCA", "MTOR")
)
prob_cas <- NetworkMotif(Gamma0_cas, Output_GTEx$GammaPst)
prob_cas
