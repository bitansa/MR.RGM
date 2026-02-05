################################################################################
## OneK1K_NetworkMotif.R
##
## What this script does
##   1) Plots 3 canonical network motifs (as clean PDFs)
##        - Feedback loop (3-cycle)
##        - Feedforward loop
##        - Cascade (3-chain)
##   2) Computes the posterior probability of each motif using
##      MR.RGM::NetworkMotif() and posterior samples Output_OneK1K$GammaPst
##
## Prerequisites (run first)
##   - OneK1K_PreProcessing.R
##   - OneK1K_ModelFitting.R
##
## Required objects created by the scripts above
##   - Output_OneK1K          : model fit object containing $GammaPst
##   - gene_names_pathway_use : character vector of gene names (final order)
################################################################################

suppressPackageStartupMessages({
  library(igraph)
  library(MR.RGM)
})

## ---- Helper: build a directed igraph from "from -> to" edges ----
make_motif_graph <- function(nodes, edges_from, edges_to) {
  stopifnot(length(edges_from) == length(edges_to))
  graph_from_data_frame(
    d = data.frame(from = edges_from, to = edges_to, stringsAsFactors = FALSE),
    directed = TRUE,
    vertices = data.frame(name = nodes, stringsAsFactors = FALSE)
  )
}

## ---- Helper: map a named layout to igraph's vertex order ----
layout_for_graph <- function(g, coords_named) {
  coords_named[V(g)$name, , drop = FALSE]
}

## ---- Helper: save motif as PDF with padding (prevents clipping) ----
save_motif_pdf <- function(g, layout_mat, file,
                           width = 4.2, height = 4.2,
                           pad = 0.35,
                           vertex_size = 30,
                           label_cex = 0.95,
                           edge_width = 2.0,
                           arrow_size = 0.6,
                           edge_curved = 0.0) {

  xlim <- range(layout_mat[, 1]) + c(-pad, pad)
  ylim <- range(layout_mat[, 2]) + c(-pad, pad)

  pdf(file, width = width, height = height, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)

  par(mar = c(1.6, 1.6, 1.6, 1.6), xpd = NA)

  plot(
    g,
    layout  = layout_mat,
    rescale = FALSE,
    xlim    = xlim,
    ylim    = ylim,
    asp     = 0,                  # avoids clipping from aspect enforcement
    vertex.size        = vertex_size,
    vertex.color       = "white",
    vertex.frame.color = "black",
    vertex.label       = V(g)$name,
    vertex.label.cex   = label_cex,
    vertex.label.color = "black",
    vertex.label.family= "serif",
    edge.width         = edge_width,
    edge.color         = "black",
    edge.arrow.size    = arrow_size,
    edge.curved        = edge_curved
  )

  box()
}

## ---- Canonical layouts ----
# Triangle (equilateral-ish): good for feedback loops
layout_triangle <- function(nodes_order) {
  coords <- matrix(
    c( 0,    1,
       -0.87, -0.5,
       0.87, -0.5),
    ncol = 2, byrow = TRUE
  )
  rownames(coords) <- nodes_order
  coords
}

# Feedforward triangle: A(left), B(top), C(right)
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

# Cascade in an L-shape: reduces overlap in small boxes
layout_cascade_L <- function(nodes_order) {
  coords <- matrix(
    c(0,  1,
      0, -0.2,
      1, -0.2),
    ncol = 2, byrow = TRUE
  )
  rownames(coords) <- nodes_order
  coords
}

## ---- Helper: build Gamma0 adjacency (binary) for NetworkMotif() ----
## MR.RGM convention: rows = to, columns = from
make_gamma0 <- function(gene_names, edges_from, edges_to) {
  p <- length(gene_names)
  Gamma0 <- matrix(0, nrow = p, ncol = p, dimnames = list(gene_names, gene_names))

  stopifnot(length(edges_from) == length(edges_to))
  for (e in seq_along(edges_from)) {
    from <- edges_from[e]
    to   <- edges_to[e]
    Gamma0[to, from] <- 1
  }
  Gamma0
}

## ---- Output filenames ----
out_onek_feedback    <- "motif_onek_feedback.pdf"
out_onek_feedforward <- "motif_onek_feedforward.pdf"
out_onek_cascade     <- "motif_onek_cascade.pdf"

################################################################################
## Motif 1: Feedback loop (3-cycle)
## PLCY2 -> RAC -> ERK1/2 -> PLCY2
################################################################################
nodes_fb <- c("PLCY2", "RAC", "ERK1/2")
from_fb  <- c("PLCY2", "RAC", "ERK1/2")
to_fb    <- c("RAC",   "ERK1/2", "PLCY2")

g_fb <- make_motif_graph(nodes = nodes_fb, edges_from = from_fb, edges_to = to_fb)
coords_fb <- layout_triangle(nodes_fb)

save_motif_pdf(
  g = g_fb,
  layout_mat = layout_for_graph(g_fb, coords_fb),
  file = out_onek_feedback,
  vertex_size = 30,
  pad = 0.38
)

Gamma0_fb <- make_gamma0(gene_names_pathway_use, from_fb, to_fb)
prob_fb   <- NetworkMotif(Gamma0_fb, Output_OneK1K$GammaPst)


################################################################################
## Motif 2: Feedforward loop
## BCR -> P70S6K, BCR -> IKK, P70S6K -> IKK
################################################################################
nodes_ff <- c("BCR", "P70S6K", "IKK")
from_ff  <- c("BCR", "BCR", "P70S6K")
to_ff    <- c("P70S6K", "IKK", "IKK")

g_ff <- make_motif_graph(nodes = nodes_ff, edges_from = from_ff, edges_to = to_ff)
coords_ff <- layout_feedforward(nodes_ff)

save_motif_pdf(
  g = g_ff,
  layout_mat = layout_for_graph(g_ff, coords_ff),
  file = out_onek_feedforward,
  vertex_size = 30,
  pad = 0.38
)

Gamma0_ff <- make_gamma0(gene_names_pathway_use, from_ff, to_ff)
prob_ff   <- NetworkMotif(Gamma0_ff, Output_OneK1K$GammaPst)


################################################################################
## Motif 3: Cascade (chain)
## LYN -> RAC -> ERK1/2
################################################################################
nodes_cas <- c("LYN", "RAC", "ERK1/2")
from_cas  <- c("LYN", "RAC")
to_cas    <- c("RAC", "ERK1/2")

g_cas <- make_motif_graph(nodes = nodes_cas, edges_from = from_cas, edges_to = to_cas)
coords_cas <- layout_cascade_L(nodes_cas)

save_motif_pdf(
  g = g_cas,
  layout_mat = layout_for_graph(g_cas, coords_cas),
  file = out_onek_cascade,
  vertex_size = 26,
  label_cex = 0.9,
  pad = 0.45
)

Gamma0_cas <- make_gamma0(gene_names_pathway_use, from_cas, to_cas)
prob_cas   <- NetworkMotif(Gamma0_cas, Output_OneK1K$GammaPst)


## ---- Print results + outputs ----
cat(
  "Motif posterior probabilities (from GammaPst):\n",
  sprintf("  Feedback loop   : %.6f\n", prob_fb),
  sprintf("  Feedforward loop: %.6f\n", prob_ff),
  sprintf("  Cascade         : %.6f\n", prob_cas),
  "\nWrote PDFs:\n",
  "  ", out_onek_feedback, "\n",
  "  ", out_onek_feedforward, "\n",
  "  ", out_onek_cascade, "\n",
  sep = ""
)
