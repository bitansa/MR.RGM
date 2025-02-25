# Function to create causal graph using AEst, BEst, zAEst, zBEst
create_causal_graph <- function(AEst, BEst, zAEst, zBEst) {

  # Extract adjacency matrices and weights
  zA <- zAEst  # Response-to-response adjacency (Y → Y)
  zB <- zBEst  # IV-to-response adjacency (X → Y)
  A_weights <- AEst  # Weights for response-response
  B_weights <- BEst  # Weights for IV-response

  p <- nrow(zA)  # Number of response variables (Y)
  k <- ncol(zB)  # Number of instrumental variables (X)

  # Create node labels
  response_nodes <- paste0("Y", 1:p)
  iv_nodes <- paste0("X", 1:k)
  all_nodes <- c(response_nodes, iv_nodes)

  edges <- c()
  edge_weights <- c()
  edge_colors <- c()
  edge_labels <- c()
  edge_curves <- c()  # Curvature to separate bidirectional edges

  ### Add Response-to-Response edges (Yj → Yi)
  for (i in 1:p) {
    for (j in 1:p) {
      if (zA[i, j] == 1) {  # If an edge exists from Yj → Yi
        edges <- c(edges, response_nodes[j], response_nodes[i])  # Correct direction
        edge_weights <- c(edge_weights, A_weights[i, j])
        edge_colors <- c(edge_colors, "blue")  # Use blue for Y → Y edges
        edge_labels <- c(edge_labels, round(A_weights[i, j], 2))

        # If the reverse edge (Yi → Yj) also exists, separate with curvature
        if (zA[j, i] == 1) {
          edge_curves <- c(edge_curves, 0.2)  # Slight curve for bidirectional edges
        } else {
          edge_curves <- c(edge_curves, 0)  # Keep unidirectional edges straight
        }
      }
    }
  }

  ### Add IV-to-Response edges (Xj → Yi)
  for (i in 1:p) {
    for (j in 1:k) {
      if (zB[i, j] != 0) {  # If zBEst[i, j] is nonzero, there is an edge Xj → Yi
        edges <- c(edges, iv_nodes[j], response_nodes[i])  # IV to response
        edge_weights <- c(edge_weights, B_weights[i, j])
        edge_colors <- c(edge_colors, "red")  # Use red for IV → Y edges
        edge_labels <- c(edge_labels, round(B_weights[i, j], 2))
        edge_curves <- c(edge_curves, 0)  # Keep IV → Y edges straight
      }
    }
  }

  # Create graph object
  graph <- igraph::graph(edges = edges, directed = TRUE)

  # Set node attributes
  igraph::V(graph)$color <- c(rep("lightblue", p), rep("lightgreen", k))  # Y in blue, X in green
  igraph::V(graph)$shape <- c(rep("circle", p), rep("square", k))  # Y as circles, X as squares
  igraph::E(graph)$color <- edge_colors
  igraph::E(graph)$width <- abs(edge_weights) * 2  # Scale edge thickness by weight
  igraph::E(graph)$label <- edge_labels  # Assign causal effect labels
  igraph::E(graph)$arrow.size <- 0.7  # Arrow size
  igraph::E(graph)$curved <- edge_curves  # Separate bidirectional edges with curves

  # Return the graph object
  return(graph)
}

# Function to create causal graph using only response-response relationships (Y → Y)
create_causal_graph_Y <- function(AEst, zAEst) {

  # Extract adjacency matrix and weights
  zA <- zAEst  # Response-to-response adjacency (indicates where edges exist)
  A_weights <- AEst  # Weights for response-response edges

  p <- nrow(zA)  # Number of response variables (Y)

  # Create node labels
  response_nodes <- paste0("Y", 1:p)

  edges <- c()
  edge_weights <- c()
  edge_colors <- c()
  edge_labels <- c()
  edge_curves <- c()  # Curve to separate bidirectional edges

  ### Add Response-to-Response edges (Yj → Yi)
  for (i in 1:p) {
    for (j in 1:p) {
      if (zA[i, j] == 1) {  # If an edge exists from Yj → Yi
        edges <- c(edges, response_nodes[j], response_nodes[i])  # Correct direction
        edge_weights <- c(edge_weights, A_weights[i, j])
        edge_colors <- c(edge_colors, "blue")  # Use blue for Y → Y edges
        edge_labels <- c(edge_labels, round(A_weights[i, j], 2))

        # If the reverse edge (Yi → Yj) also exists, separate with curvature
        if (zA[j, i] == 1) {
          edge_curves <- c(edge_curves, 0.2)  # Slight curve for bidirectional edges
        } else {
          edge_curves <- c(edge_curves, 0)  # Keep unidirectional edges straight
        }
      }
    }
  }

  # Create graph object
  graph <- igraph::graph(edges = edges, directed = TRUE)

  # Set node attributes
  igraph::V(graph)$color <- "lightblue"  # Y nodes are blue
  igraph::V(graph)$shape <- "circle"  # Y nodes are circles
  igraph::E(graph)$color <- edge_colors
  igraph::E(graph)$width <- abs(edge_weights) * 2  # Scale edge thickness by weight
  igraph::E(graph)$label <- edge_labels  # Assign causal effect labels
  igraph::E(graph)$arrow.size <- 0.7  # Arrow size
  igraph::E(graph)$curved <- edge_curves  # Separate bidirectional edges with curves

  # Return the graph object
  return(graph)
}
