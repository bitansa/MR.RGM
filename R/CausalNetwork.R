# Function to create causal graph using only response-response relationships (Y_j -> Y_i)
create_causal_graph_Y <- function(AEst, zAEst) {

  # Extract adjacency matrix and weights
  zA <- zAEst  # Response-to-response adjacency (indicates where edges exist)
  A_weights <- AEst  # Weights for response-response edges

  p <- nrow(zA)  # Number of response variables (Y)

  # Node labels: the response names (row names of AEst) or Y1, ..., Yp
  response_nodes <- if (!is.null(rownames(AEst)) && !anyDuplicated(rownames(AEst))) rownames(AEst) else paste0("Y", 1:p)

  # Initialize edge list as a data frame
  edge_list <- data.frame(
    from = character(), to = character(), color = character(), label = numeric(), curved = numeric(),
    stringsAsFactors = FALSE
  )

  ### Add Response-to-Response edges (Yj -> Yi)
  for (i in 1:p) {
    for (j in 1:p) {
      if (zA[i, j] == 1) {  # edge from Y_j to Y_i
        edge_list <- rbind(edge_list, data.frame(
          from = response_nodes[j], to = response_nodes[i],  # edge from Y_j to Y_i
          color = "grey", label = round(A_weights[i, j], 2),
          curved = ifelse(zA[j, i] == 1, 0.7, 0)  # Curvature for bidirectional edges
        ))
      }
    }
  }

  # Create graph object with all nodes included
  graph <- igraph::graph_from_data_frame(d = edge_list, directed = TRUE, vertices = data.frame(name = response_nodes))

  # Set node attributes
  igraph::V(graph)$color <- "lightyellow"  # Node colour
  igraph::V(graph)$shape <- "circle"  # Y nodes are circles
  igraph::V(graph)$size <- c(rep(25, p))  # Sizes for Y


  # Set edge attributes
  igraph::E(graph)$color <- edge_list$color
  igraph::E(graph)$width <- 2.5  # Edge width (constant)
  igraph::E(graph)$label <- edge_list$label  # Assign causal effect labels
  igraph::E(graph)$arrow.size <- 1  # Arrow size
  igraph::E(graph)$curved <- edge_list$curved  # Separate bidirectional edges with curves

  # Return the graph object
  return(graph)
}
