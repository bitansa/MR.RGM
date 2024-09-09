#' Estimating the uncertainty of a specified network
#'
#' @description The NetworkMotif function facilitates uncertainty quantification.
#'              Specifically, it determines the proportion of posterior samples that contains the given network structure. To use this function, users may use the Gamma_Pst output obtained from the RGM function.
#'
#' @param Gamma A matrix of dimension p * p that signifies a specific network structure among the response variables, where p represents the number of response variables. This matrix is the focus of uncertainty quantification.
#' @param Gamma_Pst An array of dimension p * p * n_pst, where n_pst is the number of posterior samples and p denotes the number of response variables. It comprises the posterior samples of the causal network among the response variables. This input might be obtained from the RGM function. Initially, execute the RGM function and save the resulting Gamma_Pst. Subsequently, utilize this stored Gamma_Pst as input for this function.
#'
#' @return The NetworkMotif function calculates the uncertainty quantification for the provided network structure. A value close to 1 indicates that the given network structure is frequently observed in the posterior samples, while a value close to 0 suggests that the given network structure is rarely observed in the posterior samples.
#'
#'
#'
#' @export
#'
#' @examples
#'
#'
#'#' # ---------------------------------------------------------
#'
#' # Example 1:
#' # Run NetworkMotif to do uncertainty quantification for a given network among the response variable
#'
#' # Data Generation
#' set.seed(9154)
#'
#' # Number of data points
#' n = 10000
#'
#' # Number of response variables and number of instrument variables
#' p = 3
#' k = 4
#'
#' # Initialize causal interaction matrix between response variables
#' A = matrix(sample(c(-0.1, 0.1), p^2, replace = TRUE), p, p)
#'
#' # Diagonal entries of A matrix will always be 0
#' diag(A) = 0
#'
#' # Make the network sparse
#' A[sample(which(A!=0), length(which(A!=0))/2)] = 0
#'
#' # Initialize causal interaction matrix between response and instrument variables
#' B = matrix(0, p, k)
#'
#' # Create d vector
#' d = c(2, 1, 1)
#'
#'
#' # Initialize m
#' m = 1
#'
#' # Calculate B matrix based on d vector
#' for (i in 1:p) {
#'
#'  # Update ith row of B
#'  B[i, m:(m + d[i] - 1)] = 1
#'
#'  # Update m
#'  m = m + d[i]
#'
#' }
#'
#' Sigma = diag(p)
#'
#' Mult_Mat = solve(diag(p) - A)
#'
#' Variance = Mult_Mat %*% Sigma %*% t(Mult_Mat)
#'
#' # Generate instrument data matrix
#' X = matrix(rnorm(n * k, 0, 1), nrow = n, ncol = k)
#'
#' # Initialize response data matrix
#' Y = matrix(0, nrow = n, ncol = p)
#'
#' # Generate response data matrix based on instrument data matrix
#' for (i in 1:n) {
#'
#'     Y[i, ] = MASS::mvrnorm(n = 1, Mult_Mat %*% B %*% X[i, ], Variance)
#'
#' }
#'
#'
#' # Apply RGM on individual level data for Threshold Prior
#' Output = RGM(X = X, Y = Y, d = c(2, 1, 1), prior = "Threshold")
#'
#' # Store Gamma_Pst
#' Gamma_Pst = Output$Gamma_Pst
#'
#' # Define a function to create smaller arrowheads
#' smaller_arrowheads = function(graph) {
#'     igraph::E(graph)$arrow.size = 1  # Adjust the arrow size value as needed
#'     return(graph)
#' }
#'
#' # Start with a random subgraph
#' Gamma = matrix(0, nrow = p, ncol = p)
#' Gamma[2, 1] = 1
#'
#' # Plot the subgraph to get an idea about the causal network
#' plot(smaller_arrowheads(igraph::graph_from_adjacency_matrix(Gamma,
#'          mode = "directed")), layout = igraph::layout_in_circle,
#'             main = "Subgraph")
#'
#'
#' # Do uncertainty quantification for the subgraph
#' NetworkMotif(Gamma = Gamma, Gamma_Pst = Gamma_Pst)
#'
#'
#'
#'
#'
#' @references
#' Ni, Y., Ji, Y., & Müller, P. (2018).
#' Reciprocal graphical models for integrative gene regulatory network analysis.
#' \emph{Bayesian Analysis},
#' \strong{13(4)}, 1095-1110.
#' \doi{10.1214/17-BA1087}.
NetworkMotif = function(Gamma, Gamma_Pst) {

  # Check whether Gamma_Pst is a numeric array with three dimensions
  if (!is.numeric(Gamma_Pst) || !is.array(Gamma_Pst) || length(dim(Gamma_Pst)) != 3) {

    # Print an error message
    stop("Gamma_Pst must be a numeric array with three dimensions.")

  }

  # Calculate number of rows of Gamma_Pst
  p = nrow(Gamma_Pst)

  # Check whether number of rows of Gamma_Pst is same as number of columns of Gamma_Pst
  if(ncol(Gamma_Pst) != p){

    # Print an error message
    stop("Number of rows and columns of Gamma_Pst should be equal.")

  }

  # Check whether Gamma is a numeric matrix
  if(!is.numeric(Gamma) || !is.matrix(Gamma)){

    # Print an error message
    stop("Gamma should be a numeric matrix.")

  }

  # Check whether Gamma is a square matrix with dimension p * p
  if((nrow(Gamma) != p) || (ncol(Gamma) != p)) {

    # Print an error message
    stop("Gamma should be a square matrix with number of rows and columns equal to number of rows of Gamma_Pst.")

  }

  # Return Network Motif
  return(NetworkMotif_cpp(Gamma = Gamma, Gamma_Pst = Gamma_Pst))

}
