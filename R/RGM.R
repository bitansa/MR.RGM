#' Fitting Bayesian Multivariate Bidirectional Mendelian Randomization Networks
#'
#' @description The RGM function transforms causal inference by merging Mendelian randomization and network-based methods, enabling the creation of comprehensive causal graphs within complex biological systems. RGM accommodates varied data contexts with three input options: individual-level data (X, Y matrices), summary-level data including S_YY, S_YX, and S_XX matrices, and intricate data with challenging cross-correlations, utilizing S_XX, Beta, and Sigma_Hat matrices.
#'              For the latter input, data centralization is necessary. Users can select any of these data formats to suit their needs and don’t have to specify all of them, allowing flexibility based on data availability. Crucial inputs encompass "d" (instrument count per response) and "n" (total observations, only required for summary level data), amplified by customizable parameters that refine analysis. Additionally, users can tailor the analysis by setting parameters such as "nIter" (number of MCMC iterations), "nBurnin" (number of discarded samples during burn-in for convergence), and "Thin" (thinning of posterior samples). These customizable parameters enhance the precision and relevance of the analysis.
#'              RGM provides essential causal effect/strength estimates between response variables and between response and instrument variables. Moreover, it furnishes adjacency matrices, visually mapping causal graph structures. These outputs empower researchers to untangle intricate relationships within biological networks, fostering a holistic understanding of complex systems.
#'
#' @param X A matrix of dimension n * k. In this matrix, each row signifies a distinct observation, while each column represents a specific instrument variable. The default value is set to NULL.
#' @param Y A matrix of dimension n * p. In this matrix, each row corresponds to a specific observation, and each column pertains to a particular response variable. The default value is set to NULL.
#' @param S_YY A matrix of dimensions p * p. Here, "p" signifies the count of response variables. This matrix is derived through the operation t(Y) %*% Y / n, where "Y" denotes the response data matrix and "n" stands for the total number of observations.
#' @param S_YX A matrix of dimensions p * k. Here, "p" signifies the number of response variables, and "k" represents the count of instrument variables. This matrix is calculated using the operation t(Y) %*% X / n, where "Y" is the response data matrix, "X" is the instrument data matrix and "n" is the total number of observations.
#' @param S_XX A matrix of dimensions k * k. Here, "k" signifies the count of instrument variables. This matrix is derived through the operation t(X) %*% X / n, where "X" denotes the instrument data matrix and "n" stands for the total number of observations.
#' @param Beta A matrix of dimensions p * k. In this matrix, each row corresponds to a specific response variable, and each column pertains to a distinct instrument variable. Each entry within the matrix represents the regression coefficient of the individual response variable on the specific instrument variable. To use Beta as an input, ensure you centralize each column of Y i.e. response data matrix and X i.e. instrument data matrix before calculating Beta, S_XX, and Sigma_Hat.
#' @param Sigma_Hat A matrix of dimensions p * k. In this matrix, each row corresponds to a specific response variable, and each column pertains to an individual instrument variable. Each entry in this matrix represents the mean square error associated with regressing the particular response on the specific instrument variable. To employ Sigma_Hat as an input, ensure that you centralize each column of Y i.e. response data matrix and X i.e. instrument data matrix before calculating Beta, S_XX, and Sigma_Hat.
#' @param d A vector input with a length of p i.e. number of response variables. Each element within this vector is a positive integer denoting the count of instrument variables influencing a specific response variable. The sum of all elements in the vector should be equal to the total count of instrument variables, represented as k.
#' @param n A positive integer input representing the count of data points or observations in the dataset. This input is only required when summary level data is used as input.
#' @param nIter A positive integer input representing the number of MCMC (Markov Chain Monte Carlo) sampling iterations. The default value is set to 10,000.
#' @param nBurnin A non-negative integer input representing the number of samples to be discarded during the burn-in phase of MCMC sampling. It's important that nBurnin is less than nIter. The default value is set to 2000.
#' @param Thin A positive integer input denoting the thinning factor applied to posterior samples. Thinning reduces the number of samples retained from the MCMC process for efficiency. Thin should not exceed (nIter - nBurnin). The default value is set to 1.
#' @param prior A parameter representing the prior assumption on the graph structure. It offers two options: "Threshold" or "Spike and Slab". The default value is "Spike and Slab".
#' @param a_rho A positive scalar input representing the first parameter of a Beta distribution. The default value is set to 3.
#' @param b_rho A positive scalar input representing the second parameter of a Beta distribution. The default value is set to 1.
#' @param nu_1 A positive scalar input representing the multiplication factor in the variance of the spike part in the spike and slab distribution of matrix A. The default value is set to 0.001.
#' @param a_psi A positive scalar input corresponding to the first parameter of a Beta distribution. The default value is set to 0.5.
#' @param b_psi  A positive scalar input corresponding to the second parameter of a Beta distribution. The default value is set to 0.5.
#' @param nu_2 A positive scalar input corresponding to the multiplication factor in the variance of the spike part in the spike and slab distribution of matrix B. The default value is set to 0.0001.
#' @param a_sigma A positive scalar input corresponding to the first parameter of an Inverse Gamma distribution, which is associated with the variance of the model. The default value is set to 0.01.
#' @param b_sigma A positive scalar input corresponding to the second parameter of an Inverse Gamma distribution, which is associated with the variance of the model. The default value is set to 0.01.
#' @param Prop_VarA A positive scalar input representing the variance of the normal distribution used for proposing terms within the A matrix. The default value is set to 0.01.
#' @param Prop_VarB A positive scalar input representing the variance of the normal distribution used for proposing terms within the B matrix. The default value is set to 0.01.
#'
#' @return
#'
#' \item{A_Est}{A matrix of dimensions p * p, representing the estimated causal effects or strengths between the response variables.}
#' \item{B_Est}{A matrix of dimensions p * k, representing the estimated causal effects or strengths between the response variables and the instrument variables. Each row corresponds to a specific response variable, and each column corresponds to a particular instrument variable.}
#' \item{zA_Est}{A binary adjacency matrix of dimensions p * p, indicating the graph structure between the response variables. Each entry in the matrix represents the presence (1) or absence (0) of a causal link between the corresponding response variables.}
#' \item{zB_Est}{A binary adjacency matrix of dimensions p * k, illustrating the graph structure between the response variables and the instrument variables. Each row corresponds to a specific response variable, and each column corresponds to a particular instrument variable. The presence of a causal link is denoted by 1, while the absence is denoted by 0.}
#' \item{A0_Est}{A matrix of dimensions p * p, representing the estimated causal effects or strengths between response variables before thresholding. This output is particularly relevant for cases where the "Threshold" prior assumption is utilized.}
#' \item{B0_Est}{A matrix of dimensions p * k, representing the estimated causal effects or strengths between the response variables and the instrument variables before thresholding. This output is particularly relevant for cases where the "Threshold" prior assumption is utilized. Each row corresponds to a specific response variable, and each column corresponds to a particular instrument variable.}
#' \item{Gamma_Est}{A matrix of dimensions p * p, representing the estimated probabilities of edges between response variables in the graph structure. Each entry in the matrix indicates the probability of a causal link between the corresponding response variables.}
#' \item{Tau_Est}{A matrix of dimensions p * p, representing the estimated variances of causal interactions between response variables. Each entry in the matrix corresponds to the variance of the causal effect between the corresponding response variables.}
#' \item{Phi_Est}{A matrix of dimensions p * k, representing the estimated probabilities of edges between response and instrument variables in the graph structure. Each row corresponds to a specific response variable, and each column corresponds to a particular instrument variable.}
#' \item{Eta_Est}{A matrix of dimensions p * k, representing the estimated variances of causal interactions between response and instrument variables. Each row corresponds to a specific response variable, and each column corresponds to a particular instrument variable.}
#' \item{tA_Est}{A scalar value representing the estimated thresholding value of causal interactions between response variables. This output is relevant when using the "Threshold" prior assumption.}
#' \item{tB_Est}{A scalar value representing the estimated thresholding value of causal interactions between response and instrument variables. This output is applicable when using the "Threshold" prior assumption.}
#' \item{Sigma_Est}{A vector of length p, representing the estimated variances of each response variable. Each element in the vector corresponds to the variance of a specific response variable.}
#' \item{AccptA}{The percentage of accepted entries in the A matrix, which represents the causal interactions between response variables. This metric indicates the proportion of proposed changes that were accepted during the sampling process.}
#' \item{AccptB}{The percentage of accepted entries in the B matrix, which represents the causal interactions between response and instrument variables. This metric indicates the proportion of proposed changes that were accepted during the sampling process.}
#' \item{Accpt_tA}{The percentage of accepted thresholding values for causal interactions between response variables when using the "Threshold" prior assumption. This metric indicates the proportion of proposed thresholding values that were accepted during the sampling process.}
#' \item{Accpt_tB}{The percentage of accepted thresholding values for causal interactions between response and instrument variables when using the "Threshold" prior assumption. This metric indicates the proportion of proposed thresholding values that were accepted during the sampling process.}
#' \item{LL_Pst}{A vector containing the posterior log-likelihoods of the model. Each element in the vector represents the log-likelihood of the model given the observed data and the estimated parameters.}
#' \item{Rho_Est}{A matrix of dimensions p * p, representing the estimated Bernoulli success probabilities of causal interactions between response variables when using the "Spike and Slab" prior assumption. Each entry in the matrix corresponds to the success probability of a causal interaction between the corresponding response variables.}
#' \item{Psi_Est}{A matrix of dimensions p * k, representing the estimated Bernoulli success probabilities of causal interactions between response and instrument variables when using the "Spike and Slab" prior assumption. Each row in the matrix corresponds to a specific response variable, and each column corresponds to a particular instrument variable.}
#' \item{Gamma_Pst}{An array containing the posterior samples of the network structure among the response variables.}
#'
#'
#'
#'
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
#' # Run RGM based on individual level data with Threshold prior based on the model Y = AY + BX + E
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
#' # Define a function to create smaller arrowheads
#' smaller_arrowheads = function(graph) {
#'     igraph::E(graph)$arrow.size = 0.25  # Adjust the arrow size value as needed
#'     return(graph)
#' }
#'
#' # Print true causal interaction matrices between response variables
#' # and between response and instrument variables
#' A
#' B
#'
#'
#' # Plot the true graph structure between response variables
#' plot(smaller_arrowheads(igraph::graph.adjacency(((A != 0) * 1),
#'  mode = "directed")), layout = igraph::layout_in_circle, main = "True Graph")
#'
#' # Apply RGM on individual level data for Threshold Prior
#' Output = RGM(X = X, Y = Y, d = c(2, 1, 1), prior = "Threshold")
#'
#' # Get the graph structure between response variables
#' Output$zA_Est
#'
#' # Plot the estimated graph structure between response variables
#' plot(smaller_arrowheads(igraph::graph.adjacency(Output$zA_Est,
#'  mode = "directed")), layout = igraph::layout_in_circle, main = "Estimated Graph")
#'
#' # Get the estimated causal strength matrix between response variables
#' Output$A_Est
#'
#' # Get the graph structure between response and instrument variables
#' Output$zB_Est
#'
#' # Get the estimated causal strength matrix between response and instrument variables
#' Output$B_Est
#'
#' # Plot posterior log-likelihood
#' plot(Output$LL_Pst, type = 'l', xlab = "Number of Iterations", ylab = "Log-likelihood")
#'
#'
#'
#'
#' # -----------------------------------------------------------------
#' # Example 2:
#' # Run RGM based on summary level data with Spike and Slab prior based on the model Y = AY + BX + E
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
#' # Calculate summary level data
#' S_YY = t(Y) %*% Y / n
#' S_YX = t(Y) %*% X / n
#' S_XX = t(X) %*% X / n
#'
#'
#' # Print true causal interaction matrices between response variables
#' # and between response and instrument variables
#' A
#' B
#'
#' # Plot the true graph structure between response variables
#' plot(smaller_arrowheads(igraph::graph.adjacency(((A != 0) * 1),
#'  mode = "directed")), layout = igraph::layout_in_circle, main = "True Graph")
#'
#'
#' # Apply RGM on summary level data for Spike and Slab Prior
#' Output = RGM(S_YY = S_YY, S_YX = S_YX, S_XX = S_XX,
#'           d = c(2, 1, 1), n = 10000, prior = "Spike and Slab")
#'
#' # Get the graph structure between response variables
#' Output$zA_Est
#'
#' # Plot the estimated graph structure between response variables
#' plot(smaller_arrowheads(igraph::graph.adjacency(Output$zA_Est,
#'  mode = "directed")), layout = igraph::layout_in_circle, main = "Estimated Graph")
#'
#' # Get the estimated causal strength matrix between response variables
#' Output$A_Est
#'
#' # Get the graph structure between response and instrument variables
#' Output$zB_Est
#'
#' # Get the estimated causal strength matrix between response and instrument variables
#' Output$B_Est
#'
#' # Plot posterior log-likelihood
#' plot(Output$LL_Pst, type = 'l', xlab = "Number of Iterations", ylab = "Log-likelihood")
#'
#'
#'
#'
#' # -----------------------------------------------------------------
#' # Example 3:
#' # Run RGM based on Beta and Sigma_Hat with Spike and Slab prior based on the model Y = AY + BX + E
#'
#' # Data Generation
#' set.seed(9154)
#'
#' # Number of datapoints
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
#' # Generate DNA expressions
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
#' # Centralize Data
#' Y = t(t(Y) - colMeans(Y))
#' X = t(t(X) - colMeans(X))
#'
#' # Calculate S_XX
#' S_XX = t(X) %*% X / n
#'
#' # Generate Beta matrix and Sigma_Hat
#' Beta = matrix(0, nrow = p, ncol = k)
#' Sigma_Hat = matrix(0, nrow = p, ncol = k)
#'
#' for (i in 1:p) {
#'
#'    for (j in 1:k) {
#'
#'        fit = lm(Y[, i] ~ X[, j])
#'
#'        Beta[i, j] =  fit$coefficients[2]
#'
#'        Sigma_Hat[i, j] = sum(fit$residuals^2) / n
#'
#'        }
#'
#'    }
#'
#'
#' # Print true causal interaction matrices between response variables
#' # and between response and instrument variables
#' A
#' B
#'
#'
#' # Plot the true graph structure between response variables
#' plot(smaller_arrowheads(igraph::graph.adjacency(((A != 0) * 1),
#'  mode = "directed")), layout = igraph::layout_in_circle, main = "True Graph")
#'
#'
#' # Apply RGM based on S_XX, Beta and Sigma_Hat for Spike and Slab Prior
#' Output = RGM(S_XX = S_XX, Beta = Beta, Sigma_Hat = Sigma_Hat,
#'           d = c(2, 1, 1), n = 10000, prior = "Spike and Slab")
#'
#' # Get the graph structure between response variables
#' Output$zA_Est
#'
#' # Plot the estimated graph structure between response variables
#' plot(smaller_arrowheads(igraph::graph.adjacency(Output$zA_Est,
#'  mode = "directed")), layout = igraph::layout_in_circle, main = "Estimated Graph")
#'
#' # Get the estimated causal strength matrix between response variables
#' Output$A_Est
#'
#' # Get the graph structure between response and instrument variables
#' Output$zB_Est
#'
#' # Get the estimated causal strength matrix between response and instrument variables
#' Output$B_Est
#'
#' # Plot posterior log-likelihood
#' plot(Output$LL_Pst, type = 'l', xlab = "Number of Iterations", ylab = "Log-likelihood")
#'
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
RGM = function(X = NULL, Y = NULL, S_YY = NULL, S_YX = NULL, S_XX = NULL, Beta = NULL, Sigma_Hat = NULL, d, n, nIter = 10000, nBurnin = 2000, Thin = 1, prior = c("Threshold", "Spike and Slab"), a_rho = 3, b_rho = 1, nu_1 = 0.001, a_psi = 0.5, b_psi = 0.5, nu_2 = 0.0001, a_sigma = 0.01, b_sigma = 0.01, Prop_VarA = 0.01, Prop_VarB = 0.01){

  # Check whether Y or S_YY is given as data input
  if((!is.null(Y) && is.null(X)) || (!is.null(S_YY) && is.null(S_YX) && is.null(S_XX))){

    # Check whether input Y is given
    if(!is.null(Y)){

      # Check whether Y is a numeric matrix
      if(!is.numeric(Y) || !is.matrix(Y)){

        # Print an error message
        stop("Y should be numeric matrix.")

      }

      # Calculate number of datapoints
      n = nrow(Y)

      # Calculate number of response variables from Y matrix
      p = ncol(Y)

      # Calculate S_YY
      S_YY = t(Y) %*% Y / n

    } else {

      # Check whether n is a positive integer
      if(!is.numeric(n) || n != round(n) || n <= 0){

        # Print an error message
        stop("Number of datapoints should be a positive integer.")

      }

      # Check whether S_YY is numeric matrix
      if(!is.numeric(S_YY) || !is.matrix(S_YY)){

        # Print an error message
        stop("S_YY should be a numeric matrix.")

      }

      # Calculate number of response variables from S_YY matrix
      p = ncol(S_YY)

      # Check whether number of rows of S_YY is equal to p
      if(nrow(S_YY) != p){

        # Print an error message
        stop("S_YY should be a square matrix.")

      }


    }


    # Check whether all the beta and inverse gamma parameters are positive or not
    if(!is.numeric(a_rho) || a_rho < 0 || !is.numeric(b_rho) || b_rho < 0  || !is.numeric(a_sigma) || a_sigma < 0 || !is.numeric(b_sigma) || b_sigma < 0){

      # Print an error message
      stop("All the beta and inverse gamma parameters should be positive.")

    }


    # Check whether the variance terms are positive or not
    if(!is.numeric(nu_1) || nu_1 <= 0 || !is.numeric(Prop_VarA) || Prop_VarA <= 0){

      # Print an error message
      stop("All the variance terms should be positive.")

    }


    # Check whether nIter is a postive integer
    if(!is.numeric(nIter) || nIter != round(nIter) || nIter <= 0){

      # Print an error message
      stop("Number of iterations should be a positive integer.")

    }

    # Check whether nBurnin is a non-negative integer and strictly less than nIter
    if(!is.numeric(nBurnin) || nBurnin != round(nBurnin) || nBurnin < 0 || nBurnin >= nIter){

      # Print an error message
      stop("Number of Burnin points should be a nonnegative integer strictly less than number of iterations.")

    }

    # Check whether Thin is a positive integer less than equal to (nIter - nBurnin)
    if(!is.numeric(Thin) || Thin != round(Thin) || Thin <= 0 || Thin > (nIter - nBurnin)){

      # Print an error message
      stop("Thin should be a positive integer less than or equal to the difference between number of iterations and number of burnin points.")

    }

    # Apply RGM for Spike and Slab prior and Threshold prior
    if ("Spike and Slab" %in% prior){

      # Run the algorithm for Spike and Slab prior
      Output = RGM_SpikeSlab1(S_YY, n, nIter = nIter, nBurnin = nBurnin, Thin = Thin,
                              a_rho = a_rho, b_rho = b_rho, nu_1 = nu_1,
                              a_sigma = a_sigma, b_sigma = b_sigma, Prop_VarA = Prop_VarA)



      # Get the posterior samples of gamma
      Gamma_Pst = Output$Gamma_Pst

      # Tag Gamma_Pst
      attr(Gamma_Pst, "RGM_GammaPst") = TRUE


      # Return outputs
      return(list(A_Est = Output$A_Est, zA_Est = Output$zA_Est,
                  Gamma_Est = Output$Gamma_Est, Tau_Est = Output$Tau_Est, Rho_Est = Output$Rho_Est,
                  Sigma_Est = Output$Sigma_Est,
                  AccptA = Output$AccptA,
                  LL_Pst = Output$LL_Pst, Gamma_Pst = Gamma_Pst))



    } else if ("Threshold" %in% prior){

      # Run the algorithm for Threshold prior
      Output = RGM_Threshold1(S_YY, n, nIter = nIter, nBurnin = nBurnin, Thin = Thin, nu_1 = nu_1,
                              a_sigma = a_sigma, b_sigma = b_sigma, Prop_VarA = Prop_VarA)


      # Get the posterior samples of gamma
      Gamma_Pst = Output$Gamma_Pst

      # Tag Gamma_Pst
      attr(Gamma_Pst, "RGM_GammaPst") = TRUE



      # Return outputs
      return(list(A_Est = Output$A_Est, zA_Est = Output$zA_Est,
                  A0_Est = Output$A0_Est, Gamma_Est = Output$Gamma_Est, Tau_Est = Output$Tau_Est,
                  tA_Est = Output$tA_Est,
                  Sigma_Est = Output$Sigma_Est,
                  AccptA = Output$AccptA, Accpt_tA = Output$Accpt_tA,
                  LL_Pst = Output$LL_Pst, Gamma_Pst = Gamma_Pst))





    } else {

      # Print an error message
      stop("Please specify a prior among Threshold prior and Spike and Slab prior.")


    }


  } else if((!is.null(X) && !is.null(Y)) || (!is.null(S_YY) && !is.null(S_YX) && !is.null(S_XX)) || (!is.null(S_XX) && !is.null(Beta) && !is.null(Sigma_Hat))){


    # Check whether input X and Y are given
    if(!is.null(X) && !is.null(Y)){

      # Check whether X and Y are both numeric matrices
      if(!is.numeric(X) || !is.numeric(Y) || !is.matrix(X) || !is.matrix(Y)){

        # Print an error message
        stop("X and Y should be numeric matrices.")

      }

      # Calculate number of datapoints
      n = nrow(X)

      # Check whether number of rows of X and Y are equal or not
      if(nrow(Y) != n){

        # Print an error message
        stop("Number of rows of X and Y should be equal.")

      }

      # Calculate number of response variables from Y matrix
      p = ncol(Y)

      # Calculate number of instrument variables from X matrix
      k = ncol(X)

      # Calculate S_YY, S_YX, S_XX
      S_YY = t(Y) %*% Y / n
      S_YX = t(Y) %*% X / n
      S_XX = t(X) %*% X / n

    } else if(!is.null(S_YY) && !is.null(S_YX) && !is.null(S_XX)){

      # Check whether n is a positive integer
      if(!is.numeric(n) || n != round(n) || n <= 0){

        # Print an error message
        stop("Number of datapoints should be a positive integer.")

      }

      # Check whether S_YY ,S_YX and S_XX are numeric matrices
      if(!is.numeric(S_YY) || !is.numeric(S_YX) || !is.numeric(S_XX) || !is.matrix(S_YY) || !is.matrix(S_YX) || !is.matrix(S_XX)){

        # Print an error message
        stop("S_YY, S_YX and S_XX should be numeric matrices.")

      }

      # Calculate number of response variables from S_YY matrix
      p = ncol(S_YY)

      # Calculate number of instrument variables from S_XX matrix
      k = ncol(S_XX)

      # Check whether number of rows of S_YY is equal to p
      if(nrow(S_YY) != p){

        # Print an error message
        stop("S_YY should be a square matrix.")

      }

      # Check whether number of rows of S_XX is equal to k
      if(nrow(S_XX) != k){

        # Print an error message
        stop("S_XX should be a square matrix.")

      }

      # Check whether number of rows of S_YX is equal to p and number of columns of S_YX is equal to k
      if(nrow(S_YX) != p || ncol(S_YX) != k){

        # Print an error message
        stop("Number of rows of S_YX should be equal to number of rows of S_YY and number of columns of S_YX should be equal to number of columns of S_XX.")

      }



    } else {

      # Check whether n is a positive integer
      if(!is.numeric(n) || n != round(n) || n <= 0){

        # Print an error message
        stop("Number of datapoints should be a positive integer.")

      }

      # Check whether S_XX, Beta and Sigma_Hat are numeric matrices
      if(!is.numeric(S_XX) || !is.numeric(Beta) || !is.numeric(Sigma_Hat) || !is.matrix(S_XX) || !is.matrix(Beta) || !is.matrix(Sigma_Hat)){

        # Print an error message
        stop("S_XX, Beta and SIgma_Hat should be numeric matrices.")

      }

      # Calculate number of response variables from Beta matrix
      p = nrow(Beta)

      # Calculate number of instrument variables from S_XX matrix
      k = ncol(S_XX)

      # Check whether number of rows of S_XX is equal to k
      if(nrow(S_XX) != k){

        # Print an error message
        stop("S_XX should be a square matrix.")

      }

      # Check whether number of columns of Beta is equal to k
      if(ncol(Beta) != k){

        # Print an error message
        stop("Number of cloumns of Beta should be equal to number of columns of S_XX.")

      }

      # Check whether number of rows of Sigma_Hat is equal to p and number of columns of Sigma_Hat is equal to k
      if(nrow(Sigma_Hat) != p || ncol(Sigma_Hat) != k){

        # Print an error message
        stop("Number of rows of Sigma_Hat should be equal to number of rows of Beta and number of columns of Sigma_Hat should be equal to number of columns of S_XX.")

      }



      # Calculate S_YX matrix
      S_YX = t(t(Beta) *  diag(S_XX))

      # Check whether d is a vector of non-negative integers of length p and sum od entries of d is equal to k
      if(!is.numeric(d) || sum(d != round(d)) != 0 || sum(d <= 0) != 0 || length(d) != p || sum(d) != k){

        # Print an error message
        stop("d should be a vector of positive integers of length equal to number of nodes and sum of entries should be equal to number of covariates.")

      }

      # Store indices to extract particular columns from beta matrix
      Col_Ind = c(0, cumsum(d))
      Col_Ind = Col_Ind[-length(Col_Ind)] + 1

      # Calculate Beta_New matrix by taking particular columns from Beta using Col_Ind
      Beta_New = matrix(Beta[, Col_Ind], nrow = p, ncol = p)

      # Initiate Estimate of A
      A_Est = matrix(0, nrow = p, ncol = p)

      # Update A_Est
      for (i in 1:p) {

        # Remove ith row and ith column of Beta_New
        Beta_ii = Beta_New[-i, -i]
        Beta_ii_det = det(Beta_ii)

        for (j in 1:p) {

          if(j != i){

            # Remove jth row and ith column of Beta_New
            beta_ji = Beta_New[-j, -i]

            A_Est[i, j] = det(beta_ji) / Beta_ii_det * (-1)^(i + j + 1)

          }

        }

      }

      # Calculate (I-A_Est)^(-1)
      Mult_Mat = solve(diag(p) - A_Est)

      # Calculate diagonal entries of S_YY
      S_YY_Diag = Sigma_Hat[, 1] + 2 * Beta[, 1] * S_YX[, 1] - Beta[, 1]^2 * S_XX[1, 1]

      # Calculate correlation matrix between X and between Y and X
      Cor_X = stats::cov2cor(S_XX)
      Cor_YX = t(t(Beta * (1 / sqrt(S_YY_Diag))) * sqrt(diag(S_XX)))

      # Calculate determinant of Cor_X
      Cor_X_Det = det(Cor_X)

      # Calculate sum of squares of error while fitting Y[, i] on X
      Error_Sumsq = rep(0, p)

      for (i in 1:p) {

        # Create R
        R = cbind(c(1, Cor_YX[i,]), rbind(Cor_YX[i,], Cor_X))

        # Update Error_Sumsq
        Error_Sumsq[i] = S_YY_Diag[i] * det(R) / Cor_X_Det


      }

      # calculate Sigma for the model
      Sigma = solve(Mult_Mat^2, Error_Sumsq)

      # Calculate S_YY
      S_YY = Beta %*% S_XX %*% t(Beta)  + Mult_Mat %*% (t(Mult_Mat) * Sigma)


    }


    # Check whether d is a vector of positive integers of length p and sum of entries of d is equal to k
    if(!is.numeric(d) || sum(d != round(d)) != 0 || sum(d <= 0) != 0 || length(d) != p || sum(d) != k){

      # Print an error message
      stop("d should be a vector of positive integers of length equal to number of nodes and sum of entries should be equal to number of covariates.")

    }

    # Initialize D matrix
    D = matrix(0, p, k)

    # Initialize m
    m = 1

    # Calculate D matrix based on d vector
    for (i in 1:p) {

      # Update ith row of D
      D[i, m:(m + d[i] - 1)] = 1

      # Update m
      m = m + d[i]

    }


    # Check whether all the beta and inverse gamma parameters are positive or not
    if(!is.numeric(a_rho) || a_rho < 0 || !is.numeric(b_rho) || b_rho < 0 || !is.numeric(a_psi) || a_psi < 0 || !is.numeric(b_psi) || b_psi < 0 || !is.numeric(a_sigma) || a_sigma < 0 || !is.numeric(b_sigma) || b_sigma < 0){

      # Print an error message
      stop("All the beta and inverse gamma parameters should be positive.")

    }


    # Check whether the variance terms are positive or not
    if(!is.numeric(nu_1) || nu_1 <= 0 || !is.numeric(nu_2) || nu_2 <= 0 || !is.numeric(Prop_VarA) || Prop_VarA <= 0 || !is.numeric(Prop_VarB) || Prop_VarB <= 0){

      # Print an error message
      stop("All the variance terms should be positive.")

    }


    # Check whether nIter is a postive integer
    if(!is.numeric(nIter) || nIter != round(nIter) || nIter <= 0){

      # Print an error message
      stop("Number of iterations should be a positive integer.")

    }

    # Check whether nBurnin is a non-negative integer and strictly less than nIter
    if(!is.numeric(nBurnin) || nBurnin != round(nBurnin) || nBurnin < 0 || nBurnin >= nIter){

      # Print an error message
      stop("Number of Burnin points should be a nonnegative integer strictly less than number of iterations.")

    }

    # Check whether Thin is a positive integer less than equal to (nIter - nBurnin)
    if(!is.numeric(Thin) || Thin != round(Thin) || Thin <= 0 || Thin > (nIter - nBurnin)){

      # Print an error message
      stop("Thin should be a positive integer less than or equal to the difference between number of iterations and number of burnin points.")

    }


    # Apply RGM for Spike and Slab prior and Threshold prior
    if ("Spike and Slab" %in% prior){

      # Run the algorithm for Spike and Slab prior
      Output = RGM_SpikeSlab2(S_YY, S_YX, S_XX, D, n, nIter = nIter, nBurnin = nBurnin, Thin = Thin,
                              a_rho = a_rho, b_rho = b_rho, nu_1 = nu_1, a_psi = a_psi, b_psi = b_psi,
                              nu_2 = nu_2, a_sigma = a_sigma, b_sigma = b_sigma, Prop_VarA = Prop_VarA, Prop_VarB = Prop_VarB)

      # Get the posterior samples of gamma
      Gamma_Pst = Output$Gamma_Pst

      # Tag Gamma_Pst
      attr(Gamma_Pst, "RGM_GammaPst") = TRUE

      # Return outputs
      return(list(A_Est = Output$A_Est, B_Est = Output$B_Est, zA_Est = Output$zA_Est, zB_Est = Output$zB_Est,
                  Gamma_Est = Output$Gamma_Est, Tau_Est = Output$Tau_Est, Rho_Est = Output$Rho_Est,
                  Phi_Est = Output$Phi_Est, Eta_Est = Output$Eta_Est, Psi_Est = Output$Psi_Est,
                  Sigma_Est = Output$Sigma_Est,
                  AccptA = Output$AccptA, AccptB = Output$AccptB,
                  LL_Pst = Output$LL_Pst, Gamma_Pst = Gamma_Pst))



    } else if ("Threshold" %in% prior){

      # Run the algorithm for Threshold prior
      Output = RGM_Threshold2(S_YY, S_YX, S_XX, D, n, nIter = nIter, nBurnin = nBurnin, Thin = Thin, nu_1 = nu_1, nu_2 = nu_2,
                              a_sigma = a_sigma, b_sigma = b_sigma, Prop_VarA = Prop_VarA, Prop_VarB = Prop_VarB)

      # Get the posterior samples of gamma
      Gamma_Pst = Output$Gamma_Pst

      # Tag Gamma_Pst
      attr(Gamma_Pst, "RGM_GammaPst") = TRUE


      # Return outputs
      return(list(A_Est = Output$A_Est, B_Est = Output$B_Est, zA_Est = Output$zA_Est, zB_Est = Output$zB_Est,
                  A0_Est = Output$A0_Est, B0_Est = Output$B0_Est, Gamma_Est = Output$Gamma_Est, Tau_Est = Output$Tau_Est,
                  Phi_Est = Output$Phi_Est, Eta_Est = Output$Eta_Est, tA_Est = Output$tA_Est, tB_Est = Output$tB_Est,
                  Sigma_Est = Output$Sigma_Est,
                  AccptA = Output$AccptA, AccptB = Output$AccptB, Accpt_tA = Output$Accpt_tA, Accpt_tB = Output$Accpt_tB,
                  LL_Pst = Output$LL_Pst, Gamma_Pst = Gamma_Pst))





    } else {

      # Print an error message
      stop("Please specify a prior among Threshold prior and Spike and Slab prior.")


    }


  } else {

    # Print an error message
    stop("Please give Y or S_YY or X, Y or S_YY, S_YX, S_XX or S_XX, Beta and Sigma_Hat as input.")


  }

}
