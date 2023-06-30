#' Title Fitting Reciprocal Graphical Models for Integrative Gene Regulatory Network
#'
#' @description RGM can be used to fit Reciprocal Graphical Models on gene expression data and DNA level measurements to find the relationship between different genes and the relationship between genes and DNA.
#'
#'
#' @param X n * k matrix data input, each row denotes a particular observation and each column denotes a particular DNA.
#' @param Y n * p matrix data input, each row denotes a particular observation and each column denotes a particular protein.
#' @param d Vector input of length p. Each element is a non-negative integer corresponds to number of covariates affecting a particular protein, sum of which should be k.
#' @param b_tau positive scalar input. It corresponds to the second parameter of a Inverse Gamma distribution. Default value is 0.1.
#' @param a_rho positive scalar input. It corresponds to the first parameter of a Beta distribution. Default value is 0.5.
#' @param b_rho positive scalar input. It corresponds to the first parameter of a Beta distribution. Default value is 0.5.
#' @param nu_1 positive scalar input. It corresponds to the variance of slab part of of spike and slab distribution of A. Default value is 0.0001.
#' @param a_eta positive scalar input. It corresponds to the first parameter of a Inverse Gamma distribution. Default value is 0.1.
#' @param b_eta positive scalar input. It corresponds to the second parameter of a Inverse Gamma distribution. Default value is 0.1.
#' @param a_psi positive scalar input. It corresponds to the first parameter of a Beta distribution. Default value is 0.5.
#' @param b_psi positive scalar input. It corresponds to the first parameter of a Beta distribution. Default value is 0.5.
#' @param nu_2 positive scalar input. It corresponds to the variance of slab part of of spike and slab distribution of B. Default value is 0.0001.
#' @param a_sigma positive scalar input. It corresponds to the first parameter of a Inverse Gamma distribution. Default value is 0.1.
#' @param b_sigma positive scalar input. It corresponds to the second parameter of a Inverse Gamma distribution. Default value is 0.1.
#' @param Prop_varA1 positive scalar input. It corresponds to the variance of the normal distribution of proposing A matrix terms. Default value is 0.1.
#' @param Prop_varA2 positive scalar input. It corresponds to the variance of the normal distribution of proposing A matrix terms. Default value is 5.
#' @param Prop_varB1 positive scalar input. It corresponds to the variance of the normal distribution of proposing B matrix terms. Default value is 0.1.
#' @param Prop_varB2 positive scalar input. It corresponds to the variance of the normal distribution of proposing A matrix terms. Default value is 5.
#' @param niter positive integer input. It corresponds to number of iterations the markov chain runs. Give niter as large as possible to get the convergence. Minimum value is 10,000. Default value is 10,000.
#'
#' @return A list which will hold A, B, LL:
#'
#' \item{A}{p * p matrix output of Gene-Gene interactions.}
#' \item{B}{p * k matrix output of Gene-DNA interactions.}
#' \item{LL}{niter * 1 vector output of log-likelihood values of the model for all the iterations.}
#'
#' @export
#'
#' @examples
#'
#' # -----------------------------------------------------------------
#' # Example 1:
#'
#' set.seed(123)
#'
#' # Number of datapoints
#' n = 10000
#'
#' # Number of Genes and number of DNAs
#' p = 5
#' k = 5
#'
#' # Initialize gene-gene interaction matrix
#' A = matrix(sample(c(-0.1, 0.1), p^2, replace = TRUE), p, p)
#'
#' # Diagonal entries of A matrix will always be 0
#' diag(A) = 0
#'
#' # Make the network sparse
#' A[sample(which(A!=0), length(which(A!=0))/2)] = 0
#'
#' # Initialize gene-DNA interaction matrix
#' B = diag(p)
#'
#'
#' # Indicator matrix for gene-DNA interaction
#' D = diag(p)
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
#' Y = matrix(0, nrow = n, ncol = p)
#'
#' for (i in 1:n) {
#'
#'     Y[i, ] = MASS::mvrnorm(n = 1, Mult_Mat %*% B %*% X[i, ], Variance)
#'
#' }
#'
#' # Calculate S_YY, S_YX and S_XX
#' S_YY = t(Y) %*% Y / n
#' S_YX = t(Y) %*% X / n
#' S_XX = t(X) %*% X / n
#'
#' # Apply RGM on the generated data
#' Out = RGM(S_YY, S_YX, S_XX, d = rep(1, p), n = n)
#'
#' # Get gene-gene interaction matrix, gene-DNA interaction matrix and log-likelihood
#' A_Est = Out$A_Est
#' B_Est = Out$B_Est
#' LL = Out$LL_Pst
#'
#' # Plot log-likelihood of the model at each iteration
#' plot(LL, type = 'l', xlab = "Number of Iterations", ylab = "Log-likelihood")
#'
#' @references
#' Ni, Y., Ji, Y., & Müller, P. (2018).
#' Reciprocal graphical models for integrative gene regulatory network analysis.
#' \emph{Bayesian Analysis},
#' \strong{13(4)}, 1095-1110.
#' \doi{10.1214/17-BA1087}.
RGM = function(X = NULL, Y = NULL, S_YY = NULL, S_YX = NULL, S_XX = NULL, Beta = NULL, Sigma_Hat = NULL, d, n, nIter = 10000, nBurnin = 2000, Thin = 1, prior = c("Threshold", "Spike and Slab"), a_tau = 0.01, b_tau = 0.01, a_rho = 0.5, b_rho = 0.5, nu_1 = 0.0001, a_eta = 0.01, b_eta = 0.01, a_psi = 0.5, b_psi = 0.5, nu_2 = 0.0001, a_sigma = 0.01, b_sigma = 0.01, Prop_VarA = 0.01, Prop_VarB = 0.01){

  # Check whether n is a positive integer
  if(!is.numeric(n) || n != round(n) || n <= 0){

    # Print an error message
    stop("Number of datapoints should be a positive integer.")

  }

  # Check whether input X and Y are given
  if(!is.null(X) && !is.null(Y)){

    # Check whether X and Y are both numeric matrices
    if(!is.numeric(X) || !is.numeric(Y) || !is.matrix(X) || !is.matrix(Y)){

      # Print an error message
      stop("X and Y should be numeric matrices.")

    }

    # Check whether number of rows of X and Y are n
    if(nrow(X) != n || nrow(Y) != n){

      # Print an error message
      stop("Number of rows of X and Y should be n.")

    }

    # Calculate number of nodes from Y matrix
    p = ncol(Y)

    # Calculate number of covariates from X matrix
    k = ncol(X)

    # Calculate S_YY, S_YX, S_XX
    S_YY = t(Y) %*% Y / n
    S_YX = t(Y) %*% X / n
    S_XX = t(X) %*% X / n

  } else if(!is.null(S_YY) && !is.null(S_YX) && !is.null(S_XX)){

    # Check whether S_YY ,S_YX and S_XX are numeric matrices
    if(!is.numeric(S_YY) || !is.numeric(S_YX) || !is.numeric(S_XX) || !is.matrix(S_YY) || !is.matrix(S_YX) || !is.matrix(S_XX)){

      # Print an error message
      stop("S_YY, S_YX and S_XX should be numeric matrices.")

    }

    # Calculate number of nodes from S_YY matrix
    p = ncol(S_YY)

    # Calculate number of covariates from S_XX matrix
    k = ncol(S_XX)

    # Check whether number of rows of S_YY is equal to p
    if(nrow(S_YY) != p){

      # Print an error message
      stop("S_YY should be a diagonal matrix.")

    }

    # Check whether number of rows of S_XX is equal to k
    if(nrow(S_XX) != k){

      # Print an error message
      stop("S_XX should be a diagonal matrix.")

    }

    # Check whether number of rows of S_YX is equal to p and number of columns of S_YX is equal to k
    if(nrow(S_YX) != p || ncol(S_YX) != k){

      # Print an error message
      stop("Number of rows of S_YX should be equal to number of rows of S_YY and number of columns of S_YX should be equal to number of columns of S_XX.")

    }



  } else if(!is.null(S_XX) && !is.null(Beta) && !is.null(Sigma_Hat)){

    # Check whether S_XX, Beta and Sigma_Hat are numeric matrices
    if(!is.numeric(S_XX) || !is.numeric(Beta) || !is.numeric(Sigma_Hat) || !is.matrix(S_XX) || !is.matrix(Beta) || !is.matrix(Sigma_Hat)){

      # Print an error message
      stop("S_XX, Beta and SIgma_Hat should be numeric matrices.")

    }

    # Calculate number of nodes from Beta matrix
    p = nrow(Beta)

    # Calculate number of covariates from S_XX matrix
    k = ncol(S_XX)

    # Check whether number of rows of S_XX is equal to k
    if(nrow(S_XX) != k){

      # Print an error message
      stop("S_XX should be a diagonal matrix.")

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
    Cor_X = cov2cor(S_XX)
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



  # Check whether d is a vector of non-negative integers of length p and sum of entries of d is equal to k
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


  # Check whether the inverse gamma parameters are positive or not
  if(!is.numeric(a_tau) || a_tau < 0 || !is.numeric(b_tau) || b_tau < 0 || !is.numeric(a_rho) || a_rho < 0 || !is.numeric(b_rho) || b_rho < 0 || !is.numeric(a_eta) || a_eta < 0 || !is.numeric(b_eta) || b_eta < 0 || !is.numeric(a_psi) || a_psi < 0 || !is.numeric(b_psi) || b_psi < 0 || !is.numeric(a_sigma) || a_sigma < 0 || !is.numeric(b_sigma) || b_sigma < 0){

    # Print an error message
    stop("All the inverse gamma parameters should be positive")

  }


  # Check whether the variance terms are positive or not
  if(!is.numeric(nu_1) || nu_1 <= 0 || !is.numeric(nu_2) || nu_2 <= 0 || !is.numeric(Prop_VarA) || Prop_VarA <= 0 || !is.numeric(Prop_VarB) || Prop_VarB <= 0){

    # Print an error message
    stop("All the variance terms should be positive")

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


  # Apply RGM for Threshold prior and Spike and Slab prior
  if ("Threshold" %in% prior){

    # Run the algorithm for Threshold prior
    Output = RGM_Threshold(S_YY, S_YX, S_XX, D, n, nIter = nIter, nBurnin = nBurnin, Thin = Thin, nu_1 = nu_1, nu_2 = nu_2,
                           a_sigma = a_sigma, b_sigma = b_sigma, Prop_VarA = Prop_VarA, Prop_VarB = Prop_VarB)

    # Calculate estimates from the output
    A_Est = apply(Output$A_Pst, MARGIN = c(1, 2), FUN = mean)
    B_Est = apply(Output$B_Pst, MARGIN = c(1, 2), FUN = mean)
    A0_Est = apply(Output$A0_Pst, MARGIN = c(1, 2), FUN = mean)
    B0_Est = apply(Output$B0_Pst, MARGIN = c(1, 2), FUN = mean)
    Gamma_Est = apply(Output$Gamma_Pst, MARGIN = c(1, 2), FUN = mean)
    Tau_Est = apply(Output$Tau_Pst, MARGIN = c(1, 2), FUN = mean)
    Phi_Est = apply(Output$Phi_Pst, MARGIN = c(1, 2), FUN = mean)
    Eta_Est = apply(Output$Eta_Pst, MARGIN = c(1, 2), FUN = mean)
    tA_Est = mean(Output$tA_Pst)
    tB_Est = mean(Output$tB_Pst)
    Sigma_Est = apply(Output$Sigma_Pst, MARGIN = c(1, 2), FUN = mean)
    zA_Est = (Gamma_Est > 0.5) * 1
    zB_Est = (Phi_Est > 0.5) * 1

    # Return outputs
    return(list(A_Est = A_Est, B_Est = B_Est, zA_Est = zA_Est, zB_Est = zB_Est,
                A0_Est = A0_Est, B0_Est = B0_Est, Gamma_Est = Gamma_Est, Tau_Est = Tau_Est,
                Phi_Est = Phi_Est, Eta_Est = Eta_Est, tA_Est = tA_Est, tB_Est = tB_Est,
                Sigma_Est = Sigma_Est,
                A_Pst = Output$A_Pst, B_Pst = Output$B_Pst,
                A0_Pst = Output$A0_Pst, B0_Pst = Output$B0_Pst, Gamma_Pst = Output$Gamma_Pst, Tau_Pst = Output$Tau_Pst,
                Phi_Pst = Output$Phi_Pst, Eta_Pst = Output$Eta_Pst, tA_Pst = Output$tA_Pst, tB_Pst = Output$tB_Pst,
                Sigma_Pst = Output$Sigma_Pst,
                AccptA = Output$AccptA, AccptB = Output$AccptB, Accpt_tA = Output$Accpt_tA, Accpt_tB = Output$Accpt_tB,
                LL_Pst = Output$LL_Pst))





  } else if ("Spike and Slab" %in% prior){

    # Run the algorithm for Threshold prior
    Output = RGM_SpikeSlab(S_YY, S_YX, S_XX, D, n, nIter = nIter, nBurnin = nBurnin, Thin = Thin, a_tau = a_tau, b_tau = b_tau,
                           a_rho = a_rho, b_rho = b_rho, nu_1 = nu_1, a_eta = a_eta, b_eta = b_eta, a_psi = a_psi, b_psi = b_psi,
                           nu_2 = nu_2, a_sigma = a_sigma, b_sigma = b_sigma, Prop_VarA = Prop_VarA, Prop_VarB = Prop_VarB)



    # Calculate estimates from the output
    A_Est = apply(Output$A_Pst, MARGIN = c(1, 2), FUN = mean)
    B_Est = apply(Output$B_Pst, MARGIN = c(1, 2), FUN = mean)
    Gamma_Est = apply(Output$Gamma_Pst, MARGIN = c(1, 2), FUN = mean)
    Rho_Est = apply(Output$Rho_Pst, MARGIN = c(1, 2), FUN = mean)
    Tau_Est = apply(Output$Tau_Pst, MARGIN = c(1, 2), FUN = mean)
    Phi_Est = apply(Output$Phi_Pst, MARGIN = c(1, 2), FUN = mean)
    Eta_Est = apply(Output$Eta_Pst, MARGIN = c(1, 2), FUN = mean)
    Psi_Est = apply(Output$Psi_Pst, MARGIN = c(1, 2), FUN = mean)
    Sigma_Est = apply(Output$Sigma_Pst, MARGIN = c(1, 2), FUN = mean)
    zA_Est = (Gamma_Est > 0.5) * 1
    zB_Est = (Phi_Est > 0.5) * 1

    # Return outputs
    return(list(A_Est = A_Est, B_Est = B_Est, zA_Est = zA_Est, zB_Est = zB_Est,
                Gamma_Est = Gamma_Est, Tau_Est = Tau_Est, Rho_Est = Rho_Est,
                Phi_Est = Phi_Est, Eta_Est = Eta_Est, Psi_Est = Psi_Est,
                Sigma_Est = Sigma_Est,
                A_Pst = Output$A_Pst, B_Pst = Output$B_Pst,
                Gamma_Pst = Output$Gamma_Pst, Tau_Pst = Output$Tau_Pst, Rho_Pst = Output$Rho_Pst,
                Phi_Pst = Output$Phi_Pst, Eta_Pst = Output$Eta_Pst, Psi_Pst = Output$Psi_Pst,
                Sigma_Pst = Output$Sigma_Pst,
                AccptA = Output$AccptA, AccptB = Output$AccptB,
                LL_Pst = Output$LL_Pst))



  } else{

    # Print an error message
    stop("Please specify a prior among Threshold prior and Spike and Slab prior.")


  }




}



