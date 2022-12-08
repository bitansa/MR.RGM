RGM = function(X, Y, A0 = NULL, B0 = NULL, D = NULL, a_tau = 0.1, b_tau = 0.1, a_rho = 0.5, b_rho = 0.5, nu_1 = 0.0001, a_eta = 0.1, b_eta = 0.1, a_psi = 0.5, b_psi = 0.5, nu_2 = 0.0001, a_sigma = 0.1, b_sigma = 0.1, Prop_varA = 0.1, Prop_VarB = 0.001, niter = 10000){

  # Calculate number of datapoints and number of nodes from Y matrix
  n = nrow(Y)
  p = ncol(Y)

  # Calculate number of columns of X matrix i.e. number of covariates
  k = ncol(X)

  # Check whether number of rows of X and Y are same
  if(nrow(X) != n){

    # Print an error message
    stop("Number of datapoints for both node values and covariate values should be same")

  }

  # Check all the elements of X is numeric
  if(!is.numeric(X)){

    # Print an error message
    stop("All the entries of node output should be numeric")


  }

  # Check all the elements of Y is numeric
  if(!is.numeric(Y)){

    # Print an error message
    stop("All the entries of covariate output should be numeric")


  }

  # Check whether D is given or not
  if(is.null(D)){

    # Create a D matrix with all entries 1
    D = matrix(1, nrow = p, ncol = k)

  } else {

    # Check the dimensions of D matrix
    if(nrow(D) != p){

      # Print an error message
      stop("Number of rows of the indicator matrix should be equal to the number of nodes")

    } else if(ncol(D) != k){

      # Print an error message
      stop("Number of columns of the indicator matrix should be equal to the number of covariates")

    }

   # Check whether all the entries of D matrix is either 0 or 1
   if(!is.integer(D) || max(D) > 1 || min(D) < 0){

     # Print an error message
     stop("All the entries of the indicator matrix should be either 0 or 1")

    }

  }

  # Store number of non-zero entries of D
  d = length(which(D != 0))


  # Check whether the inverse gamma parameters are positive or not
  if(!is.numeric(a_tau) || a_tau < 0 || !is.numeric(b_tau) || b_tau < 0 || !is.numeric(a_rho) || a_rho < 0 || !is.numeric(b_rho) || b_rho < 0 || !is.numeric(a_eta) || a_eta < 0 || !is.numeric(b_eta) || b_eta < 0 || !is.numeric(a_psi) || a_psi < 0 || !is.numeric(b_psi) || b_psi < 0 || !is.numeric(a_sigma) || a_sigma < 0 || !is.numeric(b_sigma) || b_sigma < 0){

    # Print an error message
    stop("All the inverse gamma parameters should be positive")

  }


  # Check whether the variance terms are positive or not
  if(!is.numeric(nu_1) || nu_1 <= 0 || !is.numeric(nu_2) || nu_2 <= 0 || !is.numeric(Prop_varA) || Prop_varA <= 0 || !is.numeric(Prop_VarB) || Prop_VarB <= 0){

    # Print an error message
    stop("All the variance terms should be positive")

  }


  # Check whether niter term is large enough or not
  if(!is.integer(niter) || niter < 10000){

    # Print an error message
    stop("Number of iterations should be large i.e. at least 10000")

  }


  # Initialize Sigma Inverse vector i.e. inverse of each Sigma_ii
  Sigma_Inv = stats::rgamma(p, a_sigma, b_sigma)

  # Initialize Bernoulli parameters Rho, Psi
  Rho = stats::rbeta(1, a_rho, b_rho)
  Psi = stats::rbeta(1, a_psi, b_psi)

  # Initialize Gamma, Tau, Phi and Eta matrices
  Gamma = matrix(stats::rbinom(p^2, 1, Rho), nrow = p, ncol = p)
  Phi = matrix(stats::rbinom(p * k, 1, Psi), nrow = p, ncol = k)
  Tau = matrix(1 / stats::rgamma(p^2, a_tau, b_tau), nrow = p, ncol = p)
  Eta = matrix(1 / stats::rgamma(p * k, a_eta, b_eta), nrow = p, ncol = k)

  # Make the diagonals of Gamma and Tau matrix to be 0
  diag(Gamma) = 0
  diag(Tau) = 0

  # Make Phi[i, j] = 0 and Eta[i, j] = 0 if D[i, j] = 0
  Phi = Phi * D
  Eta = Eta * D

  # Initialize A and B matrix
  A = matrix(0, nrow = p, ncol = p)
  B = matrix(0, nrow = p, ncol = k)

  # Check whether A0 matrix is null otherwise initialize it
  if(!is.null(A0)){

    # Check the dimensions of A0
    if(nrow(A0) != p || ncol(A0) != p){

      # Print an error message
      stop("Number of rows and columns of A0 matrix should be equal to number of nodes")

    }

    # Check all the entries of A0 is numeric and diagonal is 0
    if(!is.numeric(A0) || sum(diag(A0) != 0) != 0){

      # Print an error message
      stop("A0 should be a numeric matrix with diagonal entries all 0")

    }

    # Print a warning message
    warning("Random initialization of A matrix may take large enough number of iterations to converge")

    # Update A by A0
    A = A0

  } else {

    # Initialize A matrix
    for (i in 1:p) {

      for (j in 1:p) {

        if(j != i) {

          # Generate A[i, j]
          A[i, j] = Gamma[i, j] * stats::rnorm(1, 0, Tau[i, j]) + (1 - Gamma[i, j]) * stats::rnorm(1, 0, Tau[i, j] * nu_1)

        }

      }

    }

  }



  # Check whether B0 matrix is null otherwise initialize it
  if(!is.null(B0)){

    # Check the dimensions of B0
    if(nrow(B0) != p || ncol(B0) != k){

      # Print an error message
      stop("Number of rows of B0 matrix should be equal to number of nodes and Number of columns of B0 matrix should be equal to number of covariates")

    }

    # Check all the entries of B0 is numeric
    if(!is.numeric(B0)){

      # Print an error message
      stop("B0 should be a numeric matrix")

    }

    # Print a warning message
    warning("Random initialization of B matrix may take large enough number of iterations to converge")

    # Update B by B0
    B = B0

    # Multiply B pointwise with D to make B[i, j] = 0 if D[i, j] = 0
    B = B * D

  } else {

    # Initialize B matrix
    for (i in 1:p) {

      for (j in 1:k) {

        if(D[i, j] != 0) {

          # Generate B[i, j]
          B[i, j] = Phi[i, j] * stats::rnorm(1, 0, Eta[i, j]) + (1 - Phi[i, j]) * stats::rnorm(1, 0, Eta[i, j] * nu_2)

        }

      }

    }

  }


  # Initialize A_Update, B_Update, Gamma_update, Phi_Update matrices and Log-likelihood vector
  A_Update = matrix(0, nrow = p * p, ncol = niter)
  B_Update = matrix(0, nrow = p * k, ncol = niter)
  Gamma_Update = matrix(0, nrow = p * p, ncol = niter)
  Phi_Update = matrix(0, nrow = p * k, ncol = niter)
  LogLikelihood = rep(0, niter)


  # Run a  loop to update all the parameters
  for (i in 1:niter) {

    # Update Rho
    Rho = Generate_Rho(Gamma = Gamma, p, a_rho = a_rho, b_rho = b_rho)

    # Update Psi
    Psi = Generate_Psi(Phi = Phi, d, a_psi = a_psi, b_psi = b_psi)


    # Update Tau based on corresponding a and gamma and then Update a and gamma based on the corresponding tau
    for (j in 1:p) {

      for (l in 1:p) {

        # Don't update the diagonal entries
        if(l != j){

          # Update Tau
          Tau[j, l] = Generate_Tau(A[j, l], Gamma[j, l], a_tau = a_tau, b_tau = b_tau, nu_1 = nu_1)

          # Generate both a and gamma
          Out = Generate_Agamma(X, Y, A, j, l, Sigma_Inv = Sigma_Inv, n, p, B, gamma = Gamma[j, l], tau = Tau[j, l], rho = Rho, nu_1 = nu_1, prop_var1 = Prop_varA)

          # Update A[j, l] and Gamma[j, l]
          A[j, l] = Out$a

          Gamma[j, l] = Out$gamma


          # Update A_Update and Gamma_Update
          A_Update[p * (l - 1) + j , i] = Out$a

          Gamma_Update[p * (l - 1) + j , i] = Out$gamma


        }

      }

    }

    # Calculate (I_p - A) %*% t(Y)
    MultMat_Y = tcrossprod((diag(p) - A), Y)


    # Update Eta based on corresponding b and phi and then Update b and phi based on the corresponding eta
    for (j in 1:p) {

      for (l in 1:K) {

        # Don't update if the corresponding D entry is 0
        if(D[j, l] != 0){

          # Update Eta
          Eta[j, l] = Generate_Eta(B[j, l], Phi[j, l], a_eta = a_eta, b_eta = b_eta, nu_2 = nu_2)

          # Generate both b and phi
          Out = Generate_Bphi(X, Y, B, j, l, Sigma_Inv = Sigma_Inv,  MultMat_Y , phi = Phi[j, l], eta = Eta[j, l], psi = Psi, nu_2 = nu_2, prop_var2 = Prop_VarB)

          # Update B[j, l] and Phi[j, l]
          B[j, l] = Out$b

          Phi[j, l] = Out$phi

          # Update B_Update and Phi_Update
          B_Update[p * (l - 1) + j , i] = Out$b

          Phi_Update[p * (l - 1) + j , i] = Out$phi


        }

      }

    }


    # Calculate I_p - A
    Mult_Mat = diag(p) - A

    # Calculate Diff matrix as (I_p - A) %*% t(Y) - B %*% t(X)
    Diff = tcrossprod(Mult_Mat, Y) - tcrossprod(B, X)

    # Update Sigma_Inv vector
    for (j in 1:p) {

      Sigma_Inv[j] = 1 / Generate_Sigma(n, sum(Diff[j, ]^2), a_sigma, b_sigma)

    }


    # Update Log-likelihood
    LogLikelihood[i] = LL(A, B, X, Y, Sigma_Inv = Sigma_Inv, p, n)


  }

  # Calculate Gamma and Phi as the rowmeans of Gamma_Update and Phi_Update
  Gamma = matrix(rowMeans(Gamma_Update[, 3000:niter]), nrow = p, ncol = p)
  Phi = matrix(rowMeans(Phi_Update[, 3000:niter]), nrow = p, ncol = k)

  # Calculate A and B as the rowmeans of A_Update and B_Update
  A = matrix(rowMeans(A_Update[, 3000:niter]), nrow = p, ncol = p)
  B = matrix(rowMeans(B_Update[, 3000:niter]), nrow = p, ncol = k)

  # Take A[i, j] = A[i, j] if Gamma[i, j] > 1/2 else take 0
  A = A * (Gamma >= 1/2)

  # Take B[i, j] = B[i, j] if Phi[i, j] > 1/2 else take 0
  B = B * (Phi >= 1/2)


  # Return A, B and log-likelihood
  return(list(A = A, B = B, LL = LogLikelihood))

}
