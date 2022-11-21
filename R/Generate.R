# Generate rho
Generate_Rho = function(Gamma, p, a_rho, b_rho){

  # Calculate sum of all entries of Gamma matrix
  Gamma_sum = sum(Gamma)

  # Generate Rho from Beta distribution
  Rho = stats::rbeta(1, Gamma_sum + a_rho, p * (p - 1) - Gamma_sum + b_rho)

  # Return Rho
  return(Rho)

}


# Generate psi
Generate_Psi = function(Phi, d, a_psi, b_psi){

  # Calculate sum of all entries of Phi matrix
  Phi_sum = sum(Phi)

  # Generate psi from Beta distribution
  Psi = stats::rbeta(1, Phi_sum + a_psi, d - Phi_sum + b_psi)

  # Return Rho
  return(Psi)

}


# Generate eta
Generate_Eta = function(b, phi, a_eta, b_eta, nu_2){

  # Generate eta based on Inverse Gamma distribution
  eta = phi / stats::rgamma(1, a_eta + 1/2, b^2/2 + b_eta) + (1 - phi) / stats::rgamma(1, a_eta + 1/2, b^2/(2 * nu_2) + b_eta)

  # Return eta
  return(eta)

}


# Generate tau
Generate_Tau = function(a, gamma, a_tau, b_tau, nu_1){

  # Generate tau based on Inverse Gamma distribution
  tau = gamma / stats::rgamma(1, a_tau + 1/2, a^2/2 + b_tau) + (1 - gamma) / stats::rgamma(1, a_tau + 1/2, a^2/(2 * nu_1) + b_tau)

  # Return tau
  return(tau)

}


# Generate phi
Generate_Phi = function(b, psi, eta, nu_2){

  # Calculate C1 and C2
  C1 = psi / sqrt(eta) * exp(- b^2 / (2 * eta))

  C2 = (1 - psi) / sqrt(nu_2 * eta) * exp(- b^2 / (2 * nu_2 * eta))

  # Generate phi based on bernoulli distribution
  phi = stats::rbinom(1, 1, C1/(C1+C2))

  # Return phi
  return(phi)

}



# Generate gamma
Generate_Gamma = function(a, rho, tau, nu_1){

  # Calculate D1 and D2
  D1 = rho / sqrt(tau) * exp(- a^2 / (2 * tau))

  D2 = (1 - rho) / sqrt(nu_1 * tau) * exp(- a^2 / (2 * nu_1 * tau))

  # Generate gamma based on bernoulli distribution
  gamma = stats::rbinom(1, 1, D1/(D1+D2))

  # Return gamma
  return(gamma)

}



# Generate sigma
Generate_Sigma = function(n, z_sum, a_sigma, b_sigma){

  # Generate sigma based on Inverse gamma distribution
  sigma = 1 / stats::rgamma(1, n/2 + a_sigma, z_sum/2 + b_sigma)

  # Return sigma
  return(sigma)

}




# Generate an entry of B matrix
Generate_B = function(xz, b_vec, b, x_mat, x_vec, sigma, eta, phi, nu2){

  # Calculate numerator of mean
  mean = xz - sum(tcrossprod(b_vec, x_mat) * x_vec) / sigma + sum(x_vec^2) * b / sigma

  # Calculate variance
  variance = phi / (sum(x_vec) / sigma + 1 / eta) + (1 - phi) / (sum(x_vec) / sigma + 1 / (eta * nu2))

  # Generate b
  b = stats::rnorm(1, mean = mean * variance, variance)

  # Return beta
  return(b)

}





# Target function for A and gamma
Target_Agamma = function(X, Y, A, a, N, Sigma_Inv, p, B, gamma, tau, rho, nu_1){

  # Calculate (I_p - A)^{-1} and (I_p - A)^{-1} * B
  Mult_Mat = diag(p) - A
  Inv_MatB = solve(Mult_Mat) %*% B


  # Calculate mean matrix and variance matrix
  Mean_mat = tcrossprod(Inv_MatB, X)
  Var_mat_inv = crossprod(Mult_Mat, Sigma_Inv %*% Mult_Mat)

  # Initiate Sum
  Sum = 0

  for (i in 1:N) {

    # Calculate difference
    Diff = Y[i, ] - t(Mean_mat[, i])

    # Calculate Sum
    Sum = Sum + tcrossprod(Diff %*% Var_mat_inv, Diff)

  }

  # Target = ((det(Var_mat_inv)^(N/2)) * exp(-1/2 * Sum) * (gamma * exp(-a^2/(2 * tau)) / sqrt(tau) + (1 - gamma) * exp(- a^2/(2 * nu_1 * tau)) / sqrt(nu_1 * tau))) * (rho ^ gamma * (1-rho) ^ (1-gamma))
  Target = N/2 * determinant(Var_mat_inv, logarithm = TRUE)$modulus -1/2 * Sum + gamma * (-a^2/(2 * tau) - log(sqrt(tau))) + (1 - gamma) * (-a^2/(2 * nu_1 * tau) - log(sqrt(nu_1 * tau))) + gamma * log(rho) + (1 - gamma) * log(1 - rho)


  # Return Target
  return(Target)

}






# Generate an entry of A matrix
Generate_Agamma = function(X, Y, A, i, j, Sigma_Inv, N, p, B, gamma, tau, rho, nu_1, prop_var1){

  # Value to update
  a = A[i, j]

  # Proposed value
  a_new = stats::rnorm(1, a, prop_var1)

  # New A matrix with proposed a value
  A_new = A
  A_new[i, j] = a_new

  # Calculate r
  # r = Target_Agamma(X, Y, A_new, a_new, N, Sigma_Inv, p, B, 1-gamma, tau, rho, nu_1)  / Target_Agamma(X, Y, A, a, N, Sigma_Inv, p, B, gamma, tau, rho, nu_1)
  r = Target_Agamma(X, Y, A_new, a_new, N, Sigma_Inv, p, B, 1-gamma, tau, rho, nu_1) - Target_Agamma(X, Y, A, a, N, Sigma_Inv, p, B, gamma, tau, rho, nu_1)



  # Generate uniform u
  u = stats::runif(1, 0, 1)

  if(!is.na(r)){

    # Check whether r is big or not
    # min(1, r) >= u
    if(min(0, r) >= log(u)){

      a = a_new

      gamma = 1 - gamma

    }

  }

  return(list(a = a, gamma = gamma))

}





# Target function for A
Target_A = function(X, Y, A, a, N, Sigma_Inv, p, B, gamma, tau, nu_1){

  # Calculate (I_p - A)^{-1} and (I_p - A)^{-1} * B
  Mult_Mat = diag(p) - A

  Inv_Mat = solve(Mult_Mat)
  Inv_MatB = Inv_Mat %*% B

  # Calculate mean matrix and variance matrix
  Mean_mat = tcrossprod(Inv_MatB, X)
  Var_mat_inv = crossprod(Mult_Mat, Sigma_Inv %*% Mult_Mat)

  # Initiate Sum
  Sum = 0

  for (i in 1:N) {

    # Calculate difference
    Diff = Y[i, ] - t(Mean_mat[, i])

    # Calculate Sum
    Sum = Sum + tcrossprod(Diff %*% Var_mat_inv, Diff)

  }

  Target = (det(Var_mat_inv)^(N/2)) * exp(-1/2 * Sum) * (gamma * exp(-a^2/(2 * tau)) / sqrt(tau) + (1 - gamma) * exp(- a^2/(2 * nu_1 * tau)) / sqrt(nu_1 * tau))


  # Return Target
  return(Target)

}





# Generate an entry of A matrix
Generate_A = function(X, Y, A, i, j, Sigma_Inv, N, p, B, gamma, tau, nu_1, prop_var1){

  # Value to update
  a = A[i, j]

  # Proposed value
  a_new = stats::rnorm(1, a, prop_var1)

  # New A matrix with proposed a value
  A_new = A
  A_new[i, j] = a_new

  # Calculate r
  r = Target_A(X, Y, A_new, a_new, Sigma_Inv, N, p, B, gamma, tau, nu_1)  / Target_A(X, Y, A, a, Sigma_Inv, N, p, B, gamma, tau, nu_1)

  # Generate uniform u
  u = stats::runif(1, 0, 1)

  if(!is.na(r)){

    # Check whether r is big or not
    if(min(1, r) >= u){

      a = a_new

    }

  }

  return(a)

}





















