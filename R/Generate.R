# Generate rho
Generate_Rho = function(Gamma, p, a_rho, b_rho){

  # Calculate sum of all entries of Gamma matrix
  Gamma_sum = sum(Gamma)

  # Generate Rho from Beta distribution
  Rho = rbeta(1, Gamma_sum + a_rho, p * (p - 1) - Gamma_sum + b_rho)

  # Return Rho
  return(Rho)

}


# Generate psi
Generate_Psi = function(Phi, d, a_psi, b_psi){

  # Calculate sum of all entries of Phi matrix
  Phi_sum = sum(Phi)

  # Generate psi from Beta distribution
  Psi = rbeta(1, Phi_sum + a_psi, d - Phi_sum + b_psi)

  # Return Rho
  return(Psi)

}


# Generate eta
Generate_Eta = function(b, phi, a_eta, b_eta, nu_2){

  # Generate eta based on Inverse Gamma distribution
  eta = phi / rgamma(1, a_eta + 1/2, b^2/2 + b_eta) + (1 - phi) / rgamma(1, a_eta + 1/2, b^2/(2 * nu_2) + b_eta)

  # Return eta
  return(eta)

}


# Generate tau
Generate_Tau = function(a, gamma, a_tau, b_tau, nu_1){

  # Generate tau based on Inverse Gamma distribution
  tau = gamma / rgamma(1, a_tau + 1/2, a^2/2 + b_tau) + (1 - gamma) / rgamma(1, a_tau + 1/2, a^2/(2 * nu_1) + b_tau)

  # Return tau
  return(tau)

}


# Generate phi
Generate_Phi = function(b, psi, eta, nu_2){

  # Calculate C1 and C2
  C1 = psi / sqrt(eta) * exp(- b^2 / (2 * eta))

  C2 = (1 - psi) / sqrt(nu_2 * eta) * exp(- b^2 / (2 * nu_2 * eta))

  # Generate phi based on bernoulli distribution
  phi = rbinom(1, 1, C1/(C1+C2))

  # Return phi
  return(phi)

}



# Generate gamma
Generate_Gamma = function(a, rho, tau, nu_1){

  # Calculate D1 and D2
  D1 = rho / sqrt(tau) * exp(- a^2 / (2 * tau))

  D2 = (1 - rho) / sqrt(nu_1 * tau) * exp(- a^2 / (2 * nu_1 * tau))

  # Generate gamma based on bernoulli distribution
  gamma = rbinom(1, 1, D1/(D1+D2))

  # Return gamma
  return(gamma)

}
