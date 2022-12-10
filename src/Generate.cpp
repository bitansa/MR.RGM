#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// Generate Rho
// [[Rcpp::export]]
double Generate_Rho_c(const arma::mat& Gamma, double p, double a_rho, double b_rho){

  // Calculate sum of all entries of Gamma matrix
  double  Gamma_sum = accu(Gamma);

  // Generate Rho from Beta distribution
  NumericVector RHO = Rcpp::rbeta(1, Gamma_sum + a_rho, p * (p - 1) - Gamma_sum + b_rho);

  // Store Rho
  double Rho = RHO(0);

  // Return Rho
  return Rho;

}


// Generate psi
// [[Rcpp::export]]
double Generate_Psi_c(const arma::mat& Phi, double d, double a_psi, double b_psi){

  // Calculate sum of all entries of Phi matrix
  double Phi_sum = accu(Phi);

  // Generate Psi from Beta distribution
  NumericVector PSI = Rcpp::rbeta(1, Phi_sum + a_psi, d - Phi_sum + b_psi);

  // Store Psi
  double Psi = PSI(0);

  // Return Psi
  return(Psi);

}



// Generate eta
// [[Rcpp::export]]
double Generate_Eta_c(double b, double phi, double a_eta, double b_eta, double nu_2){

  NumericVector eta;

  if(phi == 1){

    eta = 1 / Rcpp::rgamma(1, a_eta + 1 / double(2), 1 / (b * b / 2 + b_eta));

  } else {

    eta = 1 / Rcpp::rgamma(1, a_eta + 1 / double(2), 1 / (b * b / (2 * nu_2) + b_eta));

  }

  // Store eta
  double Eta = eta(0);


  // Return eta
  return(Eta);

}



// Generate tau
// [[Rcpp::export]]
double Generate_Tau_c(double a, double gamma, double a_tau, double b_tau, double nu_1){

  NumericVector tau;

  if(gamma == 1){

    tau = 1 / Rcpp::rgamma(1, a_tau + 1 / double(2), 1 / (a * a / 2 + b_tau));

  } else {

    tau = 1 / Rcpp::rgamma(1, a_tau + 1 / double(2), 1 / (a * a / (2 * nu_1) + b_tau));

  }

  // Store Tau
  double Tau = tau(0);


  // Return tau
  return(Tau);

}

