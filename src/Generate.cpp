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
