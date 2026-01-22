#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// For Inverse Wishart distribution
#include <RcppDist.h>
// [[Rcpp::depends(RcppDist)]]


// Inverse Wishart sampler adapted from:
// Frank DiTraglia, econ722 repository
// https://github.com/fditraglia/econ722
// License: GPL-2.0



//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


// ------------------------------
// Helpers: file-local only
// ------------------------------
namespace {


// Sample Rho
double Sample_Rho(double Gamma, double a_rho, double b_rho) {

  // Sample Rho from beta distribution
  double Rho = Rcpp::rbeta(1, Gamma + a_rho, 1 - Gamma + b_rho)(0);

  // Return Rho
  return Rho;

}




// Sample Tau
double Sample_Tau(double a, double gamma, double tau, double nu_1) {

  // Sample Epsilon based on old tau
  double Epsilon = 1 / Rcpp::rgamma(1, 1, 1 / (1 + 1 / tau))(0);

  // Initialize Tau
  double Tau;

  // Check whether gamma is 0 or 1
  if (gamma == 1) {

    // Sample Tau based on a and Epsilon
    Tau = 1 / Rcpp::rgamma(1, 1, 1 / (a * a / 2 + 1 / Epsilon))(0);

  } else {

    // Sample Tau based on a, nu_1 and Epsilon
    Tau = 1 / Rcpp::rgamma(1, 1, 1 / (a * a / (2 * nu_1) + 1 / Epsilon))(0);

  }

  // Return Tau
  return Tau;

}



// Sample Gamma
double Sample_Gamma(double a, double tau, double rho, double nu_1) {

  // Calculate acceptance probability
  double p = exp(-0.5 * (a * a / tau)) * rho / (exp(-0.5 * (a * a / tau)) * rho
                                                  + 1 / sqrt(nu_1) * exp(-0.5 * (a * a / (nu_1 * tau))) * (1 - rho));

  // Check if p is NaN and assign 0.5 for it
  if (std::isnan(p)) {
    p = 0.5;  // Set p to 0.5 if NaN is encountered
  }

  // Sample Gamma from binomial distribution
  double Gamma = Rcpp::rbinom(1, 1, p)(0);

  // Return Gamma
  return Gamma;

}







// Sample from inverse wishart distribution
// From github fditraglia/econ722/RcppArmadillo/InverseWishartSampler/InverseWishart.cpp
arma::mat rinvwish(double v, arma::mat S){

  RNGScope scope;

  // Dimension of matrix S
  int p = S.n_rows;

  arma::mat L = chol(inv_sympd(S), "lower");
  arma::mat sim(p, p, arma::fill::zeros);
  arma::mat A(p,p, arma::fill::zeros);

  for(int i = 0; i < p; i++){
    int df = v - (i + 1) + 1; //zero-indexing
    A(i,i) = sqrt(R::rchisq(df));
  }

  for(int row = 1; row < p; row++){
    for(int col = 0; col < row; col++){
      A(row, col) = R::rnorm(0,1);
    }
  }

  arma::mat LA_inv = inv(trimatl(trimatl(L) * trimatl(A)));
  sim = LA_inv.t() * LA_inv;
  return sim;

}





// Calculate log-likelihood for the model when only Y is there and Sigma is full
double LL_Full_Star(const arma::mat& A, const arma::mat& S_YY, const arma::mat& Sigma_Inv, double p, double N) {

  // Calculate (I_p - A)
  const arma::mat& Mult_Mat = arma::eye(p, p) - A;

  // Calculate Sum
  double Sum = N * arma::trace(S_YY * Mult_Mat.t() * Sigma_Inv * Mult_Mat);

  // Calculate log-likelihood
  double LL = N * real(arma::log_det(Mult_Mat)) + N / 2 * real(arma::log_det(Sigma_Inv)) - Sum / 2 - N / 2 * log(2 * arma::datum::pi);

  // Return log-likelihood
  return LL;

}




// Calculate target value for a particular A for the model when only Y is there
double Target_A_Star(double a, double N, double gamma, double tau, double nu_1, double Trace3, double Trace4, double Trace5, double logdet) {

  // Calculate Sum term inside exponential in likelihood
  double Sum = Trace3 + Trace4 + Trace5;

  // Calculate Target value
  double Target = N * logdet - Sum / 2 - gamma * (a * a / (2 * tau)) - (1 - gamma) * (0.5 * log(nu_1) + a * a / (2 * nu_1 * tau));

  // Return Target
  return Target;

}



// Sample a particular entry of matrix A for the model when only Y is there
Rcpp::List Sample_A_Full_Star(const arma::mat& S_YY, const arma::mat& A, const arma::mat& A_Pseudo, double i, double j, const arma::mat& Sigma_Inv,
                              double N, double p, double gamma, double tau, double nu_1, double prop_var1, double tA, double Trace3,
                              double Trace4, double Trace5, arma::mat InvMat, double logdet) {

  // Value to update
  double a = A_Pseudo(i, j);

  // Proposed value
  double a_new = Rcpp::rnorm(1, a, sqrt(prop_var1))(0);

  // Create a copy of matrix A
  arma::mat A_new = A;

  // Modify the copy with the proposed a value
  A_new(i, j) = (fabs(a_new) > tA) * a_new;

  // Modify logdet
  double logdet_new = logdet + log(fabs(1 + ((fabs(a) > tA) * a - (fabs(a_new) > tA) * a_new) * InvMat(j, i)));


  // Calculate new trace values
  double Trace3_New = Trace3 - N * ((fabs(a_new) > tA) * a_new - (fabs(a) > tA) * a) * arma::trace(Sigma_Inv.row(i) * S_YY.col(j));
  double Trace4_New = Trace4 - N * ((fabs(a_new) > tA) * a_new - (fabs(a) > tA) * a) * arma::trace(S_YY.row(j) * Sigma_Inv.col(i));
  double Trace5_New = Trace5 + N * ((fabs(a_new) > tA) * a_new - (fabs(a) > tA) * a) * arma::trace(S_YY.row(j) * A.t() * Sigma_Inv.col(i)
                                                                                                     + S_YY.col(j).t() * A_new.t() * Sigma_Inv.row(i).t());


  // Calculate target values with a and a_new
  double Target1 = Target_A_Star(a_new, N, gamma, tau, nu_1, Trace3_New, Trace4_New, Trace5_New, logdet_new);
  double Target2 = Target_A_Star(a, N, gamma, tau, nu_1, Trace3, Trace4, Trace5, logdet);

  // Calculate r i.e. the differnce between two target values
  double r = Target1 - Target2;

  // Sample u from Uniform(0, 1)
  double u = Rcpp::runif(1, 0, 1)(0);

  // Compare u and r
  if (r >= log(u)) {

    // Update a, trace values, logdet and InvMat
    a = a_new;

    Trace3 = Trace3_New;
    Trace4 = Trace4_New;
    Trace5 = Trace5_New;

    logdet = logdet_new;

    InvMat = InvMat - (((fabs(a) > tA) * a - (fabs(a_new) > tA) * a_new) / (1 + ((fabs(a) > tA) * a -
      (fabs(a_new) > tA) * a_new) * InvMat(j, i))) * (InvMat.col(i) * InvMat.row(j));

  }

  // Return a, trace values, logdet and InvMat
  return Rcpp::List::create(Rcpp::Named("a") = a, Rcpp::Named("Trace3") = Trace3,
                            Rcpp::Named("Trace4") = Trace4, Rcpp::Named("Trace5") = Trace5,
                            Rcpp::Named("logdet") = logdet, Rcpp::Named("InvMat") = InvMat);

}


} // end anonymous namespace



//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////





// Do MCMC sampling with Spike and Slab Prior for the model when only Y is there and Sigma follows inverse wishart distribution
// [[Rcpp::export]]
Rcpp::List RGM_SpikeSlab_IW_Star(const arma::mat& S_YY, double n,
                                 int nIter, int nBurnin, int Thin,
                                 double a_rho = 3.0, double b_rho = 1.0,
                                 double nu_1 = 0.001, double Prop_VarA = 0.01){


  // Calculate number of nodes from S_YY matrix
  int p = S_YY.n_cols;

  // Take b_Tilde = p + 1 (b_Tilde should be greater than p-1; p+1 is taken to ensure mean is not infinity)
  int b_Tilde = p + 1;

  // Initialize matrix A
  arma::mat A = arma::zeros(p, p);

  // Initialize Sigma
  arma::mat Sigma = riwish(b_Tilde, arma::eye(p, p));

  // Initialize Rho, Gamma and Tau matrix
  arma::mat Rho = arma::zeros(p, p);
  arma::mat Gamma = arma::ones(p, p);
  arma::mat Tau = arma::ones(p, p);

  // Make the diagonals of Gamma and Tau matrix to be 0
  Gamma.diag().zeros();
  Tau.diag().zeros();

  // Initialize acceptance counter
  double AccptA = 0;

  // Calculate number of posterior samples
  int nPst = std::floor((nIter - nBurnin) / Thin);

  // Initialize Itr to index the posterior samples
  int Itr = 0;

  // Initialize posterior arrays and matrices
  arma::cube A_Pst = arma::zeros(p, p, nPst);
  arma::cube Gamma_Pst = arma::zeros(p, p, nPst);
  arma::cube Tau_Pst = arma::zeros(p, p, nPst);
  arma::cube Rho_Pst = arma::zeros(p, p, nPst);
  arma::cube Sigma_Pst = arma::zeros(p, p, nPst);

  // Initialize LogLikelihood vector
  arma::colvec LL_Pst = arma::zeros(nPst);

  // Calculate Sigma_Inv
  arma::mat Sigma_Inv = arma::inv(Sigma);

  // Run a loop to do MCMC sampling
  for(int i = 1; i <= nIter; i++){

    // Calculate I_p - A
    arma::mat MultMat = arma::eye(p, p) - A;

    ////////////////////
    // Update Sigma
    double param1 = b_Tilde + n;
    arma::mat param2 = n * MultMat * S_YY * MultMat.t() + arma::eye(p, p);

    // In order to make param2 completely symmetric, some entries are slightly different
    param2 = 0.5 * (param2 + param2.t());

    Sigma_Inv = arma::inv(rinvwish(param1, param2));



    ////////////////////
    // Update A
    // Calculate trace values
    double Trace3 = - n * arma::trace(S_YY * A.t() * Sigma_Inv);
    double Trace4 = - n * arma::trace(S_YY * Sigma_Inv * A);
    double Trace5 = n * arma::trace(S_YY * A.t() * Sigma_Inv * A);

    // Calculate logdet and (I - A)^(-1)
    double logdet = real(arma::log_det(MultMat));
    arma::mat InvMat = arma::inv(MultMat);

    // Update Rho, Tau, Gamma and a
    for (int j = 0; j < p; j++) {

      for (int l = 0; l < p; l++) {

        // Don't update the diagonal entries
        if (l != j) {

          // Sample Rho
          Rho(j, l) = Sample_Rho(Gamma(j, l), a_rho, b_rho);

          // Sample Tau
          Tau(j, l) = Sample_Tau(A(j, l), Gamma(j, l), Tau(j, l), nu_1);

          // Sample Gamma
          Gamma(j, l) = Sample_Gamma(A(j, l), Tau(j, l), Rho(j, l), nu_1);

          // Sample a
          Rcpp::List Output2 = Sample_A_Full_Star(S_YY, A, A, j, l, Sigma_Inv, n, p, Gamma(j, l), Tau(j, l),
                                                  nu_1, Prop_VarA, -1, Trace3, Trace4, Trace5, InvMat, logdet);


          double a = Output2[0];

          // Update acceptance counter
          if (A(j, l) != a) {

            // Increase AccptA
            AccptA = AccptA + 1;

            // Update trace values
            Trace3 = Output2[1];
            Trace4 = Output2[2];
            Trace5 = Output2[3];

            // Update logdet and (I - A)^(-1)
            logdet = Output2[4];
            InvMat = Rcpp::as<arma::mat>(Output2[5]);

          }

          // Update A
          A(j, l) = a;

        }

      }

    }



    // Store posterior samples
    if((i > nBurnin) && (i % Thin == 0)){

      A_Pst.slice(Itr) = A;
      Gamma_Pst.slice(Itr) = Gamma;
      Tau_Pst.slice(Itr) = Tau;
      Rho_Pst.slice(Itr) = Rho;
      Sigma_Pst.slice(Itr) = arma::inv(Sigma_Inv);
      LL_Pst(Itr) = LL_Full_Star(A, S_YY, Sigma_Inv, p, n);


      // Increase Itr by 1
      Itr = Itr + 1;

    }

    // Check if `i` is divisible by 100 and print progress
    //if (i % 100 == 0) {
    //Rcpp::Rcout << "Iterations " << i << " done." << std::endl;
    //}

  }

  // Calculate estimates based on posterior samples
  arma::mat A_Est = mean(A_Pst, 2);
  arma::mat Gamma_Est = mean(Gamma_Pst, 2);
  arma::mat Tau_Est = mean(Tau_Pst, 2);
  arma::mat Rho_Est = mean(Rho_Pst, 2);
  arma::mat Sigma_Est = mean(Sigma_Pst, 2);

  // Construct the graph structures
  arma::umat logicalGraph_A = (Gamma_Est > 0.5);
  arma::mat zA_Est = arma::conv_to<arma::mat>::from(logicalGraph_A);

  // Return outputs
  return Rcpp::List::create(Rcpp::Named("A_Est") = A_Est,
                            Rcpp::Named("zA_Est") = zA_Est,
                            Rcpp::Named("Gamma_Est") = Gamma_Est, Rcpp::Named("Tau_Est") = Tau_Est,
                            Rcpp::Named("Rho_Est") = Rho_Est,
                            Rcpp::Named("Sigma_Est") = Sigma_Est,
                            Rcpp::Named("AccptA") = AccptA / (p * (p - 1) * nIter) * 100,
                            Rcpp::Named("LL_Pst") = LL_Pst, Rcpp::Named("Gamma_Pst") = Gamma_Pst);


}



