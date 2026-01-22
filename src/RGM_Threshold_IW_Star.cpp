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




// Sample random number from truncated normal
double Sample_tn(double mu, double sigma, double a, double b) {

  // Calculate alpha and beta
  double alpha = (a - mu) / sigma;
  double beta = (b - mu) / sigma;

  // Calculate CDF
  double cdf_alpha = arma::normcdf(alpha, 0.0, 1.0);
  double cdf_beta = arma::normcdf(beta, 0.0, 1.0);

  // Sample from truncated normal with mean mu and sd sigma
  double u = Rcpp::runif(1, cdf_alpha, cdf_beta)(0);
  double x = R::qnorm(u, 0.0, 1.0, true, false) * sigma + mu;

  // Return x
  return x;

}


// Pdf of truncated normal distribution
double tn_pdf(double x, double mu, double sigma, double a, double b) {

  // Calculate alpha and beta
  double alpha = (a - mu) / sigma;
  double beta = (b - mu) / sigma;

  // Calculate density
  double d = arma::normpdf(x, mu, sigma) / (arma::normcdf(beta, 0.0, 1.0) - arma::normcdf(alpha, 0.0, 1.0));

  // Return d
  return d;

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




// Do MCMC sampling with threshold prior for the model when only Y is there and Sigma follows inverse wishart distribution
// [[Rcpp::export]]
Rcpp::List RGM_Threshold_IW_Star(const arma::mat& S_YY, double n,
                                 int nIter, int nBurnin, int Thin,
                                 double nu_1 = 0.0001, double Prop_VarA = 0.01){


  // Calculate number of nodes from S_YY matrix
  int p = S_YY.n_cols;

  // Take b_Tilde = p + 1 (b_Tilde should be greater than p-1; p+1 is taken to ensure mean is not infinity)
  int b_Tilde = p + 1;

  // Initialize A, A_Pseudo matrices
  arma::mat A = arma::zeros(p, p);
  arma::mat A_Pseudo = arma::zeros(p, p);

  // Initialize Sigma
  arma::mat Sigma = riwish(b_Tilde, arma::eye(p, p));

  // Initialize Gamma and Tau matrices
  arma::mat Gamma = arma::ones(p, p);
  arma::mat Tau = arma::ones(p, p);

  // Make the diagonals of Gamma and Tau matrices to be 0
  Gamma.diag().zeros();
  Tau.diag().zeros();



  // Initialize tA, t0 and t_sd
  double tA = 0;
  double t0 = 0.1;
  double t_sd = 0.1;

  // Initialize acceptance counter
  double AccptA = 0;
  double Accpt_tA = 0;

  // Calculate number of posterior samples
  int nPst = std::floor((nIter - nBurnin) / Thin);

  // Initiate Itr to index the posterior samples
  int Itr = 0;

  // Initialize posterior arrays and matrices
  arma::cube A_Pst = arma::zeros(p, p, nPst);
  arma::cube A0_Pst = arma::zeros(p, p, nPst);
  arma::cube Gamma_Pst = arma::zeros(p, p, nPst);
  arma::cube Tau_Pst = arma::zeros(p, p, nPst);
  arma::colvec tA_Pst = arma::zeros(nPst);
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

    // Calculate det(I - A) and (I - A)^(-1)
    double logdet = real(arma::log_det(MultMat));
    arma::mat InvMat = arma::inv(MultMat);

    // Update Tau based on a and then update a based on Tau
    for (int j = 0; j < p; j++) {

      for (int l = 0; l < p; l++) {

        // Don't update the diagonal entries
        if (l != j) {

          // Update Tau
          Tau(j, l) = Sample_Tau(A_Pseudo(j, l), 1, Tau(j, l), nu_1);

          // Sample a
          Rcpp::List Output2 = Sample_A_Full_Star(S_YY, A, A_Pseudo, j, l, Sigma_Inv, n, p, 1, Tau(j, l), nu_1,
                                                  Prop_VarA, tA, Trace3, Trace4, Trace5,  InvMat, logdet);
          double a = Output2[0];

          // Update Acceptance counter
          if (A_Pseudo(j, l) != a) {

            // Increase AccptA
            AccptA = AccptA + 1;

            // Update Trace values
            Trace3 = Output2[1];
            Trace4 = Output2[2];
            Trace5 = Output2[3];

            // Update logdet and (I - A)^(-1)
            logdet = Output2[4];
            InvMat = Rcpp::as<arma::mat>(Output2[5]);


          }

          // Update A_Pseudo, A and Gamma
          A_Pseudo(j, l) = a;

          A(j, l) = (std::abs(A_Pseudo(j, l)) > tA) * A_Pseudo(j, l);

          Gamma(j, l) = (std::abs(A_Pseudo(j, l)) > tA) * 1;

        }

      }

    }

    // Propose tA_new
    double tA_new = Sample_tn(tA, t_sd, 0, t0);

    // Create A_new
    arma::mat A_new = A_Pseudo % (arma::abs(A_Pseudo) > tA_new);

    // Calculate Difference
    double Diff_A = LL_Full_Star(A_new, S_YY, Sigma_Inv, p, n) - LL_Full_Star(A, S_YY, Sigma_Inv, p, n)
      + log(tn_pdf(tA, tA_new, t_sd, 0, t0)) - log(tn_pdf(tA_new, tA, t_sd, 0, t0));

    // Compare Diff with log of a random number from uniform(0, 1)
    if (Diff_A > log(Rcpp::runif(1, 0, 1)(0))) {

      // Update A, tA and Accpt_tA
      A = A_new;

      tA = tA_new;

      Accpt_tA = Accpt_tA + 1;

    }


    // Store posterior samples
    if((i > nBurnin) && (i % Thin == 0)){

      A_Pst.slice(Itr) = A;
      A0_Pst.slice(Itr) = A_Pseudo;
      Gamma_Pst.slice(Itr) = Gamma;
      Tau_Pst.slice(Itr) = Tau;
      tA_Pst(Itr) = tA;
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
  arma::mat A0_Est = mean(A0_Pst, 2);
  arma::mat Gamma_Est = mean(Gamma_Pst, 2);
  arma::mat Tau_Est = mean(Tau_Pst, 2);
  double tA_Est = mean(tA_Pst);
  arma::mat Sigma_Est = mean(Sigma_Pst, 2);

  // Construct the graph structures
  arma::umat logicalGraph_A = (Gamma_Est > 0.5);
  arma::mat zA_Est = arma::conv_to<arma::mat>::from(logicalGraph_A);


  // Return outputs
  return Rcpp::List::create(Rcpp::Named("A_Est") = A_Est,
                            Rcpp::Named("zA_Est") = zA_Est,
                            Rcpp::Named("A0_Est") = A0_Est,
                            Rcpp::Named("Gamma_Est") = Gamma_Est, Rcpp::Named("Tau_Est") = Tau_Est,
                            Rcpp::Named("tA_Est") = tA_Est,
                            Rcpp::Named("Sigma_Est") = Sigma_Est,
                            Rcpp::Named("AccptA") = AccptA / (p * (p - 1) * nIter) * 100,
                            Rcpp::Named("Accpt_tA") = Accpt_tA / (nIter) * 100,
                            Rcpp::Named("LL_Pst") = LL_Pst, Rcpp::Named("Gamma_Pst") = Gamma_Pst);



}
