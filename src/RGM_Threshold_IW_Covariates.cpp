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



// Sample Eta
double Sample_Eta(double b, double phi, double eta, double nu_2) {

  // Sample Epsilon based on old eta
  double Epsilon = 1 / Rcpp::rgamma(1, 1, 1 / (1 + 1 / eta))(0);

  // Initialize Eta
  double Eta;

  // Check whether phi is 0 or 1
  if (phi == 1) {

    // Sample Eta based on b and Epsilon
    Eta = 1 / Rcpp::rgamma(1, 1, 1 / (b * b / 2 + 1 / Epsilon))(0);

  } else {

    // Sample Eta based on b, nu_2 and Epsilon
    Eta = 1 / Rcpp::rgamma(1, 1, 1 / (b * b / (2 * nu_2) + 1 / Epsilon))(0);

  }

  // Return Eta
  return Eta;

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



// Calculate log-likelihood for the model when X, Y and U are there and Sigma is a full matrix
double LL_Full(const arma::mat& A, const arma::mat& B, const arma::mat& C, const arma::mat& S_YY, const arma::mat& S_YX,
               const arma::mat& S_XX, const arma::mat& S_YU, const arma::mat& S_XU, const arma::mat& S_UU,
               const arma::mat& Sigma_Inv, double p, double N) {

  // Calculate (I_p - A)
  const arma::mat& Mult_Mat = arma::eye(p, p) - A;

  // Calculate Sum
  double Sum = N * arma::trace(S_YY * Mult_Mat.t() * Sigma_Inv * Mult_Mat)
    - 2 * N * arma::trace(S_YX * B.t() * Sigma_Inv * Mult_Mat)
    + N * arma::trace(S_XX * B.t() * Sigma_Inv * B)
    + N * arma::trace(S_UU * C.t() * Sigma_Inv * C)
    - 2 * N * arma::trace(S_YU * C.t() * Sigma_Inv * Mult_Mat)
    + 2 * N * arma::trace(S_XU * C.t() * Sigma_Inv * B);

    // Calculate log-likelihood
    double LL = N * real(arma::log_det(Mult_Mat)) + N / 2 * real(arma::log_det(Sigma_Inv)) - Sum / 2 - N / 2 * log(2 * arma::datum::pi);

    // Return log-likelihood
    return LL;

}


// Calculate target value for a particular B
double Target_B(double b, double phi, double eta, double nu_2, double Trace1, double Trace2, double Trace3) {

  // Calculate Sum
  double Sum = Trace1 + Trace2 + Trace3;

  // Calculate Target value
  double Target = - Sum / 2 - phi * (b * b / (2 * eta)) - (1 - phi) * (0.5 * log(nu_2) + b * b / (2 * nu_2 * eta));

  // Return Target value
  return Target;

}


// Sample a particular entry of matrix B when Sigma is a full matrix
Rcpp::List Sample_B_Full(const arma::mat& S_YX, const arma::mat& S_XX, const arma::mat& S_XU, const arma::mat& B, const arma::mat& C,
                         const arma::mat& B_Pseudo, double i, double j,
                         const arma::mat& Sigma_Inv, const arma::mat& MultMat, double N, double phi, double eta, double nu_2,
                         double prop_var2, double tB, double Trace1, double Trace2, double Trace3) {

  // Value to update
  double b = B_Pseudo(i, j);

  // Proposed value
  double b_new = Rcpp::rnorm(1, b, sqrt(prop_var2))(0);

  // Create a copy of matrix B
  arma::mat B_new = B;

  // Modify the copy with the proposed b value
  B_new(i, j) = (fabs(b_new) > tB) * b_new;

  // Calculate new trace values
  double Trace1_New = Trace1 - 2 * N * arma::trace(((fabs(b_new) > tB) * b_new - (fabs(b) > tB) * b) * (Sigma_Inv.row(i) * MultMat * S_YX.col(j)));
  double Trace2_New = Trace2 + N * arma::trace(((fabs(b_new) > tB) * b_new - (fabs(b) > tB) * b) * (S_XX.row(j) * B.t() * Sigma_Inv.col(i)
                                                                                                      + S_XX.col(j).t() * B_new.t() * Sigma_Inv.row(i).t()));
  double Trace3_New = Trace3 + 2 * N * arma::trace(((fabs(b_new) > tB) * b_new - (fabs(b) > tB) * b) * (S_XU.row(j) * C.t() * Sigma_Inv.col(i)));

  // Calculate target values with b and b_new
  double Target1 = Target_B(b_new, phi, eta, nu_2, Trace1_New, Trace2_New, Trace3_New);
  double Target2 = Target_B(b, phi, eta, nu_2, Trace1, Trace2, Trace3);

  // Calculate r i.e. the difference between two target values
  double r = Target1 - Target2;

  // Sample u from Uniform(0, 1)
  double u = Rcpp::runif(1, 0, 1)(0);

  // Compare u and r
  if (r >= log(u)) {

    // Update b and trace values
    b = b_new;

    Trace1 = Trace1_New;
    Trace2 = Trace2_New;
    Trace3 = Trace3_New;

  }

  // Return b and trace values
  return Rcpp::List::create(Rcpp::Named("b") = b, Rcpp::Named("Trace1") = Trace1,
                            Rcpp::Named("Trace2") = Trace2, Rcpp::Named("Trace3") = Trace3);

}



// Calculate target value for a particular A for the model when X, Y and U are there
double Target_A(double a, double N, double gamma, double tau, double nu_1, double Trace4, double Trace5, double Trace6, double Trace7, double Trace8, double logdet) {

  // Calculate Sum term inside exponential in likelihood
  double Sum = Trace4 + Trace5 + Trace6 + Trace7 + Trace8;

  // Calculate Target value
  double Target = N * logdet - Sum / 2 - gamma * (a * a / (2 * tau)) - (1 - gamma) * (0.5 * log(nu_1) + a * a / (2 * nu_1 * tau));

  // Return Target
  return Target;

}


// Sample a particular entry of matrix A for the model when X, Y and U are there and Sigma is a full matrix
Rcpp::List Sample_A_Full(const arma::mat& S_YY, const arma::mat& S_YX, const arma::mat& S_YU, const arma::mat& A, const arma::mat& A_Pseudo, double i, double j,
                         const arma::mat& Sigma_Inv, double N, double p, const arma::mat& B, const arma::mat& C, double gamma, double tau, double nu_1,
                         double prop_var1, double tA, double Trace4, double Trace5, double Trace6, double Trace7, double Trace8, arma::mat& InvMat, double logdet) {

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
  double Trace4_New = Trace4 - N * ((fabs(a_new) > tA) * a_new - (fabs(a) > tA) * a) * arma::trace(Sigma_Inv.row(i) * S_YY.col(j));
  double Trace5_New = Trace5 - N * ((fabs(a_new) > tA) * a_new - (fabs(a) > tA) * a) * arma::trace(S_YY.row(j) * Sigma_Inv.col(i));
  double Trace6_New = Trace6 + N * ((fabs(a_new) > tA) * a_new - (fabs(a) > tA) * a) * arma::trace(S_YY.row(j) * A.t() * Sigma_Inv.col(i)
                                                                                                     + S_YY.col(j).t() * A_new.t() * Sigma_Inv.row(i).t());
  double Trace7_New = Trace7 + 2 * N * ((fabs(a_new) > tA) * a_new - (fabs(a) > tA) * a) * arma::trace(S_YX.row(j) * B.t() * Sigma_Inv.col(i));
  double Trace8_New = Trace8 + 2 * N * ((fabs(a_new) > tA) * a_new - (fabs(a) > tA) * a) * arma::trace(S_YU.row(j) * C.t() * Sigma_Inv.col(i));

  // Calculate target values with a and a_new
  double Target1 = Target_A(a_new, N, gamma, tau, nu_1, Trace4_New, Trace5_New, Trace6_New, Trace7_New, Trace8_New, logdet_new);
  double Target2 = Target_A(a, N, gamma, tau, nu_1, Trace4, Trace5, Trace6, Trace7, Trace8, logdet);

  // Calculate r i.e. the differnce between two target values
  double r = Target1 - Target2;

  // Sample u from Uniform(0, 1)
  double u = Rcpp::runif(1, 0, 1)(0);

  // Compare u and r
  if (r >= log(u)) {

    // Update a, trace values, logdet and InvMat
    a = a_new;

    Trace4 = Trace4_New;
    Trace5 = Trace5_New;
    Trace6 = Trace6_New;
    Trace7 = Trace7_New;
    Trace8 = Trace8_New;

    logdet = logdet_new;

    InvMat = InvMat - (((fabs(a) > tA) * a - (fabs(a_new) > tA) * a_new) / (1 + ((fabs(a) > tA) * a -
      (fabs(a_new) > tA) * a_new) * InvMat(j, i))) * (InvMat.col(i) * InvMat.row(j));

  }

  // Return a, trace values, logdet and InvMat
  return Rcpp::List::create(Rcpp::Named("a") = a, Rcpp::Named("Trace4") = Trace4,
                            Rcpp::Named("Trace5") = Trace5, Rcpp::Named("Trace6") = Trace6,
                            Rcpp::Named("Trace7") = Trace7, Rcpp::Named("Trace8") = Trace8,
                            Rcpp::Named("logdet") = logdet, Rcpp::Named("InvMat") = InvMat);


}


// Sample from matrix variate normal
arma::mat rmatrixvarnorm(const arma::mat& M, const arma::mat& U, const arma::mat& V) {
  Rcpp::RNGScope scope; // Ensures R's RNG is used
  int n = M.n_rows;
  int p = M.n_cols;

  // Cholesky decomposition of U (row covariance) and V (column covariance)
  arma::mat U_chol = arma::chol(U);
  arma::mat V_chol = arma::chol(V);

  // Generate matrix Z of independent standard normal random variables
  arma::mat Z = arma::randn(n, p); // Z ~ N(0, I)

  // Compute the sample from matrix-variate normal distribution: X = M + U_chol * Z * V_chol.t()
  arma::mat X = M + U_chol.t() * Z * V_chol;

  return X;
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





// Do MCMC sampling with threshold prior for the model when X, Y and U are there and Sigma is full with IW prior on it
// [[Rcpp::export]]
Rcpp::List RGM_Threshold_IW_Covariates(const arma::mat& S_YY, const arma::mat& S_YX, const arma::mat& S_XX,
                                       const arma::mat& S_YU, const arma::mat& S_XU, const arma::mat& S_UU,
                                       const arma::mat& D, double n, int nIter, int nBurnin,
                                       int Thin, double nu_1 = 0.0001, double nu_2 = 0.0001,
                                       double Prop_VarA = 0.01, double Prop_VarB = 0.01, double TAU = 2){


  // Calculate number of nodes from S_YY matrix
  int p = S_YY.n_cols;

  // Calculate number of columns of S_XX
  int k = S_XX.n_cols;

  // Calculate number of columns of S_UU
  int l = S_UU.n_cols;

  // Take b_Tilde = p + 1 (b_Tilde should be greater than p-1; p+1 is taken to ensure mean is not infinity)
  int b_Tilde = p + 1;

  // Initialize A, B, C, A_Pseudo and B_Pseudo matrices
  arma::mat A = arma::zeros(p, p);
  arma::mat B = arma::zeros(p, k);
  arma::mat C = arma::zeros(p, l);
  arma::mat A_Pseudo = arma::zeros(p, p);
  arma::mat B_Pseudo = arma::zeros(p, k);

  // Initialize Sigma
  arma::mat Sigma = riwish(b_Tilde, arma::eye(p, p));

  // Initialize Gamma, Phi, Tau and Eta matrices
  arma::mat Gamma = arma::ones(p, p);
  arma::mat Phi = arma::ones(p, k);
  arma::mat Tau = arma::ones(p, p);
  arma::mat Eta = arma::ones(p, k);


  // Make the diagonals of Gamma and Tau matrices to be 0
  Gamma.diag().zeros();
  Tau.diag().zeros();

  // Make Phi[i, j] = 0 and Eta[i, j] = 0 if D[i, j] = 0
  Phi = Phi % D;
  Eta = Eta % D;

  // Initialize tA, tB, t0 and t_sd
  double tA = 0;
  double tB = 0;
  double t0 = 0.1;
  double t_sd = 0.1;

  // Initialize acceptance counter
  double AccptA = 0;
  double AccptB = 0;
  double Accpt_tA = 0;
  double Accpt_tB = 0;

  // Calculate number of posterior samples
  int nPst = std::floor((nIter - nBurnin) / Thin);

  // Initiate Itr to index the posterior samples
  int Itr = 0;

  // Initialize posterior arrays and matrices
  arma::cube A_Pst = arma::zeros(p, p, nPst);
  arma::cube A0_Pst = arma::zeros(p, p, nPst);
  arma::cube B_Pst = arma::zeros(p, k, nPst);
  arma::cube B0_Pst = arma::zeros(p, k, nPst);
  arma::cube C_Pst = arma::zeros(p, l, nPst);
  arma::cube Gamma_Pst = arma::zeros(p, p, nPst);
  arma::cube Tau_Pst = arma::zeros(p, p, nPst);
  arma::cube Phi_Pst = arma::zeros(p, k, nPst);
  arma::cube Eta_Pst = arma::zeros(p, k, nPst);
  arma::colvec tA_Pst = arma::zeros(nPst);
  arma::colvec tB_Pst = arma::zeros(nPst);
  arma::cube Sigma_Pst = arma::zeros(p, p, nPst);

  // Initialize LogLikelihood vector
  arma::colvec LL_Pst = arma::zeros(nPst);

  // Calculate Sigma_Inv
  arma::mat Sigma_Inv = arma::inv(Sigma);

  // Caculate the sub matrix to sample C
  arma::mat SubMat_C = arma::inv(n * S_UU + arma::eye(l, l) / TAU);

  // Run a loop to do MCMC sampling
  for(int i = 1; i <= nIter; i++){

    // Update B
    // Calculate I_p - A
    arma::mat MultMat = arma::eye(p, p) - A;

    // Calculate Trace values
    double Trace1 = - 2 * n * arma::trace(S_YX * B.t() * Sigma_Inv * MultMat);
    double Trace2 = n * arma::trace(S_XX * B.t() * Sigma_Inv * B);
    double Trace3 = 2 * n * arma::trace(S_XU * C.t() * Sigma_Inv * B);


    // Update Eta based on corresponding b and then update b based on the corresponding eta
    for (int j = 0; j < p; j++) {

      for (int l = 0; l < k; l++) {

        // Don't update if the corresponding D entry is 0
        if (D(j, l) != 0) {

          // Sample Eta
          Eta(j, l) = Sample_Eta(B_Pseudo(j, l), 1, Eta(j, l), nu_2);

          // Sample b
          Rcpp::List Output1 = Sample_B_Full(S_YX, S_XX, S_XU, B, C, B_Pseudo, j, l, Sigma_Inv, MultMat, n, 1,
                                             Eta(j, l), nu_2, Prop_VarB, tB, Trace1, Trace2, Trace3);
          double b = Output1[0];


          // Update acceptance counter
          if (B_Pseudo(j, l) != b) {

            // Increase AccptB
            AccptB = AccptB + 1;

            // Update Trace Values
            Trace1 = Output1[1];
            Trace2 = Output1[2];
            Trace3 = Output1[3];

          }

          // Update B_Pseudo, B and Phi
          B_Pseudo(j, l) = b;

          B(j, l) = (std::abs(B_Pseudo(j, l)) > tB) * B_Pseudo(j, l);

          Phi(j, l) = (std::abs(B_Pseudo(j, l)) > tB) * 1;

        }

      }

    }

    // Propose tB_new
    double tB_new = Sample_tn(tB, t_sd, 0, t0);

    // Create B_new
    arma::mat B_new = B_Pseudo % (arma::abs(B_Pseudo) > tB_new);

    // Calculate difference
    double Diff = LL_Full(A, B_new, C, S_YY, S_YX, S_XX, S_YU, S_XU, S_UU, Sigma_Inv, p, n) - LL_Full(A, B, C, S_YY, S_YX, S_XX, S_YU, S_XU, S_UU, Sigma_Inv, p, n)
      + log(tn_pdf(tB, tB_new, t_sd, 0, t0)) - log(tn_pdf(tB_new, tB, t_sd, 0, t0));

    // Compare Diff with log of a random number from Uniform(0, 1)
    if (Diff > log(Rcpp::runif(1, 0, 1)(0))) {

      // Update B, tB and Accpt_tB
      B = B_new;

      tB = tB_new;

      Accpt_tB = Accpt_tB + 1;

    }


    ////////////////////
    // Update C
    arma::mat Mean_C = (n * MultMat * S_YU - n * B * S_XU) * SubMat_C;
    C = rmatrixvarnorm(Mean_C, Sigma, SubMat_C);


    ////////////////////
    // Update Sigma
    double param1 = b_Tilde + n;
    arma::mat param2 = n * MultMat * S_YY * MultMat.t() - n * MultMat * S_YX * B.t() -
      n * B * S_YX.t() * MultMat.t() + n * B * S_XX * B.t()
      + n * C * S_UU * C.t() - n * MultMat * S_YU * C.t() - n * C * S_YU.t() * MultMat.t()
      + n * B * S_XU * C.t() + n * C * S_XU.t() * B.t()
      + C * C.t() / TAU + arma::eye(p, p);

    // In order to make param2 completely symmetric, some entries are slightly different
    param2 = 0.5 * (param2 + param2.t());

    Sigma = rinvwish(param1, param2);
    Sigma_Inv = arma::inv(Sigma);



      ////////////////////
      // Update A
      // Calculate trace values
      double Trace4 = - n * arma::trace(S_YY * A.t() * Sigma_Inv);
      double Trace5 = - n * arma::trace(S_YY * Sigma_Inv * A);
      double Trace6 = n * arma::trace(S_YY * A.t() * Sigma_Inv * A);
      double Trace7 = 2 * n * arma::trace(S_YX * B.t() * Sigma_Inv * A);
      double Trace8 = 2 * n * arma::trace(S_YU * C.t() * Sigma_Inv * A);

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
            Rcpp::List Output2 = Sample_A_Full(S_YY, S_YX, S_YU, A, A_Pseudo, j, l, Sigma_Inv, n, p, B, C, 1, Tau(j, l), nu_1,
                                               Prop_VarA, tA, Trace4, Trace5, Trace6, Trace7, Trace8, InvMat, logdet);

            double a = Output2[0];


            // Update Acceptance counter
            if (A_Pseudo(j, l) != a) {

              // Increase AccptA
              AccptA = AccptA + 1;

              // Update Trace values
              Trace4 = Output2[1];
              Trace5 = Output2[2];
              Trace6 = Output2[3];
              Trace7 = Output2[4];
              Trace8 = Output2[5];

              // Update logdet and (I - A)^(-1)
              logdet = Output2[6];
              InvMat = Rcpp::as<arma::mat>(Output2[7]);


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
      double Diff_A = LL_Full(A_new, B, C, S_YY, S_YX, S_XX, S_YU, S_XU, S_UU, Sigma_Inv, p, n) - LL_Full(A, B, C, S_YY, S_YX, S_XX, S_YU, S_XU, S_UU, Sigma_Inv, p, n)
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
        B_Pst.slice(Itr) = B;
        B0_Pst.slice(Itr) = B_Pseudo;
        C_Pst.slice(Itr) = C;
        Gamma_Pst.slice(Itr) = Gamma;
        Tau_Pst.slice(Itr) = Tau;
        Phi_Pst.slice(Itr) = Phi;
        Eta_Pst.slice(Itr) = Eta;
        tA_Pst(Itr) = tA;
        tB_Pst(Itr) = tB;
        Sigma_Pst.slice(Itr) = Sigma;
        LL_Pst(Itr) = LL_Full(A, B, C, S_YY, S_YX, S_XX, S_YU, S_XU, S_UU, Sigma_Inv, p, n);

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
  arma::mat B_Est = mean(B_Pst, 2);
  arma::mat A0_Est = mean(A0_Pst, 2);
  arma::mat B0_Est = mean(B0_Pst, 2);
  arma::mat C_Est = mean(C_Pst, 2);
  arma::mat Gamma_Est = mean(Gamma_Pst, 2);
  arma::mat Tau_Est = mean(Tau_Pst, 2);
  arma::mat Phi_Est = mean(Phi_Pst, 2);
  arma::mat Eta_Est = mean(Eta_Pst, 2);
  double tA_Est = mean(tA_Pst);
  double tB_Est = mean(tB_Pst);
  arma::mat Sigma_Est = mean(Sigma_Pst, 2);

  // Construct the graph structures
  arma::umat logicalGraph_A = (Gamma_Est > 0.5);
  arma::umat logicalGraph_B = (Phi_Est > 0.5);
  arma::mat zA_Est = arma::conv_to<arma::mat>::from(logicalGraph_A);
  arma::mat zB_Est = arma::conv_to<arma::mat>::from(logicalGraph_B);


  // Return outputs
  return Rcpp::List::create(Rcpp::Named("A_Est") = A_Est, Rcpp::Named("B_Est") = B_Est,
                            Rcpp::Named("C_Est") = C_Est,
                            Rcpp::Named("zA_Est") = zA_Est, Rcpp::Named("zB_Est") = zB_Est,
                            Rcpp::Named("A0_Est") = A0_Est, Rcpp::Named("B0_Est") = B0_Est,
                            Rcpp::Named("Gamma_Est") = Gamma_Est, Rcpp::Named("Tau_Est") = Tau_Est,
                            Rcpp::Named("Phi_Est") = Phi_Est, Rcpp::Named("Eta_Est") = Eta_Est,
                            Rcpp::Named("tA_Est") = tA_Est, Rcpp::Named("tB_Est") = tB_Est,
                            Rcpp::Named("Sigma_Est") = Sigma_Est,
                            Rcpp::Named("AccptA") = AccptA / (p * (p - 1) * nIter) * 100, Rcpp::Named("AccptB") = AccptB / (arma::accu(D) * nIter) * 100,
                            Rcpp::Named("Accpt_tA") = Accpt_tA / (nIter) * 100, Rcpp::Named("Accpt_tB") = Accpt_tB / (nIter) * 100,
                            Rcpp::Named("LL_Pst") = LL_Pst, Rcpp::Named("Gamma_Pst") = Gamma_Pst);



}
