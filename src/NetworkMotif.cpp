#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// Define the function NetworkMotif_cpp with two parameters: Gamma and Gamma_Pst
// The function calculates the network motif based on the given parameters
// [[Rcpp::export]]
double NetworkMotif_cpp(const arma::mat& Gamma, const arma::cube& Gamma_Pst) {

  // Calculate the number of posterior samples
  double nPst = Gamma_Pst.n_slices;

  // Initialize Count to keep track of matching network motifs
  double Count = 0;

  // Loop through each posterior sample
  for (int i = 0; i < nPst; i++) {

    // Extract the Gamma matrix for the current posterior sample
    arma::mat Gamma_Test = Gamma_Pst.slice(i);

    // Find indices where Gamma is equal to 1
    arma::uvec indices = arma::find(Gamma == 1);

    // Calculate the absolute difference between the Gamma matrix and the Gamma_Test matrix where Gamma values are 1
    double Diff = arma::accu(arma::abs(Gamma.elem(indices) - Gamma_Test.elem(indices)));

    // If the difference is zero, it means the network motif matches
    if (Diff == 0) {

      // Increment the count of matching network motifs
      Count = Count + 1;

    }

  }

  // Return the proportion of matching network motifs to the total number of posterior samples
  return (Count / nPst);


}
