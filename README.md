
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MR.RGM

<!-- badges: start -->

[![R-CMD-check](https://github.com/bitansa/MR.RGM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bitansa/MR.RGM/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Intended use of the package

The MR.RGM R package presents a crucial advancement in Mendelian
randomization (MR) studies, providing a robust solution to a common
challenge. While MR has proven invaluable in establishing causal links
between exposures and outcomes, its traditional focus on single
exposures and specific outcomes can be limiting. Biological systems
often exhibit complexity, with interdependent outcomes influenced by
numerous factors. MR.RGM introduces a network-based approach to MR,
allowing researchers to explore the broader causal landscape.

With two available functions, RGM and NetworkMotif, the package offers
versatility in analyzing causal relationships. RGM primarily focuses on
constructing causal networks among response variables and between
responses and instrumental variables. On the other hand, NetworkMotif
specializes in quantifying uncertainty for given network structures
among response variables.

RGM accommodates both individual-level data and two types of
summary-level data, making it adaptable to various data availability
scenarios. This adaptability enhances the package’s utility across
different research contexts. The outputs of RGM include estimates of
causal effects, adjacency matrices, and other relevant parameters.
Together, these outputs contribute to a deeper understanding of the
intricate relationships within complex biological networks, thereby
enriching insights derived from MR studies.

## Installation instructions

You can install MR.RGM R package from CRAN with:

``` r

    install.packages("MR.RGM")
```

Once the MR.RGM package is installed load the library in the R
work-space.

``` r

     library("MR.RGM")
```

## Example

We offer a concise demonstration of the capabilities of the RGM function
within the package, showcasing its effectiveness in computing causal
interactions among response variables and between responses and
instrumental variables using simulated data sets. Subsequently, we
provide an example of how NetworkMotif can be applied, utilizing a
specified network structure and Gamma_Pst acquired from executing the
RGM function.

``` r

# Model: Y = AY + BX + E

# Set seed
set.seed(9154)

# Number of data points
n = 10000

# Number of response variables and number of instrument variables
p = 5
k = 6

# Create d vector
d = c(2, 1, 1, 1, 1)

# Initialize causal interaction matrix between response variables
A = matrix(sample(c(-0.1, 0.1), p^2, replace = TRUE), p, p)

# Diagonal entries of A matrix will always be 0
diag(A) = 0

# Make the network sparse
A[sample(which(A!=0), length(which(A!=0))/2)] = 0

# Initialize causal interaction matrix between response and instrument variables
B = matrix(0, p, k)

# Initialize m
m = 1

# Calculate B matrix based on d vector
for (i in 1:p) {
 
  # Update ith row of B
  B[i, m:(m + d[i] - 1)] = 1
 
  # Update m
  m = m + d[i]
 
}

# Create variance-covariance matrix
Sigma = 1 * diag(p)

Mult_Mat = solve(diag(p) - A)

Variance = Mult_Mat %*% Sigma %*% t(Mult_Mat)

# Generate instrument data matrix
X = matrix(runif(n * k, 0, 5), nrow = n, ncol = k)

# Initialize response data matrix
Y = matrix(0, nrow = n, ncol = p)

# Generate response data matrix based on instrument data matrix
for (i in 1:n) {

 Y[i, ] = MASS::mvrnorm(n = 1, Mult_Mat %*% B %*% X[i, ], Variance)

}

# Print true causal interaction matrices between response variables and between response and instrument variables
A
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]  0.0 -0.1  0.0  0.0  0.1
#> [2,]  0.1  0.0 -0.1  0.1  0.1
#> [3,]  0.0 -0.1  0.0  0.0  0.1
#> [4,]  0.0 -0.1  0.0  0.0  0.0
#> [5,]  0.0  0.1  0.0  0.0  0.0
B
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    1    0    0    0    0
#> [2,]    0    0    1    0    0    0
#> [3,]    0    0    0    1    0    0
#> [4,]    0    0    0    0    1    0
#> [5,]    0    0    0    0    0    1
```

We will now apply RGM based on individual level data, summary level data
and Beta, SigmaHat matrices to show its functionality.

``` r

# Apply RGM on individual level data with Threshold prior
Output1 = RGM(X = X, Y = Y, d = c(2, 1, 1, 1, 1), prior = "Threshold")

# Calculate summary level data
Syy = t(Y) %*% Y / n
Syx = t(Y) %*% X / n
Sxx = t(X) %*% X / n

# Apply RGM on summary level data for Spike and Slab Prior
Output2 = RGM(Syy = Syy, Syx = Syx, Sxx = Sxx,
           d = c(2, 1, 1, 1, 1), n = 10000, prior = "Spike and Slab")

# Calculate Beta and Sigma_Hat
# Centralize Data
Y = t(t(Y) - colMeans(Y))
X = t(t(X) - colMeans(X))

# Calculate Sxx
Sxx = t(X) %*% X / n

# Generate Beta matrix and SigmaHat
Beta = matrix(0, nrow = p, ncol = k)
SigmaHat = matrix(0, nrow = p, ncol = k)

for (i in 1:p) {

    for (j in 1:k) {

        fit = lm(Y[, i] ~ X[, j])

        Beta[i, j] =  fit$coefficients[2]

        SigmaHat[i, j] = sum(fit$residuals^2) / n

        }

 }


# Apply RGM on Sxx, Beta and SigmaHat for Spike and Slab Prior
Output3 = RGM(Sxx = Sxx, Beta = Beta, SigmaHat = SigmaHat,
           d = c(2, 1, 1, 1, 1), n = 10000, prior = "Spike and Slab")
```

We get the estimated causal interaction matrix between response
variables in the following way:

``` r

Output1$AEst
#>           [,1]        [,2]       [,3]      [,4]       [,5]
#> [1,] 0.0000000 -0.11032661  0.0000000 0.0000000 0.10676811
#> [2,] 0.0991208  0.00000000 -0.1104576 0.1002182 0.11012341
#> [3,] 0.0000000 -0.09370579  0.0000000 0.0000000 0.09747664
#> [4,] 0.0000000 -0.10185959  0.0000000 0.0000000 0.00000000
#> [5,] 0.0000000  0.10045256  0.0000000 0.0000000 0.00000000
Output2$AEst
#>              [,1]       [,2]          [,3]          [,4]        [,5]
#> [1,]  0.000000000 -0.1127249 -0.0006355133  0.0012237310 0.107520412
#> [2,]  0.099881480  0.0000000 -0.1079853541  0.0996589823 0.109275947
#> [3,] -0.001747246 -0.0929592  0.0000000000 -0.0006473297 0.099778267
#> [4,] -0.003863658 -0.1030056 -0.0024683725  0.0000000000 0.009191016
#> [5,]  0.001966818  0.1014136 -0.0055458038 -0.0050688662 0.000000000
Output3$AEst
#>             [,1]        [,2]         [,3]          [,4]       [,5]
#> [1,]  0.00000000 -0.08972628  0.039412442 -0.0006516754 0.09454632
#> [2,]  0.11040214  0.00000000 -0.112569739  0.0983975424 0.13676682
#> [3,]  0.01896308 -0.09780459  0.000000000  0.0100459671 0.11765744
#> [4,]  0.00000000 -0.11255315  0.001661851  0.0000000000 0.01358832
#> [5,] -0.00174675  0.12947084  0.011533919 -0.0025906139 0.00000000
```

We get the estimated causal network structure between the response
variables in the following way:

``` r

Output1$zAEst
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    1    0    0    1
#> [2,]    1    0    1    1    1
#> [3,]    0    1    0    0    1
#> [4,]    0    1    0    0    0
#> [5,]    0    1    0    0    0
Output2$zAEst
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    1    0    0    1
#> [2,]    1    0    1    1    1
#> [3,]    0    1    0    0    1
#> [4,]    0    1    0    0    0
#> [5,]    0    1    0    0    0
Output3$zAEst
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    1    0    0    1
#> [2,]    1    0    1    1    1
#> [3,]    0    1    0    0    1
#> [4,]    0    1    0    0    0
#> [5,]    0    1    0    0    0
```

We observe that the causal network structures inferred in the three
outputs mentioned are identical. To gain a clearer understanding of the
network, we compare the true network structure with the one estimated by
RGM. Since the networks derived from all three outputs are consistent,
we plot a single graph representing the estimated causal network.

``` r

# Define a function to create smaller arrowheads
smaller_arrowheads <- function(graph) {
  igraph::E(graph)$arrow.size = 0.60  # Adjust the arrow size value as needed
  return(graph)
}

# Create a layout for multiple plots
par(mfrow = c(1, 2))

# Plot the true causal network
plot(smaller_arrowheads(igraph::graph_from_adjacency_matrix((A != 0) * 1,
       mode = "directed")), layout = igraph::layout_in_circle,
          main = "True Causal Network")

# Plot the estimated causal network
plot(smaller_arrowheads(igraph::graph_from_adjacency_matrix(Output1$zAEst,
      mode = "directed")), layout = igraph::layout_in_circle,
         main = "Estimated Causal Network")
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

We get the estimated causal interaction matrix between the response and
the instrument variables from the outputs in the following way:

``` r

Output1$BEst
#>           [,1]     [,2]      [,3]      [,4]      [,5]      [,6]
#> [1,] 0.9921448 1.010301 0.0000000 0.0000000 0.0000000 0.0000000
#> [2,] 0.0000000 0.000000 0.9961312 0.0000000 0.0000000 0.0000000
#> [3,] 0.0000000 0.000000 0.0000000 0.9989298 0.0000000 0.0000000
#> [4,] 0.0000000 0.000000 0.0000000 0.0000000 0.9991987 0.0000000
#> [5,] 0.0000000 0.000000 0.0000000 0.0000000 0.0000000 0.9977224
Output2$BEst
#>           [,1]     [,2]      [,3]      [,4]     [,5]     [,6]
#> [1,] 0.9942219 1.009985 0.0000000 0.0000000 0.000000 0.000000
#> [2,] 0.0000000 0.000000 0.9931885 0.0000000 0.000000 0.000000
#> [3,] 0.0000000 0.000000 0.0000000 0.9992243 0.000000 0.000000
#> [4,] 0.0000000 0.000000 0.0000000 0.0000000 1.000336 0.000000
#> [5,] 0.0000000 0.000000 0.0000000 0.0000000 0.000000 1.001666
Output3$BEst
#>           [,1]     [,2]     [,3]      [,4]      [,5]      [,6]
#> [1,] 0.9910789 1.006024 0.000000 0.0000000 0.0000000 0.0000000
#> [2,] 0.0000000 0.000000 0.993928 0.0000000 0.0000000 0.0000000
#> [3,] 0.0000000 0.000000 0.000000 0.9987596 0.0000000 0.0000000
#> [4,] 0.0000000 0.000000 0.000000 0.0000000 0.9974289 0.0000000
#> [5,] 0.0000000 0.000000 0.000000 0.0000000 0.0000000 0.9948345
```

We get the estimated graph structure between the response and the
instrument variables from the outputs in the following way:

``` r

Output1$zBEst
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    1    0    0    0    0
#> [2,]    0    0    1    0    0    0
#> [3,]    0    0    0    1    0    0
#> [4,]    0    0    0    0    1    0
#> [5,]    0    0    0    0    0    1
Output2$zBEst
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    1    0    0    0    0
#> [2,]    0    0    1    0    0    0
#> [3,]    0    0    0    1    0    0
#> [4,]    0    0    0    0    1    0
#> [5,]    0    0    0    0    0    1
Output3$zBEst
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    1    0    0    0    0
#> [2,]    0    0    1    0    0    0
#> [3,]    0    0    0    1    0    0
#> [4,]    0    0    0    0    1    0
#> [5,]    0    0    0    0    0    1
```

We can plot the log-likelihoods from the outputs in the following way:

``` r

plot(Output1$LLPst, type = 'l', xlab = "Iterations", ylab = "Log-likelihood")
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

``` r
plot(Output2$LLPst, type = 'l', xlab = "Iterations", ylab = "Log-likelihood")
```

<img src="man/figures/README-unnamed-chunk-11-2.png" width="100%" />

``` r
plot(Output3$LLPst, type = 'l', xlab = "Iterations", ylab = "Log-likelihood")
```

<img src="man/figures/README-unnamed-chunk-11-3.png" width="100%" />

Next, we present the implementation of the NetworkMotif function. We
begin by defining a random subgraph among the response variables.
Subsequently, we collect GammaPst arrays from various outputs and
proceed to execute NetworkMotif based on these arrays.

``` r

# Start with a random subgraph
Gamma = matrix(0, nrow = p, ncol = p)
Gamma[5, 2] = Gamma[3, 5] = Gamma[2, 3] = 1

# Plot the subgraph to get an idea about the causal network
plot(smaller_arrowheads(igraph::graph_from_adjacency_matrix(Gamma,
       mode = "directed")), layout = igraph::layout_in_circle,
          main = "Subgraph")
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

``` r


# Store the GammaPst arrays from outputs
GammaPst1 = Output1$GammaPst
GammaPst2 = Output2$GammaPst
GammaPst3 = Output3$GammaPst

# Get the posterior probabilities of Gamma with these GammaPst matrices
NetworkMotif(Gamma = Gamma, GammaPst = GammaPst1)
#> [1] 1
NetworkMotif(Gamma = Gamma, GammaPst = GammaPst2)
#> [1] 0.37325
NetworkMotif(Gamma = Gamma, GammaPst = GammaPst3)
#> [1] 0.461375
```

## References

Yang Ni. Yuan Ji. Peter Müller. “Reciprocal Graphical Models for
Integrative Gene Regulatory Network Analysis.” Bayesian Anal. 13 (4)
1095 - 1110, December 2018. <https://doi.org/10.1214/17-BA1087>
