## Intended use of the package

For a given dataset of gene expressions and covariates, this R-package will estimate the gene-gene interactions and the gene-covariate interactions and then create a graph connecting them. A spike and slab prior is assumed for each interaction parameter basd on our prior knowledge to calculate posterior estimates of the interactions. The package includes various helper functions that generate gibbs samples of the parameter values for the prior distributions and interaction parameter values based on other parameter values. The package will have a main function that, using the gibbs estimates of these parameters produced by the helper functions at each iteration, will provide estimates of the values of gene-gene interactions and gene-covariate interactions.

## Installation instructions

You can install ReciprocalGraphicalModels R package from GitHub with:

    install.packages("devtools")

    devtools::install_github("bitansa/ReciprocalGraphicalModels")

Once the ReciprocalGraphicalModels library is installed load the library in the R workspace.

    library("ReciprocalGraphicalModels")

## Future works

I created all the helper functions needed to produce gibbs samplers of model parameters. Using these helper functions, I have to construct the main function that will estimate the gene-gene and gene-covariate interactions based on gene expressions and covariate data. Then, in an effort to reduce runtime, I'll check if I can enhance the computations. Finally, for large datasets, the procedure may take a while because it does a Gibbs sampling across a number of repetitions. In order to significantly reduce the runtime, I shall write the programs in RCpp.I'll complete the documentation, add some datasets and also provide some applications of the R-package on the datasets as well.
