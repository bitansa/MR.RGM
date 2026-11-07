###############################################################################
## install_packages.R  --  R packages of the simulation study that the conda environment does not provide
## (called by setup.sbatch; installs into R_LIBS_USER; packages already installed are skipped)
###############################################################################
lib <- Sys.getenv("R_LIBS_USER"); stopifnot(nzchar(lib)); dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
have <- function(p) requireNamespace(p, quietly = TRUE)

cran <- c("Rcpp", "RcppArmadillo", "GIGrvg", "igraph", "pROC",
          "rjags", "ggplot2", "scales", "microbenchmark", "future", "remotes", "testthat",
          "gmm", "ivreg", "msm", "Formula", "lmtest", "sandwich")            # the last six: dependencies of OneSampleMR
need <- cran[!vapply(cran, have, logical(1))]
if (length(need)) install.packages(need, lib = lib)
## The competing methods at the versions of the study (versions of the GitHub and r-universe packages: README.md)
for (pv in list(c("MendelianRandomization", "0.10.0"), c("mrbayes", "0.5.3"))) if (!have(pv[1]))
  remotes::install_version(pv[1], version = pv[2], lib = lib, upgrade = "never", repos = "https://cloud.r-project.org")
## OneSampleMR is archived on CRAN: version 0.1.8 (the version of the simulation study; tsri() fitted by GMM, see
## fit_onesamplemr() in sim_functions.R) from the CRAN archive, or else from the authors' GitHub repository
if (!have("OneSampleMR")) tryCatch(
  remotes::install_version("OneSampleMR", version = "0.1.8", lib = lib, upgrade = "never", repos = "https://cloud.r-project.org"),
  error = function(e) { message("CRAN archive failed (", conditionMessage(e), "); installing from GitHub")
                        remotes::install_github("remlapmot/OneSampleMR", lib = lib, upgrade = "never") })
if (!have("MR2G")) remotes::install_github("ZhaotongL/MR2G", lib = lib, upgrade = "never", dependencies = TRUE)
for (p in c("MR2", "MrDAG")) if (!have(p))
  install.packages(p, lib = lib, repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
## GrIVET: R package grivet in the authors' repository, at the commit of the study; needs cvTools and mnormt
for (p in c("cvTools", "mnormt")) if (!have(p)) install.packages(p, lib = lib)
if (!have("grivet")) remotes::install_github("chunlinli/grivet@1777dbf", subdir = "grivet", lib = lib, upgrade = "never")

all_pk <- c(cran, "MendelianRandomization", "mrbayes", "OneSampleMR", "MR2G", "MR2", "MrDAG", "cvTools", "mnormt", "grivet")
miss <- all_pk[!vapply(all_pk, have, logical(1))]
if (length(miss)) stop("not installed: ", paste(miss, collapse = ", "))
message("All R packages installed. OneSampleMR ", packageVersion("OneSampleMR"))
