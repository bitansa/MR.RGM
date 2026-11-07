###############################################################################
## OneK1K_RunAll.R  —  the complete OneK1K analysis once the number of genes is known (one R session, this folder)
##
##   data_dir <- "path/to/folder/with/donor_b_cell.rda"
##   ONEK_OUT_DIR <- "path/to/results"     # optional (default: this folder)
##   ONEK_N_GENES <- m                     # largest supported gene set (../../package_check/check_design_matched.R)
##   source("OneK1K_RunAll.R")
##
## Steps, each with its model set explicitly here (no setting carries over from an earlier step or session):
##   1. OneK1K_PreProcessing.R            (skipped if its objects are already in the session)
##   2. OneK1K_ModelFitting.R             MR.RGM, then MR.RGM+ (four chains each)
##   3. OneK1K_Competitors.R              every competing method on the data of the MR.RGM fit
##      OneK1K_PriorSensitivity.R         MR.RGM under a uniform Beta(1, 1) prior on the causal inclusion probabilities
##      OneK1K_Bootstrap.R                MR.RGM refitted to data generated from the fitted network (parametric bootstrap)
##   4. OneK1K_Results.R                  MR.RGM, then MR.RGM+ (tables, figures, comparison tables)
## The gene count is chosen before this script: ONEK_STRENGTH_ONLY <- TRUE; source("OneK1K_ModelFitting.R") writes the
## candidate genes, and ../../package_check/check_design_matched.R gives the largest supported gene set.
###############################################################################

t_start <- Sys.time()
if (!exists("ONEK_N_GENES")) stop("Set ONEK_N_GENES (the largest supported gene set from check_design_matched.R) before sourcing.")
for (v in c("ONEK_STRENGTH_ONLY", "ONEK_MODEL")) if (exists(v, envir = globalenv())) rm(list = v, envir = globalenv())
step <- function(label) message(sprintf("\n========== [%s] %s ==========", format(Sys.time(), "%H:%M:%S"), label))

prep <- c("donor", "gene.ref", "vcf.ref", "gene.promoter.ref", "RNA.log", "genotype.mat", "Subsetted_SNP_Names")
if (!all(vapply(prep, exists, logical(1)))) { step("preprocessing"); source("OneK1K_PreProcessing.R") }

for (m in c("MR.RGM", "MR.RGM+")) {
  step(paste("fit", m)); assign("ONEK_MODEL", m, envir = globalenv()); source("OneK1K_ModelFitting.R")
}
rm("ONEK_MODEL", envir = globalenv())

step("competing methods"); source("OneK1K_Competitors.R")
step("prior sensitivity"); source("OneK1K_PriorSensitivity.R")
step("parametric bootstrap"); source("OneK1K_Bootstrap.R")

for (m in c("MR.RGM", "MR.RGM+")) {
  step(paste("results", m)); assign("ONEK_MODEL", m, envir = globalenv()); source("OneK1K_Results.R")
}
rm("ONEK_MODEL", envir = globalenv())
message(sprintf("\nOneK1K analysis finished in %.0f minutes.", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
