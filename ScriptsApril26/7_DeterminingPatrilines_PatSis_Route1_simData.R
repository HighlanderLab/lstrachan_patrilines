#######################################################################################################################
#************************************* Patriline determination - PatR2_NoGE_SNP2k_Haplo1 *******************************
#######################################################################################################################
# Route 1: Using mother and father info - Haplo1
# Running 2K, noGE

rm(list = ls())
library(AlphaSimR)


args = commandArgs(trailingOnly=TRUE)
Rep = args[1]
workingDir = args[2]
trueData = as.logical(args[3])
GEType_input = args[4]
SNPArray_input = args[5]

repDir = paste0(workingDir, "/SimRep", Rep, "/")
setwd(repDir)

source(paste0(workingDir, "/ScriptsApril26/7_DeterminingPatrilines_functions.R"))
load("Pipeline/5_Haplotype_ParentAssignments.RData")

GEType <- GEType_input
SNPArray <- SNPArray_input

print("After loading data")
print(paste0("GEtype: ", GEType))
print(paste0("SNPArray: ", SNPArray))

Sim_pedigree_mat <- read.csv("Data/worker_pedigree.csv")
sister_thresholds <- c(1.0, 0.95, 0.90, 0.85, 0.75)

if (trueData == FALSE) {
  name <- paste0("PatSis_", GEType, "_SNP", SNPArray, "_Route1")
  Step5_name <- paste0("Route1_", GEType, "_SNP", SNPArray)
  recon_pedigree <- get(paste0("Rec_pedigree_", SNPArray, "_", GEType, "_filtered"))
  data_type <- paste0(GEType, "_", SNPArray)
  results_arg <- get(Step5_name)$results_flipped
} else if (trueData == TRUE) {
  name <- "PatSis_SimTrue_Route1"
  Step5_name <- "Route1_SimTrue"
  recon_pedigree <- NULL
  data_type <- "True"
  results_arg <- get(Step5_name)$real_results_flipped
}

print(paste0("Running ", name))

result <- run_sister_clustering_tests(results_arg = results_arg,
                                                    sister_thresholds = sister_thresholds,
                                                    pedigree = Sim_pedigree_mat,
                                                    recon_pedigree = recon_pedigree,
                                                    simulated = TRUE,
                                                    data_type = data_type,
                                                    haplo_assignment_type = "Used recon pedigree")

save(result, file = paste0("Pipeline/7_", name, ".Rdata"))
print(paste0("Done: ", name))
