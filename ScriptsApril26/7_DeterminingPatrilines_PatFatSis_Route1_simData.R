#######################################################################################################################
#************************************* Patriline determination - PatR2_NoGE_SNP2k_Haplo2 *******************************
#######################################################################################################################
# Route 2: Using mother info only - Haplo2
# Running 2K, NoGE

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
fathers_id_all <- Sim_pedigree_mat$father
father_id_1 <- paste(fathers_id_all, "1", sep = "_")

true_haplotypes <- pullSnpHaplo(PopMerged)

#since they're coded as identical diploids we can just use one of them
father_haplotypes_1 <- true_haplotypes[rownames(true_haplotypes) %in% father_id_1,]

#replace _1 with _paternal so that it can all be compared 
father_all_unique <- unique(fathers_id_all)
father_id_paternal <- father_all_unique[order(father_all_unique)]
father_id_paternal <- paste(father_id_paternal, "paternal", sep = "_")

father_haplotype_pat <- father_haplotypes_1
rownames(father_haplotype_pat) <- father_id_paternal

# Define the thresholds you want to test (takes a good while to run)
sister_thresholds <- c(1.0, 0.95, 0.90, 0.85, 0.75)
father_test_thresholds <- c(1,0.95, 0.9)

if (trueData == FALSE) {
  print("True data is false")
  name <- paste0("PatFatSis_", GEType, "_SNP", SNPArray, "_Route1")
  Step5_name <- paste0("Route1_", GEType, "_SNP", SNPArray)
  data_type <- paste0(GEType, "_", SNPArray)
  recon_pedigree <- get(paste0("Rec_pedigree_", SNPArray, "_", GEType, "_filtered"))
  results_arg <- get(Step5_name)$results_flipped
} else if (trueData == TRUE) {
  print("True data is true")
  name <- "PatFatSis_SimTrue_Route1"
  Step5_name <- "Route1_SimTrue"
  data_type <- "True"
  results_arg <- get(Step5_name)$real_results_flipped 
  colnames(results_arg) <- colnames(true_haplotypes)
  recon_pedigree <- NULL
}

print(paste0("Running ", name))

result <- run_paternity_tests(results_arg = results_arg,
                              sister_thresholds = sister_thresholds,
                              father_test_thresholds = father_test_thresholds,
                              pedigree = Sim_pedigree_mat,
                              recon_pedigree = recon_pedigree,
                              father_haplotypes = father_haplotype_pat,
                              data_type = data_type,
                              haplo_assignment_type = "Used recon pedigree")

save(result, file = paste0("Pipeline/7_", name, ".Rdata"))
print(paste0("Done: ", name))
