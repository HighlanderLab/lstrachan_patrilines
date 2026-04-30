#######################################################################################################################
#************************************* Patriline determination *******************************
#######################################################################################################################
# Estimates the number of patrilines by comparing the paternally assigned haplotypes workers

#Again I'm taking two routes here: 

# - ROUTE 1: Using actual drone information using the true haplotypes to test the accuracy of the threshold method
#            - checks the paternity accuracy really closely but is unrealistic unless you have father-drone information. 
#            - uses sister thresholds and father_thresholds to visualise the most accurate thresholds required

# - ROUTE 2: Use the paternally assigned haplotypes from this script and ONLY sister thresholds to determine patriline numbers
#           - again using multiple sister thresholds to determine which of the Slov real data results are the most reliable. 

################################################################################
#******* Source Functions *******************************
rm(list = ls())
# --- Libraries ---

args = commandArgs(trailingOnly=TRUE)
workingDir = args[1]
softwareDir = args[2]
pathToPlink <- softwareDir
pathToBeagle <- softwareDir


# Set working directory
setwd(paste0(workingDir, "/Real_data/"))

source(paste0(workingDir, "/ScriptsApril26/7_DeterminingPatrilines_functions.R"))
################################################################################
#******* ROUTE 1 *******************************
################################################################################

#Load the outputs of 7_Haplotype_ParentageAssignments script
load("Pipeline/5_Haplotype_ParentAssignments.RData")

sister_thresholds <- c(1.0, 0.95, 0.90, 0.85, 0.80, 0.75)
father_test_thresholds <- c(1,0.95, 0.9)

#** ••• REAL ••• **

colnames(Slov_pedigree_mat_filtered) <- c("id", "dpc", "mother")
colnames(Slov_pedigree_rec_filter) <- c("id", "dpc", "mother")

#This is route 2 for determining patrilines with haplotypes coming from route 1 of everything before
PatR2_Real_Haplo1 <- run_sister_clustering_tests(results_arg = Route1_Real$results_flipped,
                                                 sister_thresholds = sister_thresholds,
                                                 pedigree = Slov_pedigree_mat_filtered,
                                                 recon_pedigree = Slov_pedigree_rec_filter,
                                                 simulated = FALSE,
                                                 data_type = "Real",
                                                 haplo_assignment_type = "Used recon pedigree")

#This is route 2 for determining patrilines with haplotypes coming from route 2 of everything before
PatR2_Real_Haplo2 <- run_sister_clustering_tests(results_arg = Route2_Real$phased_results_flipped,
                                                 sister_thresholds = sister_thresholds,
                                                 pedigree = Slov_pedigree_mat_filtered,
                                                 recon_pedigree = NULL,
                                                 simulated = FALSE,
                                                 data_type = "Real",
                                                 haplo_assignment_type = "Used mat pedigree")



save.image(file = "Pipeline/7_DeterminingPatrilines.Rdata")



