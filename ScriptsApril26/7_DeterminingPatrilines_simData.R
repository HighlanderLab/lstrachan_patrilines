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

library(AlphaSimR)

args = commandArgs(trailingOnly=TRUE)
Rep = args[1]
workingDir = args[2]
softwareDir = args[3]
pathToPlink <- softwareDir
pathToBeagle <- softwareDir

repDir = paste0(workingDir, "/SimRep", Rep, "/")

# Set working directory
setwd(repDir)

source(paste0(workingDir, "/ScriptsApril26/7_DeterminingPatrilines_functions.R"))

################################################################################
#******* ROUTE 1 *******************************
################################################################################

#Load the outputs of 7_Haplotype_ParentageAssignments script
load("Pipeline/5_Haplotype_ParentAssignments.RData")

#••••• Simulated ••••• Can't do this with Real data since we don't have the father/drone info 

#Load your Pop object that contains the fathers and the pedigree with the fathers id in there 
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
sister_thresholds <- c(1.0, 0.95, 0.90, 0.85, 0.80, 0.75)
father_test_thresholds <- c(1,0.95, 0.9)

#Run 2 things: 
#1. Run with haplotypes that have gone through route1 of haplotype parental assignment
#2. Run with haplotypes that have gone through route2 of haplotype parental assignment 

# Haplotypes route 1 : with recon pedigree and both parents 
PatR1_SimTrue_Haplo1 <- run_paternity_tests(results_arg = Route1_SimTrue$real_results_flipped,
                                               sister_thresholds = sister_thresholds,
                                               father_test_thresholds = father_test_thresholds,
                                               pedigree = Sim_pedigree_mat,
                                               recon_pedigree = NULL,
                                               father_haplotypes = father_haplotype_pat,
                                            data_type = "True",
                                            haplo_assignment_type = "Used recon pedigree")


PatR1_NoGE_SNP2k_Haplo1 <- run_paternity_tests(results_arg = Route1_NoGE_SNP2k$results_flipped,
                                     sister_thresholds = sister_thresholds,
                                     father_test_thresholds = father_test_thresholds,
                                     pedigree = Sim_pedigree_mat,
                                     recon_pedigree = Rec_pedigree_2k_NoGE_filtered,
                                     father_haplotypes = father_haplotype_pat,
                                     data_type = "NoGE_2k",
                                     haplo_assignment_type = "Used recon pedigree")


PatR1_NoGE_SNP50k_Haplo1 <- run_paternity_tests(results_arg = Route1_NoGE_SNP50k$results_flipped,
                                               sister_thresholds = sister_thresholds,
                                               father_test_thresholds = father_test_thresholds,
                                               pedigree = Sim_pedigree_mat,
                                               recon_pedigree = Rec_pedigree_50k_NoGE_filtered,
                                               father_haplotypes = father_haplotype_pat,
                                               data_type = "NoGE_50k",
                                               haplo_assignment_type = "Used recon pedigree")

PatR1_WithGE_SNP2k_Haplo1 <- run_paternity_tests(results_arg = Route1_WithGE_SNP2k$results_flipped,
                                               sister_thresholds = sister_thresholds,
                                               father_test_thresholds = father_test_thresholds,
                                               pedigree = Sim_pedigree_mat,
                                               recon_pedigree = Rec_pedigree_2k_WithGE_filtered,
                                               father_haplotypes = father_haplotype_pat,
                                               data_type = "WithGE_SNP2k",
                                               haplo_assignment_type = "Used recon pedigree")


PatR1_WithGE_SNP50k_Haplo1 <- run_paternity_tests(results_arg = Route1_WithGE_SNP50k$results_flipped,
                                                sister_thresholds = sister_thresholds,
                                                father_test_thresholds = father_test_thresholds,
                                                pedigree = Sim_pedigree_mat,
                                                recon_pedigree = Rec_pedigree_50k_WithGE_filtered,
                                                father_haplotypes = father_haplotype_pat,
                                                data_type = "WithGE_50k",
                                                haplo_assignment_type = "Used recon pedigree")

# Haplotypes route 2 : with maternal pedigree and only maternal info 
PatR1_SimTrue_Haplo2 <- run_paternity_tests(results_arg = Route2_SimTrue$real_results_flipped,
                                            sister_thresholds = sister_thresholds,
                                            father_test_thresholds = father_test_thresholds,
                                            pedigree = Sim_pedigree_mat,
                                            recon_pedigree = NULL,
                                            father_haplotypes = father_haplotype_pat,
                                            data_type = "True",
                                            haplo_assignment_type = "Used mat pedigree")


PatR1_NoGE_SNP2k_Haplo2 <- run_paternity_tests(results_arg = Route2_NoGE_SNP2k$phased_results_flipped,
                                               sister_thresholds = sister_thresholds,
                                               father_test_thresholds = father_test_thresholds,
                                               pedigree = Sim_pedigree_mat,
                                               recon_pedigree = NULL,
                                               father_haplotypes = father_haplotype_pat,
                                               data_type = "NoGE_2k",
                                               haplo_assignment_type = "Used mat pedigree")


PatR1_NoGE_SNP50k_Haplo2 <- run_paternity_tests(results_arg = Route2_NoGE_SNP50k$phased_results_flipped,
                                                sister_thresholds = sister_thresholds,
                                                father_test_thresholds = father_test_thresholds,
                                                pedigree = Sim_pedigree_mat,
                                                recon_pedigree = NULL,
                                                father_haplotypes = father_haplotype_pat,
                                                data_type = "NoGE_50k",
                                                haplo_assignment_type = "Used mat pedigree")

PatR1_WithGE_SNP2k_Haplo2 <- run_paternity_tests(results_arg = Route2_WithGE_SNP2k$phased_results_flipped,
                                                 sister_thresholds = sister_thresholds,
                                                 father_test_thresholds = father_test_thresholds,
                                                 pedigree = Sim_pedigree_mat,
                                                 recon_pedigree = NULL,
                                                 father_haplotypes = father_haplotype_pat,
                                                 data_type = "WithGE_SNP2k",
                                                 haplo_assignment_type = "Used mat pedigree")


PatR1_WithGE_SNP50k_Haplo2 <- run_paternity_tests(results_arg = Route2_WithGE_SNP50k$phased_results_flipped,
                                                  sister_thresholds = sister_thresholds,
                                                  father_test_thresholds = father_test_thresholds,
                                                  pedigree = Sim_pedigree_mat,
                                                  recon_pedigree = NULL,
                                                  father_haplotypes = father_haplotype_pat,
                                                  data_type = "WithGE_SNP50k",
                                                  haplo_assignment_type = "Used mat pedigree")




# plot_paternity_number_grid_DRONES(PatR1_SimTrue_Haplo1)
# plot_paternity_number_grid_DRONES(PatR1_SimTrue_Haplo2)
#Do this for the others too if it works 

################################################################################
#******* ROUTE 2 *******************************
################################################################################

#••• SIMULATED ••• 

#Run 2 things: 
#1. Run with haplotypes that have gone through route1 of haplotype parental assignment
#2. Run with haplotypes that have gone through route2 of haplotype parental assignment 

# Haplotypes route 1 : with recon pedigree and both parents 

sister_thresholds <- c(1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.7)

PatR2_SimTrue_Haplo1 <- run_sister_clustering_tests(results_arg = Route1_SimTrue$real_results_flipped,
                                                       sister_thresholds = sister_thresholds,
                                                       pedigree = Sim_pedigree_mat,
                                                       recon_pedigree = NULL,
                                                       simulated = TRUE,
                                                    data_type = "True",
                                                    haplo_assignment_type = "Used recon pedigree")

PatR2_NoGE_SNP2k_Haplo1 <- run_sister_clustering_tests(results_arg = Route1_NoGE_SNP2k$results_flipped,
                                                    sister_thresholds = sister_thresholds,
                                                    pedigree = Sim_pedigree_mat,
                                                    recon_pedigree = Rec_pedigree_2k_NoGE_filtered,
                                                    simulated = TRUE,
                                                    data_type = "NoGE_SNP2k",
                                                    haplo_assignment_type = "Used recon pedigree")

PatR2_NoGE_SNP50k_Haplo1 <- run_sister_clustering_tests(results_arg = Route1_NoGE_SNP50k$results_flipped,
                                                       sister_thresholds = sister_thresholds,
                                                       pedigree = Sim_pedigree_mat,
                                                       recon_pedigree = Rec_pedigree_50k_NoGE_filtered,
                                                       simulated = TRUE,
                                                       data_type = "NoGE_SNP50k",
                                                       haplo_assignment_type = "Used recon pedigree")

PatR2_WithGE_SNP2k_Haplo1 <- run_sister_clustering_tests(results_arg = Route1_WithGE_SNP2k$results_flipped,
                                                       sister_thresholds = sister_thresholds,
                                                       pedigree = Sim_pedigree_mat,
                                                       recon_pedigree = Rec_pedigree_2k_WithGE_filtered,
                                                       simulated = TRUE,
                                                       data_type = "WithGE_SNP2k",
                                                       haplo_assignment_type = "Used recon pedigree")

PatR2_WithGE_SNP50k_Haplo1 <- run_sister_clustering_tests(results_arg = Route1_WithGE_SNP50k$results_flipped,
                                                        sister_thresholds = sister_thresholds,
                                                        pedigree = Sim_pedigree_mat,
                                                        recon_pedigree = Rec_pedigree_50k_WithGE_filtered,
                                                        simulated = TRUE,
                                                        data_type = "WithGE_SNP50k",
                                                        haplo_assignment_type = "Used recon pedigree")

# Haplotypes route 2 : with maternal pedigree and only maternal info 
PatR2_SimTrue_Haplo2 <- run_sister_clustering_tests(results_arg = Route2_SimTrue$real_results_flipped,
                                                    sister_thresholds = sister_thresholds,
                                                    pedigree = Sim_pedigree_mat,
                                                    recon_pedigree = NULL,
                                                    simulated = TRUE,
                                                    data_type = "True",
                                                    haplo_assignment_type = "Used mat pedigree")

PatR2_NoGE_SNP2k_Haplo2 <- run_sister_clustering_tests(results_arg = Route2_NoGE_SNP2k$phased_results_flipped,
                                                       sister_thresholds = sister_thresholds,
                                                       pedigree = Sim_pedigree_mat,
                                                       recon_pedigree = NULL,
                                                       simulated = TRUE,
                                                       data_type = "NoGE_2k",
                                                       haplo_assignment_type = "Used mat pedigree")

PatR2_NoGE_SNP50k_Haplo2 <- run_sister_clustering_tests(results_arg = Route2_NoGE_SNP50k$phased_results_flipped,
                                                        sister_thresholds = sister_thresholds,
                                                        pedigree = Sim_pedigree_mat,
                                                        recon_pedigree = NULL,
                                                        simulated = TRUE,
                                                        data_type = "NoGE_50k",
                                                        haplo_assignment_type = "Used mat pedigree")

PatR2_WithGE_SNP2k_Haplo2 <- run_sister_clustering_tests(results_arg = Route2_WithGE_SNP2k$phased_results_flipped,
                                                         sister_thresholds = sister_thresholds,
                                                         pedigree = Sim_pedigree_mat,
                                                         recon_pedigree = NULL,
                                                         simulated = TRUE,
                                                         data_type = "WithGE_2k",
                                                         haplo_assignment_type = "Used mat pedigree")


PatR2_WithGE_SNP50k_Haplo2 <- run_sister_clustering_tests(results_arg = Route2_WithGE_SNP50k$phased_results_flipped,
                                                          sister_thresholds = sister_thresholds,
                                                          pedigree = Sim_pedigree_mat,
                                                          recon_pedigree = NULL,
                                                          simulated = TRUE,
                                                          data_type = "WithGE_50k",
                                                          haplo_assignment_type = "Used mat pedigree")

#Do this for the others too if it works 

# plot_paternity_number_grid_SISTERONLY(
#   Dataset1 = PatR2_SimTrue_Haplo1,
#   Dataset2 = PatR2_NoGE_SNP2k_Haplo1,
#   Dataset3= PatR2_WithGE_SNP2k_Haplo1)
# 
# plot_paternity_number_grid_SISTERONLY(
#   Dataset1 = PatR2_SimTrue_Haplo2,
#   Dataset2 = PatR2_NoGE_SNP2k_Haplo2,
#   Dataset3= PatR2_WithGE_SNP2k_Haplo2)




save.image(file = "Pipeline/7_DeterminingPatrilines.Rdata")



