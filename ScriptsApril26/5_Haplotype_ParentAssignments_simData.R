#######################################################################################################################
#***************** Assign the Haplotype parent-of-origin assignments **************
#######################################################################################################################

#In this section we'll use 2 routes:

#.  - Route 1: use pedigree reconstruction and use both dam and sire pedigree ID to assign haplotype parental origins
#Use haplotypes from phasing using pedigree

#.  - Route 2: DO NOT use pedigree reconstruction and use ONLY dam pedigree ID to assign haplotype parental origins

#How does this function work? 
# - Iterates through all 16 chromosomes
# - Computes difference between the offspring's haplotypes (coded as 0/1) and each of the parents' genotypes (coded as 0/1/2).

#Offspring haplotype - Parent Genotype = difference 

#                      Parent Geno = 0. Parent Geno = 1. Parent Geno = 2
# Offspring Haplo = 0.       0                -1                -2
# Offspring Haplo = 1        1                 0                -1

#Difference of 1/-2 == allelic MISMATCH
#Difference  0/-1 == allelic MATCH

#Summarise difference across loci with a power mean of allelic mismatches Score = sum(MISMATCH)^2/length(MISMATCH)
#The higher the score the less likely the haplotype was inherited from that parent 
#The scores were checked using conditional arguments to improve the parent-of-origin assignment accuracy

###############################################################################################################
#Practical example:

#Offspring with haplotypes [0,1,0,1]
#dam with genotypes     [0,1,2,0]
#DPQ with genotypes        [1,0,1,1]

# Using Table above, we calculate mismatches for the dam:
#   0 - 0 = 0     --> MATCH
#   1 - 1 = 0     --> MATCH
#   0 - 2 = -2    --> MISMATCH
#   1 - 0 = 1     --> MISMATCH
# Resulting differences: [0, 0, -2, 1]
# Squaring these:        [0, 0,  4, 1]
# Power mean score:      (0 + 0 + 4 + 1) / 4 = 1.25
#
# Now for the father:
#   0 - 1 =  1     --> MISMATCH
#   1 - 0 = -1     --> MATCH
#   0 - 1 = -1     --> MATCH
#   1 - 1 =  0     --> MATCH
# Resulting differences: [1, -1, -1, 0]
# Squaring these:        [1,  1,  1, 0]
# Power mean score:      (1 + 1 + 1 + 0) / 4 = 0.75
#
# Conclusion:
# The lower score for the father (0.75 vs 1.25) suggests that the haplotype is paternally derived.

#######################################################################################################################
# --- Clear Workspace ---
rm(list = ls())
# --- Libraries ---
{
  library(Eagle)
  library(tidyr)
  library(AlphaSimR)
  library(SIMplyBee)
  library(readr)
  library(genio)
  library(ggplot2)
  library(dplyr)
  library(vcfR)
  library(tibble)
}


args = commandArgs(trailingOnly=TRUE)
Rep = args[1]
workingDir = args[2]
softwareDir = args[3]
pathToPlink <- softwareDir
pathToBeagle <- softwareDir

repDir = paste0(workingDir, "/SimRep", Rep, "/")

# Set working directory
setwd(repDir)

source(paste0(workingDir, "/ScriptsApril26/5_Haplotype_ParentAssignments_functions.R"))

#####################################################################################
#**** Get the true simulated haplotypes and phased haplotypes ******
#####################################################################################
print("Reading in simulated data")

setwd(repDir)

load("Data/SP_object.Rdata")
Simulated_SP <- SP
Simulated_pop <- load("Data/Pop_withFathers.Rdata")
true_haplotypes <- pullSnpHaplo(PopMerged)
true_map <- getGenMap(SP)

#Pedigree prior to reconstruction ####################################################

Worker_pedigree <- read.csv("Data/worker_pedigree.csv")

#Pedigree post reconstruction (with AlphaAssign) ########################################
Rec_pedigree_2k_NoGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_2k_NoGE.txt")
Rec_pedigree_50k_NoGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_50k_NoGE.txt")

Rec_pedigree_2k_WithGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_2k_WithGE.txt")
Rec_pedigree_50k_WithGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_50k_WithGE.txt")

cols <- c("id", "sire", "dam")
names(Rec_pedigree_2k_NoGE) <- cols
names(Rec_pedigree_50k_NoGE) <- cols
names(Rec_pedigree_2k_WithGE) <- cols
names(Rec_pedigree_50k_WithGE) <- cols


#Get haplotypes made in 6_Converting_PhasedVCF scripts
setwd(repDir)
load("Pipeline/4_Converting_PhasedVCF.Rdata")

NoGE_map_2k <- read.table("Data/Sim_NoGE/SNP_4_NoGE_QC_ACformat.map")
WithGE_map_2k <- read.table("Data/Sim_WithGE/SNP_4_WithGE_QC_ACformat.map")
NoGE_map_50k <- read.table("Data/Sim_NoGE/SNP_5_NoGE_QC_ACformat.map")
WithGE_map_50k <- read.table("Data/Sim_WithGE/SNP_5_WithGE_QC_ACformat.map")



#####################################################################################
#*** Check the different phasing versions compared to the real and prephased data just to check similarity ***
#####################################################################################
print("Checking accuracy of phasing for simulated data")
# ••• Simulated •••

print("Checking accuracy")
print("2K SNP data, rec ped, NoGE")
true_vs_NoGEphasedSNP2krecPed <- check_haplotype(true_haplotypes = true_haplotypes, results = NoGE_SNP2k_PhasedHaplotypes_recPed, pedigree = Worker_pedigree)
print("2K SNP data, mat ped, NoGE")
true_vs_NoGEphasedSNP2kmatPed <- check_haplotype(true_haplotypes = true_haplotypes, results = NoGE_SNP2k_PhasedHaplotypes_matPed, pedigree = Worker_pedigree)

print("2K SNP data, rec ped, WithGE")
true_vs_WithGEphasedSNP2krecPed <- check_haplotype(true_haplotypes = true_haplotypes, results = WithGE_SNP2k_PhasedHaplotypes_recPed, pedigree = Worker_pedigree)
print("2K SNP data, mat ped, WithGE")
true_vs_WithGEphasedSNP2kmatPed <- check_haplotype(true_haplotypes = true_haplotypes, results = WithGE_SNP2k_PhasedHaplotypes_matPed, pedigree = Worker_pedigree)

print("50K SNP data, rec ped, NoGE")
true_vs_NoGEphasedSNP50krecPed <- check_haplotype(true_haplotypes = true_haplotypes, results = NoGE_SNP50k_PhasedHaplotypes_recPed, pedigree = Worker_pedigree )
print("50K SNP data, mat ped, NoGE")
true_vs_NoGEphasedSNP50kmatPed <- check_haplotype(true_haplotypes = true_haplotypes, results = NoGE_SNP50k_PhasedHaplotypes_matPed, pedigree = Worker_pedigree)

print("50K SNP data, rec ped, WithGE")
true_vs_WithGEphasedSNP50krecPed <- check_haplotype(true_haplotypes = true_haplotypes, results = WithGE_SNP50k_PhasedHaplotypes_recPed, pedigree = Worker_pedigree )
print("50K SNP data, mat ped, WithGE")
true_vs_WithGEphasedSNP50kmatPed <- check_haplotype(true_haplotypes = true_haplotypes, results = WithGE_SNP50k_PhasedHaplotypes_matPed, pedigree = Worker_pedigree)


#####################################################################################
#**               Haplotype parent-of-origin assignments                         **

#**Route 1: use pedigree reconstruction's pedigree and use both dam and sire pedigree ID to assign haplotype parental origins **
#####################################################################################
print("Assigning haplotype PO for simulated data")
print("Route1")

Worker_pedigree <- Worker_pedigree[, c("id", "dpc", "mother")]

print("SimTrue")
#True haplotypes (should work perfectly)
Route1_SimTrue <- Route1_flipping(perfect_haplotypes = TRUE, pedigree = Worker_pedigree, method = "power_mean")


#Editing route 1 pedigrees - can't have sire's present 
Rec_pedigree_2k_NoGE_filtered <- Rec_pedigree_2k_NoGE[Rec_pedigree_2k_NoGE$sire != 0,]
Rec_pedigree_50k_NoGE_filtered <- Rec_pedigree_50k_NoGE[Rec_pedigree_50k_NoGE$sire != 0,]
Rec_pedigree_2k_WithGE_filtered <- Rec_pedigree_2k_WithGE[Rec_pedigree_2k_WithGE$sire != 0,]
Rec_pedigree_50k_WithGE_filtered <- Rec_pedigree_50k_WithGE[Rec_pedigree_50k_WithGE$sire != 0,]

print("2K SNP data")
#2k SNP
print("NoGE")
Route1_NoGE_SNP2k <- Route1_flipping(perfect_haplotypes = FALSE, pedigree = Rec_pedigree_2k_NoGE_filtered, method = "power_mean", Data_type = "NoGE_SNP2k")
print("WithGE")
Route1_WithGE_SNP2k <- Route1_flipping(perfect_haplotypes = FALSE, pedigree = Rec_pedigree_2k_WithGE_filtered, method = "power_mean", Data_type = "WithGE_SNP2k")

print("50K SNP data")
#50k SNP
print("NoGE")
Route1_NoGE_SNP50k <- Route1_flipping(perfect_haplotypes = FALSE, pedigree = Rec_pedigree_50k_NoGE_filtered, method = "power_mean", Data_type = "NoGE_SNP50k")
print("WithGE")
Route1_WithGE_SNP50k <- Route1_flipping(perfect_haplotypes = FALSE, pedigree = Rec_pedigree_50k_WithGE_filtered, method = "power_mean", Data_type = "WithGE_SNP50k")

save(Route1_SimTrue, Route1_NoGE_SNP2k, Route1_WithGE_SNP2k, Route1_NoGE_SNP50k, Route1_WithGE_SNP50k, file = paste0(repDir, "/Pipeline/5_Haplotype_ParentAssignments_Route1.RData"))

# #Checking haplotypes post flip to see if something has happened
# colnames(Route1_SimTrue$real_results_flipped) <- colnames(true_haplotypes)
# 
#Example:
# true_vs_NoGE_SNP2k_FLIPPED <- check_haplotype_postFlip(complete_haplotypes = Route1_SimTrue$real_results_flipped, results = Route1_NoGE_SNP2k$results_flipped, pedigree = Rec_pedigree_2k_NoGE)
# 

#####################################################################################
#** Route 2: DO NOT use pedigree reconstruction and use ONLY dam pedigree ID to assign haplotype parental origins **
#####################################################################################
print("Route2")

Worker_pedigree <- Worker_pedigree[, c("id", "dpc", "mother")]

print("SimTrue")
#True haplotypes (should work perfectly)
Route2_SimTrue <- Route2_flipping(perfect_haplotypes = TRUE, pedigree = Worker_pedigree, method = "power_mean", Data_type = "True")

#2k SNP
print("2K SNP data")
print("NoGE")
Route2_NoGE_SNP2k <- Route2_flipping(perfect_haplotypes = FALSE, pedigree = Worker_pedigree, method = "power_mean",Data_type = "NoGE_SNP2k")
print("WithGE")
Route2_WithGE_SNP2k <- Route2_flipping(perfect_haplotypes = FALSE, pedigree = Worker_pedigree, method = "power_mean", Data_type = "WithGE_SNP2k")

#50k SNP
print("50K SNP data")
print("NoGE")
Route2_NoGE_SNP50k <- Route2_flipping(perfect_haplotypes = FALSE, pedigree = Worker_pedigree, method = "power_mean", Data_type = "NoGE_SNP50k")
print("WithGE")
Route2_WithGE_SNP50k <- Route2_flipping(perfect_haplotypes = FALSE, pedigree = Worker_pedigree, method = "power_mean", Data_type = "WithGE_SNP50k")

print("Saving data")
save.image(file="TEST")
save(Route2_NoGE_SNP2k, Route2_WithGE_SNP2k, Route2_NoGE_SNP50k, Route2_WithGE_SNP50k, file = paste0(repDir, "/Pipeline/5_Haplotype_ParentAssignments_Route2.RData"))


print("Saving data")
save.image(paste0(repDir, "/Pipeline/5_Haplotype_ParentAssignments.RData"))

