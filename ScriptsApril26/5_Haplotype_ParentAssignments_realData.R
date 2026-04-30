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
workingDir = args[1]
softwareDir = args[2]
pathToPlink <- softwareDir 
pathToBeagle <- softwareDir

setwd(paste0(workingDir, "/Real_data/"))

load("Pipeline/4_Converting_PhasedVCF.Rdata")

#####################################################################################
#************************** FUNCTIONS ******************************************
#####################################################################################
source(paste0(workingDir, "/ScriptsApril26/5_Haplotype_ParentAssignments_functions.R"))

#####################################################################################
#**** Get the real data phased haplotypes ******
#####################################################################################
print("Reading in real data")
# Pedigree prior to ped reconstruction
cols <- c("id", "sire", "dam")
Slov_pedigree_mat <- read.table("Data/Real_Data_pedigree.txt")
colnames(Slov_pedigree_mat) <- cols

#After reconstruction
Slov_pedigree_rec <- read.table("Outputs/AlphaAssign/Alpha_pedigree_Real.txt")
colnames(Slov_pedigree_rec) <- cols

Slov_map <- read.table("Data/Slov_fM_QC_ACformat.map")

##########################################################3
print("Assigning haplotype PO for real data")


#Real data
Slov_pedigree_mat_filtered <- Slov_pedigree_mat[Slov_pedigree_mat$dam != 0,] # Remove rows with unknown mothers 
tmp <-sub("_.*", "", rownames(Slov_PhasedHaplotypes_matPed))
tmp <- unique(tmp)
Slov_pedigree_mat_filtered <- Slov_pedigree_mat_filtered[Slov_pedigree_mat_filtered$id %in% tmp,]

#Real data
Slov_pedigree_rec_filter <- Slov_pedigree_rec[Slov_pedigree_rec$sire != 0 & Slov_pedigree_rec$dam != 0,]

print("Route1")
Route1_Real <- Route1_flipping(perfect_haplotypes = FALSE, pedigree = Slov_pedigree_rec_filter, method = "power_mean", Data_type = "Real_Slov_data")
print("Route2")
Route2_Real <- Route2_flipping(perfect_haplotypes = FALSE, pedigree = Slov_pedigree_mat_filtered, method = "power_mean", Data_type = "Real_Slov_data")

print("Saving data")
save.image("Pipeline/5_Haplotype_ParentAssignments.RData")
