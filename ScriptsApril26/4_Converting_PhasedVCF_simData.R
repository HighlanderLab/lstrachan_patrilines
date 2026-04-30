#################################################################################

#******* Converting VCF files after Beagle Phasing ************

#################################################################################
#Output of Beagle is a haplotype dataframe with columns (id1_1, id1_2) and rows the SNPs coded 0/1

#Get all of the vcf files of all phased chromosomes and convert into a more manageable format 

#************ FUNCTIONS *******************

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
Rep = args[1]
workingDir = args[2]
softwareDir = args[3]
pathToPlink <- softwareDir
pathToBeagle <- softwareDir

repDir = paste0(workingDir, "/SimRep", Rep, "/")

source(paste0(workingDir, "/ScriptsApril26/4_Converting_PhasedVCF_functions.R"))
#************* SIMULATED DATA ****************************
#Pedigree prior to reconstruction 
setwd(repDir)


#Pedigree post reconstruction (with AlphaAssign)
#Alpha_pedigree_4_NoGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_2k_NoGE.txt")
#Alpha_pedigree_5_NoGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_50k_NoGE.txt") 

#Alpha_pedigree_4_WithGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_2k_WithGE.txt")
#Alpha_pedigree_5_WithGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_50k_WithGE.txt")

setwd("Outputs/Beagle_phasing")

# Combine results across all chromosomes for each dataset
for (pedPhase in c("recPedigree", "matPedigree")) {
  for (GEType in c("NoGE", "WithGE")) {
    for (n in c(4, 5)) {
      print(paste("Processing", GEType, "SNP", n, pedPhase))
      results_list <- list()
      for (chr in 1:16) { 
        print(chr)
        vcf_gz <- paste0("SNP_", n, "_", GEType, "_PHASED_", pedPhase, "_chr", chr, ".vcf.gz")
        map <- paste0("SNP_", n, "_", GEType, "_PHASED_", pedPhase, "_chr", chr, ".map")
        df <- convert_VCF(vcf_file = vcf_gz, map_file = map)
        results_list[[chr]] <- df
      }
      combined_df <- do.call(cbind, results_list)
      assign(paste0(GEType, "_SNP", n, "_", pedPhase, "_AllChrs"), combined_df)
    }
  }
}


#Reconstructed pedigree
#2 = 2K, 5 = 50K
#Get the haplotypes of the phasing with maternal pedigree
print("Arranging haplotypes for maternal pedigree")
setwd(repDir)
NoGE_SNP2k_PhasedHaplotypes_matPed <- Haplotype_using_pedigree(GenErr = "NoGE", n = "4", ped_recon = TRUE, pedigree_name = "Data/SimData_Pedigree_Full_Maternal.txt", haplo_name = "NoGE_SNP4_matPedigree_AllChrs")
WithGE_SNP2k_PhasedHaplotypes_matPed <- Haplotype_using_pedigree(GenErr = "WithGE", n = "4", ped_recon = TRUE, pedigree_name = "Data/SimData_Pedigree_Full_Maternal.txt", haplo_name = "WithGE_SNP4_matPedigree_AllChrs")
NoGE_SNP50k_PhasedHaplotypes_matPed <- Haplotype_using_pedigree(GenErr = "NoGE", n = "5", ped_recon = TRUE, pedigree_name = "Data/SimData_Pedigree_Full_Maternal.txt", haplo_name = "NoGE_SNP5_matPedigree_AllChrs")
WithGE_SNP50k_PhasedHaplotypes_matPed <- Haplotype_using_pedigree(GenErr = "WithGE", n = "5", ped_recon = TRUE, pedigree_name = "Data/SimData_Pedigree_Full_Maternal.txt", haplo_name = "WithGE_SNP5_matPedigree_AllChrs")

print("Arranging haplotypes for reconstructed pedigree")
#Get the haplotypes phased with reconstructed pedugree
NoGE_SNP2k_PhasedHaplotypes_recPed <- Haplotype_using_pedigree(GenErr = "NoGE", n = "4", ped_recon = TRUE, pedigree_name = "Outputs/AlphaAssign/Alpha_pedigree_2k_NoGE.txt", haplo_name = "NoGE_SNP4_recPedigree_AllChrs")
WithGE_SNP2k_PhasedHaplotypes_recPed <- Haplotype_using_pedigree(GenErr = "WithGE", n = "4", ped_recon = TRUE, pedigree_name = "Outputs/AlphaAssign/Alpha_pedigree_2k_WithGE.txt", haplo_name = "WithGE_SNP4_recPedigree_AllChrs")
NoGE_SNP50k_PhasedHaplotypes_recPed <- Haplotype_using_pedigree(GenErr = "NoGE", n = "5", ped_recon = TRUE, pedigree_name = "Outputs/AlphaAssign/Alpha_pedigree_50k_NoGE.txt", haplo_name = "NoGE_SNP5_recPedigree_AllChrs")
WithGE_SNP50k_PhasedHaplotypes_recPed <- Haplotype_using_pedigree(GenErr = "WithGE", n = "5", ped_recon = TRUE, pedigree_name = "Outputs/AlphaAssign/Alpha_pedigree_50k_WithGE.txt", haplo_name = "WithGE_SNP5_recPedigree_AllChrs")




save.image(file = paste0(repDir, "/Pipeline/4_Converting_PhasedVCF.Rdata"))