#################################################################################

#******* Converting VCF files after Beagle Phasing ************

#################################################################################
#Output of Beagle is a haplotype dataframe with columns (id1_1, id1_2) and rows the SNPs coded 0/1

#Get all of the vcf files of all phased chromosomes and convert into a more manageable format 

#************ FUNCTIONS *******************
rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
workingDir = args[1]
softwareDir = args[2]
pathToPlink <- softwareDir 
pathToBeagle <- softwareDir

source(paste0(workingDir, "/ScriptsApril26/4_Converting_PhasedVCF_functions.R"))
######################################################################################3
#********* REAL DATA ******************
print("Processing real data")
setwd(paste0(workingDir, "/Real_data"))
setwd("Outputs/Beagle_phasing")

#Phased with pedigree
results_list <- list()
for (pedPhase in c("recPedigree", "matPedigree")) {
  for (chr in 1:16){
    vcf_gz <- paste0("Slov_PHASED_", pedPhase, "_chr",chr,".vcf.gz")
    map <- paste0("Slov_PHASED_", pedPhase, "_chr",chr,".map")
    df <- convert_VCF(vcf_file = vcf_gz, map_file = map)
    results_list[[chr]] <- df
  }
  combined_df <- do.call(cbind, results_list)
  assign(paste0("Slov_PhasedHaplotypes_", pedPhase, "_AllChrs"), combined_df)
}

setwd(paste0(workingDir, "/Real_data"))
#Prephased file (used for comparison in the next script )
vcf_file = "Data/Slov_fM_QC_ACformat_sorted_biallelic_noDupPos.vcf.gz"
base <- sub("\\.vcf\\.gz$", "", vcf_file)
plink_cmd <- paste0(pathToPlink, "/plink --vcf ", vcf_file, " --recode --double-id --allow-extra-chr --out ", base)
system(plink_cmd)
map_file = "Data/Slov_fM_QC_ACformat_sorted_biallelic_noDupPos.map"
Slov_prephased_haplotypes <- convert_VCF(vcf_file = vcf_file, map_file = map_file)


real_data_pedigree = read.csv("Data/Real_Data_pedigree.csv", header = TRUE)
write.table(real_data_pedigree, file = "Data/Real_Data_pedigree.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

Slov_PhasedHaplotypes_matPed <- Haplotype_using_pedigree_realData(pedigree_name = "Data/Real_Data_pedigree.txt", haplo_name = "Slov_PhasedHaplotypes_matPedigree_AllChrs") #Slov_PhasedHaplotypes_NOPed
Slov_PhasedHaplotypes_recPed <- Haplotype_using_pedigree_realData(pedigree_name = "Outputs/AlphaAssign/Alpha_pedigree_Real.txt", haplo_name = "Slov_PhasedHaplotypes_recPedigree_AllChrs") #Slov_PhasedHaplotypes_WithPed

print("Saving data")
save(list = c("Slov_PhasedHaplotypes_matPed", "Slov_PhasedHaplotypes_recPed"),
 file = "Outputs/Beagle_phasing/All_Phased_Haplotypes.RData")


save.image(file = "Pipeline/4_Converting_PhasedVCF.Rdata")
