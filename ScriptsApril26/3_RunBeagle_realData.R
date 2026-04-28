#######################################################################################################################
#########   BEAGLE 4.0 Phasing          ########################################################
#######################################################################################################################

#Everything has been generated, pedigree files have been generated, quality control has been done. Now to phase!

#We will phase with and without pedigree information to see what difference this makes to phasing
#We are phasing by chromosome:ensures biological accuracy, improves computational efficiency, 
#                             and avoids confusing the phasing algorithm with unrelated recombination patterns across different chromosomes.

rm(list=ls())


args = commandArgs(trailingOnly=TRUE)
workingDir = args[2]
softwareDir = args[3]
pathToPlink <- softwareDir 
pathToBeagle <- softwareDir


# Function to convert allele codes
convert_alleles <- function(allele) {
  if (allele == "1") {
    return("A")
  } else if (allele == "2") {
    return("C")
  } else if (allele == "9"){
    return("0")
    
  } else {
    return(allele)
} }



convert_ped_genotypes <- function(ped_data, GE = NULL) {
  # Applying the conversion function to each allele in the genotype columns
  genotype_data <- ped_data[, 7:ncol(ped_data)]

  # Convert each allele using the sapply function
  corrected_genotype_data <- apply(genotype_data, 2, function(column) sapply(column, convert_alleles))

  
  # Replace the original genotype data with the converted data
  ped_data[, 7:ncol(ped_data)] <- corrected_genotype_data
  
  return(ped_data)
}





setwd(workingDir)
setwd("Real_data")
dir.create("Beagle_phasing", showWarnings = FALSE)

#******** REAL DATA *********************************************************************************************


#Need to modify the ped files for BEAGLE (genotypes need to being the 0ACGT format not 0,1,2,9) 
print("Converting to AC format for Beagle")
realPed = read.table("Slov_fM_QC.ped", header = F, stringsAsFactors = F)
realPed_AC <- convert_ped_genotypes(ped_data = realPed)

write.table(realPed_AC, file = paste0("Slov_fM_QC_ACformat.ped"), quote = FALSE, sep = " ", row.names = FALSE, col.names = F)
system(paste0("cp Slov_fM_QC.map Slov_fM_QC_ACformat.map")) 

system(paste0(pathToPlink,"/plink --file Slov_fM_QC_ACformat  --recode vcf --out Slov_fM_QC_ACformat"))   #make a vcf of the QC for phasing
system("bcftools sort Slov_fM_QC_ACformat.vcf -Oz -o Slov_fM_QC_ACformat_sorted.vcf.gz")
system("tabix -p vcf Slov_fM_QC_ACformat_sorted.vcf.gz")

# # Step 3.5: Check for duplicate positions
system("bcftools query -f '%CHROM\\t%POS\\n' Slov_fM_QC_ACformat_sorted.vcf.gz | sort | uniq -d")

# # Step 3.6: Normalize and remove duplicate records
system("bcftools norm -m -any Slov_fM_QC_ACformat_sorted.vcf.gz -Oz -o Slov_fM_QC_ACformat_sorted_biallelic.vcf.gz")
system("tabix -p vcf Slov_fM_QC_ACformat_sorted_biallelic.vcf.gz")
system("bcftools norm -d all Slov_fM_QC_ACformat_sorted_biallelic.vcf.gz -Oz -o Slov_fM_QC_ACformat_sorted_biallelic_noDupPos.vcf.gz")
system("tabix -p vcf Slov_fM_QC_ACformat_sorted_biallelic_noDupPos.vcf.gz")



#PHASING With Mat PEDIGREE
for (i in 1:16) {

  system("awk -F, 'NR>2 {print $1,$1,$2,$3}' Real_Data_pedigree.csv > Real_Data_pedigree_Beagle.csv")
  
  beagle_cmd <- paste0(
    "java -jar ", pathToBeagle, "/beagle.4.0.jar gt=Slov_fM_QC_ACformat_sorted_biallelic_noDupPos.vcf.gz ped=Real_Data_pedigree_Beagle.csv out=Beagle_phasing/Slov_PHASED_matPedigree_chr", 
    i, " phase-its=20 impute-its=20 burnin-its=20 chrom=", i
  )
  #window=200 overlap=50 
  system(beagle_cmd)
}

#TODO: Need to make the Slov_Alpha_pedigree using the chosen sires from AlphaAssign <----------------------------------
#PHASING WITH full reconstructed PEDIGREE
system("awk '{print $1,$1,$2,$3}' AlphaAssign/Alpha_pedigree_Real.txt > AlphaAssign/Alpha_pedigree_Real_Beagle.txt")

for (i in 1:16) {
  beagle_cmd <- paste0(
    "java -jar ", pathToBeagle, "/beagle.4.0.jar gt=Slov_fM_QC_ACformat_sorted_biallelic_noDupPos.vcf.gz ped=AlphaAssign/Alpha_pedigree_Real_Beagle.txt out=Beagle_phasing/Slov_PHASED_recPedigree_chr", 
    i, " phase-its=20 impute-its=20 burnin-its=20 chrom=", i
  )
  #window=200 overlap=50 
  system(beagle_cmd)
}


# Use PLINK to convert phased VCFs to .map files
setwd("Beagle_phasing")
vcf_files <- list.files(pattern = "Slov_PHASED.*\\.vcf\\.gz$")
for (f in vcf_files) {
 base <- sub("\\.vcf\\.gz$", "", f)
 plink_cmd <- paste0(pathToPlink, "/plink --vcf ", f, " --recode --double-id --allow-extra-chr --out ", base)
 system(plink_cmd)
}##

#Compress and index all uncompressed VCFs
vcf_plain <- list.files(pattern = "\\.vcf$")
for (f in vcf_plain) {
 if (!file.exists(paste0(f, ".gz"))) {
   system(paste("bgzip", f))
 }
}
vcf_gz <- list.files(pattern = "\\.vcf\\.gz$")
for (f in vcf_gz) {
 if (!file.exists(paste0(f, ".tbi"))) {
   system(paste("tabix -p vcf", f))
 }
}


