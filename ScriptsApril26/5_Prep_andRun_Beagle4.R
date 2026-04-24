#######################################################################################################################
#########   BEAGLE 4.0 Phasing          ########################################################
#######################################################################################################################

#Everything has been generated, pedigree files have been generated, quality control has been done. Now to phase!

#We will phase with and without pedigree information to see what difference this makes to phasing
#We are phasing by chromosome:ensures biological accuracy, improves computational efficiency, 
#                             and avoids confusing the phasing algorithm with unrelated recombination patterns across different chromosomes.

rm(list=ls())
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




workingDir <- "/home/jana/github/lstrachan_patrilines/"
pathToBeagle <- "/home/jana/bin/"
pathToPlink <- "/home/jana/bin/"
setwd(workingDir)
dir.create("Outputs/Beagle_phasing", showWarnings = FALSE)

nSNP_array <- rbind(c(3125, 5L), c(108, 4L), c(50, 3L), c(10, 2L), c(1,1L))

# First, convert the PED after QC to AC format
for (n in 1:nrow(nSNP_array)){
  #Need to modify the ped files for BEAGLE (genotypes need to being the 0ACGT format not 0,1,2,9) 
  print("Converting to AC format for Beagle")
  # NoGE
  ped_NoGE = read.table(paste0("Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE_QC.ped"), header = F, stringsAsFactors = F)
  ped_NoGE_AC <- convert_ped_genotypes(ped_data = ped_NoGE)

  ped_WithGE = read.table(paste0("Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE_QC.ped"), header = F, stringsAsFactors = F)
  ped_WithGE_AC <- convert_ped_genotypes(ped_data = ped_WithGE)

  write.table(ped_NoGE_AC, file = paste0("Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE_QC_ACformat.ped"), quote = FALSE, sep = " ", row.names = FALSE, col.names = F)
  system(paste0("cp Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE_QC.map Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE_QC_ACformat.map")) 

  write.table(ped_WithGE_AC, file = paste0("Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE_QC_ACformat.ped"), quote = FALSE, sep = " ", row.names = FALSE, col.names = F)
  system(paste0("cp Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE_QC.map Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE_QC_ACformat.map"))


  system(paste0(pathToPlink,"/plink --file Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE_QC_ACformat  --recode vcf --out Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE_QC_ACformat"))   #make a vcf of the QC for phasing
  system(paste0(pathToPlink,"/plink --file Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE_QC_ACformat --recode vcf --out Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE_QC_ACformat"))   #make a vcf of the QC for phasing
}

#******** REAL DATA *********************************************************************************************


setwd(workingDir)

#Need to modify the ped files for BEAGLE (genotypes need to being the 0ACGT format not 0,1,2,9) 
print("Converting to AC format for Beagle")
realPed = read.table("Data/Real_data/Slov_fM_QC.ped", header = F, stringsAsFactors = F)
realPed_AC <- convert_ped_genotypes(ped_data = realPed)

write.table(realPed_AC, file = paste0("Data/Real_data/Slov_fM_QC_ACformat.ped"), quote = FALSE, sep = " ", row.names = FALSE, col.names = F)
system(paste0("cp Data/Real_data/Slov_fM_QC.map Data/Real_data/Slov_fM_QC_ACformat.map")) 

system(paste0(pathToPlink,"/plink --file Data/Real_data/Slov_fM_QC_ACformat  --recode vcf --out Data/Real_data/Slov_fM_QC_ACformat"))   #make a vcf of the QC for phasing
system("bcftools sort Data/Real_data/Slov_fM_QC_ACformat.vcf -Oz -o Data/Real_data/Slov_fM_QC_ACformat_sorted.vcf.gz")
system("tabix -p vcf Data/Real_data/Slov_fM_QC_ACformat_sorted.vcf.gz")

# # Step 3.5: Check for duplicate positions
system("bcftools query -f '%CHROM\\t%POS\\n' Data/Real_data/Slov_fM_QC_ACformat_sorted.vcf.gz | sort | uniq -d")

# # Step 3.6: Normalize and remove duplicate records
system("bcftools norm -m -any Data/Real_data/Slov_fM_QC_ACformat_sorted.vcf.gz -Oz -o Data/Real_data/Slov_fM_QC_ACformat_sorted_biallelic.vcf.gz")
system("tabix -p vcf Data/Real_data/Slov_fM_QC_ACformat_sorted_biallelic.vcf.gz")
system("bcftools norm -d all Data/Real_data/Slov_fM_QC_ACformat_sorted_biallelic.vcf.gz -Oz -o Data/Real_data/Slov_fM_QC_ACformat_sorted_biallelic_noDupPos.vcf.gz")
system("tabix -p vcf Data/Real_data/Slov_fM_QC_ACformat_sorted_biallelic_noDupPos.vcf.gz")



#PHASING With Mat PEDIGREE
for (i in 1:16) {

  system("awk -F, 'NR>2 {print $1,$1,$2,$3}' Data/Real_data/Real_Data_pedigree.csv > Data/Real_data/Real_Data_pedigree_Beagle.csv")
  
  beagle_cmd <- paste0(
    "java -jar ", pathToBeagle, "/beagle.4.0.jar gt=Data/Real_data/Slov_fM_QC_ACformat_sorted_biallelic_noDupPos.vcf.gz ped=Data/Real_data/Real_Data_pedigree_Beagle.csv out=Outputs/Beagle_phasing/Slov_PHASED_matPedigree_chr", 
    i, " phase-its=20 impute-its=20 burnin-its=20 chrom=", i
  )
  #window=200 overlap=50 
  system(beagle_cmd)
}

#TODO: Need to make the Slov_Alpha_pedigree using the chosen sires from AlphaAssign <----------------------------------
#PHASING WITH full reconstructed PEDIGREE
system("awk '{print $1,$1,$2,$3}' Outputs/AlphaAssign/Alpha_pedigree_Real.txt > Outputs/AlphaAssign/Alpha_pedigree_Real_Beagle.txt")

for (i in 1:16) {
  beagle_cmd <- paste0(
    "java -jar ", pathToBeagle, "/beagle.4.0.jar gt=Data/Real_data/Slov_fM_QC_ACformat_sorted_biallelic_noDupPos.vcf.gz ped=Outputs/AlphaAssign/Alpha_pedigree_Real_Beagle.txt out=Outputs/Beagle_phasing/Slov_PHASED_recPedigree_chr", 
    i, " phase-its=20 impute-its=20 burnin-its=20 chrom=", i
  )
  #window=200 overlap=50 
  system(beagle_cmd)
}


# Use PLINK to convert phased VCFs to .map files
setwd("Outputs/Beagle_phasing")
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


#******** SIMULATED DATA *********************************************************************************************


#******** No Genotyping Errors *************************************************************************************
#TODO: set working directory- not sure where this will be stored yet 

nSNP_array <- c(4,5)
setwd(workingDir)

#PHASING With maternal true PEDIGREE
for (n in nSNP_array){
  for (i in 1:16) {

    system("awk  '{print \"1_\"$1,\"1_\"$1,\"1_\"$2,\"1_\"$3}' Data/SimData_Pedigree_Full_Maternal.csv > Data/SimData_Pedigree_Full_Maternal_Beagle.txt")
    system(paste0("sed -i 's/1_0/0/g' Data/SimData_Pedigree_Full_Maternal_Beagle.csv"))
    beagle_cmd_noGE <- paste0(
      "java -jar ", pathToBeagle, "/beagle.4.0.jar gt=Data/Sim_NoGE/SNP_", n, "_NoGE_QC_ACformat.vcf ped=Data/SimData_Pedigree_Full_Maternal_Beagle.txt out=Outputs/Beagle_phasing/SNP_", n, "_NoGE_PHASED_matPedigree_chr", 
      i, " phase-its=20 impute-its=20 burnin-its=20 chrom=", i
    )
    #window=200 overlap=50 
    system(beagle_cmd_noGE)

    
    beagle_cmd_withGE <- paste0(
      "java -jar ", pathToBeagle, "/beagle.4.0.jar gt=Data/Sim_WithGE/SNP_", n, "_WithGE_QC_ACformat.vcf ped=Data/SimData_Pedigree_Full_Maternal_Beagle.txt out=Outputs/Beagle_phasing/SNP_", n, "_WithGE_PHASED_matPedigree_chr", 
      i, " phase-its=20 impute-its=20 burnin-its=20 chrom=", i
    )
    #window=200 overlap=50 
    system(beagle_cmd_withGE)
  }
}


#PHASING WITH PEDIGREE
for (n in nSNP_array){
  for (i in 1:16) {
    snp_size = c("4" = "2k", "5" = "50k")

    system(paste0("awk   '{print \"1_\"$1,\"1_\"$1,\"1_\"$2,\"1_\"$3}' Outputs/AlphaAssign/Alpha_pedigree_", snp_size[as.character(n)], "_NoGE.txt > Outputs/AlphaAssign/Alpha_pedigree_", snp_size[as.character(n)], "_NoGE_Beagle.txt"))
    system(paste0("sed -i 's/1_0/0/g' Outputs/AlphaAssign/Alpha_pedigree_", snp_size[as.character(n)], "_NoGE_Beagle.txt"))

    beagle_cmd_noGE <- paste0(
      "java -jar ", pathToBeagle, "/beagle.4.0.jar gt=Data/Sim_NoGE/SNP_", n, "_NoGE_QC_ACformat.vcf ped=Outputs/AlphaAssign/Alpha_pedigree_", snp_size[as.character(n)], "_NoGE_Beagle.txt out=Outputs/Beagle_phasing/SNP_", n, "_NoGE_PHASED_recPedigree_chr", 
      i, " phase-its=20 impute-its=20 burnin-its=20 chrom=", i
    )
    system(beagle_cmd_noGE)



    system(paste0("awk   '{print \"1_\"$1,\"1_\"$1,\"1_\"$2,\"1_\"$3}' Outputs/AlphaAssign/Alpha_pedigree_", snp_size[as.character(n)], "_WithGE.txt > Outputs/AlphaAssign/Alpha_pedigree_", snp_size[as.character(n)], "_WithGE_Beagle.txt"))
    system(paste0("sed -i 's/1_0/0/g' Outputs/AlphaAssign/Alpha_pedigree_", snp_size[as.character(n)], "_WithGE_Beagle.txt"))

    beagle_cmd_withGE <- paste0(
      "java -jar ", pathToBeagle, "/beagle.4.0.jar gt=Data/Sim_WithGE/SNP_", n, "_WithGE_QC_ACformat.vcf ped=Outputs/AlphaAssign/Alpha_pedigree_", snp_size[as.character(n)], "_WithGE_Beagle.txt out=Outputs/Beagle_phasing/SNP_", n, "_WithGE_PHASED_recPedigree_chr", 
      i, " phase-its=20 impute-its=20 burnin-its=20 chrom=", i
    )
    system(beagle_cmd_withGE)
  }
}



# Use PLINK to convert phased VCFs to .map files
setwd("Outputs/Beagle_phasing")
vcf_files <- list.files(pattern = "SNP_.*\\GE_PHASED_.*\\chr.*\\.vcf\\.gz$")
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



  