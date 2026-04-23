#######################################################################################################################
#########   BEAGLE 4.0 Phasing          ########################################################
#######################################################################################################################

#Everything has been generated, pedigree files have been generated, quality control has been done. Now to phase!

#We will phase with and without pedigree information to see what difference this makes to phasing
#We are phasing by chromosome:ensures biological accuracy, improves computational efficiency, 
#                             and avoids confusing the phasing algorithm with unrelated recombination patterns across different chromosomes.

workingDir <- "/home/jana/github/lstrachan_patrilines/"
pathToBeagle <- "/home/jana/bin/"

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
dir.create("Outputs/Beagle_phasing", showWarnings = FALSE)

#PHASING WITHOUT PEDIGREE
for (i in 1:16) {
  beagle_cmd <- paste0(
    "java -jar ", pathToBeagle, "/beagle.4.0.jar gt=Data/Real_data/Slov_fM_AC_QC_noDupPos.vcf.gz out=Outputs/Beagle_phasing/Slov_PHASED_chr", 
    i, " window=200 overlap=50 phase-its=10 impute-its=10 burnin-its=10 chrom=", i
  )
  system(beagle_cmd)
}

#TODO: Need to make the Slov_Alpha_pedigree using the chosen sires from AlphaAssign <----------------------------------
#PHASING WITH PEDIGREE
realPed = read.table("Outputs/AlphaAssign/Alpha_pedigree_Real.txt", header = F, stringsAsFactors = F)
realPedBeagle = realPed[, c(1, 1, 2, 3)]
write.table(realPedBeagle, file = "Outputs/AlphaAssign/Alpha_pedigree_Real_Beagle.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = F)
for (i in 1:16) {
  beagle_cmd <- paste0(
    "java -jar ", pathToBeagle, "/beagle.4.0.jar gt=Data/Real_data/Slov_fM_AC_QC_noDupPos.vcf.gz ped=Outputs/AlphaAssign/Alpha_pedigree_Real_Beagle.txt out=Outputs/Beagle_phasing/Slov_PHASED_pedigree_chr", 
    i, " window=200 overlap=50 phase-its=10 impute-its=10 burnin-its=10 chrom=", i
  )
  system(beagle_cmd)
}


########### STOPPED HERE
# Use PLINK to convert phased VCFs to .map files
vcf_files <- list.files(pattern = "Slov_PHASED_chr.*\\.vcf\\.gz$")
for (f in vcf_files) {
  base <- sub("\\.vcf\\.gz$", "", f)
  plink_cmd <- paste0("./plink --vcf ", f, " --recode --double-id --allow-extra-chr --out ", base)
  system(plink_cmd)
}

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

for (n in 1:nSNP_array){
#PHASING WITHOUT PEDIGREE - NoGE
for (i in 1:16) {
  beagle_cmd <- paste0(
    "java -jar ", pathToBeagle, "/beagle.4.0.jar gt=Data/Sim_NoGE/SNP_", n, "_NoGE_QC_ACformat.vcf out=Outputs/Beagle_phasing/SNP_", n, "_NoGE_PHASED_chr", 
    i, " window=200 overlap=50 phase-its=10 impute-its=10 burnin-its=10 chrom=", i
  )
  system(beagle_cmd)
}
}

#TODO: Need to make the Slov_Alpha_pedigree using the chosen sires from AlphaAssign <----------------------------------
#PHASING WITH PEDIGREE
for (i in 1:16) {
  beagle_cmd <- paste0(
    "java -jar beagle.4.0.jar gt=Slov_fM_AC_QC_noDupPos.vcf.gz ped=Slov_Alpha_pedigree.txt out=Slov_PHASED_chr", 
    i, " window=200 overlap=50 phase-its=200 impute-its=1000 burnin-its=200 chrom=", i
  )
  system(beagle_cmd)
}

# Use PLINK to convert phased VCFs to .map files
vcf_files <- list.files(pattern = "Slov_PHASED_chr.*\\.vcf\\.gz$")
for (f in vcf_files) {
  base <- sub("\\.vcf\\.gz$", "", f)
  plink_cmd <- paste0("./plink --vcf ", f, " --recode --double-id --allow-extra-chr --out ", base)
  system(plink_cmd)
}

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




  