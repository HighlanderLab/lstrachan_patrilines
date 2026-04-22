#######################################################################################################################
#########   BEAGLE 4.0 Phasing          ########################################################
#######################################################################################################################

#Everything has been generated, pedigree files have been generated, quality control has been done. Now to phase!

#We will phase with and without pedigree information to see what difference this makes to phasing
#We are phasing by chromosome:ensures biological accuracy, improves computational efficiency, 
#                             and avoids confusing the phasing algorithm with unrelated recombination patterns across different chromosomes.


#******** REAL DATA *********************************************************************************************

#TODO: set working directory- not sure where this will be stored yet 

#PHASING WITHOUT PEDIGREE
for (i in 1:16) {
  beagle_cmd <- paste0(
    "java -jar beagle.4.0.jar gt=Slov_fM_AC_QC_noDupPos.vcf.gz out=Slov_PHASED_chr", 
    i, " window=200 overlap=50 phase-its=200 impute-its=1000 burnin-its=200 chrom=", i
  )
  system(beagle_cmd)
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



#******** SIMULATED DATA *********************************************************************************************


#******** No Genotyping Errors *************************************************************************************
#TODO: set working directory- not sure where this will be stored yet 

nSPN_arrays <- 5

for (n in 1:nrow(nSNP_array)){
#PHASING WITHOUT PEDIGREE
for (i in 1:16) {
  beagle_cmd <- paste0(
    "java -jar beagle.4.0.jar gt=Slov_fM_AC_QC_noDupPos.vcf.gz out=Slov_PHASED_chr", 
    i, " window=200 overlap=50 phase-its=200 impute-its=1000 burnin-its=200 chrom=", i
  )
  system(beagle_cmd)
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




