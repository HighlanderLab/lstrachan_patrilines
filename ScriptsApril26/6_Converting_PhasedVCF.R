#################################################################################

#******* Converting VCF files after Beagle Phasing ************

#################################################################################
#Output of Beagle is a haplotype dataframe with columns (id1_1, id1_2) and rows the SNPs coded 0/1

#Get all of the vcf files of all phased chromosomes and convert into a more manageable format 

#************ FUNCTIONS *******************

workingDir <- "/home/jana/github/lstrachan_patrilines/"
pathToBeagle <- "/home/jana/bin/"
pathToPlink <- "/home/jana/bin/"

get_out_haplotypes <- function(ped_matrix, ind_id_1, ind_id_2, realData = FALSE) {
  
  # Initialize matrices to hold haplotype data
  haplotype_matrix_0 <- matrix(nrow = nrow(ped_matrix), ncol = ncol(ped_matrix))
  haplotype_matrix_1 <- matrix(nrow = nrow(ped_matrix), ncol = ncol(ped_matrix))
  
  # Fill the matrices with haplotype data
  for (i in 1:nrow(ped_matrix)) {
    for (haplotype in 0:1) {
      if (haplotype == 0) {
        haplotype_matrix_0[i, ] <- sapply(ped_matrix[i, ], function(gt) {
          haplotypes <- unlist(strsplit(as.character(gt), " "))
          return(haplotypes[haplotype + 1])
        })
      } else {
        haplotype_matrix_1[i, ] <- sapply(ped_matrix[i, ], function(gt) {
          haplotypes <- unlist(strsplit(as.character(gt), " "))
          return(haplotypes[haplotype + 1])
        })
      }
    }
  }
  
  # Convert matrices to dataframes and set row and column names
  colnames(haplotype_matrix_0) <- colnames(ped_matrix)
  rownames(haplotype_matrix_0) <- ind_id_1
  
  colnames(haplotype_matrix_1) <- colnames(ped_matrix)
  rownames(haplotype_matrix_1) <- ind_id_2
  
  
  # Combine haplotype matrices into a list
  haplotype_matrices <- rbind(haplotype_matrix_0, haplotype_matrix_1)

  #order them
  if (realData) {
    haplotype_matrices <- order_by_prefix_realData(haplotype_matrices)
  } else {
    haplotype_matrices <- order_by_prefix(haplotype_matrices)
  }
}

convert_genotypes <- function(genotypes) {
  #genotypes[genotypes == '1'] <- 0
  genotypes[genotypes == '0'] <- NA
  #genotypes[genotypes == '2'] <- 1
  genotypes[genotypes == 'A'] <- 0
  genotypes[genotypes == 'C'] <- 1
  
  return(as.numeric(genotypes))
}

convert_VCF <- function(vcf_file = NULL, map_file = NULL){
  
  chr_map <- read.table(map_file)
  phasedVCF=vcf_file
  
  system(paste0('bcftools query -f "%REF %ALT [ %GT]\n" ', phasedVCF, ' > PhasedGT.txt'))
  system(paste0("bcftools query -l ", phasedVCF, " > SampleIDs.txt"))
  
  df = read.table("PhasedGT.txt")
  colnames(df) = c("REF", "ALT", paste0("Ind", 1:(ncol(df)-2)))
  
  refAlt=df[,1:2]     
  gt=df[3:ncol(df)]
  
  gt_t = t(gt)
  
  for (col in 1:ncol(gt_t)) {
    ref = refAlt$REF[col]
    alt = refAlt$ALT[col]
    gt_t[,col] = gsub("0", ref, gt_t[,col])
    gt_t[,col] = gsub("1", alt, gt_t[,col])
    gt_t[,col] = gsub("\\|", " ", gt_t[,col])
  }
  
  #ped order is always mum/dpc/workers
  ids = read.table("SampleIDs.txt")
  ids$V1 <- gsub("AMEL_", "", ids$V1)
  rownames(gt_t) <- ids$V1
  
  #get true haplotypes for colnames
  colnames(gt_t) <- chr_map$V2
  
  return(gt_t)
}


Haplotype_using_pedigree <- function(GenErr = NULL, n = NULL, ped_recon = NULL, pedigree_name = NULL, haplo_name = NULL) {
  
  # Validate inputs early
 
 if (is.null(GenErr) || is.null(n) || is.null(ped_recon)) {
    stop("GenErr, n, and ped_recon must all be provided.")
  }
  
  # Get pedigree
  pedigree <- read.table(pedigree_name)
  colnames(pedigree) <- c("id", "dpc", "mother")
  pedigree = pedigree[pedigree$mother != 0,]
  
  # Get haplotypes
  if (!exists(haplo_name)) {
    stop(paste("Object", haplo_name, "not found"))
  }
  All_chroms <- get(haplo_name)
  rownames(All_chroms) <- gsub("1_", "", rownames(All_chroms))
  
  # Subset haplotypes
  Sim_mother_haplo  <- All_chroms[rownames(All_chroms) %in% pedigree$mother, , drop = FALSE]
  Sim_dpc_haplo     <- All_chroms[rownames(All_chroms) %in% pedigree$dpc, , drop = FALSE]
  Sim_workers_haplo <- All_chroms[rownames(All_chroms) %in% pedigree$id, , drop = FALSE]
  
  Sim_all_haplo <- rbind(Sim_mother_haplo, Sim_dpc_haplo, Sim_workers_haplo)
  
  # IDs
  mother_ids <- unique(pedigree$mother)
  dpc_ids    <- setdiff(unique(pedigree$dpc), "0")
  worker_ids <- unique(pedigree$id)
  
  all_ids <- unique(c(mother_ids, dpc_ids, worker_ids))
  
  ind_id_1 <- paste0(all_ids, "_1")
  ind_id_2 <- paste0(all_ids, "_2")
  
  # Process haplotypes
  tmp  <- get_out_haplotypes(Sim_all_haplo, ind_id_1 = ind_id_1, ind_id_2 = ind_id_2, realData = FALSE)
  tmp2 <- apply(tmp, 2, convert_genotypes)
  
  rownames(tmp2) <- rownames(tmp)
  
  return(tmp2)
}

Haplotype_using_pedigree_realData <- function(pedigree_name = NULL, haplo_name = NULL) {
  
  # Get pedigree
  pedigree <- read.table(pedigree_name)
  colnames(pedigree) <- c("id", "dpc", "mother")
  pedigree = pedigree[pedigree$mother != 0,]
  
  # Get haplotypes
  if (!exists(haplo_name)) {
    stop(paste("Object", haplo_name, "not found"))
  }
  All_chroms <- get(haplo_name)
  rownames(All_chroms) <- gsub("1_", "", rownames(All_chroms))
  pedigree = pedigree[pedigree$id %in% rownames(All_chroms),]
  
  # Subset haplotypes
  Sim_mother_haplo  <- All_chroms[rownames(All_chroms) %in% pedigree$mother, , drop = FALSE]
  Sim_dpc_haplo     <- All_chroms[rownames(All_chroms) %in% pedigree$dpc, , drop = FALSE]
  Sim_workers_haplo <- All_chroms[rownames(All_chroms) %in% pedigree$id, , drop = FALSE]
  
  Sim_all_haplo <- rbind(Sim_mother_haplo, Sim_dpc_haplo, Sim_workers_haplo)
  
  # IDs
  mother_ids <- unique(pedigree$mother)
  dpc_ids    <- setdiff(unique(pedigree$dpc), "0")
  worker_ids <- unique(pedigree$id)
  
  all_ids <- unique(c(mother_ids, dpc_ids, worker_ids))
  
  ind_id_1 <- paste0(all_ids, "_1")
  ind_id_2 <- paste0(all_ids, "_2")
  
  # Process haplotypes
  tmp  <- get_out_haplotypes(Sim_all_haplo, ind_id_1 = ind_id_1, ind_id_2 = ind_id_2, realData = TRUE)
  tmp2 <- apply(tmp, 2, convert_genotypes)
  
  rownames(tmp2) <- rownames(tmp)
  
  return(tmp2)
}

order_by_prefix <- function(df) {
  # Extract the part before the underscore
  prefix <- as.numeric(sapply(rownames(df), function(x) strsplit(x, "_")[[1]][1]))
  
  # Order the data frame based on the extracted prefix
  df_ordered <- df[order(prefix, rownames(df)), , drop = FALSE]
  
  return(df_ordered)
}

order_by_prefix_realData <- function(df) {
  # Extract the part before the underscore
  prefix <- sapply(rownames(df), function(x) strsplit(x, "_")[[1]][1])
  
  # Order the data frame based on the extracted prefix
  df_ordered <- df[order(prefix, rownames(df)), , drop = FALSE]
  
  return(df_ordered)
}

#************* SIMULATED DATA ****************************
#Pedigree prior to reconstruction 
setwd(workingDir)


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
      for (chr in 1:1) { #16
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
setwd(workingDir)
NoGE_SNP2k_PhasedHaplotypes_matPed <- Haplotype_using_pedigree(GenErr = "NoGE", n = "4", ped_recon = TRUE, pedigree_name = "Data/SimData_Pedigree_Full_Maternal.txt", haplo_name = "NoGE_SNP4_matPedigree_AllChrs")
WithGE_SNP2k_PhasedHaplotypes_matPed <- Haplotype_using_pedigree(GenErr = "WithGE", n = "4", ped_recon = TRUE, pedigree_name = "Data/SimData_Pedigree_Full_Maternal.txt", haplo_name = "WithGE_SNP4_matPedigree_AllChrs")
NoGE_SNP50k_PhasedHaplotypes_matPed <- Haplotype_using_pedigree(GenErr = "NoGE", n = "5", ped_recon = TRUE, pedigree_name = "Data/SimData_Pedigree_Full_Maternal.txt", haplo_name = "NoGE_SNP5_matPedigree_AllChrs")
WithGE_SNP50k_PhasedHaplotypes_matPed <- Haplotype_using_pedigree(GenErr = "WithGE", n = "5", ped_recon = TRUE, pedigree_name = "Data/SimData_Pedigree_Full_Maternal.txt", haplo_name = "WithGE_SNP5_matPedigree_AllChrs")

#Get the haplotypes phased with reconstructed pedugree
NoGE_SNP2k_PhasedHaplotypes_recPed <- Haplotype_using_pedigree(GenErr = "NoGE", n = "4", ped_recon = TRUE, pedigree_name = "Outputs/AlphaAssign/Alpha_pedigree_2k_NoGE.txt", haplo_name = "NoGE_SNP4_recPedigree_AllChrs")
WithGE_SNP2k_PhasedHaplotypes_recPed <- Haplotype_using_pedigree(GenErr = "WithGE", n = "4", ped_recon = TRUE, pedigree_name = "Outputs/AlphaAssign/Alpha_pedigree_2k_WithGE.txt", haplo_name = "WithGE_SNP4_recPedigree_AllChrs")
NoGE_SNP50k_PhasedHaplotypes_recPed <- Haplotype_using_pedigree(GenErr = "NoGE", n = "5", ped_recon = TRUE, pedigree_name = "Outputs/AlphaAssign/Alpha_pedigree_50k_NoGE.txt", haplo_name = "NoGE_SNP5_recPedigree_AllChrs")
WithGE_SNP50k_PhasedHaplotypes_recPed <- Haplotype_using_pedigree(GenErr = "WithGE", n = "5", ped_recon = TRUE, pedigree_name = "Outputs/AlphaAssign/Alpha_pedigree_50k_WithGE.txt", haplo_name = "WithGE_SNP5_recPedigree_AllChrs")



######################################################################################3
#********* REAL DATA ******************

setwd("Outputs/Beagle_phasing")

#Phased with pedigree
results_list <- list()
for (pedPhase in c("recPedigree", "matPedigree")) {
  for (chr in 1:1){
    vcf_gz <- paste0("Slov_PHASED_", pedPhase, "_chr",chr,".vcf.gz")
    map <- paste0("Slov_PHASED_", pedPhase, "_chr",chr,".map")
    df <- convert_VCF(vcf_file = vcf_gz, map_file = map)
    results_list[[chr]] <- df
  }
  combined_df <- do.call(cbind, results_list)
  assign(paste0("Slov_PhasedHaplotypes_", pedPhase, "_AllChrs"), combined_df)
}

setwd(workingDir)
#Prephased file (used for comparison in the next script )
vcf_file = "Data/Real_data/Slov_fM_QC_ACformat_sorted_biallelic_noDupPos.vcf.gz"
base <- sub("\\.vcf\\.gz$", "", vcf_file)
plink_cmd <- paste0(pathToPlink, "/plink --vcf ", vcf_file, " --recode --double-id --allow-extra-chr --out ", base)
system(plink_cmd)
map_file = "Data/Real_data/Slov_fM_QC_ACformat_sorted_biallelic_noDupPos.map"
Slov_prephased_haplotypes <- convert_VCF(vcf_file = vcf_file, map_file = map_file)


Slov_PhasedHaplotypes_matPed <- Haplotype_using_pedigree_realData(pedigree_name = "Data/Real_data/Real_Data_pedigree.txt", haplo_name = "Slov_PhasedHaplotypes_matPedigree_AllChrs") #Slov_PhasedHaplotypes_NOPed
Slov_PhasedHaplotypes_recPed <- Haplotype_using_pedigree_realData(pedigree_name = "Outputs/AlphaAssign/Alpha_pedigree_Real.txt", haplo_name = "Slov_PhasedHaplotypes_recPedigree_AllChrs") #Slov_PhasedHaplotypes_WithPed


save(list = c("Slov_PhasedHaplotypes_matPed", "Slov_PhasedHaplotypes_recPed",
              "NoGE_SNP2k_PhasedHaplotypes_matPed", "WithGE_SNP2k_PhasedHaplotypes_matPed", "NoGE_SNP50k_PhasedHaplotypes_matPed", "WithGE_SNP50k_PhasedHaplotypes_matPed", 
              "NoGE_SNP2k_PhasedHaplotypes_recPed", "WithGE_SNP2k_PhasedHaplotypes_recPed", "NoGE_SNP50k_PhasedHaplotypes_recPed", "WithGE_SNP50k_PhasedHaplotypes_recPed"),
            file = "Outputs/Beagle_phasing/All_Phased_Haplotypes.RData")

save.image(file = paste0(workingDir, "/Data/Pipeline/6_Converting_PhasedVCF.Rdata"))