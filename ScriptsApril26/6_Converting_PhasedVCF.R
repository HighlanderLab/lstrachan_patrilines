#################################################################################

#******* Converting VCF files after Beagle Phasing ************

#################################################################################
#Output of Beagle is a haplotype dataframe with columns (id1_1, id1_2) and rows the SNPs coded 0/1

#Get all of the vcf files of all phased chromosomes and convert into a more manageable format 

#************ FUNCTIONS *******************
convert_VCF <- function(vcf_file = NULL, map_file = NULL){
  
  chr_map <- read.table(map_file)
  phasedVCF=vcf_file
  
  system(paste0('bcftools query -f "%REF %ALT [ %GT]\n" ', phasedVCF, ' > PhasedGT.txt'))
  
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
  rownames(gt_t) <- ids_all
  
  #get true haplotypes for colnames
  colnames(gt_t) <- chr_map$V2
  
  return(gt_t)
}
get_out_haplotypes <- function(ped_matrix, ind_id_1, ind_id_2) {
  
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
  haplotype_matrices <- order_by_prefix(haplotype_matrices)
  
  return(haplotype_matrices)
}
convert_genotypes <- function(genotypes) {
  genotypes[genotypes == '1'] <- 0
  genotypes[genotypes == '0'] <- NA
  genotypes[genotypes == '2'] <- 1
  genotypes[genotypes == 'A'] <- 0
  genotypes[genotypes == 'C'] <- 1
  
  return(as.numeric(genotypes))
}
convert_VCF_Slov <- function(vcf_file = NULL, map_file = NULL){
  
  chr_map <- read.table(map_file)
  phasedVCF=vcf_file
  
  system(paste0('bcftools query -f "%REF %ALT [ %GT]\n" ', phasedVCF, ' > PhasedGT.txt'))
  
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
  rownames(gt_t) <- Slov_ids_preRecon
  
  #get true haplotypes for colnames
  colnames(gt_t) <- chr_map$V2
  
  return(gt_t)
}
Haplotype_using_pedigree <- function(GenErr = NULL, n = NULL, ped_recon = NULL) {
  
  # Validate inputs early
  if (is.null(GenErr) || is.null(n) || is.null(ped_recon)) {
    stop("GenErr, n, and ped_recon must all be provided.")
  }
  
  # Get pedigree
  if (ped_recon) {
    pedigree_name <- paste0("Alpha_pedigree_", n, "_", GenErr)
    if (!exists(pedigree_name)) {
      stop(paste("Object", pedigree_name, "not found"))
    }
    pedigree <- get(pedigree_name)
  } else {
    if (!exists("Sim_pedigree_pre")) {
      stop("Sim_pedigree_pre not found")
    }
    pedigree <- Sim_pedigree_pre
  }
  
  # Get haplotypes
  haplo_name <- paste0(GenErr, "_SNP", n, "_NoPED_Allchroms")
  if (!exists(haplo_name)) {
    stop(paste("Object", haplo_name, "not found"))
  }
  All_chroms <- get(haplo_name)
  
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
  tmp  <- get_out_haplotypes(Sim_all_haplo, ind_id_1 = ind_id_1, ind_id_2 = ind_id_2)
  tmp2 <- apply(tmp, 2, convert_genotypes)
  
  rownames(tmp2) <- rownames(tmp)
  
  return(tmp2)
}



#************* SIMULATED DATA ****************************
#Pedigree prior to reconstruction 
Sim_pedigree_pre <- read.csv("/Data/worker_pedigree.csv")

#Pedigree post reconstruction (with AlphaAssign)
Alpha_pedigree_2k_NoGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_2k_NoGE.txt")
Alpha_pedigree_50k_NoGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_50k_NoGE.txt")

Alpha_pedigree_2k_WithGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_2k_WithGE.txt")
Alpha_pedigree_50k_WithGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_50k_WithGE.txt")

#Chromosomes phased with pedigree
{
#nGE_SNPn_WithPED_ChrN loop
results_list <- list()
for (n in 1:16){
  vcf_gz <- load("nGE_SNP2k_WithPED_Chr",n,"PHASED.vcf.gz")
  map <- load("nGE_SNP2k_WithPED_Chr",n,"PHASED.map")
  df <- convert_VCF(vcf_file = vcf_gz, map_file = map)
  results_list[[n]] <- df
}
NoGE_SNP2k_WithPED_Allchroms <- do.call(cbind, results_list)

#WithGE_SNP2k_WithPED_ChrN loop
results_list <- list()
for (n in 1:16){
  vcf_gz <- load("WithGE_SNP2k_WithPED_Chr",n,"PHASED.vcf.gz")
  map <- load("WithGE_SNP2k_WithPED_Chr",n,"PHASED.map")
  df <- convert_VCF(vcf_file = vcf_gz, map_file = map)
  results_list[[n]] <- df
}
WithGE_SNP2k_WithPED_Allchroms <- do.call(cbind, results_list)
#NoGE_SNP50k_WithPED_ChrN loop

results_list <- list()
for (n in 1:16){
  vcf_gz <- load("nGE_SNP50k_WithPED_Chr",n,"PHASED.vcf.gz")
  map <- load("nGE_SNP50k_WithPED_Chr",n,"PHASED.map")
  df <- convert_VCF(vcf_file = vcf_gz, map_file = map)
  results_list[[n]] <- df
}
NoGE_SNP50k_WithPED_Allchroms <- do.call(cbind, results_list)

#WithGE_SNP50k_WithPED_ChrN loop
results_list <- list()
for (n in 1:16){
  vcf_gz <- load("WithGE_SNP50k_WithPED_Chr",n,"PHASED.vcf.gz")
  map <- load("WithGE_SNP50k_WithPED_Chr",n,"PHASED.map")
  df <- convert_VCF(vcf_file = vcf_gz, map_file = map)
  results_list[[n]] <- df
}
WithGE_SNP50k_WithPED_Allchroms <- do.call(cbind, results_list)
}

#Chromosomes phased with NO pedigree
{
  #nGE_SNPn_NoPED_ChrN loop
  results_list <- list()
  for (n in 1:16){
    vcf_gz <- load("nGE_SNP2k_NoPED_Chr",n,"PHASED.vcf.gz")
    map <- load("nGE_SNP2k_NoPED_Chr",n,"PHASED.map")
    df <- convert_VCF(vcf_file = vcf_gz, map_file = map)
    results_list[[n]] <- df
  }
  NoGE_SNP2k_NoPED_Allchroms <- do.call(cbind, results_list)
  
  #WithGE_SNP2k_NoPED_ChrN loop
  results_list <- list()
  for (n in 1:16){
    vcf_gz <- load("WithGE_SNP2k_NoPED_Chr",n,"PHASED.vcf.gz")
    map <- load("WithGE_SNP2k_NoPED_Chr",n,"PHASED.map")
    df <- convert_VCF(vcf_file = vcf_gz, map_file = map)
    results_list[[n]] <- df
  }
  WithGE_SNP2k_NoPED_Allchroms <- do.call(cbind, results_list)
  #NoGE_SNP50k_NoPED_ChrN loop
  
  results_list <- list()
  for (n in 1:16){
    vcf_gz <- load("nGE_SNP50k_NoPED_Chr",n,"PHASED.vcf.gz")
    map <- load("nGE_SNP50k_NoPED_Chr",n,"PHASED.map")
    df <- convert_VCF(vcf_file = vcf_gz, map_file = map)
    results_list[[n]] <- df
  }
  NoGE_SNP50k_NoPED_Allchroms <- do.call(cbind, results_list)
  
  #WithGE_SNP50k_NoPED_ChrN loop
  results_list <- list()
  for (n in 1:16){
    vcf_gz <- load("WithGE_SNP50k_NoPED_Chr",n,"PHASED.vcf.gz")
    map <- load("WithGE_SNP50k_NoPED_Chr",n,"PHASED.map")
    df <- convert_VCF(vcf_file = vcf_gz, map_file = map)
    results_list[[n]] <- df
  }
  WithGE_SNP50k_NoPED_Allchroms <- do.call(cbind, results_list)
}


#Get the haplotypes of the pedigree made by the pedigree reconstruction
NoGE_SNP2k_PhasedHaplotypes_WithPed <- Haplotype_using_pedigree(GenErr = "NoGE", n = "2k", ped_recon = TRUE)
WithGE_SNP2k_PhasedHaplotypes_WithPed <- Haplotype_using_pedigree(GenErr = "WithGE", n = "2k", ped_recon = TRUE)
NoGE_SNP50k_PhasedHaplotypes_WithPed <- Haplotype_using_pedigree(GenErr = "NoGE", n = "50k", ped_recon = TRUE)
WithGE_SNP50k_PhasedHaplotypes_WithPed <- Haplotype_using_pedigree(GenErr = "WithGE", n = "50k", ped_recon = TRUE)

#Get the haplotypes of the pedigree prior to pedigree reconstruction
NoGE_SNP2k_PhasedHaplotypes_NoPed <- Haplotype_using_pedigree(GenErr = "NoGE", n = "2k", ped_recon = FALSE)
WithGE_SNP2k_PhasedHaplotypes_NoPed <- Haplotype_using_pedigree(GenErr = "WithGE", n = "2k", ped_recon = FALSE)
NoGE_SNP50k_PhasedHaplotypes_NoPed <- Haplotype_using_pedigree(GenErr = "NoGE", n = "50k", ped_recon = FALSE)
WithGE_SNP50k_PhasedHaplotypes_NoPed <- Haplotype_using_pedigree(GenErr = "WithGE", n = "50k", ped_recon = FALSE)

#********* REAL DATA ******************

  setwd("PLACE WHERE THE REAL PHASED DATA IS STORED")

  # Pedigree prior to ped reconstruction
  Slov_pedigree_pre <- read.table("~/Desktop/lstrachan_patrilines/Data/Real_data/AlphaAssign/Pedigree.txt")
  colnames(Slov_pedigree_pre) <- c("id","dpc","mother")
  
  #After reconstruction
  Slov_pedigree_post <- read.table("Outputs/AlphaAssign/Alpha_pedigree_Real.txt")
  colnames(Slov_pedigree_post) <- c("id","dpc","mother")
  
  #Phased with pedigree
  results_list <- list()
  for (n in 1:16){
    vcf_gz <- load("Slov_PHASED_noPED_Chr",n,".vcf.gz")
    map <- load("Slov_PHASED_noPED_Chr",n,".map")
    df <- convert_VCF_Slov(vcf_file = vcf_gz, map_file = map)
    results_list[[n]] <- df
  }
  Slov_allchroms_phased_NoPed <- do.call(cbind, results_list)
  
  #Phased with no pedigree
  results_list <- list()
  for (n in 1:16){
    vcf_gz <- load("Slov_PHASED_WithPED_Chr",n,".vcf.gz")
    map <- load("Slov_PHASED_WithPED_Chr",n,".map")
    df <- convert_VCF_Slov(vcf_file = vcf_gz, map_file = map)
    results_list[[n]] <- df
  }
  Slov_allchroms_phased_WithPed <- do.call(cbind, results_list)
  
  
  
  #Prephased file (used for comparison in the next script )
  Slov_prephased_haplotypes <- convert_VCF_Slov(vcf_file = "/Data/Real_data/Slov_fM_AC_QC_filtered.vcf.gz", map_file = "/Data/Real_data/Slov_fM_AC_QC.map")

  
  #Get the haplotypes of the pedigree made by the pedigree reconstruction
  {
    Slov_mother_haplo <- Slov_allchroms_phasedwithPed[rownames(Slov_allchroms_phasedwithPed) %in% Slov_pedigree_post$mother,]
    Slov_dpc_haplo <- Slov_allchroms_phasedwithPed[rownames(Slov_allchroms_phasedwithPed) %in% Slov_pedigree_post$dpc,]
    Slov_workers_haplo <- Slov_allchroms_phasedwithPed[rownames(Slov_allchroms_phasedwithPed) %in% Slov_pedigree_post$id,]
    
    Slov_all_haplo_phasedwithPed <-  rbind(Slov_mother_haplo, Slov_dpc_haplo, Slov_workers_haplo)
    
    mother_ids <- unique(Slov_pedigree_post$mother)
    dpc_ids <- unique(Slov_pedigree_post$dpc)
    dpc_ids <- dpc_ids[dpc_ids != "0"]
    worker_ids <- unique(Slov_pedigree_post$id)
    Slov_all_ids <- c(mother_ids, dpc_ids, worker_ids)
    Slov_ind_id_1 <- paste(Slov_all_ids, "1", sep = "_")
    Slov_ind_id_2 <- paste(Slov_all_ids, "2", sep = "_")
    
    Slov_haplotypes_phasedwithPed <- get_out_haplotypes(ped_matrix = Slov_all_haplo_phasedwithPed, ind_id_1 = Slov_ind_id_1, ind_id_2 = Slov_ind_id_2)
    Slov_PhasedHaplotypes_WithPed <- apply(Slov_haplotypes_phasedwithPed, 2, convert_genotypes)
    rownames(Slov_PhasedHaplotypes_WithPed) <- rownames(Slov_haplotypes_phasedwithPed)
  }
  
  #without pedigree
  {
    Slov_mother_haplo_noped <- Slov_allchroms_phasedNoPed[rownames(Slov_allchroms_phasedNoPed) %in% Slov_pedigree_pre$mother,]
    Slov_dpc_haplo_noped <- Slov_allchroms_phasedNoPed[rownames(Slov_allchroms_phasedNoPed) %in% Slov_pedigree_pre$dpc,]
    Slov_workers_haplo_noped <- Slov_allchroms_phasedNoPed[rownames(Slov_allchroms_phasedNoPed) %in% Slov_pedigree_pre$id,]
    
    Slov_all_haplo_phasedNoPed <-  rbind(Slov_mother_haplo_noped, Slov_dpc_haplo_noped, Slov_workers_haplo_noped)
    
    mother_ids <- unique(Slov_pedigree_pre$mother)
    dpc_ids <- unique(Slov_pedigree_pre$dpc)
    dpc_ids <- dpc_ids[dpc_ids != "0"]
    worker_ids <- unique(Slov_pedigree_pre$id)
    Slov_all_ids <- c(mother_ids, dpc_ids, worker_ids)
    Slov_ind_id_1 <- paste(Slov_all_ids, "1", sep = "_")
    Slov_ind_id_2 <- paste(Slov_all_ids, "2", sep = "_")
    
    Slov_haplotypes_phasedNoPed <- get_out_haplotypes(ped_matrix = Slov_all_haplo_phasedNoPed, ind_id_1 = Slov_ind_id_1, ind_id_2 = Slov_ind_id_2)
    Slov_PhasedHaplotypes_NoPed <- apply(Slov_haplotypes_phasedNoPed, 2, convert_genotypes)
    rownames(Slov_PhasedHaplotypes_NoPed) <- rownames(Slov_haplotypes_phasedNoPed)
  }

  #Prephased - used in the next script for comparison 
  {
  Slov_mother_haplo_noped <- Slov_prephased_haplotypes[rownames(Slov_prephased_haplotypes) %in% Slov_pedigree_pre$mother,]
  Slov_dpc_haplo_noped <- Slov_prephased_haplotypes[rownames(Slov_prephased_haplotypes) %in% Slov_pedigree_pre$dpc,]
  Slov_workers_haplo_noped <- Slov_prephased_haplotypes[rownames(Slov_prephased_haplotypes) %in% Slov_pedigree_pre$id,]
  
  Slov_all__prephasedHaplotypes_NoPed <-  rbind(Slov_mother_haplo_noped, Slov_dpc_haplo_noped, Slov_workers_haplo_noped)
  
  mother_ids <- unique(Slov_pedigree_pre$mother)
  dpc_ids <- unique(Slov_pedigree_pre$dpc)
  dpc_ids <- dpc_ids[dpc_ids != "0"]
  worker_ids <- unique(Slov_pedigree_pre$id)
  Slov_all_ids <- c(mother_ids, dpc_ids, worker_ids)
  Slov_ind_id_1 <- paste(Slov_all_ids, "1", sep = "_")
  Slov_ind_id_2 <- paste(Slov_all_ids, "2", sep = "_")
  
  Slov_haplotypes_Prephased <- get_out_haplotypes(ped_matrix = Slov_all__prephasedHaplotypes_NoPed, ind_id_1 = Slov_ind_id_1, ind_id_2 = Slov_ind_id_2)
  Slov_haplotypes_Prephased <- apply(Slov_haplotypes_Prephased, 2, convert_genotypes)
  rownames(Slov_haplotypes_Prephased) <- rownames(Slov_haplotypes_Prephased)
  }
