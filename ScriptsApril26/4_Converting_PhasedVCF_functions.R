
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
