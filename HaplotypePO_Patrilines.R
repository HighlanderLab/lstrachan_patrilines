rm(list = ls())

library(Eagle)
library(tidyr)
library(readr)
library(genio)
library(AlphaSimR)
library(ggplot2)
library(dplyr)
library(VariantAnnotation)

############################################################################
########   GENERIC FUNCTIONS
############################################################################

name_genotypes <- function(map_file) {
  # Initialize an empty list to store the transformed rows
  transformed_list <- list()
  
  # Loop through each row in the dataframe
  for (i in 1:nrow(map_file)) {
    # Get the current row
    current_row <- map_file[i, ]
    
    # Create two copies of the current row with modified marker values
    row_copy_1 <- current_row
    row_copy_2 <- current_row
    row_copy_1$V2 <- paste(current_row$V2, "1", sep = "_")
    row_copy_2$V2 <- paste(current_row$V2, "2", sep = "_")
    
    # Add the modified rows to the list
    transformed_list[[length(transformed_list) + 1]] <- row_copy_1
    transformed_list[[length(transformed_list) + 1]] <- row_copy_2
  }
  
  # Convert the list back to a dataframe
  Long_map_file <- do.call(rbind, transformed_list)
  
  # Return the column V2
  return(Long_map_file$V2)
}

convert_genotypes <- function(genotypes) {
  genotypes[genotypes == '1'] <- 0
  genotypes[genotypes == '0'] <- NA
  genotypes[genotypes == '2'] <- 1
  genotypes[genotypes == 'A'] <- 0
  genotypes[genotypes == 'C'] <- 1
  
  return(as.numeric(genotypes))
}

transform_genotype <- function(genotype_matrix) {
  # Initialize an empty list to store the new rows
  new_rows <- list()
  
  # Get the row names of the genotype matrix
  individual_ids <- rownames(genotype_matrix)
  
  # Get the unique markers by removing the last two characters (_1 or _2) from the column names
  markers <- unique(sub("_.$", "", colnames(genotype_matrix)))
  
  # Loop through each row of the genotype matrix
  for (i in seq_len(nrow(genotype_matrix))) {
    # Get the current individual ID
    ind_id <- individual_ids[i]
    
    # Loop through each marker
    for (marker in markers) {
      # Create new row names for _1 and _2
      ind_id_1 <- paste(ind_id, "1", sep = "_")
      ind_id_2 <- paste(ind_id, "2", sep = "_")
      
      # Extract the values for _1 and _2 markers
      value_1 <- genotype_matrix[i, paste(marker, "1", sep = "_")]
      value_2 <- genotype_matrix[i, paste(marker, "2", sep = "_")]
      
      # Create new rows and add them to the list
      new_rows[[length(new_rows) + 1]] <- data.frame(ID = ind_id_1, Marker = marker, Value = value_1)
      new_rows[[length(new_rows) + 1]] <- data.frame(ID = ind_id_2, Marker = marker, Value = value_2)
    }
  }
  
  # Combine all the new rows into a single data frame
  result <- do.call(rbind, new_rows)
  
  # Reshape the result to have individuals as row names and markers as column names
  result_wide <- reshape(result, timevar = "Marker", idvar = "ID", direction = "wide")
  
  # Clean up the column names
  colnames(result_wide) <- sub("Value\\.", "", colnames(result_wide))
  
  # Set the row names to the IDs
  rownames(result_wide) <- result_wide$ID
  
  # Remove the ID column
  result_wide$ID <- NULL
  
  return(result_wide)
}

get_genotypes <- function(ped, id) {
  row <- ped %>% filter(IID == id)
  genotypes <- row[,7:ncol(row)]
  return(genotypes)
}

Hap1 <- function(df, map) {
  # Get the column names that end in _2
  columns_ending_in_1 <- grep("_1$", colnames(df), value = TRUE)
  
  # Subset the dataframe with the selected columns
  new_df <- df[, columns_ending_in_1, drop = FALSE]
  
  colnames(new_df) <- map[,2]
  return(new_df)
}

Hap2 <- function(df, map) {
  # Get the column names that end in _2
  columns_ending_in_2 <- grep("_2$", colnames(df), value = TRUE)
  
  # Subset the dataframe with the selected columns
  new_df <- df[, columns_ending_in_2, drop = FALSE]
  
  colnames(new_df) <- map[,2]
  return(new_df)
}

calcGeno <- function(Hap1, Hap2, map){
  Hap1 <- t(Hap1)
  Hap2 <- t(Hap2)
  Geno <- matrix(data = Hap1 + Hap2, nrow = 1)
  colnames(Geno) <- map[,2]
  
  return(Geno)
}

assign_parent_haplo <- function(df, maternal = NULL, paternal = NULL, offspring_id = NULL) {
  
  if(!is.null(maternal)){
    row_copy <- paste(offspring_id, "maternal", sep = "_")
  }
  if(!is.null(paternal)){
    row_copy <- paste(offspring_id, "paternal", sep = "_")
  }
  
  rownames(df) <- row_copy
  # Return the column V2
  return(df)
}

get_phased_haplotypes <- function(Geno_Error = NULL){
  # Set the working directory and read the pedigree file
  setwd("~/Desktop/Slovenia data/Attempt2/Nested/General Data")
  pedigree <- read_csv("worker_pedigree.csv")
  
  if (Geno_Error == "GE"){
    # Set the working directory and read the phased data
    #Geno Error 
    setwd("~/Desktop/Slovenia data/Attempt2/Nested/SNP4/Geno Err/phased")
    Nested_ped_phased <- read.table("phased_SNP4_ped_withGenoError_ACformat_QC.ped", header = FALSE)
    map <- read.table("phased_SNP4_ped_withGenoError_ACformat_QC.map", header = FALSE)
    
  } else if (Geno_Error == "nGE") {
    #No GenoError 
    setwd("~/Desktop/Slovenia data/Attempt2/Nested/SNP4/No_GE/phased")
    Nested_ped_phased <- read.table("SNP4Nested_nGE_phasingTest.ped", header = FALSE)
    map <- read.table("SNP4Nested_nGE_phasingTest.map", header = FALSE)
  } else{
    stop("No Geno_Error provided, must be GE or nGE")
  }
  
  
  # Extract genotype columns and convert them to 0/1 format
  genotypes_columns_phased <- Nested_ped_phased[, 7:ncol(Nested_ped_phased)]
  genotypes_01_phased <- apply(genotypes_columns_phased, 2, convert_genotypes)
  
  # Generate IDs for the pedigree
  IDs_forPed <- name_genotypes(map_file = map)
  rownames(genotypes_01_phased) <- Nested_ped_phased$V2
  colnames(genotypes_01_phased) <- IDs_forPed
  
  # Transform genotypes
  phased_haplotypes <- transform_genotype(genotypes_01_phased)
  
  return(phased_haplotypes)
}

order_by_prefix <- function(df) {
  # Extract the part before the underscore
  prefix <- as.numeric(sapply(rownames(df), function(x) strsplit(x, "_")[[1]][1]))
  
  # Order the data frame based on the extracted prefix
  df_ordered <- df[order(prefix, rownames(df)), , drop = FALSE]
  
  return(df_ordered)
}


####################################################################################
###### Manually convert vcf to ped file - not using PLINK
#####################################################################################

setwd("~/Desktop/Slovenia data/Test_simulated")
pedigree <- read.csv("worker_pedigree.csv")
#ped order is always mum/dpc/workers
ids_all <- c(unique(pedigree$mother), unique(pedigree$dpc), unique(pedigree$id))
ind_id_1 <- paste(ids_all, "1", sep = "_")
ind_id_2 <- paste(ids_all, "2", sep = "_")

setwd("~/Desktop/Slovenia data/PLINK/Phasing/QC/")

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

nGE_chr1 <- convert_VCF(vcf_file = "nGE_QC_PHASED_chr1.vcf", map_file = "nGE_QC_PHASED_chr1.map")
nGE_chr2 <- convert_VCF(vcf_file = "nGE_QC_PHASED_chr2.vcf", map_file = "nGE_QC_PHASED_chr2.map")
nGE_chr3 <- convert_VCF(vcf_file = "nGE_QC_PHASED_chr3.vcf", map_file = "nGE_QC_PHASED_chr3.map")
nGE_chr4 <- convert_VCF(vcf_file = "nGE_QC_PHASED_chr4.vcf", map_file = "nGE_QC_PHASED_chr4.map")
nGE_chr5 <- convert_VCF(vcf_file = "nGE_QC_PHASED_chr5.vcf", map_file = "nGE_QC_PHASED_chr5.map")
nGE_chr6 <- convert_VCF(vcf_file = "nGE_QC_PHASED_chr6.vcf", map_file = "nGE_QC_PHASED_chr6.map")
nGE_chr7 <- convert_VCF(vcf_file = "nGE_QC_PHASED_chr7.vcf", map_file = "nGE_QC_PHASED_chr7.map")
nGE_chr8 <- convert_VCF(vcf_file = "nGE_QC_PHASED_chr8.vcf", map_file = "nGE_QC_PHASED_chr8.map")
nGE_chr9 <- convert_VCF(vcf_file = "nGE_QC_PHASED_chr9.vcf", map_file = "nGE_QC_PHASED_chr9.map")
nGE_chr10 <- convert_VCF(vcf_file = "nGE_QC_PHASED_chr10.vcf", map_file = "nGE_QC_PHASED_chr10.map")
nGE_chr11 <- convert_VCF(vcf_file = "nGE_QC_PHASED_chr11.vcf", map_file = "nGE_QC_PHASED_chr11.map")
nGE_chr12 <- convert_VCF(vcf_file = "nGE_QC_PHASED_chr12.vcf", map_file = "nGE_QC_PHASED_chr12.map")
nGE_chr13 <- convert_VCF(vcf_file = "nGE_QC_PHASED_chr13.vcf", map_file = "nGE_QC_PHASED_chr13.map")
nGE_chr14 <- convert_VCF(vcf_file = "nGE_QC_PHASED_chr14.vcf", map_file = "nGE_QC_PHASED_chr14.map")
nGE_chr15 <- convert_VCF(vcf_file = "nGE_QC_PHASED_chr15.vcf", map_file = "nGE_QC_PHASED_chr15.map")
nGE_chr16 <- convert_VCF(vcf_file = "nGE_QC_PHASED_chr16.vcf", map_file = "nGE_QC_PHASED_chr16.map")

nGE_all <- cbind(nGE_chr1,nGE_chr2,nGE_chr3,nGE_chr4,nGE_chr5,nGE_chr6,nGE_chr7,nGE_chr8,
                 nGE_chr9,nGE_chr10,nGE_chr11,nGE_chr12,nGE_chr13,nGE_chr14,nGE_chr15,nGE_chr16)

rm(nGE_chr1,nGE_chr2,nGE_chr3,nGE_chr4,nGE_chr5,nGE_chr6,nGE_chr7,nGE_chr8,
      nGE_chr9,nGE_chr10,nGE_chr11,nGE_chr12,nGE_chr13,nGE_chr14,nGE_chr15,nGE_chr16)


GE_chr1 <- convert_VCF(vcf_file = "GE_QC_PHASED_chr1.vcf", map_file = "GE_QC_PHASED_chr1.map")
GE_chr2 <- convert_VCF(vcf_file = "GE_QC_PHASED_chr2.vcf", map_file = "GE_QC_PHASED_chr2.map")
GE_chr3 <- convert_VCF(vcf_file = "GE_QC_PHASED_chr3.vcf", map_file = "GE_QC_PHASED_chr3.map")
GE_chr4 <- convert_VCF(vcf_file = "GE_QC_PHASED_chr4.vcf", map_file = "GE_QC_PHASED_chr4.map")
GE_chr5 <- convert_VCF(vcf_file = "GE_QC_PHASED_chr5.vcf", map_file = "GE_QC_PHASED_chr5.map")
GE_chr6 <- convert_VCF(vcf_file = "GE_QC_PHASED_chr6.vcf", map_file = "GE_QC_PHASED_chr6.map")
GE_chr7 <- convert_VCF(vcf_file = "GE_QC_PHASED_chr7.vcf", map_file = "GE_QC_PHASED_chr7.map")
GE_chr8 <- convert_VCF(vcf_file = "GE_QC_PHASED_chr8.vcf", map_file = "GE_QC_PHASED_chr8.map")
GE_chr9 <- convert_VCF(vcf_file = "GE_QC_PHASED_chr9.vcf", map_file = "GE_QC_PHASED_chr9.map")
GE_chr10 <- convert_VCF(vcf_file = "GE_QC_PHASED_chr10.vcf", map_file = "GE_QC_PHASED_chr10.map")
GE_chr11 <- convert_VCF(vcf_file = "GE_QC_PHASED_chr11.vcf", map_file = "GE_QC_PHASED_chr11.map")
GE_chr12 <- convert_VCF(vcf_file = "GE_QC_PHASED_chr12.vcf", map_file = "GE_QC_PHASED_chr12.map")
GE_chr13 <- convert_VCF(vcf_file = "GE_QC_PHASED_chr13.vcf", map_file = "GE_QC_PHASED_chr13.map")
GE_chr14 <- convert_VCF(vcf_file = "GE_QC_PHASED_chr14.vcf", map_file = "GE_QC_PHASED_chr14.map")
GE_chr15 <- convert_VCF(vcf_file = "GE_QC_PHASED_chr15.vcf", map_file = "GE_QC_PHASED_chr15.map")
GE_chr16 <- convert_VCF(vcf_file = "GE_QC_PHASED_chr16.vcf", map_file = "GE_QC_PHASED_chr16.map")

GE_all <- cbind(GE_chr1,GE_chr2,GE_chr3,GE_chr4,GE_chr5,GE_chr6,GE_chr7,GE_chr8,
                 GE_chr9,GE_chr10,GE_chr11,GE_chr12,GE_chr13,GE_chr14,GE_chr15,GE_chr16)

rm(GE_chr1,GE_chr2,GE_chr3,GE_chr4,GE_chr5,GE_chr6,GE_chr7,GE_chr8,
      GE_chr9,GE_chr10,GE_chr11,GE_chr12,GE_chr13,GE_chr14,GE_chr15,GE_chr16)



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
nGE_haplotypes <- get_out_haplotypes(ped_matrix = nGE_all, ind_id_1 = ind_id_1, ind_id_2 = ind_id_2)
GE_haplotypes <- get_out_haplotypes(ped_matrix = GE_all, ind_id_1 = ind_id_1, ind_id_2 = ind_id_2)
#convert to 0/1 format 
nGE4.0_phasedhaplotypes <- apply(nGE_haplotypes, 2, convert_genotypes)
rownames(nGE4.0_phasedhaplotypes) <- rownames(nGE_haplotypes)

GE4.0_phasedhaplotypes <- apply(GE_haplotypes, 2, convert_genotypes)
rownames(GE4.0_phasedhaplotypes) <- rownames(GE_haplotypes)

####################################################################################
# Get the true haplotypes
#####################################################################################
setwd("~/Desktop/Slovenia data/Attempt2/Nested/General Data")
pedigree <- read.csv("worker_pedigree.csv")
load("~/Desktop/Slovenia data/Attempt2/Nested/General Data/SP_object.Rdata")
load("~/Desktop/Slovenia data/Attempt2/Nested/General Data/Pop_withNoFathers.Rdata")
##### ALPHASIMR FORMAT Real haplotypes #####
real_haplotypes <- pullSnpHaplo(PopMerged_noFathers)
real_map <- getGenMap(SP)


#####################################################################################
#### Check the different phasing versions compared to the real and prephased data just to check similarity 
#####################################################################################

#function for checking  
check_haplotype <- function(real_haplotypes = NULL, results = NULL, pedigree = NULL) {
  comparison_results <- data.frame()
  
  cols <- colnames(results)
  real_haplotypes_cols <- real_haplotypes[, cols, drop = FALSE]
  
  offspring_ids <- pedigree$id
  
  for (i in seq_along(offspring_ids)) {
    offspring_id <- offspring_ids[i]
    
    # Pull rows starting with offspring_id from the real haplotypes and results
    real_haplo_maternal <- rownames(real_haplotypes_cols)[grep(paste0("^", offspring_id, "_1"), rownames(real_haplotypes_cols))]
    real_haplo_maternal <- real_haplotypes_cols[real_haplo_maternal, , drop = FALSE]
    
    real_haplo_paternal <- rownames(real_haplotypes_cols)[grep(paste0("^", offspring_id, "_2"), rownames(real_haplotypes_cols))]
    real_haplo_paternal <- real_haplotypes_cols[real_haplo_paternal, , drop = FALSE]
    
    results_maternal <- rownames(results)[grep(paste0("^", offspring_id, "_1"), rownames(results))]
    results_maternal <- results[results_maternal, , drop = FALSE]
    
    results_paternal <- rownames(results)[grep(paste0("^", offspring_id, "_2"), rownames(results))]
    results_paternal <- results[results_paternal, , drop = FALSE]
    
    # Compare haplotypes
    maternal_comparison <- as.data.frame(matrix(NA, nrow = 1, ncol = 2))
    rownames(maternal_comparison) <- rownames(results_maternal)
    colnames(maternal_comparison) <- c("_1", "_2")
    
    total_elements <- length(results_maternal)
    
    maternal_comparison$`_1` <- (sum(results_maternal == real_haplo_maternal, na.rm = TRUE)/total_elements) * 100
    maternal_comparison$`_2` <- (sum(results_maternal == real_haplo_paternal, na.rm = TRUE)/total_elements) * 100
    
    paternal_comparison <- as.data.frame(matrix(NA, nrow = 1, ncol = 2))
    rownames(paternal_comparison) <- rownames(results_paternal)
    colnames(paternal_comparison) <- c("_1", "_2")
    
    paternal_comparison$`_1` <- (sum(results_paternal == real_haplo_maternal)/total_elements) * 100
    paternal_comparison$`_2` <- (sum(results_paternal == real_haplo_paternal)/total_elements) * 100
    
    # Combine results for current offspring_id
    current_comparison <- rbind(maternal_comparison, paternal_comparison)
    
    # Append current_comparison to comparison_results
    comparison_results <- rbind(comparison_results, current_comparison)
  }
  # Check conditions and print messages
  comparison_results_maternal <- rownames(comparison_results)[grep(paste0("_1"), rownames(comparison_results))]
  comparison_results_maternal <- comparison_results[comparison_results_maternal, , drop = FALSE]
  
  comparison_results_paternal <- rownames(comparison_results)[grep(paste0("_2"), rownames(comparison_results))]
  comparison_results_paternal <- comparison_results[comparison_results_paternal, , drop = FALSE]
  
  maternal_correct <- all(comparison_results_maternal[, "_1"] > comparison_results_maternal[, "_2"])
  paternal_correct <- all(comparison_results_paternal[, "_1"] < comparison_results_paternal[, "_2"])
  
  print(mean(comparison_results_maternal$`_1`))
  print(mean(comparison_results_maternal$`_2`))
  print(mean(comparison_results_paternal$`_1`))
  print(mean(comparison_results_paternal$`_2`))
  
  return(comparison_results)
}

real_vc_nGE_phased <- check_haplotype(real_haplotypes = real_haplotypes, results = nGE4.0_phasedhaplotypes, pedigree = pedigree)
real_vc_GE_phased <- check_haplotype(real_haplotypes = real_haplotypes, results = GE4.0_phasedhaplotypes, pedigree = pedigree)
nGE_phased_vc_GE_phased <- check_haplotype(real_haplotypes = nGE4.0_phasedhaplotypes, results = GE4.0_phasedhaplotypes, pedigree = pedigree)

real_vs_real <- check_haplotype(real_haplotypes = real_haplotypes, results = real_haplotypes, pedigree = pedigree)


#####################################################################################
####                  Haplotype parent-of-origin assignments 
#####################################################################################

#Need to make a flipping function first - haplotype 1 and 2 are not assigned parents 

####  using TRUE/FALSE to determine parent ##############################

#functions for the flipping functions
weighing_haplotypes <- function(off_hap, par_geno, method = NULL, result = NULL, j){
  
  # Description of the parameters:
  # off_hap: Vector of offspring haplotypes.
  # par_geno: Vector of parent genotypes.
  # method: The method used for scoring; options include mean, weighted_mean, etc.
  # j: Chromosome identifier.
  
  # Calculate the logical vector
  logical_vec <- (off_hap - par_geno) %in% c(-2, 1)
  
  # Run length encoding
  rle_result <- rle(logical_vec)
  
  # Separate the runs of FALSE and TRUE
  false_runs <- rle_result$lengths[rle_result$values == FALSE]
  true_runs <- rle_result$lengths[rle_result$values == TRUE]
  
  # Define weights as the lengths of the runs
  false_weights <- false_runs
  true_weights <- true_runs
  
  # Calculate power parameter p as the mean of run lengths if not empty
  p <- if(length(false_runs) > 0) {
    mean(false_runs)  # just make it the mean for now
  } else {
    1  # default for p if there are no false_runs
  }
  
  
  # Handle scoring methods
  if (method == "mean") {  
    false_score <- mean(false_runs)  
    true_score <- mean(true_runs)  
    overall_score <- false_score - true_score  
  } else if (method == "sum"){  # Total Sum
    false_score <- sum(false_runs) 
    true_score <- sum(true_runs)  
    overall_score <- false_score - true_score 
    
  } else if (method == "quadratic") {  # Quadratic Mean
    false_score <- sqrt(sum(false_runs^2) / length(false_runs))  # High =  more mismatches & less likely to be parent
    true_score <- sqrt(sum(true_runs^2) / length(false_runs))  # High = longer runs of matches & more likely to be parent
    overall_score <- false_score - true_score  # High =  more mismatches & less likely to be parent
    
  } else if (method == "logarithmic") {  # Logarithmic Mean
    false_score <- sum(log(false_runs + 1)) / length(false_runs)  # High =  more mismatches & less likely to be parent
    true_score <- sum(log(true_runs + 1)) / length(false_runs)  # High = longer runs of matches & more likely to be parent
    overall_score <- false_score - true_score  # High =  more mismatches & less likely to be parent
    
  } else if (method == "geometric") {  # Geometric Mean
    false_score <- prod(false_runs)^(1/length(false_runs))  # High =  more mismatches & less likely to be parent
    true_score <- prod(true_runs)^(1/length(true_runs))  # High = longer runs of matches & more likely to be parent
    overall_score <- false_score - true_score  # High =  more mismatches & less likely to be parent
    
  } else if (method == "composite"){  # Composite Score
    false_score <- (sum(false_runs)^2) / length(false_runs)  # High =  more mismatches & less likely to be parent
    true_score <- (sum(true_runs)^2) / length(true_runs)  # High = longer runs of matches & more likely to be parent
    overall_score <- false_score - true_score  # High =  more mismatches & less likely to be parent
    
  } else if (method == "exponential"){  # Exponential Mean
    false_score <- sum(2^false_runs) / length(false_runs)  # High =  more mismatches & less likely to be parent
    true_score <- sum(2^true_runs) / length(false_runs)  # High = longer runs of matches & more likely to be parent
    overall_score <- false_score - true_score  # High =  more mismatches & less likely to be parent
    
  } else if (method == "harmonic"){  # Harmonic Mean
    false_score <- length(false_runs) / sum(1 / false_runs)  # High =  more mismatches & less likely to be parent
    true_score <- length(true_runs) / sum(1 / true_runs)  # High = longer runs of matches & more likely to be parent
    overall_score <- false_score - true_score  # High =  more mismatches & less likely to be parent
    
  } else if (method == "power_mean" && !is.null(p)){  # Power Mean
    false_score <- (sum(false_runs^p) / length(false_runs))^(1/p)  # High =  more mismatches & less likely to be parent
    true_score <- (sum(true_runs^p) / length(true_runs))^(1/p)  # High = longer runs of matches & more likely to be parent
    overall_score <- false_score - true_score  # High =  more mismatches & less likely to be parent
    
  } else {
    stop("Invalid method specified")
  }
  
  weighed_haplotypes <- data.frame(
    false_score = false_score,
    true_score = true_score,
    overall_score = overall_score,
    method = method,
    chr = j
  )
  
  return(weighed_haplotypes)
}
#method options :  c("mean","sum","quadratic","logarithmic" ,"geometric","composite", "exponential", "harmonic","power_mean") 

check_dist <- function(off_hap, par_geno, j){
  
  logical_vec <- (off_hap - par_geno) %in% c(-2, 1)
  
  # Run length encoding
  rle_result <- rle(logical_vec)
  
  # Separate the runs of FALSE and TRUE
  false_runs <- rle_result$lengths[rle_result$values == FALSE]
  true_runs <- rle_result$lengths[rle_result$values == TRUE]
  
  if(length(false_runs) != 0){
    false_runsdf <- data.frame(length = false_runs, type = "FALSE")
  } else {
    false_runsdf <- data.frame(length = NA, type = "FALSE")
  }
  if (length(true_runs) != 0) {
    true_runsdf <- data.frame(length = true_runs, type = "TRUE")
  } else {
    true_runsdf <- data.frame(length = NA, type = "TRUE")
  }
  
  runs_df <- rbind(false_runsdf, true_runsdf)
  runs_df$chr <- j
  
  return(runs_df)
}

TRUE_FALSE_flipping <- function(Beagle_version = NULL, Data_type = NULL, pedigree = NULL, Real_haplotypes = NULL, method = NULL){
  
  
  if (Real_haplotypes == TRUE & Data_type == "Real" & Beagle_version == "Non"){
     real_haplotypes <- real_haplotypes
     map <- real_map
    
    real_results_methods <- list()
    real_results_lengths <- list()
    flipped_haplotypes <- list()
    
    
    for (i in 1:nrow(pedigree)) {
      offspring_id <- pedigree$id[i]
      mother_id <- pedigree$mother[i]
      dpc_id <- pedigree$dpc[i]
      
      offspring_row <- rownames(real_haplotypes)[sapply(rownames(real_haplotypes), function(x) strsplit(x,'_')[[1]][1]) %in% offspring_id]
      offspring_haplotypes_i <- real_haplotypes[offspring_row,]
      Offspring_Hap1 <- t(as.data.frame(offspring_haplotypes_i[1,]))
      Offspring_Hap2 <- t(as.data.frame(offspring_haplotypes_i[2,]))
      
      mother_row <- rownames(real_haplotypes)[sapply(rownames(real_haplotypes), function(x) strsplit(x,'_')[[1]][1]) %in% mother_id]
      mother_haplotypes_i <- real_haplotypes[mother_row,]
      Maternal_Hap1 <- t(as.data.frame(mother_haplotypes_i[1,]))
      rownames(Maternal_Hap1) <- mother_id
      Maternal_Hap2 <- t(as.data.frame(mother_haplotypes_i[2,]))
      rownames(Maternal_Hap2) <- mother_id
      Maternal_Geno <- matrix(data = Maternal_Hap1 + Maternal_Hap2, nrow = 1)
      colnames(Maternal_Geno) <- colnames(Maternal_Hap1)
      rownames(Maternal_Geno) <- mother_id
      
      dpc_row <- rownames(real_haplotypes)[sapply(rownames(real_haplotypes), function(x) strsplit(x,'_')[[1]][1]) %in% dpc_id]
      dpc_haplotypes_i <- real_haplotypes[dpc_row,]
      Dpc_Hap1 <- t(as.data.frame(dpc_haplotypes_i[1,]))
      rownames(Dpc_Hap1) <- dpc_id
      Dpc_Hap2 <- t(as.data.frame(dpc_haplotypes_i[2,]))
      rownames(Dpc_Hap2) <- dpc_id
      Dpc_Geno <- matrix(data = Dpc_Hap1 + Dpc_Hap2, nrow = 1)
      colnames(Dpc_Geno) <- colnames(Dpc_Hap1)
      rownames(Dpc_Geno) <- dpc_id
      
      method_results <- list()
      length_results <- list()
      pat_results <- list()
      mat_results <- list() 
      
      
      for (j in unique(map[, 2])) {
        print(x = c(i,j))
        Chr_map <- map[map$chr == j, ]
        Chr_markerNames <- Chr_map[, 1]
        
        
        Offspring_Hap1_chr <- Offspring_Hap1[, colnames(Offspring_Hap1) %in% Chr_markerNames, drop = FALSE]
        Offspring_Hap2_chr <- Offspring_Hap2[, colnames(Offspring_Hap2) %in% Chr_markerNames, drop = FALSE]
        Maternal_Geno_chr <- Maternal_Geno[, colnames(Maternal_Geno) %in% Chr_markerNames, drop = FALSE]
        Dpc_Geno_chr <- Dpc_Geno[, colnames(Dpc_Geno) %in% Chr_markerNames, drop = FALSE]
        
        scores_h1_p <- weighing_haplotypes(Offspring_Hap1_chr, Dpc_Geno_chr, method = method, j =j)
        scores_h1_m <- weighing_haplotypes(Offspring_Hap1_chr, Maternal_Geno_chr, method = method, j =j)
        scores_h2_p <- weighing_haplotypes(Offspring_Hap2_chr, Dpc_Geno_chr, method = method, j =j)
        scores_h2_m <- weighing_haplotypes(Offspring_Hap2_chr, Maternal_Geno_chr, method = method, j =j)
        
        # Combine scores for this method
        scores <- rbind(scores_h1_p, scores_h1_m, scores_h2_p, scores_h2_m)
        scores$Parent <- c("h1_p", "h1_m", "h2_p", "h2_m")
        
        #do check
        score_1 <- scores_h1_p$false_score
        score_2 <- scores_h1_m$false_score
        score_3 <- scores_h2_p$false_score
        score_4 <- scores_h2_m$false_score
        
        abs_diff <- abs(c(score_1, score_2, score_3, score_4))
        max_diff <- max(abs_diff)
        
        
        if (max_diff == abs(score_1) & max_diff != abs(score_2) & max_diff != abs(score_3) & max_diff != abs(score_4)) {
          message("Hap1 is paternal - assume Hap2 is maternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, paternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, maternal = TRUE, offspring_id = offspring_id)
        } 
        if (max_diff == abs(score_2) & max_diff != abs(score_1) & max_diff != abs(score_3) & max_diff != abs(score_4)) {
          message("Hap1 is maternal - assume Hap2 is paternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, maternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, paternal = TRUE, offspring_id = offspring_id)
        } 
        if (max_diff == abs(score_3) & max_diff != abs(score_2) & max_diff != abs(score_1) & max_diff != abs(score_4)) {
          message("Hap2 is paternal - assume Hap1 is maternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, maternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, paternal = TRUE, offspring_id = offspring_id)
        } 
        if (max_diff == abs(score_4) & max_diff != abs(score_2) & max_diff != abs(score_1) & max_diff != abs(score_3))  {
          message("Hap2 is maternal - assume Hap1 is paternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, paternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, maternal = TRUE, offspring_id = offspring_id)
        }
        if (max_diff == abs(score_4) & max_diff == abs(score_3) &  abs(score_1) > abs(score_2))  {
          message("Hap2 is maternal - assume Hap1 is paternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, paternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, maternal = TRUE, offspring_id = offspring_id)
        }
        if (max_diff == abs(score_4) & max_diff == abs(score_3) &  abs(score_1) < abs(score_2))  {
          message("Hap2 is maternal - assume Hap1 is paternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, maternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, paternal = TRUE, offspring_id = offspring_id)
        }
        if (max_diff == abs(score_1) & max_diff == abs(score_2) &  abs(score_3) < abs(score_4))  {
          message("Hap2 is maternal - assume Hap1 is paternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, paternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, maternal = TRUE, offspring_id = offspring_id)
        }
        if (max_diff == abs(score_1) & max_diff == abs(score_2) &  abs(score_3) > abs(score_4))  {
          message("Hap2 is maternal - assume Hap1 is paternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, maternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, paternal = TRUE, offspring_id = offspring_id)
        }
        if (max_diff == abs(score_1) & max_diff == abs(score_2) &  max_diff == abs(score_3) & max_diff == abs(score_4))  {
          stop("All values are equal")
          
        }
        if (max_diff == abs(score_1) & max_diff == abs(score_4) &  max_diff != abs(score_2) & max_diff != abs(score_3))  {
          message("Hap1 is paternal and Hap2 is maternal - both are max")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, paternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, maternal = TRUE, offspring_id = offspring_id)
        }
        if (max_diff == abs(score_2) & max_diff == abs(score_3) &  max_diff != abs(score_1) & max_diff != abs(score_4))  {
          message("Hap1 is maternal and Hap2 is paternal - both are max")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, maternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, paternal = TRUE, offspring_id = offspring_id)
          
        }
        
        merged_haps <- rbind(Offspring_Hap1_chr, Offspring_Hap2_chr)
        identified_dpc_haplo <- grep("_paternal$", rownames(merged_haps), value = TRUE)
        identified_maternal_haplo <- grep("_maternal$", rownames(merged_haps), value = TRUE)
        
        dpc_assigned <- merged_haps[identified_dpc_haplo, , drop = FALSE]
        maternal_assigned <- merged_haps[identified_maternal_haplo, , drop = FALSE]
        
        pat_results[[j]] <- as.data.frame(dpc_assigned)
        mat_results[[j]] <- as.data.frame(maternal_assigned)
        
        # Combine all scores into one data frame
        testing_methods <- do.call(rbind, scores)
        
        
        lengths_h1_p <- check_dist(Offspring_Hap1_chr, Dpc_Geno_chr, j =j)
        lengths_h1_p$Parent <- rep("h1_p")
        lengths_h1_m <- check_dist(Offspring_Hap1_chr, Maternal_Geno_chr, j =j)
        lengths_h1_m$Parent <- rep("h1_m")
        lengths_h2_p <- check_dist(Offspring_Hap2_chr, Dpc_Geno_chr, j =j)
        lengths_h2_p$Parent <- rep("h2_p")
        lengths_h2_m <- check_dist(Offspring_Hap2_chr, Maternal_Geno_chr, j =j)
        lengths_h2_m$Parent <- rep("h2_m")
        
        lengths_all <- rbind(lengths_h1_p, lengths_h1_m, lengths_h2_p, lengths_h2_m)
        
        
        
        
        method_results[[j]] <- as.data.frame(testing_methods)
        length_results[[j]] <- as.data.frame(lengths_all)
      }
      
      real_results_methods[[i]] <- do.call(rbind,method_results)
      real_results_lengths[[i]] <- do.call(rbind, length_results)
      
      pat_results <- do.call(cbind, pat_results)
      mat_results <- do.call(cbind, mat_results)
      results <- rbind(mat_results, pat_results)
      flipped_haplotypes[[i]] <- as.data.frame(results)
      
    }  
    real_results_methods_all <- do.call(rbind, real_results_methods)
    real_results_lengths_all <- do.call(rbind, real_results_lengths)
    real_results_flipped <- do.call(rbind, flipped_haplotypes)
    
    results <- list(real_results_methods_all = real_results_methods_all, real_results_lengths_all = real_results_lengths_all,
                    real_results_flipped = real_results_flipped)
    
    
  }
  
  if (Real_haplotypes == FALSE){
    
    phased_results_methods <- list()
    phased_results_lengths <- list()
    flipped_haplotypes <- list()
    
    if(Data_type == "Pre_nGE_nQC" & Beagle_version == "Non"){
      haplotypes <- nGE_nQC_nPhased_haplotypes
      map <- nGE_nQC_nPhased_map
    }else if (Data_type == "Pre_GE_nQC" & Beagle_version == "Non"){
      haplotypes <- GE4.0
      map <- GE_nQC_nPhased_map
    } else if (Data_type == "Phased_nGE" & Beagle_version == "4.0"){
      haplotypes <- nGE4.0_phasedhaplotypes
      map <- nGE_phased_map
    } else if (Data_type == "Phased_GE" & Beagle_version == "4.0"){
      haplotypes <- GE4.0_phasedhaplotypes
      map <- GE_phased_map
    } else if (Data_type == "Phased_nGE" & Beagle_version == "4.1"){
      haplotypes <- nGE4.1_phasedhaplotypes
      map <- nGE_4.1_map
    } else if (Data_type == "Phased_GE" & Beagle_version == "4.1"){
      haplotypes <- GE4.1_phasedhaplotypes
      map <- GE_4.1_map
     } else if (Data_type == "Phased_GE" & Beagle_version == "Slov4.0"){
        haplotypes <- Slov4.0_phasedhaplotype
        map <- Slov_map
    } else (stop("Arguments not met"))
    
    for (i in 1:nrow(pedigree)) {
      offspring_id <- pedigree$id[i]
      mother_id <- pedigree$mother[i]
      dpc_id <- pedigree$dpc[i]
      
      offspring_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x, '_')[[1]][1]) %in% offspring_id]
      offspring_haplotypes_i <- haplotypes[offspring_row, ]
      Offspring_Hap1 <- as.data.frame(offspring_haplotypes_i[1, , drop = FALSE])
      Offspring_Hap2 <- as.data.frame(offspring_haplotypes_i[2, , drop = FALSE])
      
      mother_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x, '_')[[1]][1]) %in% mother_id]
      mother_haplotypes_i <- haplotypes[mother_row, ]
      Maternal_Hap1 <- as.data.frame(mother_haplotypes_i[1, , drop = FALSE])
      rownames(Maternal_Hap1) <- mother_id
      Maternal_Hap2 <- as.data.frame(mother_haplotypes_i[2, , drop = FALSE])
      rownames(Maternal_Hap2) <- mother_id
      Maternal_Geno <- matrix(data = Maternal_Hap1 + Maternal_Hap2, nrow = 1)
      colnames(Maternal_Geno) <- colnames(Maternal_Hap1)
      rownames(Maternal_Geno) <- mother_id
      
      dpc_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x, '_')[[1]][1]) %in% dpc_id]
      dpc_haplotypes_i <- haplotypes[dpc_row, ]
      Dpc_Hap1 <- as.data.frame(dpc_haplotypes_i[1, , drop = FALSE])
      rownames(Dpc_Hap1) <- dpc_id
      Dpc_Hap2 <- as.data.frame(dpc_haplotypes_i[2, , drop = FALSE])
      rownames(Dpc_Hap2) <- dpc_id
      Dpc_Geno <- matrix(data = Dpc_Hap1 + Dpc_Hap2, nrow = 1)
      colnames(Dpc_Geno) <- colnames(Dpc_Hap1)
      rownames(Dpc_Geno) <- dpc_id
      
      method_results <- list()
      length_results <- list()
      pat_results <- list()
      mat_results <- list() 
      
      for (j in unique(map[, 1])) {
        print(x = c(i,j))
        Chr_map <- map[map$V1 == j, ]
        Chr_markerNames <- Chr_map[, 2]
        
        Offspring_Hap1_chr <- Offspring_Hap1[, colnames(Offspring_Hap1) %in% Chr_markerNames, drop = FALSE]
        Offspring_Hap2_chr <- Offspring_Hap2[, colnames(Offspring_Hap2) %in% Chr_markerNames, drop = FALSE]
        Maternal_Geno_chr <- Maternal_Geno[, colnames(Maternal_Geno) %in% Chr_markerNames, drop = FALSE]
        Dpc_Geno_chr <- Dpc_Geno[, colnames(Dpc_Geno) %in% Chr_markerNames, drop = FALSE]
        
        scores_h1_p <- weighing_haplotypes(Offspring_Hap1_chr, Dpc_Geno_chr, method = method, j =j)
        scores_h1_m <- weighing_haplotypes(Offspring_Hap1_chr, Maternal_Geno_chr, method = method, j =j)
        scores_h2_p <- weighing_haplotypes(Offspring_Hap2_chr, Dpc_Geno_chr, method = method, j =j)
        scores_h2_m <- weighing_haplotypes(Offspring_Hap2_chr, Maternal_Geno_chr, method = method, j =j)
        
        # Combine scores for this method
        scores <- rbind(scores_h1_p, scores_h1_m, scores_h2_p, scores_h2_m)
        scores$Parent <- c("h1_p", "h1_m", "h2_p", "h2_m")
        
        #do check
        score_1 <- scores_h1_p$false_score
        score_2 <- scores_h1_m$false_score
        score_3 <- scores_h2_p$false_score
        score_4 <- scores_h2_m$false_score
        
        if (is.nan(score_1)) score_1 <- 0
        if (is.nan(score_2)) score_2 <- 0
        if (is.nan(score_3)) score_3 <- 0
        if (is.nan(score_4)) score_4 <- 0
        
        
        abs_diff <- abs(c(score_1, score_2, score_3, score_4))
        
        max_diff <- max(abs_diff, na.rm = TRUE)
        
        # if  (max_diff == abs(score_4) | max_diff == abs(score_1)) {
        #   print(c(i,j))
        #   print(c(score_1, score_2, score_3, score_4))
        #   stop("Something is funky here")
        # 
        # } 
        
        #Only 1 is the absolute clear parent 
        if (max_diff == abs(score_1) & max_diff != abs(score_2) & max_diff != abs(score_3) & max_diff != abs(score_4)) {
          message("Hap1 is paternal - assume Hap2 is maternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, paternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, maternal = TRUE, offspring_id = offspring_id)
        } 
        if (max_diff == abs(score_2) & max_diff != abs(score_1) & max_diff != abs(score_3) & max_diff != abs(score_4)) {
          message("Hap1 is maternal - assume Hap2 is paternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, maternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, paternal = TRUE, offspring_id = offspring_id)
        } 
        if (max_diff == abs(score_3) & max_diff != abs(score_2) & max_diff != abs(score_1) & max_diff != abs(score_4)) {
          message("Hap2 is paternal - assume Hap1 is maternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, maternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, paternal = TRUE, offspring_id = offspring_id)
        } 
        if (max_diff == abs(score_4) & max_diff != abs(score_2) & max_diff != abs(score_1) & max_diff != abs(score_3))  {
          message("Hap2 is maternal - assume Hap1 is paternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, paternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, maternal = TRUE, offspring_id = offspring_id)
        }
        
        #Two on the same haplotype have equally high scores, use the other haplotype to work it out 
        if (max_diff == abs(score_4) & max_diff == abs(score_3) &  abs(score_1) > abs(score_2))  {
          message("Hap2 is maternal - assume Hap1 is paternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, paternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, maternal = TRUE, offspring_id = offspring_id)
        }
        if (max_diff == abs(score_4) & max_diff == abs(score_3) &  abs(score_1) < abs(score_2))  {
          message("Hap2 is maternal - assume Hap1 is paternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, maternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, paternal = TRUE, offspring_id = offspring_id)
        }
        if (max_diff == abs(score_1) & max_diff == abs(score_2) &  abs(score_3) < abs(score_4))  {
          message("Hap2 is maternal - assume Hap1 is paternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, paternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, maternal = TRUE, offspring_id = offspring_id)
        }
        if (max_diff == abs(score_1) & max_diff == abs(score_2) &  abs(score_3) > abs(score_4))  {
          message("Hap2 is maternal - assume Hap1 is paternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, maternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, paternal = TRUE, offspring_id = offspring_id)
        }
        if (max_diff == abs(score_1) & max_diff == abs(score_2) &  max_diff == abs(score_3) & max_diff == abs(score_4))  {
          stop("All values are equal")
          
        }
        if (max_diff == abs(score_1) & max_diff == abs(score_4) &  max_diff != abs(score_2) & max_diff != abs(score_3))  {
          message("Hap1 is paternal and Hap2 is maternal - both are max")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, paternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, maternal = TRUE, offspring_id = offspring_id)
        }
        if (max_diff == abs(score_2) & max_diff == abs(score_3) &  max_diff != abs(score_1) & max_diff != abs(score_4))  {
          message("Hap1 is maternal and Hap2 is paternal - both are max")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, maternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, paternal = TRUE, offspring_id = offspring_id)
          
        }
        
        
        merged_haps <- rbind(Offspring_Hap1_chr, Offspring_Hap2_chr)
        identified_dpc_haplo <- grep("_paternal$", rownames(merged_haps), value = TRUE)
        identified_maternal_haplo <- grep("_maternal$", rownames(merged_haps), value = TRUE)
        
        dpc_assigned <- merged_haps[identified_dpc_haplo, , drop = FALSE]
        maternal_assigned <- merged_haps[identified_maternal_haplo, , drop = FALSE]
        
        pat_results[[j]] <- as.data.frame(dpc_assigned)
        mat_results[[j]] <- as.data.frame(maternal_assigned)
        
        # Combine all scores into one data frame
        testing_methods <- do.call(rbind, scores)
        
        
        lengths_h1_p <- check_dist(Offspring_Hap1_chr, Dpc_Geno_chr, j =j)
        lengths_h1_p$Parent <- rep("h1_p")
        lengths_h1_m <- check_dist(Offspring_Hap1_chr, Maternal_Geno_chr, j =j)
        lengths_h1_m$Parent <- rep("h1_m")
        lengths_h2_p <- check_dist(Offspring_Hap2_chr, Dpc_Geno_chr, j =j)
        lengths_h2_p$Parent <- rep("h2_p")
        lengths_h2_m <- check_dist(Offspring_Hap2_chr, Maternal_Geno_chr, j =j)
        lengths_h2_m$Parent <- rep("h2_m")
        
        lengths_all <- rbind(lengths_h1_p, lengths_h1_m, lengths_h2_p, lengths_h2_m)
        
        
        
        
        method_results[[j]] <- as.data.frame(testing_methods)
        length_results[[j]] <- as.data.frame(lengths_all)
      }
      
      phased_results_methods[[i]] <- do.call(rbind,method_results)
      phased_results_lengths[[i]] <- do.call(rbind, length_results)
      
      pat_results <- do.call(cbind, pat_results)
      mat_results <- do.call(cbind, mat_results)
      results <- rbind(mat_results, pat_results)
      flipped_haplotypes[[i]] <- as.data.frame(results)
      
    }  
    results_methods_all <- do.call(rbind, phased_results_methods)
    results_lengths_all <- do.call(rbind, phased_results_lengths)
    results_flipped <- do.call(rbind, flipped_haplotypes)
    results <- list(results_methods_all = results_methods_all, results_lengths_all = results_lengths_all,
                    results_flipped = results_flipped)
  }
  
  return(results)
}

Real_FLIP <- TRUE_FALSE_flipping(Beagle_version = "Non", Real_haplotypes = TRUE, Data_type = "Real", pedigree = pedigree, method = "power_mean")
GE_phased_FLIP <- TRUE_FALSE_flipping(Beagle_version = "4.0", Real_haplotypes = FALSE, Data_type = "Phased_GE", pedigree = pedigree, method = "power_mean")
nGE_phased_FLIP <- TRUE_FALSE_flipping(Beagle_version = "4.0", Real_haplotypes = FALSE, Data_type = "Phased_nGE", pedigree = pedigree, method = "power_mean")
Slov_FLIP <- TRUE_FALSE_flipping(Beagle_version = "Slov4.0", Real_haplotypes = FALSE, Data_type = "Phased_GE", pedigree = Workers_with_fathers_id, method = "power_mean")

#function for checking  
check_haplotype_postFlip <- function(complete_haplotypes = NULL, results = NULL, pedigree = NULL) {
  comparison_results <- data.frame()
  
  cols <- colnames(results)
  real_haplotypes_cols <- complete_haplotypes[, cols, drop = FALSE]
  
  offspring_ids <- pedigree$id
  
  for (i in seq_along(offspring_ids)) {
    offspring_id <- offspring_ids[i]
    
    # Pull rows starting with offspring_id from the real haplotypes and results
    real_haplo_maternal <- rownames(real_haplotypes_cols)[grep(paste0("^", offspring_id, "_maternal"), rownames(real_haplotypes_cols))]
    real_haplo_maternal <- real_haplotypes_cols[real_haplo_maternal, , drop = FALSE]
    
    real_haplo_paternal <- rownames(real_haplotypes_cols)[grep(paste0("^", offspring_id, "_paternal"), rownames(real_haplotypes_cols))]
    real_haplo_paternal <- real_haplotypes_cols[real_haplo_paternal, , drop = FALSE]
    
    results_maternal <- rownames(results)[grep(paste0("^", offspring_id, "_maternal"), rownames(results))]
    results_maternal <- results[results_maternal, , drop = FALSE]
    
    results_paternal <- rownames(results)[grep(paste0("^", offspring_id, "_paternal"), rownames(results))]
    results_paternal <- results[results_paternal, , drop = FALSE]
    
    # Compare haplotypes

    total_elements <- length(results_maternal)
    
    maternal_comparison <- (sum(results_maternal == real_haplo_maternal, na.rm = TRUE)/total_elements) * 100
    paternal_comparison <- (sum(results_paternal == real_haplo_paternal, na.rm = TRUE)/total_elements) * 100
    
    # Combine results for current offspring_id
    current_comparison <- rbind(maternal_comparison, paternal_comparison)
    
    # Append current_comparison to comparison_results
    comparison_results <- rbind(comparison_results, current_comparison)
  }
  # Check conditions and print messages
  comparison_results_maternal <- rownames(comparison_results)[grep(paste0("maternal_"), rownames(comparison_results))]
  comparison_results_maternal <- comparison_results[comparison_results_maternal, , drop = FALSE]
  comparison_results_maternal <- t(comparison_results_maternal)
  
  comparison_results_paternal <- rownames(comparison_results)[grep(paste0("paternal_"), rownames(comparison_results))]
  comparison_results_paternal <- comparison_results[comparison_results_paternal, , drop = FALSE]
  comparison_results_paternal <- t(comparison_results_paternal)
  
  print("maternal mean")
  print(mean(comparison_results_maternal))
  print("maternal range")
  print(range(comparison_results_maternal))
  print("paternal mean")
  print(mean(comparison_results_paternal))
  print("paternal range")
  print(range(comparison_results_paternal))

  return(comparison_results)
}

colnames(Real_FLIP$real_results_flipped) <- colnames(real_haplotypes)

real_vc_nGE_FLIP <- check_haplotype_postFlip(complete_haplotypes = Real_FLIP$real_results_flipped, results = nGE_phased_FLIP$results_flipped, pedigree = pedigree)
real_vc_GE_FLIP <- check_haplotype_postFlip(complete_haplotypes = Real_FLIP$real_results_flipped, results = GE_phased_FLIP$results_flipped, pedigree = pedigree)


#####################################################################################
####                  Gametic/Haplotype Mendelian sampling 
# gi,1 = 1/2gm(i),1 + 1/2gm(i),2 + ri,1
# gi,2 = 1/2gf(i),1 + 1/2gf(i),2 + ri,2
#####################################################################################


##### Main function to calculate gametic relatedness ##################################
calculate_gametic_relatedness <- function(sorted_offspring_haplotypes, all_haplotypes, pedigree) {
  # Initialize an empty data frame to store results
  gametic_results <- data.frame(
    Offspring_ID = character(),
    Mother_ID = character(),
    DPC_ID = character(),
    Maternal_Mendelian = numeric(),
    Dpc_Mendelian = numeric(),
    stringsAsFactors = FALSE
  )
  
  cols <- colnames(sorted_offspring_haplotypes)
  all_haplotypes_cols <- all_haplotypes[, cols, drop = FALSE]
  
  # Loop through each offspring in the pedigree file
  for (i in 1:nrow(pedigree)) {
    offspring_id <- pedigree$id[i]
    mother_id <- pedigree$mother[i]
    dpc_id <- pedigree$dpc[i]
    
    
    #pull out the offspring maternal and paternal genotypes 
    offspring_maternal <- rownames(sorted_offspring_haplotypes)[grep(paste0("^", offspring_id, "_maternal"), rownames(sorted_offspring_haplotypes))]
    offspring_maternal <- sorted_offspring_haplotypes[offspring_maternal, , drop = FALSE]
    
    offspring_paternal <- rownames(sorted_offspring_haplotypes)[grep(paste0("^", offspring_id, "_paternal"), rownames(sorted_offspring_haplotypes))]
    offspring_paternal <- sorted_offspring_haplotypes[offspring_paternal, , drop = FALSE]
    
    # Separate parents into Haplotype 1 and Haplotype 2
    mother_row <- rownames(all_haplotypes_cols)[sapply(rownames(all_haplotypes_cols), function(x) strsplit(x,'_')[[1]][1]) %in% mother_id]
    mother_haplotypes_i <- all_haplotypes_cols[mother_row,]
    Maternal_Hap1 <- as.numeric(mother_haplotypes_i[1,])
    Maternal_Hap2 <- as.numeric(mother_haplotypes_i[2,])
    
    dpc_row <- rownames(all_haplotypes_cols)[sapply(rownames(all_haplotypes_cols), function(x) strsplit(x,'_')[[1]][1]) %in% dpc_id]
    dpc_haplotypes_i <- all_haplotypes_cols[dpc_row,]
    Dpc_Hap1 <- as.numeric(dpc_haplotypes_i[1,])
    Dpc_Hap2 <- as.numeric(dpc_haplotypes_i[2,])
    
    # Calculate parental average and Mendelian sampling for each haplotype
    ma_hap <- 0.5 * (Maternal_Hap1 + Maternal_Hap2)
    ri_hap1 <- offspring_maternal - ma_hap
    ri_hap1_sum <- sum(ri_hap1)
    
    pa_hap <- 0.5 * (Dpc_Hap1 + Dpc_Hap2)
    ri_hap2 <- offspring_paternal - pa_hap
    ri_hap2_sum <- sum(ri_hap2)
    
    # Add the values to the results data frame
    gametic_results <- rbind(gametic_results, data.frame(
      Offspring_ID = offspring_id,
      Mother_ID = mother_id,
      DPC_ID = dpc_id,
      Maternal_Mendelian = ri_hap1_sum,
      Dpc_Mendelian = ri_hap2_sum,
      stringsAsFactors = FALSE
    ))
  }
  
  return(gametic_results)
}


Gametic_real <- calculate_gametic_relatedness(sorted_offspring_haplotypes = Real_FLIP$real_results_flipped , all_haplotypes = real_haplotypes, pedigree = pedigree)
Gametic_real$Phasing <- "True"

Gametic_nGE_phased <- calculate_gametic_relatedness(sorted_offspring_haplotypes = nGE_phased_FLIP$results_flipped, all_haplotypes = nGE4.0_phasedhaplotypes, pedigree = pedigree)
Gametic_nGE_phased$Phasing <- "Phased_nGE"

Gametic_GE_phased <- calculate_gametic_relatedness(sorted_offspring_haplotypes = GE_phased_FLIP$results_flipped, all_haplotypes = GE4.0_phasedhaplotypes, pedigree = pedigree)
Gametic_GE_phased$Phasing <- "Phased_GE"

Gametic_Slov <- calculate_gametic_relatedness(sorted_offspring_haplotypes = Slov_FLIP$results_flipped, all_haplotypes = Slov4.0_phasedhaplotypes, pedigree = Slov_pedigree)
Gametic_Slov$Phasing <- "Real"

Gametic_df <- rbind(Gametic_real, Gametic_nGE_phased,  Gametic_GE_phased)

library(ggplot2)
library(gridExtra)
library(cowplot)

Plotting_Gametic <- function(df, phased_type, plotting_styles) {
  # Ensure the phased_type is valid
  if (!all(phased_type %in% c("Phased_GE", "Phased_nGE",  "Pre_Phased_nGE", "Pre_Phased_GE","True", "Real"))) {
    stop("Invalid phased_type. Choose either 'Phased_GE', 'Phased_nGE',  'Pre_Phased_nGE', 'Pre_Phased_GE', 'True' or 'Real'.")
  }
  
  # Ensure the plotting_styles are valid
  valid_styles <- c("histogram", "density", "scatter")
  if (!all(plotting_styles %in% valid_styles)) {
    stop("Invalid plotting_styles. Choose either 'histogram', 'density', or 'scatter'.")
  }
  
  # Define color and fill scales for consistency
  color_scale <- scale_color_manual(name = "Assigned parent haplotype",
                                    values = c("Maternal" = "blue", "Paternal" = "red"),
                                    labels = c("Maternal", "Paternal"))
  fill_scale <- scale_fill_manual(name = "Assigned parent haplotype",
                                  values = c("Maternal" = "blue", "Paternal" = "red"),
                                  labels = c("Maternal", "Paternal"))
  
  # Initialize an empty list to store plots
  plots <- list()
  
  # Generate plots for each phased_type
  for (phase in phased_type) {
    # Initialize list to store plots for current phased_type
    phase_plots <- list()
    
    # Generate plots for each plotting_style
    for (style in plotting_styles) {
      if (style == "histogram") {
        p <- ggplot(df[df$Phasing == phase, ], aes(x = Maternal_Mendelian, fill = "Maternal")) +
          geom_histogram(aes(y = ..count..), binwidth = 1, alpha = 0.7, position = 'identity') +
          geom_histogram(aes(x = Dpc_Mendelian, y = ..count.., fill = "Paternal"), 
                         binwidth = 1, alpha = 0.7, position = 'identity') +
          labs(title = paste(phase),
               x = "Mendelian Sampling", y = "Count") +
          theme_minimal() +
          theme(strip.text = element_blank(),
                plot.title = element_text(size = 20),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text = element_text(size = 14),
                legend.position = "none") +  # Remove legend from individual plots
          fill_scale +  # Apply fill scale
          facet_wrap(~ Phasing, ncol = 1)  # Separate plots by Phasing type
        
      } else if (style == "density") {
        p <- ggplot(df[df$Phasing == phase, ], aes(x = Maternal_Mendelian, color = "Maternal")) +
          geom_density(aes(y = ..scaled..), alpha = 0.7) +
          geom_density(aes(x = Dpc_Mendelian, y = ..scaled.., color = "Paternal"), alpha = 0.7) +
          labs(title = paste(phase),
               x = "Mendelian Sampling", y = "Density") +
          theme_minimal() +
          theme(strip.text = element_blank(),
                plot.title = element_text(size = 20),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text = element_text(size = 14),
                legend.position = "none") +  # Remove legend from individual plots
          color_scale +  # Apply color scale
          facet_wrap(~ Phasing, ncol = 1)  # Separate plots by Phasing type
        
      }
      # Add plot to current phased_type plots list
      phase_plots[[style]] <- p
    }
    
    # Combine plots for the current phased_type into a single row
    phase_plots_combined <- do.call(cowplot::plot_grid, c(phase_plots, list(nrow = 1)))
    
    # Add combined plots to the main plots list
    plots[[phase]] <- phase_plots_combined
  }
  
  # Combine all phased_type plots into a single row
  combined_plots <- do.call(cowplot::plot_grid, c(plots, list(nrow = 1)))
  
  # Add legend only if "Phased_GE" is in the phased_type
  if ("Phased_GE" %in% phased_type) {
    # Extract legend from one of the "Phased_GE" plots
    example_plot <- ggplot(df[df$Phasing == "Phased_GE", ], aes(x = Maternal_Mendelian, fill = "Maternal")) +
      geom_histogram(aes(y = ..count..), binwidth = 1, alpha = 0.7, position = 'identity') +
      geom_histogram(aes(x = Dpc_Mendelian, y = ..count.., fill = "Paternal"), 
                     binwidth = 1, alpha = 0.7, position = 'identity') +
      fill_scale +  # Apply fill scale
      theme_minimal()
    
    # Extract legend using get_legend from cowplot
    legend <- get_legend(
      example_plot + 
        theme(legend.position = "right",
              legend.title = element_text(size = 16),   # Increase legend title font size
              legend.text = element_text(size = 14))   # Increase legend text font size
    )
    
    # Combine plots and legend
    combined_plots_with_legend <- plot_grid(
      combined_plots,
      legend,
      ncol = 2,
      rel_widths = c(4, 1)  # Adjust widths as needed
    )
  } else {
    # Print combined plots without a legend
    combined_plots_with_legend <- combined_plots
  }
  
  # Print the combined plots with or without the legend
  print(combined_plots_with_legend)
}

Plotting_Gametic(df = Gametic_df, phased_type = c("True", "Phased_nGE", "Phased_GE"), plotting_styles = c("density"))
Plotting_Gametic(df = Gametic_Slov, phased_type = c("Real"), plotting_styles = c("density", "histogram"))
Plotting_Gametic(df = Gametic_nGE_phased, phased_type = c("Phased_nGE"), plotting_styles = c("density", "histogram"))



#####################################################################################
####                  Determine the number of fathers per queen 
#####################################################################################

#threshold since there could be genotyping errors 

calc_nPaternity_perQueen <- function(results, threshold, pedigree) {
  
  # Extract paternal haplotypes
  results_paternal <- rownames(results)[grep(paste0("_paternal"), rownames(results))]
  results_paternal <- results[results_paternal, , drop = FALSE]
  
  # Create a loop grouping together sisters using pedigree 
  queen_ids <- unique(pedigree$mother)
  
  # Initialize a data frame to store the results
  nPaternity <- data.frame(queen_id = queen_ids, num_fathers = rep(0, length(queen_ids)), stringsAsFactors = FALSE)
  
  for (i in 1:length(queen_ids)) {
    print(i)
    queen_id <- queen_ids[i]
    sister_worker_ids <- pedigree[pedigree$mother == queen_id, "id"]
    
    # Pull out the results_paternal with the worker IDs
    pattern <- paste0("^(", paste(sister_worker_ids, collapse="|"), ")_paternal$")
    matching_indices <- grep(pattern, rownames(results_paternal))
    sister_paternal <- results_paternal[matching_indices, , drop = FALSE]
    
    if (nrow(sister_paternal) == 0) next
    
    # Create similarity matrix 
    num_haplotypes <- nrow(sister_paternal)
    similarity_matrix <- matrix(0, nrow = num_haplotypes, ncol = num_haplotypes)
    rownames(similarity_matrix) <- rownames(sister_paternal)[1:num_haplotypes]
    colnames(similarity_matrix) <- rownames(sister_paternal)[1:num_haplotypes]
    
    # Function to calculate match percentage between two haplotypes
    calculate_similarity <- function(hap1, hap2) {
      sum(hap1 == hap2) / length(hap1)
    }
    
    # Initialize list to store father groups
    father_groups <- list()
    
    # Loop through each haplotype and create father groups based on similarity
    for (j in 1:num_haplotypes) {
      print(j)
      haplotype_j <- sister_paternal[j, ]
      
      # Check if this haplotype is already in a father group
      already_grouped <- any(sapply(father_groups, function(group) rownames(sister_paternal)[j] %in% group))
      if (already_grouped) next
      
      # Find matching fathers
      matching_fathers <- c()
      for (k in 1:num_haplotypes) {
        if (j != k) {
          haplotype_k <- sister_paternal[k, ]
          match_percentage <- sum(haplotype_j == haplotype_k) / length(haplotype_j)
          similarity_matrix[j, k] <- match_percentage
          if (match_percentage >= threshold) {
            matching_fathers <- c(matching_fathers, rownames(sister_paternal)[k])
          }
        }
      }
      
      # Create a new father group
      if (length(matching_fathers) > 0) {
        matching_fathers <- c(rownames(sister_paternal)[j], matching_fathers)
        father_groups <- append(father_groups, list(unique(matching_fathers)))
      } else {
        father_groups <- append(father_groups, list(rownames(sister_paternal)[j]))
      }
    }
    
    # Remove duplicate father groups
    unique_father_groups <- list()
    for (group in father_groups) {
      group_names <- unlist(group)
      is_duplicate <- any(sapply(unique_father_groups, function(ug) all(group_names %in% ug)))
      if (!is_duplicate) {
        unique_father_groups <- append(unique_father_groups, list(group_names))
      }
    }
    
    # Store the number of unique father groups in the data frame
    nPaternity$num_fathers[i] <- length(unique_father_groups)
  }
  
  return(nPaternity)
}

Real_nFathers <- calc_nPaternity_perQueen(results = Real_FLIP$real_results_flipped, threshold = 0.95, pedigree = pedigree)
nGE_nFathers <- calc_nPaternity_perQueen(results = nGE_phased_FLIP$results_flipped, threshold = 0.95, pedigree = pedigree)
GE_nFathers <- calc_nPaternity_perQueen(results = GE_phased_FLIP$results_flipped, threshold = 0.85, pedigree = pedigree)
Slov_nFathers85 <- calc_nPaternity_perQueen(results = Slov_FLIP$results_flipped, threshold = 0.85, pedigree = pedigree)

#################################
#check the accuracy of the assignments 
################################
#Load your Pop object that contains the fathers and the pedigree with the fathers id in there 
father_id_all <- pedigree$father
father_id_1 <- paste(father_id_all, "1", sep = "_")
father_haplotypes <- pullSnpHaplo(PopMerged, simParam = SP)

#since they're coded as identical diploids we can just use one of them
father_haplotypes_1 <- father_haplotypes[rownames(father_haplotypes) %in% father_id_1,]

#replace _1 with _paternal so that it can all be compared 
father_id_paternal <- unique(father_id_paternal)
father_all_unique <- unique(father_id_all)
father_id_paternal <- father_all_unique[order(father_all_unique)]
father_id_paternal <- paste(father_id_paternal, "paternal", sep = "_")

father_haplotype_pat <- father_haplotypes_1
rownames(father_haplotype_pat) <- father_id_paternal


#function for determining the accuracy of the determined number of patrilines
calc_nPaternity_Accuracy <- function(results, sister_threshold, pedigree, father_haplotypes, father_test_threshold) {
  
  # Extract paternal haplotypes
  results_paternal <- rownames(results)[grep(paste0("_paternal"), rownames(results))]
  results_paternal <- results[results_paternal, , drop = FALSE]
  
  # Create a loop grouping together sisters using pedigree 
  queen_ids <- unique(pedigree$mother)
  
  # Initialize a data frame to store the results
  # Initialize a data frame to store the results
  nPaternity <- data.frame(queen_id = queen_ids, 
                           num_fathers_actual = rep(0, length(queen_ids)),
                           num_workers = rep(0, length(queen_ids)),
                           num_fathers_estimated = rep(0, length(queen_ids)),
                           num_fathers_correct = rep(0, length(queen_ids)),
                           father_accuracy_threshold = father_test_threshold,
                           sister_threshold = sister_threshold,
                           stringsAsFactors = FALSE)
  
  for (i in 1:length(queen_ids)) {
    print(i)
    queen_id <- queen_ids[i]
    sister_worker_ids <- t(pedigree[pedigree$mother == queen_id, "id"])
    
    # Pull out the results_paternal with the worker IDs
    pattern <- paste0("^(", paste(sister_worker_ids, collapse = "|"), ")_paternal$")
    matching_indices <- grep(pattern, rownames(results_paternal))
    sister_paternal <- results_paternal[matching_indices, , drop = FALSE]
    
    nPaternity$num_workers[i] <- nrow(sister_paternal)
    
    # Pull out the father IDs for this queen from the pedigree
    father_ids <- t(unique(pedigree[pedigree$mother == queen_id, "father"]))
    
    # Get the actual father haplotypes for this queen
    pattern_fathers <- paste0("^(", paste(father_ids, collapse="|"), ")_paternal$")
    matching_fathers_indices <- grep(pattern_fathers, rownames(father_haplotypes))
    actual_fathers_haplotypes <- father_haplotypes[matching_fathers_indices, , drop = FALSE]
    
    # Ensure column names (alleles) of haplotypes match before comparison
    common_columns <- intersect(colnames(sister_paternal), colnames(actual_fathers_haplotypes))
    
    sister_paternal <- sister_paternal[, common_columns, drop = FALSE]
    actual_fathers_haplotypes <- actual_fathers_haplotypes[, common_columns, drop = FALSE]
    
    # Store the actual number of fathers
    nPaternity$num_fathers_actual[i] <- length(father_ids)
    
    # Create similarity matrix 
    num_worker_haplotypes <- nrow(sister_paternal)
    num_father_haplotypes <- nrow(actual_fathers_haplotypes)
    similarity_matrix <- matrix(0, nrow = num_father_haplotypes, ncol = num_worker_haplotypes)
    rownames(similarity_matrix) <- rownames(actual_fathers_haplotypes)[1:num_father_haplotypes]
    colnames(similarity_matrix) <- rownames(sister_paternal)[1:num_worker_haplotypes]
    
    # Updated calculate_similarity function with NA handling
    calculate_similarity <- function(hap1, hap2) {
      # Remove positions where either haplotype has an NA value
      valid_indices <- !is.na(hap1) & !is.na(hap2)
      
      if (sum(valid_indices) == 0) {
        return(NA)  # Return NA if no valid comparisons can be made
      }
      
      # Calculate similarity only on valid indices
      sum(hap1[valid_indices] == hap2[valid_indices]) / length(hap1[valid_indices])
    }
    
    # Calculate the similarity matrix
    for (j in 1:num_father_haplotypes) {
      for (k in 1:num_worker_haplotypes) {
        haplotype_j <- actual_fathers_haplotypes[j, ]
        haplotype_k <- sister_paternal[k, ]
        match_percentage <- calculate_similarity(haplotype_j, haplotype_k)
        similarity_matrix[j, k] <- match_percentage
      }
    }
    
    # Step 1: Assign workers to father groups based on similarity
    father_groups <- list()
    unmatched_workers <- c()
    
    for (k in 1:num_worker_haplotypes) {
      matching_fathers <- c()
      
      for (j in 1:num_father_haplotypes) {
        if (similarity_matrix[j, k] >= father_test_threshold) {
          matching_fathers <- c(matching_fathers, rownames(actual_fathers_haplotypes)[j])
        }
      }
      
      if (length(matching_fathers) > 1) {
        stop("Threshold is too relaxed, multiple fathers match a single worker.")
      }
      
      if (length(matching_fathers) == 1) {
        father <- matching_fathers[1]
        if (!father %in% names(father_groups)) {
          father_groups[[father]] <- c()
        }
        father_groups[[father]] <- c(father_groups[[father]], rownames(sister_paternal)[k])
      } else {
        unmatched_workers <- c(unmatched_workers, rownames(sister_paternal)[k])
      }
    }
    
    # Step 2: Compare unmatched workers with matched workers in father groups
    still_unmatched <- c()
    
    for (worker in unmatched_workers) {
      matched_to_group <- FALSE
      
      for (father in names(father_groups)) {
        for (sister_worker in father_groups[[father]]) {
          similarity_score <- calculate_similarity(sister_paternal[worker, ], sister_paternal[sister_worker, ])
          if (!is.na(similarity_score) && similarity_score >= sister_threshold) {
            father_groups[[father]] <- c(father_groups[[father]], worker)
            matched_to_group <- TRUE
            break
          }
        }
        if (matched_to_group) {
          break
        }
      }
      
      if (!matched_to_group) {
        still_unmatched <- c(still_unmatched, worker)
      }
    }
    
    # Step 3: Compare still unmatched workers with each other
    num_still_unmatched <- length(still_unmatched)
    new_groups <- list()
    
    if (num_still_unmatched > 1) {
      for (a in 1:(num_still_unmatched - 1)) {
        for (b in (a + 1):num_still_unmatched) {
          worker_a <- still_unmatched[a]
          worker_b <- still_unmatched[b]
          
          similarity_score <- calculate_similarity(sister_paternal[worker_a, ], sister_paternal[worker_b, ])
          if (!is.na(similarity_score) && similarity_score >= sister_threshold) {
            new_groups[[paste0("group_", length(father_groups) + length(new_groups) + 1)]] <- c(worker_a, worker_b)
            still_unmatched <- still_unmatched[!(still_unmatched %in% c(worker_a, worker_b))]
          }
        }
      }
    }
    
    # Assign remaining unmatched workers to their own groups
    for (worker in still_unmatched) {
      new_groups[[paste0("group_", length(father_groups) + length(new_groups) + 1)]] <- c(worker)
    }
    
    # Combine father groups and new groups
    father_groups <- c(father_groups, new_groups)
    
    # Store the number of estimated fathers
    nPaternity$num_fathers_estimated[i] <- length(father_groups)
    
    # Count correct fathers (those that match the actual fathers)
    correct_fathers <- intersect(names(father_groups), rownames(actual_fathers_haplotypes))
    nPaternity$num_fathers_correct[i] <- length(correct_fathers)
  }
  
  return(nPaternity)
}


# Define the function to run the test with dynamic results and thresholds, running them separately
run_paternity_tests <- function(results_arg, sister_thresholds, father_test_thresholds, father_haplotypes) {
  
  # Initialize a list to store the results
  results_list <- list()
  
  # Loop over each combination of sister_threshold and father_test_threshold
  for (sister_threshold in sister_thresholds) {
    for (father_test_threshold in father_test_thresholds) {
      
      # Create a descriptive name for the result
      result_name <- paste0("sis", sister_threshold, "_fat", father_test_threshold)
      
      # Print to ensure thresholds are being passed correctly
      print(paste("Running for sister_threshold:", sister_threshold, "and father_test_threshold:", father_test_threshold))
      
      # Run the function separately for each combination of thresholds
      result <- calc_nPaternity_Accuracy(
        results = results_arg,  # Pass the dynamic results argument here
        sister_threshold = sister_threshold,  # sister threshold
        pedigree = pedigree,
        father_haplotypes = father_haplotypes,
        father_test_threshold = father_test_threshold  # father test threshold
      )
      
      # Store the result with the unique name
      results_list[[result_name]] <- result
      
      # Print to confirm each result has been stored
      print(paste("Finished:", result_name))
    }
  }
  
  # Return the list of results
  return(results_list)
}

# Define the thresholds you want to test
sister_thresholds <- c(1, 0.95, 0.85, 0.75)
father_test_thresholds <- c(0.9, 0.95, 1)

# Run the function and specify the results dataset to use
Real_FatherTest_Results <- run_paternity_tests(
  results_arg = Real_FLIP$real_results_flipped, 
  sister_thresholds = 1,
  father_test_thresholds = 1,
  father_haplotypes = father_haplotype_pat
)

nGE_Father_Test_Results <- run_paternity_tests(
  results_arg = nGE_phased_FLIP$results_flipped, 
  sister_thresholds = sister_thresholds,
  father_test_thresholds = father_test_thresholds,
  father_haplotypes = father_haplotype_pat
)

GE_Father_Test_Results <- run_paternity_tests(
  results_arg = GE_phased_FLIP$results_flipped,  
  sister_thresholds = sister_thresholds,
  father_test_thresholds = father_test_thresholds,
  father_haplotypes = father_haplotype_pat
)

Slov_FatherTest_results <- run_paternity_tests_Slov(
  results_arg = Slov_FLIP$results_flipped,  
  sister_thresholds = sister_thresholds,
  father_test_thresholds = father_test_thresholds,
  father_haplotypes = Slov_Dpc_Haplo
)


# Define the plotting function
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # Ensure this is loaded

# Define the plotting function
plot_paternity_accuracy <- function(Results) {
  
  # Ensure Results has the correct column names and data types
  Results <- as.data.frame(Results)
  
  # Add columns for `sister_threshold` and `father_test_threshold` by extracting them from rownames if they are not already present
  if (!("sister_threshold" %in% colnames(Results))) {
    Results <- Results %>%
      rownames_to_column("threshold_combination") %>%
      separate(threshold_combination, into = c("prefix", "sis", "fat"), sep = "_", remove = FALSE) %>%
      mutate(sister_threshold = as.numeric(gsub("sis", "", sis)),
             father_test_threshold = as.numeric(gsub("fat", "", fat))) %>%
      select(-prefix, -sis, -fat, -threshold_combination)
  }
  
  # Calculate the accuracy as num_fathers_correct / num_fathers_actual
  Results <- Results %>%
    mutate(accuracy = (num_fathers_correct/ num_fathers_estimated) *100,
           accuracy = ifelse(is.nan(accuracy), 0, accuracy))  # Handle NaNs by replacing them with 0
  
  summarised_results <- Results %>%
    group_by(sister_threshold, father_accuracy_threshold) %>%
    summarize(average_accuracy = mean(accuracy, na.rm = TRUE), .groups = 'drop')
  
  
  # Plot 2: Accuracy of determined sires
  plot2 <- ggplot(summarised_results, aes(x = as.factor(sister_threshold), y = average_accuracy, fill = as.factor(father_accuracy_threshold))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Sister Threshold", y = "Accuracy of the estimated patrilines (%)", fill = "Father Threshold") +
    theme_minimal()+
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.position = "right"  # Ensure this is inside the theme() function
    )
  
  # Return the combined plot
  return(plot2)
}

plot_paternity_accuracy(nGE_Father_Test_Results)
plot_paternity_accuracy(GE_Father_Test_Results)
plot_paternity_accuracy(Real_FatherTest_Results)




# The real data is a bit different in structure and names so needs a different function
calc_nPaternity_Accuracy_SLOV <- function(results, sister_threshold, pedigree, father_haplotypes, father_test_threshold) {
  
  # Extract paternal haplotypes
  results_paternal <- rownames(results)[grep(paste0("_paternal"), rownames(results))]
  results_paternal <- results[results_paternal, , drop = FALSE]
  
  # Create a loop grouping together sisters using pedigree 
  queen_ids <- unique(pedigree$mother)
  
  # Initialize a data frame to store the results
  nPaternity <- data.frame(queen_id = queen_ids, 
                           num_fathers_actual = rep(0, length(queen_ids)),
                           num_workers = rep(0, length(queen_ids)),
                           num_fathers_estimated = rep(0, length(queen_ids)),
                           num_fathers_correct = rep(0, length(queen_ids)),
                           father_accuracy_threshold = father_test_threshold,
                           sister_threshold = sister_threshold,
                           stringsAsFactors = FALSE)
  
  for (i in 1:length(queen_ids)) {
    print(i)
    queen_id <- queen_ids[i]
    sister_worker_ids <- t(pedigree[pedigree$mother == queen_id, "id"])
    
    # Pull out the results_paternal with the worker IDs
    pattern <- paste0("^(", paste(sister_worker_ids, collapse = "|"), ")_paternal$")
    matching_indices <- grep(pattern, rownames(results_paternal))
    sister_paternal <- results_paternal[matching_indices, , drop = FALSE]
    
    nPaternity$num_workers[i] <- nrow(sister_paternal)
    
    actual_fathers_haplotypes <- father_haplotypes
    
    # Ensure column names (alleles) of haplotypes match before comparison
    common_columns <- intersect(colnames(sister_paternal), colnames(actual_fathers_haplotypes))
    
    sister_paternal <- sister_paternal[, common_columns, drop = FALSE]
    actual_fathers_haplotypes <- actual_fathers_haplotypes[, common_columns, drop = FALSE]
    
    # Store the actual number of fathers
    father_ids <- rownames(actual_fathers_haplotypes)
    nPaternity$num_fathers_actual[i] <- length(father_ids)
    
    # Create similarity matrix 
    num_worker_haplotypes <- nrow(sister_paternal)
    num_father_haplotypes <- nrow(actual_fathers_haplotypes)
    similarity_matrix <- matrix(0, nrow = num_father_haplotypes, ncol = num_worker_haplotypes)
    rownames(similarity_matrix) <- rownames(actual_fathers_haplotypes)
    colnames(similarity_matrix) <- rownames(sister_paternal)
    
    # Updated calculate_similarity function with NA handling
    calculate_similarity <- function(hap1, hap2) {
      # Remove positions where either haplotype has an NA value
      valid_indices <- !is.na(hap1) & !is.na(hap2)
      
      if (sum(valid_indices) == 0) {
        return(NA)  # Return NA if no valid comparisons can be made
      }
      
      # Calculate similarity only on valid indices
      sum(hap1[valid_indices] == hap2[valid_indices]) / length(hap1[valid_indices])
    }
    
    # Calculate the similarity matrix
    for (j in 1:num_father_haplotypes) {
      for (k in 1:num_worker_haplotypes) {
        haplotype_j <- actual_fathers_haplotypes[j, ]
        haplotype_k <- sister_paternal[k, ]
        match_percentage <- calculate_similarity(haplotype_j, haplotype_k)
        similarity_matrix[j, k] <- match_percentage
      }
    }
    
    # Step 1: Assign workers to father groups based on similarity
    father_groups <- list()
    unmatched_workers <- c()
    
    for (k in 1:num_worker_haplotypes) {
      matching_fathers <- c()
      
      for (j in 1:num_father_haplotypes) {
        if (similarity_matrix[j, k] >= father_test_threshold) {
          matching_fathers <- c(matching_fathers, rownames(actual_fathers_haplotypes)[j])
        }
      }
      
      if (length(matching_fathers) > 1) {
        stop("Threshold is too relaxed, multiple fathers match a single worker.")
      }
      
      if (length(matching_fathers) == 1) {
        father <- matching_fathers[1]
        if (!father %in% names(father_groups)) {
          father_groups[[father]] <- c()
        }
        father_groups[[father]] <- c(father_groups[[father]], rownames(sister_paternal)[k])
      } else {
        unmatched_workers <- c(unmatched_workers, rownames(sister_paternal)[k])
      }
    }
    
    # Step 2: Compare unmatched workers with matched workers in father groups
    still_unmatched <- c()
    
    for (worker in unmatched_workers) {
      matched_to_group <- FALSE
      
      for (father in names(father_groups)) {
        for (sister_worker in father_groups[[father]]) {
          similarity_score <- calculate_similarity(sister_paternal[worker, ], sister_paternal[sister_worker, ])
          if (!is.na(similarity_score) && similarity_score >= sister_threshold) {
            father_groups[[father]] <- c(father_groups[[father]], worker)
            matched_to_group <- TRUE
            break
          }
        }
        if (matched_to_group) {
          break
        }
      }
      
      if (!matched_to_group) {
        still_unmatched <- c(still_unmatched, worker)
      }
    }
    
    # Step 3: Compare still unmatched workers with each other
    num_still_unmatched <- length(still_unmatched)
    new_groups <- list()
    
    if (num_still_unmatched > 1) {
      for (a in 1:(num_still_unmatched - 1)) {
        for (b in (a + 1):num_still_unmatched) {
          worker_a <- still_unmatched[a]
          worker_b <- still_unmatched[b]
          
          similarity_score <- calculate_similarity(sister_paternal[worker_a, ], sister_paternal[worker_b, ])
          if (!is.na(similarity_score) && similarity_score >= sister_threshold) {
            new_groups[[paste0("group_", length(father_groups) + length(new_groups) + 1)]] <- c(worker_a, worker_b)
            still_unmatched <- still_unmatched[!(still_unmatched %in% c(worker_a, worker_b))]
          }
        }
      }
    }
    
    # Assign remaining unmatched workers to their own groups
    for (worker in still_unmatched) {
      new_groups[[paste0("group_", length(father_groups) + length(new_groups) + 1)]] <- c(worker)
    }
    
    # Combine father groups and new groups
    father_groups <- c(father_groups, new_groups)
    
    # Store the number of estimated fathers
    nPaternity$num_fathers_estimated[i] <- length(father_groups)
    
    # Count correct fathers (those that match the actual fathers)
    correct_fathers <- intersect(names(father_groups), rownames(actual_fathers_haplotypes))
    nPaternity$num_fathers_correct[i] <- length(correct_fathers)
  }
  
  return(nPaternity)
}
run_paternity_tests_Slov <- function(results_arg, sister_thresholds, father_test_thresholds, father_haplotypes) {
  
  # Initialize a list to store the results
  results_list <- list()
  
  # Loop over each combination of sister_threshold and father_test_threshold
  for (sister_threshold in sister_thresholds) {
    for (father_test_threshold in father_test_thresholds) {
      
      # Create a descriptive name for the result
      result_name <- paste0("sis", sister_threshold, "_fat", father_test_threshold)
      
      # Print to ensure thresholds are being passed correctly
      print(paste("Running for sister_threshold:", sister_threshold, "and father_test_threshold:", father_test_threshold))
      
      # Run the function separately for each combination of thresholds
      result <- calc_nPaternity_Accuracy_SLOV(
        results = results_arg,  # Pass the dynamic results argument here
        sister_threshold = sister_threshold,  # sister threshold
        pedigree = pedigree,
        father_haplotypes = father_haplotypes,
        father_test_threshold = father_test_threshold  # father test threshold
      )
      
      # Store the result with the unique name
      results_list[[result_name]] <- result
      
      # Print to confirm each result has been stored
      print(paste("Finished:", result_name))
    }
  }
  
  # Return the list of results
  return(results_list)
}



# a plot for the paternity 
plot_paternity_number <- function(Results) {
  
  # Ensure Results has the correct column names and data types
  Results <- as.data.frame(Results)
  
  # Add columns for `sister_threshold` and `father_accuracy_threshold` by extracting them from rownames if they are not already present
  if (!("sister_threshold" %in% colnames(Results)) || !("father_accuracy_threshold" %in% colnames(Results))) {
    Results <- Results %>%
      rownames_to_column("threshold_combination") %>%
      separate(threshold_combination, into = c("prefix", "sis", "fat"), sep = "_", remove = FALSE) %>%
      mutate(sister_threshold = as.numeric(gsub("sis", "", sis)),
             father_accuracy_threshold = as.numeric(gsub("fat", "", fat))) %>%
      select(-prefix, -sis, -fat, -threshold_combination)
  }
  
  # Convert sister_threshold to factor and handle worker_label
  Results <- Results %>%
    mutate(sister_threshold = as.factor(sister_threshold),
           worker_label = "nWorkers",
           known_label = "nKnown Patrilines")  # Add a label for the legend
  
  # Calculate the accuracy as num_fathers_correct / num_fathers_actual
  Results <- Results %>%
    mutate(accuracy = (num_fathers_correct / num_fathers_actual) * 100,
           accuracy = ifelse(is.nan(accuracy), 0, accuracy))  # Handle NaNs by replacing them with 0
  
  # Create a unique shape mapping for the father accuracy and worker labels
  shape_values <- c(16, 17, 18, 19)  # Adjust these values based on your needs
  
  # Plot 1: Scatter plot of number of determined sires (for each queen on X-axis)
  plot1 <- ggplot(Results, aes(x = as.factor(queen_id), y = num_fathers_estimated, 
                               color = sister_threshold, shape = as.factor(father_accuracy_threshold))) +
    geom_point(size = 7, alpha = 0.7, position = position_dodge(width = 0.5)) +  # Points for number of determined sires
    geom_point(aes(y = num_workers, shape = worker_label), color = "black", size = 7, stroke = 1, position = position_dodge(width = 0.5)) +  # Large black crosses for num_workers
    geom_point(aes(y = num_fathers_actual, shape = known_label), color = "red", size = 7, stroke = 2, position = position_dodge(width = 0.5)) + 
    scale_color_discrete(name = "Sister Threshold") +  # Scale for color legend
    scale_shape_manual(values = c(16, 17,4, 4), name = "Father Threshold") +  # Manual shapes for father accuracy and worker labels
    labs(x = "Colony ID", y = "Number of Determined Patrilines") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.position = "right"  # Ensure this is inside the theme() function
    )
  
  # Return the combined plot
  return(plot1)
}

plot_paternity_number(Slov_FatherTest_results)
plot_paternity_number(nGE_Father_Test_Results)
plot_paternity_number(GE_Father_Test_Results)




