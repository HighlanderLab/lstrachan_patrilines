#######################################################################################################################
#***************** Assign the Haplotype parent-of-origin assignments **************
#######################################################################################################################

#In this section we'll use 2 routes:

#.  - Route 1: use pedigree reconstruction and use both dam and sire pedigree ID to assign haplotype parental origins
#Use haplotypes from phasing using pedigree

#.  - Route 2: DO NOT use pedigree reconstruction and use ONLY dam pedigree ID to assign haplotype parental origins

#How does this function work? 
# - Iterates through all 16 chromosomes
# - Computes difference between the offspring's haplotypes (coded as 0/1) and each of the parents' genotypes (coded as 0/1/2).

#Offspring haplotype - Parent Genotype = difference 

#                      Parent Geno = 0. Parent Geno = 1. Parent Geno = 2
# Offspring Haplo = 0.       0                -1                -2
# Offspring Haplo = 1        1                 0                -1

#Difference of 1/-2 == allelic MISMATCH
#Difference  0/-1 == allelic MATCH

#Summarise difference across loci with a power mean of allelic mismatches Score = sum(MISMATCH)^2/length(MISMATCH)
#The higher the score the less likely the haplotype was inherited from that parent 
#The scores were checked using conditional arguments to improve the parent-of-origin assignment accuracy

###############################################################################################################
#Practical example:

#Offspring with haplotypes [0,1,0,1]
#dam with genotypes     [0,1,2,0]
#DPQ with genotypes        [1,0,1,1]

# Using Table above, we calculate mismatches for the dam:
#   0 - 0 = 0     --> MATCH
#   1 - 1 = 0     --> MATCH
#   0 - 2 = -2    --> MISMATCH
#   1 - 0 = 1     --> MISMATCH
# Resulting differences: [0, 0, -2, 1]
# Squaring these:        [0, 0,  4, 1]
# Power mean score:      (0 + 0 + 4 + 1) / 4 = 1.25
#
# Now for the father:
#   0 - 1 =  1     --> MISMATCH
#   1 - 0 = -1     --> MATCH
#   0 - 1 = -1     --> MATCH
#   1 - 1 =  0     --> MATCH
# Resulting differences: [1, -1, -1, 0]
# Squaring these:        [1,  1,  1, 0]
# Power mean score:      (1 + 1 + 1 + 0) / 4 = 0.75
#
# Conclusion:
# The lower score for the father (0.75 vs 1.25) suggests that the haplotype is paternally derived.

#######################################################################################################################
# --- Clear Workspace ---
rm(list = ls())
# --- Libraries ---
{
  library(Eagle)
  library(tidyr)
  library(AlphaSimR)
  library(SIMplyBee)
  library(readr)
  library(genio)
  library(ggplot2)
  library(dplyr)
  library(vcfR)
  library(tibble)
}



workingDir = "~/Desktop/lstrachan_patrilines"
setwd(workingDir)

#####################################################################################
#************************** FUNCTIONS ******************************************
#####################################################################################
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

check_haplotype <- function(true_haplotypes = NULL, results = NULL, pedigree = NULL) {
  
  comparison_results <- data.frame()
  
  # Align columns between true and predicted data
  cols <- colnames(results)
  true_haplotypes_cols <- true_haplotypes[, cols, drop = FALSE]
  
  offspring_ids <- pedigree$id
  
  for (i in seq_along(offspring_ids)) {
    offspring_id <- offspring_ids[i]
    
    # 1. Extract TRUE haplotypes
    real_haplo_maternal <- true_haplotypes_cols[grep(paste0("^", offspring_id, "_1"), 
                                                     rownames(true_haplotypes_cols)), , drop = FALSE]
    real_haplo_paternal <- true_haplotypes_cols[grep(paste0("^", offspring_id, "_2"), 
                                                     rownames(true_haplotypes_cols)), , drop = FALSE]
    
    # 2. Extract PREDICTED haplotypes
    results_maternal <- results[grep(paste0("^", offspring_id, "_1"), 
                                     rownames(results)), , drop = FALSE]
    results_paternal <- results[grep(paste0("^", offspring_id, "_2"), 
                                     rownames(results)), , drop = FALSE]
    
    total_elements <- length(results_maternal)
    
    # 3. Compare maternal predictions
    maternal_comparison <- as.data.frame(matrix(NA, nrow = 1, ncol = 2))
    colnames(maternal_comparison) <- c("_1", "_2")
    # CRITICAL: Re-add rownames so grep works later
    rownames(maternal_comparison) <- rownames(results_maternal) 
    
    maternal_comparison$`_1` <- (sum(results_maternal == real_haplo_maternal, na.rm = TRUE) / total_elements) * 100
    maternal_comparison$`_2` <- (sum(results_maternal == real_haplo_paternal, na.rm = TRUE) / total_elements) * 100
    
    # 4. Compare paternal predictions
    paternal_comparison <- as.data.frame(matrix(NA, nrow = 1, ncol = 2))
    colnames(paternal_comparison) <- c("_1", "_2")
    # CRITICAL: Re-add rownames so grep works later
    rownames(paternal_comparison) <- rownames(results_paternal)
    
    paternal_comparison$`_1` <- (sum(results_paternal == real_haplo_maternal, na.rm = TRUE) / total_elements) * 100
    paternal_comparison$`_2` <- (sum(results_paternal == real_haplo_paternal, na.rm = TRUE) / total_elements) * 100
    
    # 5. Store results
    current_comparison <- rbind(maternal_comparison, paternal_comparison)
    comparison_results <- rbind(comparison_results, current_comparison)
  }
  
  # Split results into maternal and paternal prediction sets
  # This relies on the rownames we assigned inside the loop
  comparison_results_maternal <- comparison_results[grep("_1$", rownames(comparison_results)), , drop = FALSE]
  comparison_results_paternal <- comparison_results[grep("_2$", rownames(comparison_results)), , drop = FALSE]
  
  # Create summary table
  summary_stats <- data.frame(
    comparison = c(
      "true_maternal vs predicted_maternal",
      "true_paternal vs predicted_maternal",
      "true_maternal vs predicted_paternal",
      "true_paternal vs predicted_paternal"
    ),
    percent = c(
      mean(comparison_results_maternal$`_1`, na.rm = TRUE),
      mean(comparison_results_maternal$`_2`, na.rm = TRUE),
      mean(comparison_results_paternal$`_1`, na.rm = TRUE),
      mean(comparison_results_paternal$`_2`, na.rm = TRUE)
    )
  )
  
  print(summary_stats)
  
  return(list(
    comparison_table = comparison_results,
    summary = summary_stats
  ))
}
# ------------------------------------------------------------
# KEY STEPS IN `check_haplotype()`:
# ------------------------------------------------------------
# 1. Align true haplotypes with result columns.
#
# 2. For each offspring:
#     - Extract true and predicted maternal and paternal haplotypes.
#     - Compare all true vs predicted combinations.
#
# 3. Compute MATCHING PERCENTAGE:
#     - % of identical SNP/haplotype entries between two matrices.
#     - (matches / total elements) × 100 → direct concordance.
#
# 4. Store per-offspring maternal and paternal comparison results.
#
# 5. Summarise mean % matching across all individuals.
#
# OUTPUT:
# - Per-offspring haplotype concordance (% matching)
# - Summary table of mean true vs predicted matching percentages
# ------------------------------------------------------------

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

Route1_flipping <- function(Data_type = NULL, pedigree = NULL, perfect_haplotypes = NULL, method = NULL){
  
  if (perfect_haplotypes == TRUE){
    haplotypes <- true_haplotypes
    map <- true_map
    
    real_results_methods <- list()
    real_results_lengths <- list()
    flipped_haplotypes <- list()
    
    
    for (i in 1:nrow(pedigree)) {
      offspring_id <- as.character(pedigree$id[i])
      dam_id <- as.character(pedigree$mother[i])
      sire_id <- as.character(pedigree$dpc[i])
      
      offspring_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x,'_')[[1]][1]) %in% offspring_id]
      offspring_haplotypes_i <- haplotypes[offspring_row,]
      Offspring_Hap1 <- t(as.data.frame(offspring_haplotypes_i[1,]))
      Offspring_Hap2 <- t(as.data.frame(offspring_haplotypes_i[2,]))
      
      dam_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x,'_')[[1]][1]) %in% dam_id]
      dam_haplotypes_i <- haplotypes[dam_row,]
      Maternal_Hap1 <- t(as.data.frame(dam_haplotypes_i[1,]))
      rownames(Maternal_Hap1) <- dam_id
      Maternal_Hap2 <- t(as.data.frame(dam_haplotypes_i[2,]))
      rownames(Maternal_Hap2) <- dam_id
      Maternal_Geno <- matrix(data = Maternal_Hap1 + Maternal_Hap2, nrow = 1)
      colnames(Maternal_Geno) <- colnames(Maternal_Hap1)
      rownames(Maternal_Geno) <- dam_id
      
      sire_row <- rownames(true_haplotypes)[sapply(rownames(true_haplotypes), function(x) strsplit(x,'_')[[1]][1]) %in% sire_id]
      sire_haplotypes_i <- true_haplotypes[sire_row,]
      sire_Hap1 <- t(as.data.frame(sire_haplotypes_i[1,]))
      rownames(sire_Hap1) <- sire_id
      sire_Hap2 <- t(as.data.frame(sire_haplotypes_i[2,]))
      rownames(sire_Hap2) <- sire_id
      sire_Geno <- matrix(data = sire_Hap1 + sire_Hap2, nrow = 1)
      colnames(sire_Geno) <- colnames(sire_Hap1)
      rownames(sire_Geno) <- sire_id
      
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
        sire_Geno_chr <- sire_Geno[, colnames(sire_Geno) %in% Chr_markerNames, drop = FALSE]
        
        scores_h1_p <- weighing_haplotypes(Offspring_Hap1_chr, sire_Geno_chr, method = method, j =j)
        scores_h1_m <- weighing_haplotypes(Offspring_Hap1_chr, Maternal_Geno_chr, method = method, j =j)
        scores_h2_p <- weighing_haplotypes(Offspring_Hap2_chr, sire_Geno_chr, method = method, j =j)
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
        identified_sire_haplo <- grep("_paternal$", rownames(merged_haps), value = TRUE)
        identified_maternal_haplo <- grep("_maternal$", rownames(merged_haps), value = TRUE)
        
        sire_assigned <- merged_haps[identified_sire_haplo, , drop = FALSE]
        maternal_assigned <- merged_haps[identified_maternal_haplo, , drop = FALSE]
        
        pat_results[[j]] <- as.data.frame(sire_assigned)
        mat_results[[j]] <- as.data.frame(maternal_assigned)
        
        # Combine all scores into one data frame
        testing_methods <- do.call(rbind, scores)
        
        
        lengths_h1_p <- check_dist(Offspring_Hap1_chr, sire_Geno_chr, j =j)
        lengths_h1_p$Parent <- rep("h1_p")
        lengths_h1_m <- check_dist(Offspring_Hap1_chr, Maternal_Geno_chr, j =j)
        lengths_h1_m$Parent <- rep("h1_m")
        lengths_h2_p <- check_dist(Offspring_Hap2_chr, sire_Geno_chr, j =j)
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
  
  if (perfect_haplotypes == FALSE){
    
    phased_results_methods <- list()
    phased_results_lengths <- list()
    flipped_haplotypes <- list()
    
    if (Data_type == "NoGE_SNP2k"){
      haplotypes <- NoGE_SNP2k_PhasedHaplotypes_recPed
      map <- NoGE_map_2k
    } else if (Data_type == "WithGE_SNP2k"){
      haplotypes <- WithGE_SNP2k_PhasedHaplotypes_recPed
      map <- WithGE_map_2k
    } else if (Data_type == "NoGE_SNP50k"){
      haplotypes <- NoGE_SNP50k_PhasedHaplotypes_recPed
      map <- NoGE_map_50k
    }else if (Data_type == "WithGE_SNP50k"){
      haplotypes <- WithGE_SNP50k_PhasedHaplotypes_recPed
      map <- WithGE_map_50k
    }else if (Data_type == "Real_Slov_data"){
      haplotypes <- Slov_phasedhaplotypes_phasedrecPed
      map <- Slov_map
    } else (stop("Arguments not met"))
    
    for (i in 1:nrow(pedigree)) {
      offspring_id <- pedigree$id[i]
      dam_id <- pedigree$dam[i]
      sire_id <- pedigree$sire[i]
      
      offspring_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x, '_')[[1]][1]) %in% offspring_id]
      offspring_haplotypes_i <- haplotypes[offspring_row, ]
      Offspring_Hap1 <- as.data.frame(offspring_haplotypes_i[1, , drop = FALSE])
      Offspring_Hap2 <- as.data.frame(offspring_haplotypes_i[2, , drop = FALSE])
      
      dam_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x, '_')[[1]][1]) %in% dam_id]
      dam_haplotypes_i <- haplotypes[dam_row, ]
      Maternal_Hap1 <- as.data.frame(dam_haplotypes_i[1, , drop = FALSE])
      rownames(Maternal_Hap1) <- dam_id
      Maternal_Hap2 <- as.data.frame(dam_haplotypes_i[2, , drop = FALSE])
      rownames(Maternal_Hap2) <- dam_id
      Maternal_Geno <- matrix(data = Maternal_Hap1 + Maternal_Hap2, nrow = 1)
      colnames(Maternal_Geno) <- colnames(Maternal_Hap1)
      rownames(Maternal_Geno) <- dam_id
      
      sire_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x, '_')[[1]][1]) %in% sire_id]
      sire_haplotypes_i <- haplotypes[sire_row, ]
      sire_Hap1 <- as.data.frame(sire_haplotypes_i[1, , drop = FALSE])
      rownames(sire_Hap1) <- sire_id
      sire_Hap2 <- as.data.frame(sire_haplotypes_i[2, , drop = FALSE])
      rownames(sire_Hap2) <- sire_id
      sire_Geno <- matrix(data = sire_Hap1 + sire_Hap2, nrow = 1)
      colnames(sire_Geno) <- colnames(sire_Hap1)
      rownames(sire_Geno) <- sire_id
      
      method_results <- list()
      length_results <- list()
      pat_results <- list()
      mat_results <- list() 
      
      for (j in 1){#unique(map[, 1])) {
        print(x = c(i,j))
        Chr_map <- map[map$V1 == j, ]
        Chr_markerNames <- Chr_map[, 2]
        
        Offspring_Hap1_chr <- Offspring_Hap1[, colnames(Offspring_Hap1) %in% Chr_markerNames, drop = FALSE]
        Offspring_Hap2_chr <- Offspring_Hap2[, colnames(Offspring_Hap2) %in% Chr_markerNames, drop = FALSE]
        Maternal_Geno_chr <- Maternal_Geno[, colnames(Maternal_Geno) %in% Chr_markerNames, drop = FALSE]
        sire_Geno_chr <- sire_Geno[, colnames(sire_Geno) %in% Chr_markerNames, drop = FALSE]
        
        scores_h1_p <- weighing_haplotypes(Offspring_Hap1_chr, sire_Geno_chr, method = method, j =j)
        scores_h1_m <- weighing_haplotypes(Offspring_Hap1_chr, Maternal_Geno_chr, method = method, j =j)
        scores_h2_p <- weighing_haplotypes(Offspring_Hap2_chr, sire_Geno_chr, method = method, j =j)
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
        identified_sire_haplo <- grep("_paternal$", rownames(merged_haps), value = TRUE)
        identified_maternal_haplo <- grep("_maternal$", rownames(merged_haps), value = TRUE)
        
        sire_assigned <- merged_haps[identified_sire_haplo, , drop = FALSE]
        maternal_assigned <- merged_haps[identified_maternal_haplo, , drop = FALSE]
        
        pat_results[[j]] <- as.data.frame(sire_assigned)
        mat_results[[j]] <- as.data.frame(maternal_assigned)
        
        # Combine all scores into one data frame
        testing_methods <- do.call(rbind, scores)
        
        
        lengths_h1_p <- check_dist(Offspring_Hap1_chr, sire_Geno_chr, j =j)
        lengths_h1_p$Parent <- rep("h1_p")
        lengths_h1_m <- check_dist(Offspring_Hap1_chr, Maternal_Geno_chr, j =j)
        lengths_h1_m$Parent <- rep("h1_m")
        lengths_h2_p <- check_dist(Offspring_Hap2_chr, sire_Geno_chr, j =j)
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

Route2_flipping <- function(Data_type = NULL, pedigree = NULL, perfect_haplotypes = NULL, method = NULL) {
  
  # First block: Real haplotypes case
  if (perfect_haplotypes == TRUE & Data_type == "True") {
    haplotypes <- true_haplotypes
    map <- true_map
    
    real_results_methods <- list()
    real_results_lengths <- list()
    flipped_haplotypes <- list()
    
    for (i in 1:nrow(pedigree)) {
      offspring_id <- pedigree$id[i]
      dam_id <- pedigree$dam[i]
      
      offspring_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x,'_')[[1]][1]) %in% offspring_id]
      offspring_haplotypes_i <- haplotypes[offspring_row,]
      Offspring_Hap1 <- t(as.data.frame(offspring_haplotypes_i[1,]))
      Offspring_Hap2 <- t(as.data.frame(offspring_haplotypes_i[2,]))
      
      dam_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x,'_')[[1]][1]) %in% dam_id]
      dam_haplotypes_i <- haplotypes[dam_row,]
      Maternal_Hap1 <- t(as.data.frame(dam_haplotypes_i[1,]))
      rownames(Maternal_Hap1) <- dam_id
      Maternal_Hap2 <- t(as.data.frame(dam_haplotypes_i[2,]))
      rownames(Maternal_Hap2) <- dam_id
      Maternal_Geno <- matrix(data = Maternal_Hap1 + Maternal_Hap2, nrow = 1)
      colnames(Maternal_Geno) <- colnames(Maternal_Hap1)
      rownames(Maternal_Geno) <- dam_id
      
      method_results <- list()
      length_results <- list()
      pat_results <- list()
      mat_results <- list()
      
      for (j in unique(map[, 2])) {
        print(x = c(i, j))
        Chr_map <- map[map$chr == j, ]
        Chr_markerNames <- Chr_map[, 1]
        
        Offspring_Hap1_chr <- Offspring_Hap1[, colnames(Offspring_Hap1) %in% Chr_markerNames, drop = FALSE]
        Offspring_Hap2_chr <- Offspring_Hap2[, colnames(Offspring_Hap2) %in% Chr_markerNames, drop = FALSE]
        Maternal_Geno_chr <- Maternal_Geno[, colnames(Maternal_Geno) %in% Chr_markerNames, drop = FALSE]
        
        scores_h1_m <- weighing_haplotypes(Offspring_Hap1_chr, Maternal_Geno_chr, method = method, j = j)
        scores_h2_m <- weighing_haplotypes(Offspring_Hap2_chr, Maternal_Geno_chr, method = method, j = j)
        
        # Combine scores for this method
        scores <- rbind(scores_h1_m, scores_h2_m)
        scores$Parent <- c("h1_m", "h2_m")
        
        # Check logic
        score_1 <- scores_h1_m$false_score
        score_2 <- scores_h2_m$false_score
        
        abs_diff <- abs(c(score_1, score_2))
        max_diff <- max(abs_diff)
        
        if (max_diff == abs(score_1) & max_diff != abs(score_2)) {
          message("Hap1 is maternal - assume Hap2 is paternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, maternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, paternal = TRUE, offspring_id = offspring_id)
        }
        if (max_diff == abs(score_2) & max_diff != abs(score_1)) {
          message("Hap2 is maternal - assume Hap1 is paternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, paternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, maternal = TRUE, offspring_id = offspring_id)
        }
        if (max_diff == abs(score_1) & max_diff == abs(score_2)) {
          stop("Both haplotypes are assumed to be maternal")
        }
        
        merged_haps <- rbind(Offspring_Hap1_chr, Offspring_Hap2_chr)
        identified_sire_haplo <- grep("_paternal$", rownames(merged_haps), value = TRUE)
        identified_maternal_haplo <- grep("_maternal$", rownames(merged_haps), value = TRUE)
        
        sire_assigned <- merged_haps[identified_sire_haplo, , drop = FALSE]
        maternal_assigned <- merged_haps[identified_maternal_haplo, , drop = FALSE]
        
        pat_results[[j]] <- as.data.frame(sire_assigned)
        mat_results[[j]] <- as.data.frame(maternal_assigned)
        
        # Combine all scores into one data frame
        testing_methods <- do.call(rbind, scores)
        
        lengths_h1_m <- check_dist(Offspring_Hap1_chr, Maternal_Geno_chr, j = j)
        lengths_h1_m$Parent <- rep("h1_m")
        lengths_h2_m <- check_dist(Offspring_Hap2_chr, Maternal_Geno_chr, j = j)
        lengths_h2_m$Parent <- rep("h2_m")
        
        lengths_all <- rbind(lengths_h1_m, lengths_h2_m)
        
        method_results[[j]] <- as.data.frame(testing_methods)
        length_results[[j]] <- as.data.frame(lengths_all)
      }
      
      real_results_methods[[i]] <- do.call(rbind, method_results)
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
  
  # Second block: Phased haplotypes case
  if (perfect_haplotypes == FALSE) {
    phased_results_methods <- list()
    phased_results_lengths <- list()
    flipped_haplotypes <- list()
    
    if (Data_type == "NoGE_SNP2k"){
      haplotypes <- NoGE_SNP2k_PhasedHaplotypes_matPed
      map <- NoGE_map
    } else if (Data_type == "WithGE_SNP2k"){
      haplotypes <- WithGE_SNP2k_PhasedHaplotypes_matPed
      map <- WithGE_map
    } else if (Data_type == "NoGE_SNP50k"){
      haplotypes <- NoGE_SNP50k_PhasedHaplotypes_matPed
      map <- NoGE_map
    }else if (Data_type == "WithGE_SNP50k"){
      haplotypes <- WithGE_SNP50k_PhasedHaplotypes_matPed
      map <- WithGE_map
    }else if (Data_type == "Real_Slov_data"){
      haplotypes <- Slov_phasedhaplotypes_phasedmatPed
      map <- Slov_map
    } else (stop("Arguments not met"))
    
    for (i in 1:nrow(pedigree)) {
      offspring_id <- pedigree$id[i]
      dam_id <- pedigree$dam[i]
      
      offspring_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x, '_')[[1]][1]) %in% offspring_id]
      offspring_haplotypes_i <- haplotypes[offspring_row, ]
      Offspring_Hap1 <- as.data.frame(offspring_haplotypes_i[1, , drop = FALSE])
      Offspring_Hap2 <- as.data.frame(offspring_haplotypes_i[2, , drop = FALSE])
      
      dam_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x, '_')[[1]][1]) %in% dam_id]
      dam_haplotypes_i <- haplotypes[dam_row, ]
      Maternal_Hap1 <- as.data.frame(dam_haplotypes_i[1, , drop = FALSE])
      Maternal_Hap2 <- as.data.frame(dam_haplotypes_i[2, , drop = FALSE])
      Maternal_Geno <- matrix(data = Maternal_Hap1 + Maternal_Hap2, nrow = 1)
      colnames(Maternal_Geno) <- colnames(Maternal_Hap1)
      rownames(Maternal_Geno) <- dam_id  
      
      
      method_results <- list()
      length_results <- list()
      pat_results <- list()
      mat_results <- list()
      
      for (j in unique(map[, 1])) {
        print(x = c(i, j))
        Chr_map <- map[map[,1] == j, ]
        Chr_markerNames <- Chr_map[, 2]
        
        Offspring_Hap1_chr <- Offspring_Hap1[, colnames(Offspring_Hap1) %in% Chr_markerNames, drop = FALSE]
        Offspring_Hap2_chr <- Offspring_Hap2[, colnames(Offspring_Hap2) %in% Chr_markerNames, drop = FALSE]
        Maternal_Geno_chr <- Maternal_Geno[, colnames(Maternal_Geno) %in% Chr_markerNames, drop = FALSE]
        
        scores_h1_m <- weighing_haplotypes(Offspring_Hap1_chr, Maternal_Geno_chr, method = method, j = j)
        scores_h2_m <- weighing_haplotypes(Offspring_Hap2_chr, Maternal_Geno_chr, method = method, j = j)
        
        scores <- rbind(scores_h1_m, scores_h2_m)
        scores$Parent <- c("h1_m", "h2_m")
        
        score_1 <- scores_h1_m$false_score
        score_2 <- scores_h2_m$false_score
        
        abs_diff <- abs(c(score_1, score_2))
        max_diff <- max(abs_diff)
        
        if (max_diff == abs(score_1) & max_diff != abs(score_2)) {
          message("Hap1 is maternal - assume Hap2 is paternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, maternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, paternal = TRUE, offspring_id = offspring_id)
        }
        if (max_diff == abs(score_2) & max_diff != abs(score_1)) {
          message("Hap2 is maternal - assume Hap1 is paternal")
          Offspring_Hap1_chr <- assign_parent_haplo(Offspring_Hap1_chr, paternal = TRUE, offspring_id = offspring_id)
          Offspring_Hap2_chr <- assign_parent_haplo(Offspring_Hap2_chr, maternal = TRUE, offspring_id = offspring_id)
        }
        if (max_diff == abs(score_1) & max_diff == abs(score_2)) {
          stop("Both haplotypes are assumed to be maternal")
        }
        
        merged_haps <- rbind(Offspring_Hap1_chr, Offspring_Hap2_chr)
        identified_sire_haplo <- grep("_paternal$", rownames(merged_haps), value = TRUE)
        identified_maternal_haplo <- grep("_maternal$", rownames(merged_haps), value = TRUE)
        
        sire_assigned <- merged_haps[identified_sire_haplo, , drop = FALSE]
        maternal_assigned <- merged_haps[identified_maternal_haplo, , drop = FALSE]
        
        pat_results[[j]] <- as.data.frame(sire_assigned)
        mat_results[[j]] <- as.data.frame(maternal_assigned)
        
        testing_methods <- do.call(rbind, scores)
        
        lengths_h1_m <- check_dist(Offspring_Hap1_chr, Maternal_Geno_chr, j = j)
        lengths_h1_m$Parent <- rep("h1_m")
        lengths_h2_m <- check_dist(Offspring_Hap2_chr, Maternal_Geno_chr, j = j)
        lengths_h2_m$Parent <- rep("h2_m")
        
        lengths_all <- rbind(lengths_h1_m, lengths_h2_m)
        
        method_results[[j]] <- as.data.frame(testing_methods)
        length_results[[j]] <- as.data.frame(lengths_all)
      }
      
      phased_results_methods[[i]] <- do.call(rbind, method_results)
      phased_results_lengths[[i]] <- do.call(rbind, length_results)
      
      pat_results <- do.call(cbind, pat_results)
      mat_results <- do.call(cbind, mat_results)
      results <- rbind(mat_results, pat_results)
      flipped_haplotypes[[i]] <- as.data.frame(results)
    }
    
    phased_results_methods_all <- do.call(rbind, phased_results_methods)
    phased_results_lengths_all <- do.call(rbind, phased_results_lengths)
    phased_results_flipped <- do.call(rbind, flipped_haplotypes)
    
    results <- list(phased_results_methods_all = phased_results_methods_all,
                    phased_results_lengths_all = phased_results_lengths_all,
                    phased_results_flipped = phased_results_flipped)
  }
  
  return(results)
}

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


#####################################################################################
#**** Get the true simulated haplotypes and phased haplotypes ******
#####################################################################################

load(paste0(workingDir,"/Data/SP_object.Rdata"))
Simulated_SP <- SP
Simulated_pop <- load(paste0(workingDir,"/Data/Pop_withFathers.Rdata"))
true_haplotypes <- pullSnpHaplo(PopMerged)
true_map <- getGenMap(SP)

#Pedigree prior to reconstruction ####################################################

Worker_pedigree <- read.csv("Data/worker_pedigree.csv")

#Pedigree post reconstruction (with AlphaAssign) ########################################
Rec_pedigree_2k_NoGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_2k_NoGE.txt")
Rec_pedigree_50k_NoGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_50k_NoGE.txt")

Rec_pedigree_2k_WithGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_2k_WithGE.txt")
Rec_pedigree_50k_WithGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_50k_WithGE.txt")

cols <- c("id", "sire", "dam")
names(Rec_pedigree_2k_NoGE) <- cols
names(Rec_pedigree_50k_NoGE) <- cols
names(Rec_pedigree_2k_WithGE) <- cols
names(Rec_pedigree_50k_WithGE) <- cols


#Get haplotypes made in 6_Converting_PhasedVCF scripts
NoGE_SNP2k_PhasedHaplotypes_recPed <- read.table("FILELOCATI|ON/NoGE_SNP2k_PhasedHaplotypes_recPed.txt")
WithGE_SNP2k_PhasedHaplotypes_recPed <- read.table("FILELOCATI|ON/WithGE_SNP2k_PhasedHaplotypes_recPed.txt")
NoGE_SNP50k_PhasedHaplotypes_recPed <- read.table("FILELOCATI|ON/NoGE_SNP50k_PhasedHaplotypes_recPed.txt")
WithGE_SNP50k_PhasedHaplotypes_recPed <- read.table("FILELOCATI|ON/WithGE_SNP50k_PhasedHaplotypes_recPed.txt")

NoGE_SNP2k_PhasedHaplotypes_matPed <- read.table("FILELOCATI|ON/NoGE_SNP2k_PhasedHaplotypes_matPed.txt")
WithGE_SNP2k_PhasedHaplotypes_matPed <- read.table("FILELOCATI|ON/WithGE_SNP2k_PhasedHaplotypes_matPed.txt")
NoGE_SNP50k_PhasedHaplotypes_matPed <- read.table("FILELOCATI|ON/NoGE_SNP50k_PhasedHaplotypes_matPed.txt")
WithGE_SNP50k_PhasedHaplotypes_matPed <- read.table("FILELOCATI|ON/WithGE_SNP50k_PhasedHaplotypes_matPed.txt")

NoGE_map_2k <- read.table("Data/Sim_NoGE/SNP_4_NoGE_QC_ACformat.map")
WithGE_map_2k <- read.table("Data/Sim_WithGE/SNP_4_WithGE_QC_ACformat.map")
NoGE_map_50k <- read.table("Data/Sim_NoGE/SNP_5_NoGE_QC_ACformat.map")
WithGE_map_50k <- read.table("Data/Sim_WithGE/SNP_5_WithGE_QC_ACformat.map")

#####################################################################################
#**** Get the real data phased haplotypes ******
#####################################################################################
setwd(workingDir)
# Pedigree prior to ped reconstruction
Slov_pedigree_pre <- read.table("Data/Real_data/Real_Data_pedigree.txt")
colnames(Slov_pedigree_pre) <- cols

#After reconstruction
Slov_pedigree_post <- read.table("Outputs/AlphaAssign/Alpha_pedigree_Real.txt")
colnames(Slov_pedigree_post) <- cols

Slov_PhasedHaplotypes_recPed <- read.table("FILELOCATION/ Slov_PhasedHaplotypes_recPed.txt")
Slov_PhasedHaplotypes_matPed <- read.table("FILELOCATION/ Slov_PhasedHaplotypes_matPed.txt")

#get the pre-phased data for comparison
Slov_haplotypes_Prephased <- read.table("FILELOCATION/ Slov_haplotypes_Prephased.txt")

Slov_map <- read.table("/Data/Real_data/Slov_fM_AC_QC.map")

#####################################################################################
#*** Check the different phasing versions compared to the real and prephased data just to check similarity ***
#####################################################################################

# ••• Simulated •••
true_vs_NoGEphasedSNP2krecPed <- check_haplotype(true_haplotypes = true_haplotypes, results = NoGE_SNP2k_PhasedHaplotypes_recPed, pedigree = Worker_pedigree)
true_vs_NoGEphasedSNP2kmatPed <- check_haplotype(true_haplotypes = true_haplotypes, results = NoGE_SNP2k_PhasedHaplotypes_matPed, pedigree = Worker_pedigree)

true_vs_WithGEphasedSNP2krecPed <- check_haplotype(true_haplotypes = true_haplotypes, results = WithGE_SNP2k_PhasedHaplotypes_recPed, pedigree = Worker_pedigree)
true_vs_WithGEphasedSNP2kmatPed <- check_haplotype(true_haplotypes = true_haplotypes, results = WithGE_SNP2k_PhasedHaplotypes_matPed, pedigree = Worker_pedigree)

true_vs_NoGEphasedSNP50krecPed <- check_haplotype(true_haplotypes = true_haplotypes, results = NoGE_SNP50k_PhasedHaplotypes_recPed, pedigree = Rec_pedigree_50k_NoGE )
true_vs_NoGEphasedSNP50kmatPed <- check_haplotype(true_haplotypes = true_haplotypes, results = NoGE_SNP50k_PhasedHaplotypes_matPed, pedigree = Worker_pedigree)

true_vs_WithGEphasedSNP50krecPed <- check_haplotype(true_haplotypes = true_haplotypes, results = WithGE_SNP50k_PhasedHaplotypes_recPed, pedigree = Rec_pedigree_50k_WithGE )
true_vs_WithGEphasedSNP50kmatPed <- check_haplotype(true_haplotypes = true_haplotypes, results = WithGE_SNP50k_PhasedHaplotypes_matPed, pedigree = Worker_pedigree)


#####################################################################################
#**               Haplotype parent-of-origin assignments                         **

#**Route 1: use pedigree reconstruction's pedigree and use both dam and sire pedigree ID to assign haplotype parental origins **
#####################################################################################

#True haplotypes (should work perfectly)
Route1_SimTrue <- Route1_flipping(perfect_haplotypes = TRUE, pedigree = Worker_pedigree, method = "power_mean")
dir.create("Data/Haplo_Assignment")
save(Route1_SimTrue, file = "Data/Haplo_Assignment/Route1_SimTrue.Rdata")

#Editing route 1 pedigrees - can't have sire's present 
Rec_pedigree_2k_NoGE_filtered <- Rec_pedigree_2k_NoGE[Rec_pedigree_2k_NoGE$sire != 0,]
Rec_pedigree_50k_NoGE_filtered <- Rec_pedigree_50k_NoGE[Rec_pedigree_50k_NoGE$sire != 0,]
Rec_pedigree_2k_WithGE_filtered <- Rec_pedigree_2k_WithGE[Rec_pedigree_2k_WithGE$sire != 0,]
Rec_pedigree_50k_WithGE_filtered <- Rec_pedigree_50k_WithGE[Rec_pedigree_50k_WithGE$sire != 0,]


#2k SNP
Route1_NoGE_SNP2k <- Route1_flipping(perfect_haplotypes = FALSE, pedigree = Rec_pedigree_2k_NoGE_filtered, method = "power_mean", Data_type = "NoGE_SNP2k")
Route1_WithGE_SNP2k <- Route1_flipping(perfect_haplotypes = FALSE, pedigree = Rec_pedigree_2k_WithGE_filtered, method = "power_mean", Data_type = "WithGE_SNP2k")
save(Route1_NoGE_SNP2k, file ="/Data/Haplo_Assignment/Route1_NoGE_SNP2k.Rdata")
save(Route1_WithGE_SNP2k, file = "/Data/Haplo_Assignment/Route1_WithGE_SNP2k.Rdata")

#50k SNP
Route1_NoGE_SNP50k <- Route1_flipping(perfect_haplotypes = FALSE, pedigree = Rec_pedigree_50k_NoGE_filtered, method = "power_mean", Data_type = "NoGE_SNP50k")
Route1_WithGE_SNP50k <- Route1_flipping(perfect_haplotypes = FALSE, pedigree = Rec_pedigree_50k_WithGE_filtered, method = "power_mean", Data_type = "WithGE_SNP50k")
save(Route1_NoGE_SNP50k, file = "/Data/Haplo_Assignment/Route1_NoGE_SNP50k.Rdata")
save(Route1_WithGE_SNP50k, file = "/Data/Haplo_Assignment/Route1_WithGE_SNP50k.Rdata")

#Real data
Route1_Real <- Route1_flipping(perfect_haplotypes = FALSE, pedigree = Slov_pedigree_post, methods = "power_mean", Data_type = "Real_Slov_data")
write.csv("/Data/Haplo_Assignment/Route1_Real.csv")


# #Checking haplotypes post flip to see if something has happened
# colnames(Route1_SimTrue$real_results_flipped) <- colnames(true_haplotypes)
# 
#Example:
# true_vs_NoGE_SNP2k_FLIPPED <- check_haplotype_postFlip(complete_haplotypes = Route1_SimTrue$real_results_flipped, results = Route1_NoGE_SNP2k$results_flipped, pedigree = Rec_pedigree_2k_NoGE)
# 

#####################################################################################
#** Route 2: DO NOT use pedigree reconstruction and use ONLY dam pedigree ID to assign haplotype parental origins **
#####################################################################################

#True haplotypes (should work perfectly)
Route2_SimTrue <- Route2_flipping(perfect_haplotypes = TRUE, pedigree = Mat_pedigree, method = "power_mean")
write.csv("/Data/Haplo_Assignment/Route2_SimTrue.csv")
#2k SNP
Route2_NoGE_SNP2k <- Route2_flipping(perfect_haplotypes = FALSE, pedigree = Rec_pedigree_2k_NoGE, method = "power_mean",Data_type = "NoGE_SNP2k")
Route2_WithGE_SNP2k <- Route2_flipping(perfect_haplotypes = FALSE, pedigree = Rec_pedigree_2k_WithGE, method = "power_mean", Data_type = "WithGE_SNP2k")
write.csv("/Data/Haplo_Assignment/Route2_NoGE_SNP2k.csv")
write.csv("/Data/Haplo_Assignment/Route2_WithGE_SNP2k.csv")
#50k SNP
Route2_NoGE_SNP50k <- Route2_flipping(perfect_haplotypes = FALSE, pedigree = Rec_pedigree_50k_NoGE, method = "power_mean", Data_type = "NoGE_SNP50k")
Route2_WithGE_SNP50k <- Route2_flipping(perfect_haplotypes = FALSE, pedigree = Rec_pedigree_50k_WithGE, method = "power_mean", Data_type = "WithGE_SNP2k")
write.csv("/Data/Haplo_Assignment/Route2_NoGE_SNP50k.csv")
write.csv("/Data/Haplo_Assignment/Route2_WithGE_SNP50k.csv")

#Real data
Route2_Real <- Route2_flipping(perfect_haplotypes = FALSE, pedigree = Slov_pedigree_post, methods = "power_mean", Data_type = "Real_Slov_data")
write.csv("/Data/Haplo_Assignment/Route2_Real.csv")






