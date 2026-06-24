Process_AlphaAssign_output <- function(GE = NULL, True_pedigree = NULL, nOffspring = NULL) {
  
  if (isTRUE(GE)) {

    test <- "WithGE"
    filetype <- "_WithGE"
  } else if (isFALSE(GE)) {

    test <- "NoGE"
    filetype <- "_NoGE"
  } else {
    stop("GE must be TRUE or FALSE")
  }
  
  results_list <- list()
  
  for (n in 1:5){
    
    # FIX: proper file name construction
    file_name <- paste0("AlphaGeno_SNP", n, filetype, ".sires")
    
    Alpha_output <- read.table(file_name, header = TRUE)
    
    Sires_assigned <- Alpha_output[Alpha_output$chosen == 1, ]
    nSires_assigned <- nrow(Sires_assigned)
    
    Offspring_and_candidateParent <- Sires_assigned[, c(1,2)]
    colnames(Offspring_and_candidateParent) <- c("id", "sire")
    
    compare_sires <- merge(
      Offspring_and_candidateParent,
      True_pedigree,
      by = "id",
      suffixes = c("_assigned", "_known")
    )
    
    nCorrect_sires <- sum(
      compare_sires$sire_assigned == compare_sires$sire_known,
      na.rm = TRUE
    )
    
    df <- data.frame(
      Test = test,
      nOffspring = length(unique(Alpha_output$id)),
      nDpc = length(unique(Alpha_output$candidate)),
      SNP_group = n,
      nSires_assigned = nSires_assigned,
      nCorrect_sires = nCorrect_sires,
      Software = "AlphaAssign"
    )
    
    results_list[[n]] <- df
  }
  
  # bind all iterations into one dataframe
  final_df <- do.call(rbind, results_list)
  
  return(final_df)
}

ped_to_raw <- function(ped_file, map_file, output_file) {
  # Read the .ped file
  ped_data <- read.table(ped_file, header = FALSE, stringsAsFactors = FALSE)
  map_data <- read.table(map_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Extract the first 6 columns (family ID, individual ID, paternal ID, maternal ID, sex, phenotype)
  ped_info <- ped_data[, 1:6] 
  colnames(ped_info) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
  ped_snps <- ped_data[, -c(1:6)] # Extract the SNP genotype data
  



  # Create new dataframe with paired pasted columns
  ped_snps_raw <- data.frame(
    lapply(seq(1, ncol(ped_snps), by = 2), function(i) {
      paste(ped_snps[[i]], ped_snps[[i+1]], sep = "")
    })
  )
  colnames(ped_snps_raw) <- map_data$V2

  map <- c("11" = 0, "22" = 2, "12" = 1, "21" = 1, "00" = 9)

  ped_snps_raw[] <- lapply(ped_snps_raw, function(col) {
    map[col]
  })

  ped_snps_raw_plink <- cbind(ped_info, ped_snps_raw)
  write.table(ped_snps_raw_plink, file = output_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
}