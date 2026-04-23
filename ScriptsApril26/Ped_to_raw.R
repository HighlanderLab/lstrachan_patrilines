ped_to_raw <- function(ped_file, output_file) {
  # Read the .ped file
  ped_data <- read.table(ped_file, header = FALSE, stringsAsFactors = FALSE)
  
  # Extract the first 6 columns (family ID, individual ID, paternal ID, maternal ID, sex, phenotype)
  ped_info <- ped_data[, 1:6] 
  ped_snps <- ped_data[, -c(1:6)] # Extract the SNP genotype data


  # Create new dataframe with paired pasted columns
  ped_snps_raw <- data.frame(
    lapply(seq(1, ncol(ped_snps), by = 2), function(i) {
      paste(ped_snps[[i]], ped_snps[[i+1]], sep = "")
    })
  )
  colnames(ped_snps_raw) <- NULL

  map <- c("CC" = 2, "AA" = 0, "AC" = 1, "CA" = 1, "00" = 9)

  ped_snps_raw[] <- lapply(ped_snps_raw, function(col) {
    map[col]
  })

  ped_snps_raw_plink <- cbind(ped_info, ped_snps_raw)
  write.table(ped_snps_raw_plink, file = output_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
}
