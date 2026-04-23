AB_to_12 <- function(ped_file, map_ref_file, output_file_prefix) {  
  ped_data <- read.table(ped_file, header = FALSE, stringsAsFactors = FALSE)
  map_ref_data <- read.table(map_ref_file, header = TRUE, stringsAsFactors = FALSE)
  map_ref_data[map_ref_data$RefAllele == "C", "RefAllele"] <- "A"
  map_ref_data[map_ref_data$RefAllele == "D", "RefAllele"] <- "B"
  
  # Extract the first 6 columns (family ID, individual ID, paternal ID, maternal ID, sex, phenotype)
  ped_info <- ped_data[, 1:6] 
  colnames(ped_info) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
  ped_snps <- ped_data[, -c(1:6)] # Extract the SNP genotype data
  ped_snps_AB = ped_snps
  
  # Write a function that loops through columns of ped_snps and if the value in the RefAllele column is A, then all the As are turned into 1s and all the Bs are turned into 2s
  #If the RefAllele is B, then all the As are turned into 2s and all the Bs are turned into 1s
  # The value of 0 stays 0
  for (row in seq(from = 1, to = ncol(ped_snps), by=2)) {
    ref_allele = map_ref_data$RefAllele[(row+1)/2]
    if (ref_allele == "A") {
      ped_snps[, row][ped_snps[, row] == "A"] <- "1"
      ped_snps[, row][ped_snps[, row] == "B"] <- "2"
      ped_snps[, row+1][ped_snps[, row+1] == "A"] <- "1"
      ped_snps[, row+1][ped_snps[, row+1] == "B"] <- "2"
    } else if (ref_allele == "B") {
      ped_snps[, row][ped_snps[, row] == "A"] <- "2"
      ped_snps[, row][ped_snps[, row] == "B"] <- "1"
      ped_snps[, row+1][ped_snps[, row+1] == "A"] <- "2"
      ped_snps[, row+1][ped_snps[, row+1] == "B"] <- "1"
    } 
  }
  colnames(ped_snps) <- paste0(map_ref_data$Chromosome, "_", map_ref_data$Position)

  ped_snps_plink <- cbind(ped_info, ped_snps)  
  new_map = data.frame(Chromosome = map_ref_data$Chromosome, ID = paste0(map_ref_data$Chromosome, "_", map_ref_data$Position), Position1 = 0, Position2 = map_ref_data$Position)
  
  write.table(ped_snps_plink, file = paste0(output_file_prefix, ".ped"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(new_map, file = paste0(output_file_prefix, ".map"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")
}
