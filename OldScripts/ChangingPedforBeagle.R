#CHANGING PED. SO THAT THE VCF FILE WILL WORK IN BEAGLE 


# Define file paths
setwd("~/Desktop/Slovenia data/Simulated /nested/160SNP/data")

# Function to convert allele codes
convert_alleles <- function(allele) {
  if (allele == "A") {
    return("A")
  } else if (allele == "B") {
    return("C")
  } else {
    return(allele)
  }
}

# Read the original PED file
ped_data <- read.table("newPat.ped", header = FALSE, stringsAsFactors = FALSE)

# Convert alleles
corrected_genotype_data <- apply(ped_data[, 7:ncol(ped_data)], 2, function(allele) sapply(allele, convert_alleles))

# Combine the non-genotype and corrected genotype data
corrected_ped_data <- cbind(ped_data[, 1:6], corrected_genotype_data)

# Write the corrected PED file
write.table(corrected_ped_data, file = "corrected_pedSLOV.ped", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)




setwd("~/Desktop/Slovenia data/Beagle ")
# Read the original .map file
map_data <- read.table("corrected_ped_file16.map", header = FALSE, stringsAsFactors = FALSE)

# Check for duplicated genetic positions
if(any(duplicated(map_data$V3))) {
  cat("Duplicated genetic positions found. Creating unique positions...\n")
  
  # Assign unique genetic positions
  unique_positions <- seq(0, by = 0.01, length.out = nrow(map_data))
  map_data$V3 <- unique_positions
}

# Write the corrected map file
write.table(map_data, file = corrected_map_file_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

cat("Corrected MAP file has been saved to", corrected_map_file_path, "\n")

