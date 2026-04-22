###### Phasing prep #########

rm(list = ls())

library(Eagle)
library(tidyr)
library(dplyr)

pathToPlink <- "/home/jana/bin/"
workingDir = "/home/jana/github/lstrachan_patrilines/"
setwd(workingDir)

Slov_raw_ped <- read.table("Data/Real_data/newPat.ped")
Slov_raw_map <- read.table("Data/Real_data/newPat.map")
#Slov map is in the wrong format so we need to fix it first
Fix_map <- Slov_raw_map[,2]

# Function to extract the desired parts and create a data frame
extract_parts <- function(x) {
  parts <- unlist(strsplit(x, "[._]"))
  chromosome <- paste(parts[2], paste(parts[3], parts[4], sep="."), sep="_")
  position <- as.numeric(parts[5])
  data.frame(Chromosome = chromosome, Position = position, stringsAsFactors = FALSE)
}

# Organise the map (we don't have the Genetic distance colummn yet and will add that at the end)
Slov_map_organised <- do.call(rbind, lapply(Fix_map, extract_parts))
Slov_map_organised$markerID <- Fix_map

#we need to add the column names to the pedigree so we get rid of the right columns 
transformed_list <- list()

# Loop through each row in the dataframe
for (i in 1:nrow(Slov_map_organised)) {
  # Get the current row
  current_row <- Slov_map_organised[i, ]
  
  # Create two copies of the current row with modified Chromosome values
  row_copy_1 <- current_row
  row_copy_2 <- current_row
  row_copy_1$markerID <- paste(current_row$markerID, "1", sep = "_")
  row_copy_2$markerID <- paste(current_row$markerID, "2", sep = "_")
  
  # Add the modified rows to the list
  transformed_list[[length(transformed_list) + 1]] <- row_copy_1
  transformed_list[[length(transformed_list) + 1]] <- row_copy_2
}

# Convert the list back to a dataframe
Long_map_file <- do.call(rbind, transformed_list)

Chr_IDs_forPed <- Long_map_file$markerID
#Make these the colnames for the ped file 
colnames(Slov_raw_ped)[7:ncol(Slov_raw_ped)] <- Chr_IDs_forPed


# Next we only want chromosomes beginning with NC (NW are those with unknown chromomsomes)
filter_chromosome_nc <- function(df) {
  # Filter the dataframe to keep only rows where the Chromosome column starts with "NC"
  filtered_df <- df[grep("^NC", df$Chromosome), ]
  return(filtered_df)
}

Slov_map_filtered <- filter_chromosome_nc(Slov_map_organised)
nrow(Slov_map_filtered)

filtered_list <- list()

# Loop through each row in the dataframe
for (i in 1:nrow(Slov_map_filtered)) {
  # Get the current row
  current_row <- Slov_map_filtered[i, ]
  
  # Create two copies of the current row with modified Chromosome values
  row_copy_1 <- current_row
  row_copy_2 <- current_row
  row_copy_1$markerID <- paste(current_row$markerID, "1", sep = "_")
  row_copy_2$markerID <- paste(current_row$markerID, "2", sep = "_")
  
  # Add the modified rows to the list
  filtered_list[[length(filtered_list) + 1]] <- row_copy_1
  filtered_list[[length(filtered_list) + 1]] <- row_copy_2
}

# Convert the list back to a dataframe
Long_file <- do.call(rbind, filtered_list)
Chr_IDs <- Long_file$markerID

#Do the same for the ped file with the column names 
filter_columns_nc <- function(df) {
  
  # Filter the dataframe to keep only those columns
  filtered_df <- df[, Chr_IDs]
  
  return(filtered_df)
}

Slov_ped_filtered <- filter_columns_nc(Slov_raw_ped)
ncol(Slov_ped_filtered)
#We just got rid of the start of the ped file so lets pop that back in
Slov_ped_filtered <- cbind(Slov_raw_ped[,c(1:6)], Slov_ped_filtered)
ncol(Slov_ped_filtered) #just check the six were added ok

#These filtered files SHOULD now be alright to save and put into PLINK for QC and phasing 
#Add in the genetic distance column with 0 for missing 
Slov_map_filtered$GenDis <- rep(0)
Slov_map_filtered <- Slov_map_filtered[,c(1,3,4,2)] #put it in the right order for plink 


# ---  Map Chromosomes numerically ---

Slov_map_filtered_chrom <- Slov_map_filtered %>%
  dplyr::mutate(Chromosome = case_when(
    Chromosome == "NC_007070.3" ~ 1,
    Chromosome == "NC_007071.3" ~ 2,
    Chromosome == "NC_007072.3" ~ 3,
    Chromosome == "NC_007073.3" ~ 4,
    Chromosome == "NC_007074.3" ~ 5,
    Chromosome == "NC_007075.3" ~ 6,
    Chromosome == "NC_007076.3" ~ 7,
    Chromosome == "NC_007077.3" ~ 8,
    Chromosome == "NC_007078.3" ~ 9,
    Chromosome == "NC_007079.3" ~ 10,
    Chromosome == "NC_007080.3" ~ 11,
    Chromosome == "NC_007081.3" ~ 12,
    Chromosome == "NC_007082.3" ~ 13,
    Chromosome == "NC_007083.3" ~ 14,
    Chromosome == "NC_007084.3" ~ 15,
    Chromosome == "NC_007085.3" ~ 16,
    TRUE ~ 0
  ))


Slov_map_filtered_chrom <- Slov_map_filtered_chrom[, c("Chromosome", "GenDis", "markerID", "Position")]

# --- Save Final Map File --- (must be the same name as the ped file to work in PLINK)
write.table(Slov_map_filtered_chrom, file = "Data/Real_data/Slov_fM_AC.map", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")



## ped file is in 0/A/B format which can't be read by Beagle. Needs to be changed to 0/A/C

# Function to convert allele codes
convert_ped_AB_to_AC <- function(ped_df) {
  ped_df[, 7:ncol(ped_df)] <- apply(
    ped_df[, 7:ncol(ped_df)],
    2,
    function(col) {
      col[col == "B"] <- "C"
      col
    }
  )
  return(ped_df)
}

Slov_ped_filtered[2:10, 7:20]
Slov_ped_filtered_AC <- convert_ped_AB_to_AC(Slov_ped_filtered)
Slov_ped_filtered_AC[2:10, 7:20]


write.table(Slov_ped_filtered_AC, file = "Data/Real_data/Slov_fM_AC.ped", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")


#### ---- PLINK quality control ###################################################################################################
setwd("Data/Real_data")
{# Step 1: Quality control with PLINK
  system(paste0(pathToPlink, "plink --file Slov_fM_AC --make-bed --geno 0.1 --mind 0.1 --maf 0.01 --out Slov_fM_AC_QC"))

  # Step 2: Convert to VCF
  system(paste0(pathToPlink, "plink --bfile Slov_fM_AC_QC --recode vcf --out Slov_fM_AC_QC"))
  
  # Step 3: Sort and index the VCF
  system("bcftools sort Slov_fM_AC_QC.vcf -Oz -o Slov_fM_AC_QC_sorted.vcf.gz")
  system("tabix -p vcf Slov_fM_AC_QC_sorted.vcf.gz")
  
  # Step 3.5: Check for duplicate positions
  system("bcftools query -f '%CHROM\\t%POS\\n' Slov_fM_AC_QC_sorted.vcf.gz | sort | uniq -d")
  
  # Step 3.6: Normalize and remove duplicate records
  system("bcftools norm -m -any Slov_fM_AC_QC_sorted.vcf.gz -Oz -o Slov_fM_AC_QC_biallelic.vcf.gz")
  system("tabix -p vcf Slov_fM_AC_QC_biallelic.vcf.gz")
  system("bcftools norm -d all Slov_fM_AC_QC_biallelic.vcf.gz -Oz -o Slov_fM_AC_QC_noDupPos.vcf.gz")
  system("tabix -p vcf Slov_fM_AC_QC_noDupPos.vcf.gz")
  #TODO: Add the manual transformtion from VCF to Ped and Map
  
}  





#******* Don't know what this is yet - May make sense during phasing *************
#Haploid genome data not working. Try formatting haplotypes as 0|/ 1| 

# Define the input and output file paths
vcf_file <- "Genome_SloDrones.vcf"
output_file <- "Genome_SloDrones_haploid_fixed.vcf"

# Open the input and output files
infile <- file(vcf_file, "r")
outfile <- file(output_file, "w")

# Process the file line by line
while(TRUE) {
  line <- readLines(infile, n = 1)
  
  # Stop if we reach the end of the file
  if(length(line) == 0) {
    break
  }
  
  # Write header lines directly to the output file
  if(grepl("^#", line)) {
    writeLines(line, outfile)
  } else {
    fields <- strsplit(line, "\t")[[1]]
    
    # Check if the format field is GT (Genotype)
    if(fields[9] == "GT") {
      genotypes <- fields[10:length(fields)]
      
      # Add a separator to each haploid genotype
      fixed_genotypes <- sapply(genotypes, function(g) {
        if(g %in% c("0", "1")) {
          return(paste0(g, "|"))
        } else {
          return(g)
        }
      })
      
      # Write the modified line to the output file
      writeLines(paste(c(fields[1:9], fixed_genotypes), collapse = "\t"), outfile)
    } else {
      writeLines(line, outfile)
    }
  }
}

# Close the input and output files
close(infile)
close(outfile)

cat("Fixed VCF file saved to", output_file, "\n")

phased_Slov_vcf <- read.vcfR("phased_Slov_filtered_Beagle4.vcf")

phased_Slov_haplo <- extract.gt(phased_Slov_vcf)



