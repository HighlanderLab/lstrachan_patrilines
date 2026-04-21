###### Phasing prep #########

rm(list = ls())

library(Eagle)
library(tidyr)
setwd("~/Desktop/Slovenia data/Attempt2/Slov/Not-phased/AfterQC-plink-files/AC format ")
#Get the raw files Jana gave for the Slovenian data
Slov_ped <- read.table("Slov_QC_AC.ped")
Slov_map <- read.table("Slov_QC_AC.map")

setwd("~/Desktop/Slovenia data/Attempt2/Slov/Not-phased/Raw data from Jana")
Slov_raw_ped <- read.table("newPat.ped")
Slov_raw_map <- read.table("newPat.map")
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


#Woo thats worked!
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


original_chromosomes <- unique(Slov_map_filtered$Chromosome)
new_chromosomes <- 1:16
chromosome_map <- setNames(new_chromosomes, original_chromosomes)

# Replace the chromosome names in the dataframe with the new numbers
Slov_map_filtered$Chromosome <- chromosome_map[Slov_map_filtered$Chromosome]

#These filtered files SHOULD now be alright to save and put into PLINK for QC and phasing 
#Add in the genetic distance column with 0 for missing 
Slov_map_filtered$GenDis <- rep(0)
Slov_map_filtered <- Slov_map_filtered[,c(1,3,4,2)] #put it in the right order for plink 


#These are duplicated positions but all have the same call rate so are true duplicates. Just get rid of the first one
# CHR	POS	ALLELES	IDS
# 1	10570868	A,C	WP3DKS_NC_007070.3_10570868_A WP3ES_NC_007070.3_10570868_A
# 1	17742947	A,C	WP3DKS_NC_007070.3_17742947_A WP3ES_NC_007070.3_17742947_A
# 5	7287901	A,C	WP3DKS_NC_007074.3_7287901_A WP3ES_NC_007074.3_7287901_A
# 5	7368521	A,C	WP3DKS_NC_007074.3_7368521_A WP3ES_NC_007074.3_7368521_A
# 5	11335903	A,C	WP3DKS_NC_007074.3_11335903_A WP3ES_NC_007074.3_11335903_A
# 8	679379	A,C	WP3DKS_NC_007077.3_679379_A WP3ES_NC_007077.3_679379_A
# 8	11497122	A,C	WP3DKS_NC_007077.3_11497122_B WP3ES_NC_007077.3_11497122_A
# 8	11502065	A,C	WP3DKL_NC_007077.3_11502065_A WP3DKS_NC_007077.3_11502065_A
# 8	11502167	A,C	WP3DKL_NC_007077.3_11502167_A WP3DKS_NC_007077.3_11502167_A
# 10	12473010	A,C	WP3DKL_NC_007079.3_12473010_A WP3DKS_NC_007079.3_12473010_A
# 11	14304532	A,C	WP3DKL_NC_007080.3_14304532_A WP3DKS_NC_007080.3_14304532_A
# 11	14305861	A,C	WP3DKL_NC_007080.3_14305861_A WP3DKS_NC_007080.3_14305861_A
# 12	3272406	A,C	WP3DKL_NC_007081.3_3272406_A WP3DKS_NC_007081.3_3272406_A
# 12	5162719	A,C	WP3DKS_NC_007081.3_5162719_A WP3ES_NC_007081.3_5162719_A
# 13	8293069	A,C	WP3DKS_NC_007082.3_8293069_A WP3ES_NC_007082.3_8293069_A
# 14	4662725	A,C	WP3DKS_NC_007083.3_4662725_A WP3ES_NC_007083.3_4662725_A

#get the ids using 
#./plink --file Slov_filtered_QC --list-duplicate-vars ids-only suppress-first --out Slov_filtered_QC
Dup_Ids <- read.table("Slov_filtered_QC.dupvar")
Dup_Ids <- Dup_Ids$V1

Slov_map_filtered_noDup <- Slov_map_filtered[!Slov_map_filtered$markerID %in% Dup_Ids$V1, ]

#Woops we need to change the Chromsome names to numbers! 

Slov_map_filtered_chrom <- Slov_map_filtered_noDup %>%
  mutate(Chromosome = ifelse(Chromosome == "NC_007070.3", 1, 
                      ifelse(Chromosome == "NC_007071.3", 2,
                      ifelse(Chromosome == "NC_007072.3", 3,
                      ifelse(Chromosome == "NC_007073.3", 4,
                      ifelse(Chromosome == "NC_007074.3", 5,
                      ifelse(Chromosome == "NC_007075.3", 6,
                      ifelse(Chromosome == "NC_007076.3", 7,
                      ifelse(Chromosome == "NC_007077.3", 8,
                      ifelse(Chromosome == "NC_007078.3", 9,
                      ifelse(Chromosome == "NC_007079.3", 10,
                      ifelse(Chromosome == "NC_007080.3", 11,
                      ifelse(Chromosome == "NC_007081.3", 12,
                      ifelse(Chromosome == "NC_007082.3", 13,
                      ifelse(Chromosome == "NC_007083.3", 14,
                      ifelse(Chromosome == "NC_007084.3", 15,
                      ifelse(Chromosome == "NC_007085.3", 16, Chromosome)))))))))))))))))

unique(Slov_map_filtered_chrom$Chromosome)


write.table(Slov_map_filtered_chrom, file = "Slov_filtered_QC.map", quote = F, col.names = F, row.names = F, sep = " ")
write.table(Slov_ped_filtered, file = "Slov_filtered_QC.ped", quote = F, col.names = F, row.names = F, sep = " ")

Slov_ped_edit <- read.table("Slov_filtered_QC.ped")
Slov_ped_edit[,c(1:6)] <- Slov_ped[,c(1:6)]
write.table(Slov_ped_edit, file = "Slov_filtered_QC_test.ped", quote = F, col.names = F, row.names = F, sep = " ")

Slov_map_filtered_chrom_sorted <- Slov_map_filtered_chrom[order(Slov_map_filtered_chrom$Position),]
Slov_map_filtered_chrom_sorted$GenDis <- rep("0")
write.table(Slov_map_filtered_chrom_sorted, file = "Slov_filtered_ordered_QC.map", quote = F, col.names = F, row.names = F, sep = " ")


#Lets check the percentages of A,B,0

calculate_genotype_percentages <- function(df) {
  # Select only the genotype columns (columns after the first 6)
  genotype_df <- df[, -c(1:6)]
  
  # Function to calculate percentages of A, B, and 0 for a single row
  calculate_row_percentages <- function(row) {
    total <- length(row)
    count_A <- sum(row == "A")
    count_B <- sum(row == "B")
    count_0 <- sum(row == "0")
    
    percent_A <- (count_A / total) * 100
    percent_B <- (count_B / total) * 100
    percent_0 <- (count_0 / total) * 100
    
    return(c(percent_A, percent_B, percent_0))
  }
  
  # Apply the function to each row and store the results in a dataframe
  percentages <- t(apply(genotype_df, 1, calculate_row_percentages))
  colnames(percentages) <- c("percent_A", "percent_B", "percent_0")
  
  # Combine the original dataframe with the calculated percentages
  result_df <- cbind(df[, 1:6], percentages)
  
  return(result_df)
}

percentages_df_filter <- calculate_genotype_percentages(Slov_ped_filtered)

percetage_df_prefiltered <- calculate_genotype_percentages(Slov_ped)

#Now to make the pedigree file for Beagle 4 to phased more accurately 
# Family ID: A unique identifier for the family.
# Individual ID: A unique identifier for the individual within the family.
# Paternal ID: The Individual ID of the father (0 if father is unknown).
# Maternal ID: The Individual ID of the mother (0 if mother is unknown).
# Sex: Sex of the individual coded as 1 = male, 2 = female, and 0 = unknown.
# Phenotype: The phenotype for the individual (commonly used for affection status, with 1 = unaffected, 2 = affected, and 0/-9 = unknown).

Phased_pedigree_file <- Slov_ped[,c(1:6)]
colnames(Phased_pedigree_file) <- c("FamilyID", "IndID", "PaternalID","MaternalID", "Sex", "Phenotype")
write.table(Phased_pedigree_file, file = "Phased_pedigree_file.txt", col.names = F, row.names = F, quote = F, sep = " ")


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




