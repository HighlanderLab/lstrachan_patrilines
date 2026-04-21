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
library(VariantAnnotation)
library(vcfR)
library(tibble)
}

# --- Clear Workspace ---
rm(list = ls())


############################################################################
########   STORE GENERIC FUNCTIONS
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



#######################################################################################################################
#########   Step 1: Slovenian Real Data Setup          ########################################################
#######################################################################################################################

## ----- Sort out the map file of the raw data ----
# Currently map looks like this: 0 WP3DKL_NC_007070.3_10930297_A 0 0
# We need it to look like this:  1 WP3ES_NC_007070.3_159948_A 0 159948

setwd("~/Desktop/Slovenia data/Attempt2/Slov/Not-phased/Raw data from Jana")
Slov_raw_map <- read.table("newPat.map", stringsAsFactors = FALSE)

# -- Take out the column with all the info in it --
Map_info <- Slov_raw_map[, 2]

# --- Function to Extract Chromosome and Position ---
extract_parts <- function(x) {
  parts <- unlist(strsplit(x, "[._]"))
  chromosome <- paste(parts[2], paste(parts[3], parts[4], sep = "."), sep = "_")
  position <- as.numeric(parts[5])
  data.frame(Chromosome = chromosome, Position = position, stringsAsFactors = FALSE)
}

# --- Organise and Format Map File ---
Slov_map_organised <- do.call(rbind, lapply(Map_info, extract_parts))
Slov_map_organised$markerID <- Map_info
View(Slov_map_organised)

# --- Filter for Chromosomes Starting with "NC" ---
Slov_map_filtered <- Slov_map_organised[grep("^NC", Slov_map_organised$Chromosome), ]

View(Slov_map_filtered)

# ---  Map Chromosomes numerically ---

Slov_map_filtered_chrom <- Slov_map_filtered %>%
  mutate(Chromosome = case_when(
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

View(Slov_map_filtered_chrom)

# --- Add Genetic Distance Column and Format ---
Slov_map_filtered_chrom$GenDis <- 0
Slov_map_filtered_chrom <- Slov_map_filtered_chrom[, c("Chromosome", "GenDis", "markerID", "Position")]


# --- Save Final Map File --- (must be the same name as the ped file to work in PLINK)
setwd("~/Desktop/Slovenia data/TestMay25")
write.table(Slov_map_filtered_chrom, file = "Slov_fM_AC.map", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")


## ---- Now the map file is formatted lets look at the ped file  ##################################################################
# -- ped file is in 0/A/B format which can't be read by Beagle. Needs to be changed to 0/A/C
setwd("~/Desktop/Slovenia data/Attempt2/Slov/Not-phased/Raw data from Jana")
Slov_raw_ped <- read.table("newPat.ped", stringsAsFactors = FALSE)

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
Slov_ped_AC <- convert_ped_AB_to_AC(Slov_raw_ped)

#--- lets save the final ped file (must be the same name as the map file to work in PLINK)
setwd("~/Desktop/Slovenia data/TestMay25")
write.table(Slov_ped_filtered, file = "Slov_fM_AC.ped", quote = F, col.names = F, row.names = F, sep = " ")


#### ---- PLINK quality control ###################################################################################################
{# Step 1: Quality control with PLINK
system("./plink --file Slov_fM_AC --make-bed --geno 0.1 --mind 0.1 --maf 0.01 --out Slov_fM_AC_QC")

# Step 2: Remove non-1:16 chromosomes
system("./plink --bfile Slov_fM_AC_QC --chr 1-16 --make-bed --out Slov_fM_AC_QC_noChr0")

# Step 3: Convert to VCF
system("./plink --bfile Slov_fM_AC_QC_noChr0 --recode vcf --out Slov_fM_AC_QC_noChr0")

# Step 4: Sort and index the VCF
system("bcftools sort Slov_fM_AC_QC_noChr0.vcf -Oz -o Slov_fM_AC_QC_noChr0_sorted.vcf.gz")
system("tabix -p vcf Slov_fM_AC_QC_noChr0_sorted.vcf.gz")

# Step 4.5: Check for duplicate positions
system("bcftools query -f '%CHROM\\t%POS\\n' Slov_fM_AC_QC_noChr0_sorted.vcf.gz | sort | uniq -d")

# Step 4.6: Normalize and remove duplicate records
system("bcftools norm -m -any Slov_fM_AC_QC_noChr0_sorted.vcf.gz -Oz -o Slov_fM_AC_QC_biallelic.vcf.gz")
system("tabix -p vcf Slov_fM_AC_QC_biallelic.vcf.gz")
system("bcftools norm -d all Slov_fM_AC_QC_biallelic.vcf.gz -Oz -o Slov_fM_AC_QC_noDupPos.vcf.gz")
system("tabix -p vcf Slov_fM_AC_QC_noDupPos.vcf.gz")

#Alternatively you could just make a bashscript to do all of this and run the script. 
#This script can be used for the real data and the simulated data to speed things up
{
  bash_script <- "
#!/bin/bash

# Step 1: Quality control
./plink --file Slov_fM_AC --make-bed --geno 0.1 --mind 0.1 --maf 0.01 --out Slov_fM_AC_QC

# Step 2: Remove non-1:16 chromosomes
./plink --bfile Slov_fM_AC_QC --chr 1-16 --make-bed --out Slov_fM_AC_QC_noChr0

# Step 3: Recode to VCF
./plink --bfile Slov_fM_AC_QC_noChr0 --recode vcf --out Slov_fM_AC_QC_noChr0

# Step 4: Sort and index
bcftools sort Slov_fM_AC_QC_noChr0.vcf -Oz -o Slov_fM_AC_QC_noChr0_sorted.vcf.gz
tabix -p vcf Slov_fM_AC_QC_noChr0_sorted.vcf.gz

# Step 4.5: Check duplicate positions
bcftools query -f '%CHROM\\t%POS\\n' Slov_fM_AC_QC_noChr0_sorted.vcf.gz | sort | uniq -d

# Step 4.6: Normalize and remove duplicates
bcftools norm -m -any Slov_fM_AC_QC_noChr0_sorted.vcf.gz -Oz -o Slov_fM_AC_QC_biallelic.vcf.gz
tabix -p vcf Slov_fM_AC_QC_biallelic.vcf.gz
bcftools norm -d all Slov_fM_AC_QC_biallelic.vcf.gz -Oz -o Slov_fM_AC_QC_noDupPos.vcf.gz
tabix -p vcf Slov_fM_AC_QC_noDupPos.vcf.gz

# Step 5: Beagle phasing for each chromosome
for chr in {1..16}; do
  java -jar beagle.4.0.jar gt=Slov_fM_AC_QC_noDupPos.vcf.gz out=Slov_PHASED_chr$chr window=200 overlap=50 phase-its=200 impute-its=1000 burnin-its=200 chrom=$chr
done

# Step 6: Convert phased VCFs to PLINK format
for f in Slov_PHASED_chr*.vcf.gz; do
  base=\"${f%.vcf.gz}\"
  ./plink --vcf \"$f\" --recode --double-id --allow-extra-chr --out \"$base\"
done

# Step 7: Compress and index any leftover .vcf files
for f in *.vcf; do
  if [[ ! -f \"$f.gz\" ]]; then
    bgzip \"$f\"
  fi
done

for f in *.vcf.gz; do
  if [[ ! -f \"$f.tbi\" ]]; then
    tabix -p vcf \"$f\"
  fi
done

echo 'All steps complete.'
"
  
  # Step 1: Write the script to a file
  writeLines(bash_script, "run_pipeline.sh")
  system("chmod +x run_pipeline.sh")
  
  # Step 2: Run the script
  system("./run_pipeline.sh")
}
}



#######################################################################################################################
#########   Step 2: Simulated Data Setup          ########################################################
#######################################################################################################################

#Here we make 4 SNP chip sizes with an SNP array (1728 (close to Slovenian data), 800, 160,16)
#Doing smaller SNPs to test the accuracy of pedigree reconstruction (in a different script)

# Step 1: Create the founder genomes from SIMplyBee
{
founderGenomes <- simulateHoneyBeeGenomes(nCar = 10,
                                           nChr = 16,
                                           nSegSites = 1000)
}

#Step 2: Set up the simulation parameters
{
SP <- SimParamBee$new(founderGenomes, csdChr = ifelse(16 >= 3, 3, 1), nCsdAlleles = 128)
SP$nWorkers <- 30
# Track the pedigree
SP$setTrackPed(TRUE)
# Track the recombination
SP$setTrackRec(TRUE)
# Add a SNP chip (Audrey's) 
createArray(array_name = 'BigArray', array_number = 1, nChr = 16, segSites = rep(1000, 16), nSNPPerChr = rep(261, 16), pop = founderGenomes)
#save array
array = colnames(pullSnpGeno(pop = founderGenomes, snpChip = 1, simParam = SP))
write.table(array, "Bigarray.txt", sep = " ", na = "NA", quote = F, row.names = FALSE, col.names = FALSE)

# define csd chromomsome
csdChr <- SP$csdChr
# Add traits - taken from the QuantGen vignettte (find on simplybee.info)
mean <- c(20, 0)
varA <- c(1, 1 / SP$nWorkers)
corA <- matrix(data = c( 1.0, -0.5,
                         -0.5,  1.0), nrow = 2, byrow = TRUE)
SP$addTraitA(nQtlPerChr = 100, mean = mean, var = varA, corA = corA,
             name = c("queenTrait", "workersTrait"))

varE <- c(3, 3 / SP$nWorkers)
corE <- matrix(data = c(1.0, 0.3,
                        0.3, 1.0), nrow = 2, byrow = TRUE)
SP$setVarE(varE = varE, corE = corE)
}

#Step 3: Establish each of the SNP chip sizes 
{
#Create vector for array size 1728, 800, 160,16
nSNP_array <- rbind(c(108, 4L), c(50, 3L), c(10, 2L), c(1,1L))

#nest array 
for (n in 2:5){ #we create the first array (larger one) outside of the loop and subset it to create the others
  tmp1 = vector()
  for (chr in 1:16){
    tmp2 = sample(colnames(pullSnpGeno(pop = founderGenomes, snpChip = (n-1), chr = chr, simParam = SP)), size = nSNP_array[(n-1),1], replace = FALSE)
    
    tmp1 = c(tmp1, tmp2)
  }
  SP$addSnpChipByName(tmp1, name=paste0('SNP_', nSNP_array[n-1,2]))
  write.table(tmp1, paste0("SNP_", nSNP_array[(n-1),2], "_array.txt"), sep = " ", na = "NA", quote = F, row.names = FALSE, col.names = FALSE)
}

}

#Step 4: Step up the simulation
{
#Create base population
basePop <- createVirginQueens(founderGenomes)
#Create drones to create dpc (THAT ARE RELATED BECAUSE MATING STATION!!!)
pFathers <- nFathersPoisson
nDronesPerQueen = 50
drones <- createDrones(basePop[1], nInd = nDronesPerQueen)
fathers <- pullDroneGroupsFromDCA(drones, n = 1, nDrones = pFathers)

#Create a colony and cross it
colony1 <- createColony(x = basePop[2])
colony1 <- SIMplyBee::cross(colony1, drones = fathers[[1]])


#create Mating station Drones where the dpcs are
DPQs <- createVirginQueens(colony1, nInd = 4)
matingStation_drones <- createDrones(DPQs, nInd = nDronesPerQueen)

# Introduce 8 visiting virgin queens
virgin_queens <- basePop[3:10]

# cross mating station DCA  with virgin queens
queens <- SIMplyBee::cross(x = virgin_queens, drones = pullDroneGroupsFromDCA(DCA = matingStation_drones, n = 8, nDrones = pFathers))
queen_colonies <- createMultiColony(queens)
#create 30 workers per queen
queen_colonies <- buildUp(queen_colonies, nWorkers = 30, exact = TRUE)

#Collect the SNP genotypes of the individuals we want.
workers <- mergePops(getWorkers(queen_colonies))
queens1 <- mergePops(getQueen(queen_colonies))
workers_fathers <- workers@father
fathers <- matingStation_drones[matingStation_drones@id %in% workers_fathers, ]


PopList <- list(queens1, DPQs, workers, fathers)
PopMerged <- mergePops(PopList)

Pop2List <- list(queens1, DPQs, workers)
PopMerged_noFathers <- mergePops(Pop2List)

col <- as.data.frame(cbind(workers@id, workers@mother, workers@father))
colnames(col) <- c("id", "mother", "father")

# Create a dataframe with father IDs and corresponding mother IDs
father_mother_df <- data.frame(fathers@id, fathers@mother)

# Create an empty vector to store the values
dpc_values <- character(nrow(col))

# Loop through each row of col
for (i in 1:nrow(col)) {
  # Check if the father in col matches any id in father_mother_df
  if (col$father[i] %in% father_mother_df$fathers.id) {
    # If there's a match, assign the corresponding mother value to dpc_values
    dpc_values[i] <- father_mother_df$fathers.mother[which(father_mother_df$fathers.id == col$father[i])]
  } else {
    # If there's no match, assign NA
    dpc_values[i] <- NA
  }
}
col$dpc <- dpc_values

#Save the pedigree
write_csv(col, file = "worker_pedigree.csv", col_names = TRUE)
}

#Step 4: Add errors (Genotyping and sampling errors)
{
error  = 0.0001
error2 = 0.00001
missing = 0.005 #using this as it matches the missingness in the real Slovenian Data 
sampling_error = 0.005

# Genotyping error function
generateGenoErr <- function(geno, error, error2, sampling_error, missing) {
  nLoci <- ncol(geno)
  nInd <- nrow(geno)
  # error for sampling error
  for (ind in 1:nInd) {
    for (locus in 1:nLoci) {
      if (geno[ind, locus] != 9 && rbinom(1, 1, sampling_error) == 1) {
        geno[ind, locus] <- 9
      }
    }
  }
  # error for genotyping error
  for (ind in 1:nInd) {
    for (locus in 1:nLoci) {
      if (geno[ind, locus] != 9 && rbinom(1, 1, error) == 1) {
        if (geno[ind, locus] == 0) {
          geno[ind, locus] <- sample(c(1, 2, 9), 1, prob = c(1 - error2 - missing, error2, missing))
        } else if (geno[ind, locus] == 1) {
          geno[ind, locus] <- sample(c(0, 2, 9), 1, prob = c(1 - error2 - missing, 1 - error2 - missing, missing))
        } else if (geno[ind, locus] == 2) {
          geno[ind, locus] <- sample(c(0, 1, 9), 1, prob = c(error2, 1 - error2 - missing, missing))
        }
      }
    }
  }
  return(geno)
}

# Function to read and process a PLINK PED file
process_ped_file <- function(ped_data, error, error2, sampling_error, missing) {
  
  # Extract the genotype data (columns 7 onwards in PED format)
  geno_data <- ped_data[, 7:ncol(ped_data)]
  
  # Convert genotype data from character to numeric
  geno_data <- as.data.frame(apply(geno_data, 2, as.numeric))
  
  # Generate genotyping errors
  geno_data <- generateGenoErr(geno_data, error, error2, sampling_error, missing)
  
  # Combine the original metadata with the modified genotype data
  ped_data[, 7:ncol(ped_data)] <- geno_data
  
  return(ped_data)
}

SNP4_GE_nQC_ped01 <- process_ped_file(ped_data = SNP4_nGE_nQC_ped01, error = error, error2 = error2, sampling_error = sampling_error, missing = missing)

switchToNA <- function(genoErr) {
  # Check if genoErr is a matrix
  if (!is.matrix(genoErr)) {
    stop("Input must be a matrix.")
  }
  
  # Replace 9s with NA
  genoErr[genoErr == 9] <- NA
  genoErr[genoErr == -9] <- NA
  return(genoErr)
}


#You cant have a 9 in the ped file for the vcf it needs to be  0 0 so lets change it 
# Function to check and correct half-missing genotypes
correct_half_missing <- function(genotype1, genotype2) {
  if (genotype1 == 9) genotype1 <- 0
  if (genotype2 == 9) genotype2 <- 0
  if ((genotype1 == 0 && genotype2 != 0) || (genotype1 != 0 && genotype2 == 0)) {
    return(c(0, 0)) # Replace half-missing with full-missing
  } else {
    return(c(genotype1, genotype2))
  }
}

# Apply the function to each genotype pair in the dataframe to preQC dataframe 
for (i in seq(7, ncol(GE_nQC_ped), by = 2)) { # Start from the 7th column (first genotype) and step by 2
  GE_nQC_ped[, c(i, i+1)] <- t(apply(GE_nQC_ped[, c(i, i+1)], 1, function(x) correct_half_missing(x[1], x[2])))
}

#Now we make a loop for all of the different arrays, going through and saving the needed info for each one. 
#We are saving the genotype matrix (both full version and one with genotypic errors), along with the map files of the SNP chips and the haplotypes for the full genotypic data
for (n in 1:nrow(nSNP_array)){
  SNP_array = as.vector(read_table(paste0('SNP_', nSNP_array[n,2], '_array.txt'), show_col_types = FALSE, col_names = FALSE))
  SP$addSnpChipByName(SNP_array$X1, name=paste0('array_', nSNP_array[n,2]))
  
  #Using PopMerged_noFathers since the genotypic info of the drones is not available in the real Slov data 
  geno_noError = pullSnpGeno(PopMerged_noFathers, snpChip= paste0('array_', nSNP_array[n,2]), chr = NULL, asRaw = FALSE, simParam=SP)
  haplo_noError = pullSnpHaplo(PopMerged_noFathers, snpChip= paste0('array_', nSNP_array[n,2]), chr = NULL, asRaw = FALSE, simParam=SP)
  SNPmap_file = getSnpMap(snpChip= paste0('array_', nSNP_array[n,2]), simParam=SP)
  geno_WithError = generateGenoErr(geno = geno_noError, error = error, error2 = error2, sampling_error = sampling_error, missing = missing)
  
  write.table(SNPmap_file[c("chr", "id", "site")], paste0("SNPMapFile_", nSNP_array[n,2] , ".txt"), sep = " ", na = "0", quote = F, row.names = TRUE, col.names = TRUE)
  write.table(geno_noError, paste0("SNP", nSNP_array[n,2], "_geno_noError.txt"), sep = " ", na = "0", quote = F, row.names = TRUE, col.names = TRUE)
  write.table(haplo_noError, paste0("SNP", nSNP_array[n,2], "_haplo_noError.txt"), sep = " ", na = "0", quote = F, row.names = TRUE, col.names = TRUE)
  write.table(geno_WithError, paste0("SNP", nSNP_array[n,2], "_geno_WithError.txt"), sep = " ", na = "0", quote = F, row.names = TRUE, col.names = TRUE)
  
  #write plink files (ped and map) for each SNP array size
  #we use the population without fathers to mimic the real data since we would not have the drone data
  writePlink(PopMerged_noFathers, baseName = paste0('SNP_', nSNP_array[n,2], '_arrayNoGenoError'), snpChip = paste0('array_', nSNP_array[n,2]), simParam = SP)
  
  rm(geno_noError, haplo_noError, geno_WithError, SNPmap_file, SNP_array)
}

#we have the genotyping errors, but we need to make them into a plink ped file and that has to be done manually 
#first we take the true haplotypes so that we can use them as a reference 
true_haplotype1 <- haplo_noError[grep("_1$", rownames(haplo_noError)),]
true_haplotype2 <- haplo_noError[grep("_2$", rownames(haplo_noError)),]

#this function will go through the genotypes with errors and use the true haplotypes as a reference for the genotype with value 1. 
genotypes_to_haplotypes <- function(genotypes, ref_hap1, ref_hap2) {
  inferred_hap1 <- matrix(0, nrow = nrow(genotypes), ncol = ncol(genotypes))
  inferred_hap2 <- matrix(0, nrow = nrow(genotypes), ncol = ncol(genotypes))
  
  for (i in 1:nrow(genotypes)) {
    for (j in 1:ncol(genotypes)) {
      if (genotypes[i, j] == 9) {
        inferred_hap1[i, j] <- 9
        inferred_hap2[i, j] <- 9
      } else if (genotypes[i, j] == ref_hap1[i, j] + ref_hap2[i, j]) {
        inferred_hap1[i, j] <- ref_hap1[i, j]
        inferred_hap2[i, j] <- ref_hap2[i, j]
      } else if (genotypes[i, j] == 2) {
        inferred_hap1[i, j] <- 1
        inferred_hap2[i, j] <- 1
      } else if (genotypes[i, j] == 0) {
        inferred_hap1[i, j] <- 0
        inferred_hap2[i, j] <- 0
      } else if (genotypes[i, j] == 1) {
        inferred_hap1[i, j] <- ref_hap1[i, j]
        inferred_hap2[i, j] <- 1 - ref_hap1[i, j]
      }
    }
  }
  
  list(haplotype_1 = inferred_hap1, haplotype_2 = inferred_hap2)
}

genoError_haplotypes <- genotypes_to_haplotypes(genotypes = geno_WithError, ref_hap1 = true_haplotype1, ref_hap2 = true_haplotype2)
genoError_hap1 <- genoError_haplotypes[["haplotype_1"]]
rownames(genoError_hap1) <- rownames(true_haplotype1)
genoError_hap2 <- genoError_haplotypes[["haplotype_2"]]
rownames(genoError_hap2) <- rownames(true_haplotype2)

#now we have to order it correctly 
combined_matrix <- rbind(genoError_hap1, genoError_hap2)
colnames(combined_matrix) <- colnames(haplo_noError)
# Extract the numeric part and whether it's _1 or _2
numeric_part <- as.numeric(gsub("_.*", "", rownames(combined_matrix)))
suffix_part <- gsub(".*_", "", rownames(combined_matrix))
# Create a data frame to order the rows
order_df <- data.frame(numeric_part, suffix_part)
# Order first by the numeric part, then by the suffix (_1 before _2)
ordered_indices <- order(order_df$numeric_part, order_df$suffix_part)
# Reorder the combined matrix
ordered_matrix <- combined_matrix[ordered_indices, ]

combine_haplotypes <- function(haplotype_matrix) {
  # Extract unique IDs
  unique_ids <- unique(gsub("_.*", "", rownames(haplotype_matrix)))
  
  # Create new column names with the correct order
  original_colnames <- colnames(haplotype_matrix)
  new_colnames <- unlist(lapply(original_colnames, function(col) {
    paste(col, c("_1", "_2"), sep = "")
  }))
  
  # Initialize new matrix
  new_matrix <- matrix(NA, nrow = length(unique_ids), ncol = length(new_colnames))
  rownames(new_matrix) <- unique_ids
  colnames(new_matrix) <- new_colnames
  
  # Fill the new matrix
  for (id in unique_ids) {
    rows_to_combine <- haplotype_matrix[grep(paste0("^", id, "_"), rownames(haplotype_matrix)), ]
    # Fill in columns in the correct order
    for (i in seq_along(original_colnames)) {
      new_matrix[id, (2 * (i - 1) + 1):(2 * i)] <- rows_to_combine[, i]
    }
  }
  
  return(new_matrix)
}

# Apply the function
GE_ped_format <- combine_haplotypes(ordered_matrix)

GE_ped <- cbind(nGE_ped[,c(1:6)], GE_ped_format)
#the map is the same so just use this ped file and the map file to go and do the rest of the stuff.
}

#We come out with three versions: True haplotypes (perfect, fully phased), GE haplotypes (with the genotyping errors), 
#also make a nGE haplotyes (also true haplotypes) at this stage but will go through phasing to see if phasing changes results.

#These then undergo the same quality controls as the slovenian data :) 




#######################################################################################################################
#########   Step 3: Pedigree reconstruction          ########################################################
#######################################################################################################################

#Run the genotypes of all data through different pedigree recontruction software (sequoia, KING, AlphaAssign, Colony)
#Use the simulated data with the known pedigree to assess the accuracy of the pedigree reconstruction 
#At 1.7k SNP all of the software were 100% accurate in the simulated. 
# AlphaAssign gave us the highest number of "sires" so we use that as the Slovenian "pedigree". 
#We know the Slovenian mothers but not sires (DPQs) so this creates a most accurate pedigree




#######################################################################################################################
#########   Step 4: Phasing          ########################################################
#######################################################################################################################

#Everything has been generated, pedigree files have been generated, quality control has been done. Now to phase!

#We will phase with and without pedigree information to see what difference this makes to phasing
#We are phasing by chromosome:ensures biological accuracy, improves computational efficiency, 
#                             and avoids confusing the phasing algorithm with unrelated recombination patterns across different chromosomes.

#PHASING WITH PEDIGREE
for (i in 1:16) {
  beagle_cmd <- paste0(
    "java -jar beagle.4.0.jar gt=Slov_fM_AC_QC_noDupPos.vcf.gz out=Slov_PHASED_chr", 
    i, " window=200 overlap=50 phase-its=200 impute-its=1000 burnin-its=200 chrom=", i
  )
  system(beagle_cmd)
}

#PHASING WITHOUT PEDIGREE
for (i in 1:16) {
  beagle_cmd <- paste0(
    "java -jar beagle.4.0.jar gt=Slov_fM_AC_QC_noDupPos.vcf.gz ped=Slov_pedigree.txt out=Slov_PHASED_chr", 
    i, " window=200 overlap=50 phase-its=200 impute-its=1000 burnin-its=200 chrom=", i
  )
  system(beagle_cmd)
}

# Use PLINK to convert phased VCFs to .map files
vcf_files <- list.files(pattern = "Slov_PHASED_chr.*\\.vcf\\.gz$")
for (f in vcf_files) {
  base <- sub("\\.vcf\\.gz$", "", f)
  plink_cmd <- paste0("./plink --vcf ", f, " --recode --double-id --allow-extra-chr --out ", base)
  system(plink_cmd)
}

#Compress and index all uncompressed VCFs
vcf_plain <- list.files(pattern = "\\.vcf$")
for (f in vcf_plain) {
  if (!file.exists(paste0(f, ".gz"))) {
    system(paste("bgzip", f))
  }
}
vcf_gz <- list.files(pattern = "\\.vcf\\.gz$")
for (f in vcf_gz) {
  if (!file.exists(paste0(f, ".tbi"))) {
    system(paste("tabix -p vcf", f))
  }
}




#######################################################################################################################
#########   Step 5: Manually convert the vcf files in ped files (something weird with the PLINK conversion) ###########
#######################################################################################################################

#Output is a haplotype dataframe with columns (id1_1, id1_2) and rows the SNPs coded 0/1
#E.g We will get Slovenian haplotypes from phasing with pedigree AND phased without pedigree (This differentiation is important for the next calculations)

#1. Process the phased simulated data 
{
setwd("~/Desktop/Slovenia data/Test_simulated")
simulated_pedigree <- read.csv("worker_pedigree.csv")
#ped order is always mum/dpc/workers
ids_all <- c(unique(simulated_pedigree$mother), unique(simulated_pedigree$dpc), unique(simulated_pedigree$id))
ind_id_1 <- paste(ids_all, "1", sep = "_")
ind_id_2 <- paste(ids_all, "2", sep = "_")

#Get all of the vcf files of all phased chromosomes and convert into a more manageable format 
#Function to do this:  
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

#Phased with Pedigree from pedigree reconstruction
#nGE = no genotyping errors
{
  setwd("~/Desktop/Slovenia data/May25/Phasing/phased/simulated/nGE/withPed")
  
  nGE_chr1 <- convert_VCF(vcf_file = "nGE_QC_PHASEDwithPed_chr1.vcf.gz", map_file = "nGE_QC_PHASEDwithPed_chr1.map")
  nGE_chr2 <- convert_VCF(vcf_file = "nGE_QC_PHASEDwithPed_chr2.vcf.gz", map_file = "nGE_QC_PHASEDwithPed_chr2.map")
  nGE_chr3 <- convert_VCF(vcf_file = "nGE_QC_PHASEDwithPed_chr3.vcf.gz", map_file = "nGE_QC_PHASEDwithPed_chr3.map")
  nGE_chr4 <- convert_VCF(vcf_file = "nGE_QC_PHASEDwithPed_chr4.vcf.gz", map_file = "nGE_QC_PHASEDwithPed_chr4.map")
  nGE_chr5 <- convert_VCF(vcf_file = "nGE_QC_PHASEDwithPed_chr5.vcf.gz", map_file = "nGE_QC_PHASEDwithPed_chr5.map")
  nGE_chr6 <- convert_VCF(vcf_file = "nGE_QC_PHASEDwithPed_chr6.vcf.gz", map_file = "nGE_QC_PHASEDwithPed_chr6.map")
  nGE_chr7 <- convert_VCF(vcf_file = "nGE_QC_PHASEDwithPed_chr7.vcf.gz", map_file = "nGE_QC_PHASEDwithPed_chr7.map")
  nGE_chr8 <- convert_VCF(vcf_file = "nGE_QC_PHASEDwithPed_chr8.vcf.gz", map_file = "nGE_QC_PHASEDwithPed_chr8.map")
  nGE_chr9 <- convert_VCF(vcf_file = "nGE_QC_PHASEDwithPed_chr9.vcf.gz", map_file = "nGE_QC_PHASEDwithPed_chr9.map")
  nGE_chr10 <- convert_VCF(vcf_file = "nGE_QC_PHASEDwithPed_chr10.vcf.gz", map_file = "nGE_QC_PHASEDwithPed_chr10.map")
  nGE_chr11 <- convert_VCF(vcf_file = "nGE_QC_PHASEDwithPed_chr11.vcf.gz", map_file = "nGE_QC_PHASEDwithPed_chr11.map")
  nGE_chr12 <- convert_VCF(vcf_file = "nGE_QC_PHASEDwithPed_chr12.vcf.gz", map_file = "nGE_QC_PHASEDwithPed_chr12.map")
  nGE_chr13 <- convert_VCF(vcf_file = "nGE_QC_PHASEDwithPed_chr13.vcf.gz", map_file = "nGE_QC_PHASEDwithPed_chr13.map")
  nGE_chr14 <- convert_VCF(vcf_file = "nGE_QC_PHASEDwithPed_chr14.vcf.gz", map_file = "nGE_QC_PHASEDwithPed_chr14.map")
  nGE_chr15 <- convert_VCF(vcf_file = "nGE_QC_PHASEDwithPed_chr15.vcf.gz", map_file = "nGE_QC_PHASEDwithPed_chr15.map")
  nGE_chr16 <- convert_VCF(vcf_file = "nGE_QC_PHASEDwithPed_chr16.vcf.gz", map_file = "nGE_QC_PHASEDwithPed_chr16.map")
  
  nGE_allchroms_phasedwithPed <- cbind(nGE_chr1,nGE_chr2,nGE_chr3,nGE_chr4,nGE_chr5,nGE_chr6,nGE_chr7,nGE_chr8,
                                       nGE_chr9,nGE_chr10,nGE_chr11,nGE_chr12,nGE_chr13,nGE_chr14,nGE_chr15,nGE_chr16)
  
  rm(nGE_chr1,nGE_chr2,nGE_chr3,nGE_chr4,nGE_chr5,nGE_chr6,nGE_chr7,nGE_chr8,
     nGE_chr9,nGE_chr10,nGE_chr11,nGE_chr12,nGE_chr13,nGE_chr14,nGE_chr15,nGE_chr16)
}
#GE with genotyping errors
{
  setwd("~/Desktop/Slovenia data/May25/Phasing/phased/simulated/GE/WithPed")
  
  GE_chr1 <- convert_VCF(vcf_file = "GE_QC_PHASEDwithPed_chr1.vcf.gz", map_file = "GE_QC_PHASEDwithPed_chr1.map")
  GE_chr2 <- convert_VCF(vcf_file = "GE_QC_PHASEDwithPed_chr2.vcf.gz", map_file = "GE_QC_PHASEDwithPed_chr2.map")
  GE_chr3 <- convert_VCF(vcf_file = "GE_QC_PHASEDwithPed_chr3.vcf.gz", map_file = "GE_QC_PHASEDwithPed_chr3.map")
  GE_chr4 <- convert_VCF(vcf_file = "GE_QC_PHASEDwithPed_chr4.vcf.gz", map_file = "GE_QC_PHASEDwithPed_chr4.map")
  GE_chr5 <- convert_VCF(vcf_file = "GE_QC_PHASEDwithPed_chr5.vcf.gz", map_file = "GE_QC_PHASEDwithPed_chr5.map")
  GE_chr6 <- convert_VCF(vcf_file = "GE_QC_PHASEDwithPed_chr6.vcf.gz", map_file = "GE_QC_PHASEDwithPed_chr6.map")
  GE_chr7 <- convert_VCF(vcf_file = "GE_QC_PHASEDwithPed_chr7.vcf.gz", map_file = "GE_QC_PHASEDwithPed_chr7.map")
  GE_chr8 <- convert_VCF(vcf_file = "GE_QC_PHASEDwithPed_chr8.vcf.gz", map_file = "GE_QC_PHASEDwithPed_chr8.map")
  GE_chr9 <- convert_VCF(vcf_file = "GE_QC_PHASEDwithPed_chr9.vcf.gz", map_file = "GE_QC_PHASEDwithPed_chr9.map")
  GE_chr10 <- convert_VCF(vcf_file = "GE_QC_PHASEDwithPed_chr10.vcf.gz", map_file = "GE_QC_PHASEDwithPed_chr10.map")
  GE_chr11 <- convert_VCF(vcf_file = "GE_QC_PHASEDwithPed_chr11.vcf.gz", map_file = "GE_QC_PHASEDwithPed_chr11.map")
  GE_chr12 <- convert_VCF(vcf_file = "GE_QC_PHASEDwithPed_chr12.vcf.gz", map_file = "GE_QC_PHASEDwithPed_chr12.map")
  GE_chr13 <- convert_VCF(vcf_file = "GE_QC_PHASEDwithPed_chr13.vcf.gz", map_file = "GE_QC_PHASEDwithPed_chr13.map")
  GE_chr14 <- convert_VCF(vcf_file = "GE_QC_PHASEDwithPed_chr14.vcf.gz", map_file = "GE_QC_PHASEDwithPed_chr14.map")
  GE_chr15 <- convert_VCF(vcf_file = "GE_QC_PHASEDwithPed_chr15.vcf.gz", map_file = "GE_QC_PHASEDwithPed_chr15.map")
  GE_chr16 <- convert_VCF(vcf_file = "GE_QC_PHASEDwithPed_chr16.vcf.gz", map_file = "GE_QC_PHASEDwithPed_chr16.map")
  
  GE_allchroms_phasedwithPed <- cbind(GE_chr1,GE_chr2,GE_chr3,GE_chr4,GE_chr5,GE_chr6,GE_chr7,GE_chr8,
                                      GE_chr9,GE_chr10,GE_chr11,GE_chr12,GE_chr13,GE_chr14,GE_chr15,GE_chr16)
  
  rm(GE_chr1,GE_chr2,GE_chr3,GE_chr4,GE_chr5,GE_chr6,GE_chr7,GE_chr8,
     GE_chr9,GE_chr10,GE_chr11,GE_chr12,GE_chr13,GE_chr14,GE_chr15,GE_chr16)
}

#Phased without a pedigree file
#nGE
{
  setwd("~/Desktop/Slovenia data/May25/Phasing/phased/simulated/nGE/NoPed")
  
  nGE_chr1 <- convert_VCF(vcf_file = "nGE_QC_PHASED-NoPed_chr1.vcf.gz", map_file = "nGE_QC_PHASED-NoPed_chr1.map")
  nGE_chr2 <- convert_VCF(vcf_file = "nGE_QC_PHASED-NoPed_chr2.vcf.gz", map_file = "nGE_QC_PHASED-NoPed_chr2.map")
  nGE_chr3 <- convert_VCF(vcf_file = "nGE_QC_PHASED-NoPed_chr3.vcf.gz", map_file = "nGE_QC_PHASED-NoPed_chr3.map")
  nGE_chr4 <- convert_VCF(vcf_file = "nGE_QC_PHASED-NoPed_chr4.vcf.gz", map_file = "nGE_QC_PHASED-NoPed_chr4.map")
  nGE_chr5 <- convert_VCF(vcf_file = "nGE_QC_PHASED-NoPed_chr5.vcf.gz", map_file = "nGE_QC_PHASED-NoPed_chr5.map")
  nGE_chr6 <- convert_VCF(vcf_file = "nGE_QC_PHASED-NoPed_chr6.vcf.gz", map_file = "nGE_QC_PHASED-NoPed_chr6.map")
  nGE_chr7 <- convert_VCF(vcf_file = "nGE_QC_PHASED-NoPed_chr7.vcf.gz", map_file = "nGE_QC_PHASED-NoPed_chr7.map")
  nGE_chr8 <- convert_VCF(vcf_file = "nGE_QC_PHASED-NoPed_chr8.vcf.gz", map_file = "nGE_QC_PHASED-NoPed_chr8.map")
  nGE_chr9 <- convert_VCF(vcf_file = "nGE_QC_PHASED-NoPed_chr9.vcf.gz", map_file = "nGE_QC_PHASED-NoPed_chr9.map")
  nGE_chr10 <- convert_VCF(vcf_file = "nGE_QC_PHASED-NoPed_chr10.vcf.gz", map_file = "nGE_QC_PHASED-NoPed_chr10.map")
  nGE_chr11 <- convert_VCF(vcf_file = "nGE_QC_PHASED-NoPed_chr11.vcf.gz", map_file = "nGE_QC_PHASED-NoPed_chr11.map")
  nGE_chr12 <- convert_VCF(vcf_file = "nGE_QC_PHASED-NoPed_chr12.vcf.gz", map_file = "nGE_QC_PHASED-NoPed_chr12.map")
  nGE_chr13 <- convert_VCF(vcf_file = "nGE_QC_PHASED-NoPed_chr13.vcf.gz", map_file = "nGE_QC_PHASED-NoPed_chr13.map")
  nGE_chr14 <- convert_VCF(vcf_file = "nGE_QC_PHASED-NoPed_chr14.vcf.gz", map_file = "nGE_QC_PHASED-NoPed_chr14.map")
  nGE_chr15 <- convert_VCF(vcf_file = "nGE_QC_PHASED-NoPed_chr15.vcf.gz", map_file = "nGE_QC_PHASED-NoPed_chr15.map")
  nGE_chr16 <- convert_VCF(vcf_file = "nGE_QC_PHASED-NoPed_chr16.vcf.gz", map_file = "nGE_QC_PHASED-NoPed_chr16.map")
  
  nGE_allchroms_phasedNoPed <- cbind(nGE_chr1,nGE_chr2,nGE_chr3,nGE_chr4,nGE_chr5,nGE_chr6,nGE_chr7,nGE_chr8,
                                     nGE_chr9,nGE_chr10,nGE_chr11,nGE_chr12,nGE_chr13,nGE_chr14,nGE_chr15,nGE_chr16)
  
  rm(nGE_chr1,nGE_chr2,nGE_chr3,nGE_chr4,nGE_chr5,nGE_chr6,nGE_chr7,nGE_chr8,
     nGE_chr9,nGE_chr10,nGE_chr11,nGE_chr12,nGE_chr13,nGE_chr14,nGE_chr15,nGE_chr16)
}
#GE
{
  setwd("~/Desktop/Slovenia data/May25/Phasing/phased/simulated/GE/NoPed")
  
  GE_chr1 <- convert_VCF(vcf_file = "GE_QC_PHASED-NoPed_chr1.vcf.gz", map_file = "GE_QC_PHASED-NoPed_chr1.map")
  GE_chr2 <- convert_VCF(vcf_file = "GE_QC_PHASED-NoPed_chr2.vcf.gz", map_file = "GE_QC_PHASED-NoPed_chr2.map")
  GE_chr3 <- convert_VCF(vcf_file = "GE_QC_PHASED-NoPed_chr3.vcf.gz", map_file = "GE_QC_PHASED-NoPed_chr3.map")
  GE_chr4 <- convert_VCF(vcf_file = "GE_QC_PHASED-NoPed_chr4.vcf.gz", map_file = "GE_QC_PHASED-NoPed_chr4.map")
  GE_chr5 <- convert_VCF(vcf_file = "GE_QC_PHASED-NoPed_chr5.vcf.gz", map_file = "GE_QC_PHASED-NoPed_chr5.map")
  GE_chr6 <- convert_VCF(vcf_file = "GE_QC_PHASED-NoPed_chr6.vcf.gz", map_file = "GE_QC_PHASED-NoPed_chr6.map")
  GE_chr7 <- convert_VCF(vcf_file = "GE_QC_PHASED-NoPed_chr7.vcf.gz", map_file = "GE_QC_PHASED-NoPed_chr7.map")
  GE_chr8 <- convert_VCF(vcf_file = "GE_QC_PHASED-NoPed_chr8.vcf.gz", map_file = "GE_QC_PHASED-NoPed_chr8.map")
  GE_chr9 <- convert_VCF(vcf_file = "GE_QC_PHASED-NoPed_chr9.vcf.gz", map_file = "GE_QC_PHASED-NoPed_chr9.map")
  GE_chr10 <- convert_VCF(vcf_file = "GE_QC_PHASED-NoPed_chr10.vcf.gz", map_file = "GE_QC_PHASED-NoPed_chr10.map")
  GE_chr11 <- convert_VCF(vcf_file = "GE_QC_PHASED-NoPed_chr11.vcf.gz", map_file = "GE_QC_PHASED-NoPed_chr11.map")
  GE_chr12 <- convert_VCF(vcf_file = "GE_QC_PHASED-NoPed_chr12.vcf.gz", map_file = "GE_QC_PHASED-NoPed_chr12.map")
  GE_chr13 <- convert_VCF(vcf_file = "GE_QC_PHASED-NoPed_chr13.vcf.gz", map_file = "GE_QC_PHASED-NoPed_chr13.map")
  GE_chr14 <- convert_VCF(vcf_file = "GE_QC_PHASED-NoPed_chr14.vcf.gz", map_file = "GE_QC_PHASED-NoPed_chr14.map")
  GE_chr15 <- convert_VCF(vcf_file = "GE_QC_PHASED-NoPed_chr15.vcf.gz", map_file = "GE_QC_PHASED-NoPed_chr15.map")
  GE_chr16 <- convert_VCF(vcf_file = "GE_QC_PHASED-NoPed_chr16.vcf.gz", map_file = "GE_QC_PHASED-NoPed_chr16.map")
  
  GE_allchroms_phasedNoPed <- cbind(GE_chr1,GE_chr2,GE_chr3,GE_chr4,GE_chr5,GE_chr6,GE_chr7,GE_chr8,
                                    GE_chr9,GE_chr10,GE_chr11,GE_chr12,GE_chr13,GE_chr14,GE_chr15,GE_chr16)
  
  rm(GE_chr1,GE_chr2,GE_chr3,GE_chr4,GE_chr5,GE_chr6,GE_chr7,GE_chr8,
     GE_chr9,GE_chr10,GE_chr11,GE_chr12,GE_chr13,GE_chr14,GE_chr15,GE_chr16)
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
nGE_haplotypes_phasedwithPed <- get_out_haplotypes(ped_matrix = nGE_allchroms_phasedwithPed  , ind_id_1 = ind_id_1, ind_id_2 = ind_id_2)
GE_haplotypes_phasedwithPed <- get_out_haplotypes(ped_matrix = GE_allchroms_phasedwithPed, ind_id_1 = ind_id_1, ind_id_2 = ind_id_2)

nGE_haplotypes_phasedNoPed <- get_out_haplotypes(ped_matrix = nGE_allchroms_phasedNoPed  , ind_id_1 = ind_id_1, ind_id_2 = ind_id_2)
GE_haplotypes_phasedNoPed <- get_out_haplotypes(ped_matrix = GE_allchroms_phasedNoPed, ind_id_1 = ind_id_1, ind_id_2 = ind_id_2)


#convert to 0/1 format 
nGE_phasedhaplotypes_phasedwithPed <- apply(nGE_haplotypes_phasedwithPed, 2, convert_genotypes)
rownames(nGE_phasedhaplotypes_phasedwithPed) <- rownames(nGE_haplotypes_phasedwithPed)

GE_phasedhaplotypes_phasedwithPed <- apply(GE_haplotypes_phasedwithPed, 2, convert_genotypes)
rownames(GE_phasedhaplotypes_phasedwithPed) <- rownames(GE_haplotypes_phasedwithPed)

nGE_phasedhaplotypes_phasedNoPed <- apply(nGE_haplotypes_phasedNoPed, 2, convert_genotypes)
rownames(nGE_phasedhaplotypes_phasedNoPed) <- rownames(nGE_haplotypes_phasedNoPed)

GE_phasedhaplotypes_phasedNoPed <- apply(GE_haplotypes_phasedNoPed, 2, convert_genotypes)
rownames(GE_phasedhaplotypes_phasedNoPed) <- rownames(GE_haplotypes_phasedNoPed)

setwd("~/Desktop/Slovenia data/May25/Phasing/pre-phase/simulated")
GE_map <- read.table("GE_QC_prePhase.map")
nGE_map <- read.table("nGE_QC_prePhase.map")
}

#2. Get the true haplotypes
{
load("~/Desktop/Slovenia data/Test_simulated/SP_object.Rdata")
load("~/Desktop/Slovenia data/Test_simulated/Pop_withNoFathers.Rdata")
##### ALPHASIMR FORMAT Real haplotypes #####
true_haplotypes <- pullSnpHaplo(PopMerged_noFathers)
true_map <- getGenMap(SP)
}

#3. Process the Slovenian real data
{
setwd("~/Desktop/Slovenia data/PLINK/Phased_by_chromosome/Slov/All")
Slov_map <-  read.table("Slov_all_phased_chr.map")

Slov_pedigree <- read.table("Slov_full_pedigree.txt")
colnames(Slov_pedigree) <- c("id","dpc","mother")

#slov Pedigree prior to ped reconstruction 
setwd("~/Desktop/Slovenia data/May25/Phasing/pre-phase/real")
Slov_ids_preRecon <- read.table("Slov_filtered_QC_noDup.ped")
Slov_ids_preRecon <- Slov_ids_preRecon$V2


#After reconstruction
Slov_mothers <- read.table("~/Desktop/Slovenia data/Attempt2/Slov/Not-phased/Colony/data/KnownMothers.txt")
Slov_sire <- read.table("~/Desktop/Slovenia data/PLINK/Slov/Slov AlphaAssign/output/Slov_AlphaAssignOutput.sires")

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

#Phased with pedigree
{
  setwd("~/Desktop/Slovenia data/May25/Phasing/phased/REAL Redo/WithPed")
  
  Slov_chr1 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr1.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr1.map")
  Slov_chr2 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr2.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr2.map")
  Slov_chr3 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr3.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr3.map")
  Slov_chr4 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr4.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr4.map")
  Slov_chr5 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr5.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr5.map")
  Slov_chr6 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr6.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr6.map")
  Slov_chr7 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr7.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr7.map")
  Slov_chr8 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr8.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr8.map")
  Slov_chr9 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr9.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr9.map")
  Slov_chr10 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr10.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr10.map")
  Slov_chr11 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr11.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr11.map")
  Slov_chr12 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr12.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr12.map")
  Slov_chr13 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr13.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr13.map")
  Slov_chr14 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr14.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr14.map")
  Slov_chr15 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr15.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr15.map")
  Slov_chr16 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr16.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_withPed_chr16.map")
  
  
  Slov_allchroms_phasedwithPed <- cbind(Slov_chr1,Slov_chr2,Slov_chr3,Slov_chr4,Slov_chr5,Slov_chr6,Slov_chr7,Slov_chr8,
                                        Slov_chr9,Slov_chr10,Slov_chr11,Slov_chr12,Slov_chr13,Slov_chr14,Slov_chr15,Slov_chr16)
  
  rm(Slov_chr1,Slov_chr2,Slov_chr3,Slov_chr4,Slov_chr5,Slov_chr6,Slov_chr7,Slov_chr8,
     Slov_chr9,Slov_chr10,Slov_chr11,Slov_chr12,Slov_chr13,Slov_chr14,Slov_chr15,Slov_chr16)
}

#Phased with no pedigree
{
  setwd("~/Desktop/Slovenia data/May25/Phasing/phased/REAL Redo/NoPed")
  
  Slov_chr1 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr1.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr1.map")
  Slov_chr2 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr2.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr2.map")
  Slov_chr3 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr3.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr3.map")
  Slov_chr4 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr4.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr4.map")
  Slov_chr5 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr5.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr5.map")
  Slov_chr6 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr6.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr6.map")
  Slov_chr7 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr7.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr7.map")
  Slov_chr8 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr8.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr8.map")
  Slov_chr9 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr9.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr9.map")
  Slov_chr10 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr10.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr10.map")
  Slov_chr11 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr11.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr11.map")
  Slov_chr12 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr12.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr12.map")
  Slov_chr13 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr13.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr13.map")
  Slov_chr14 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr14.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr14.map")
  Slov_chr15 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr15.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr15.map")
  Slov_chr16 <- convert_VCF_Slov(vcf_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr16.vcf.gz", map_file = "Slov_fM_AC_QC_nChr0_nDupPos_PHASED_NoPed_chr16.map")
  
  
  Slov_allchroms_phasedNoPed <- cbind(Slov_chr1,Slov_chr2,Slov_chr3,Slov_chr4,Slov_chr5,Slov_chr6,Slov_chr7,Slov_chr8,
                                      Slov_chr9,Slov_chr10,Slov_chr11,Slov_chr12,Slov_chr13,Slov_chr14,Slov_chr15,Slov_chr16)
  
  rm(Slov_chr1,Slov_chr2,Slov_chr3,Slov_chr4,Slov_chr5,Slov_chr6,Slov_chr7,Slov_chr8,
     Slov_chr9,Slov_chr10,Slov_chr11,Slov_chr12,Slov_chr13,Slov_chr14,Slov_chr15,Slov_chr16)
}


#Get the haplotypes of the pedigree made by the pedigree reconstruction
{
Slov_mother_haplo <- Slov_allchroms_phasedwithPed[rownames(Slov_allchroms_phasedwithPed) %in% Slov_pedigree$mother,]
Slov_dpc_haplo <- Slov_allchroms_phasedwithPed[rownames(Slov_allchroms_phasedwithPed) %in% Slov_pedigree$dpc,]
Slov_workers_haplo <- Slov_allchroms_phasedwithPed[rownames(Slov_allchroms_phasedwithPed) %in% Slov_pedigree$id,]

Slov_all_haplo_phasedwithPed <-  rbind(Slov_mother_haplo, Slov_dpc_haplo, Slov_workers_haplo)

mother_ids <- unique(Slov_pedigree$mother)
dpc_ids <- unique(Slov_pedigree$dpc)
dpc_ids <- dpc_ids[dpc_ids != "0"]
worker_ids <- unique(Slov_pedigree$id)
Slov_all_ids <- c(mother_ids, dpc_ids, worker_ids)
Slov_ind_id_1 <- paste(Slov_all_ids, "1", sep = "_")
Slov_ind_id_2 <- paste(Slov_all_ids, "2", sep = "_")

Slov_haplotypes_phasedwithPed <- get_out_haplotypes(ped_matrix = Slov_all_haplo_phasedwithPed, ind_id_1 = Slov_ind_id_1, ind_id_2 = Slov_ind_id_2)
Slov_phasedhaplotypes_phasedwithPed <- apply(Slov_haplotypes_phasedwithPed, 2, convert_genotypes)
rownames(Slov_phasedhaplotypes_phasedwithPed) <- rownames(Slov_haplotypes_phasedwithPed)
}

#without pedigree
{
Slov_mother_haplo_noped <- Slov_allchroms_phasedNoPed[rownames(Slov_allchroms_phasedNoPed) %in% Slov_pedigree$mother,]
Slov_dpc_haplo_noped <- Slov_allchroms_phasedNoPed[rownames(Slov_allchroms_phasedNoPed) %in% Slov_pedigree$dpc,]
Slov_workers_haplo_noped <- Slov_allchroms_phasedNoPed[rownames(Slov_allchroms_phasedNoPed) %in% Slov_pedigree$id,]

Slov_all_haplo_phasedNoPed <-  rbind(Slov_mother_haplo_noped, Slov_dpc_haplo_noped, Slov_workers_haplo_noped)


Slov_haplotypes_phasedNoPed <- get_out_haplotypes(ped_matrix = Slov_all_haplo_phasedNoPed, ind_id_1 = Slov_ind_id_1, ind_id_2 = Slov_ind_id_2)
Slov_phasedhaplotypes_phasedNoPed <- apply(Slov_haplotypes_phasedNoPed, 2, convert_genotypes)
rownames(Slov_phasedhaplotypes_phasedNoPed) <- rownames(Slov_haplotypes_phasedNoPed)
}
}





#######################################################################################################################
#########   Step 6: Assign the Haplotype parent-of-origin assignments ################################################
#######################################################################################################################

#In this section I'll be using 2 routes:

#.  - Route 1: use pedigree reconstruction and use both dam and dpc pedigree ID to assign haplotype parental origins
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
#Mother with genotypes     [0,1,2,0]
#DPQ with genotypes        [1,0,1,1]

# Using Table above, we calculate mismatches for the mother:
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

#ROUTE 1 of Haplotype parental origin assignments : Using the pedigree
#Slov real data only has 117 offspring individuals assigned for this route since only 117 were identified for the pedigree using pedigree reconstruction
{
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

Route1_flipping <- function(Data_type = NULL, pedigree = NULL, perfect_haplotypes = NULL, method = NULL){
  
  if (perfect_haplotypes == TRUE){
    haplotypes <- true_haplotypes
    map <- true_map
    
    real_results_methods <- list()
    real_results_lengths <- list()
    flipped_haplotypes <- list()
    
    
    for (i in 1:nrow(pedigree)) {
      offspring_id <- pedigree$id[i]
      mother_id <- pedigree$mother[i]
      dpc_id <- pedigree$dpc[i]
      
      offspring_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x,'_')[[1]][1]) %in% offspring_id]
      offspring_haplotypes_i <- haplotypes[offspring_row,]
      Offspring_Hap1 <- t(as.data.frame(offspring_haplotypes_i[1,]))
      Offspring_Hap2 <- t(as.data.frame(offspring_haplotypes_i[2,]))
      
      mother_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x,'_')[[1]][1]) %in% mother_id]
      mother_haplotypes_i <- haplotypes[mother_row,]
      Maternal_Hap1 <- t(as.data.frame(mother_haplotypes_i[1,]))
      rownames(Maternal_Hap1) <- mother_id
      Maternal_Hap2 <- t(as.data.frame(mother_haplotypes_i[2,]))
      rownames(Maternal_Hap2) <- mother_id
      Maternal_Geno <- matrix(data = Maternal_Hap1 + Maternal_Hap2, nrow = 1)
      colnames(Maternal_Geno) <- colnames(Maternal_Hap1)
      rownames(Maternal_Geno) <- mother_id
      
      dpc_row <- rownames(true_haplotypes)[sapply(rownames(true_haplotypes), function(x) strsplit(x,'_')[[1]][1]) %in% dpc_id]
      dpc_haplotypes_i <- true_haplotypes[dpc_row,]
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
  
  if (perfect_haplotypes == FALSE){
    
    phased_results_methods <- list()
    phased_results_lengths <- list()
    flipped_haplotypes <- list()
    
    if (Data_type == "Phased_nGE"){
      haplotypes <- nGE_phasedhaplotypes_phasedwithPed
      map <- nGE_map
    } else if (Data_type == "Phased_GE"){
      haplotypes <- GE_phasedhaplotypes_phasedwithPed
      map <- GE_map
    } else if (Data_type == "Real_Slov_data"){
      haplotypes <- Slov_phasedhaplotypes_phasedwithPed
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

Route1_True_FLIP <- Route1_flipping(perfect_haplotypes = TRUE, pedigree = simulated_pedigree, method = "power_mean")
Route1_GE_phased_FLIP <- Route1_flipping(perfect_haplotypes = FALSE, Data_type = "Phased_GE", pedigree = simulated_pedigree, method = "power_mean")
Route1_nGE_phased_FLIP <- Route1_flipping(perfect_haplotypes = FALSE, Data_type = "Phased_nGE", pedigree = simulated_pedigree, method = "power_mean")

#Since it uses the pedigree reconstruction info, you cant have a pedigree with a DPC that is unidentified. 
Slov_pedigree_withDPC <- Slov_pedigree[Slov_pedigree$dpc != "0",]
Route1_Slov_FLIP <- Route1_flipping(perfect_haplotypes = FALSE, Data_type = "Real_Slov_data", pedigree = Slov_pedigree_withDPC, method = "power_mean")
}

#ROUTE2 of Haplotype parental origin assignments : ONLY maternal information used (no pedigree from pedigree reconstruction)
#Slov has 235 individuals assigned this time since we know the maternal information and sire information isn't required for this method
{
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
        mother_id <- pedigree$mother[i]
        
        offspring_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x,'_')[[1]][1]) %in% offspring_id]
        offspring_haplotypes_i <- haplotypes[offspring_row,]
        Offspring_Hap1 <- t(as.data.frame(offspring_haplotypes_i[1,]))
        Offspring_Hap2 <- t(as.data.frame(offspring_haplotypes_i[2,]))
        
        mother_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x,'_')[[1]][1]) %in% mother_id]
        mother_haplotypes_i <- haplotypes[mother_row,]
        Maternal_Hap1 <- t(as.data.frame(mother_haplotypes_i[1,]))
        rownames(Maternal_Hap1) <- mother_id
        Maternal_Hap2 <- t(as.data.frame(mother_haplotypes_i[2,]))
        rownames(Maternal_Hap2) <- mother_id
        Maternal_Geno <- matrix(data = Maternal_Hap1 + Maternal_Hap2, nrow = 1)
        colnames(Maternal_Geno) <- colnames(Maternal_Hap1)
        rownames(Maternal_Geno) <- mother_id
        
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
          identified_dpc_haplo <- grep("_paternal$", rownames(merged_haps), value = TRUE)
          identified_maternal_haplo <- grep("_maternal$", rownames(merged_haps), value = TRUE)
          
          dpc_assigned <- merged_haps[identified_dpc_haplo, , drop = FALSE]
          maternal_assigned <- merged_haps[identified_maternal_haplo, , drop = FALSE]
          
          pat_results[[j]] <- as.data.frame(dpc_assigned)
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
      
      if (Data_type == "Phased_nGE") {
        haplotypes <- nGE_phasedhaplotypes_phasedNoPed
        map <- nGE_map
      } else if (Data_type == "Phased_GE" ) {
        haplotypes <- GE_phasedhaplotypes_phasedNoPed
        map <- GE_map
      } else if (Data_type == "Slov_real_data"){
        haplotypes <- Slov_phasedhaplotypes_phasedNoPed
        map <- Slov_map
      }
      
      for (i in 1:nrow(pedigree)) {
        offspring_id <- pedigree$id[i]
        mother_id <- pedigree$mother[i]
        
        offspring_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x, '_')[[1]][1]) %in% offspring_id]
        offspring_haplotypes_i <- haplotypes[offspring_row, ]
        Offspring_Hap1 <- as.data.frame(offspring_haplotypes_i[1, , drop = FALSE])
        Offspring_Hap2 <- as.data.frame(offspring_haplotypes_i[2, , drop = FALSE])
        
        mother_row <- rownames(haplotypes)[sapply(rownames(haplotypes), function(x) strsplit(x, '_')[[1]][1]) %in% mother_id]
        mother_haplotypes_i <- haplotypes[mother_row, ]
        Maternal_Hap1 <- as.data.frame(mother_haplotypes_i[1, , drop = FALSE])
        Maternal_Hap2 <- as.data.frame(mother_haplotypes_i[2, , drop = FALSE])
        Maternal_Geno <- matrix(data = Maternal_Hap1 + Maternal_Hap2, nrow = 1)
        colnames(Maternal_Geno) <- colnames(Maternal_Hap1)
        rownames(Maternal_Geno) <- mother_id  
        
        
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
          identified_dpc_haplo <- grep("_paternal$", rownames(merged_haps), value = TRUE)
          identified_maternal_haplo <- grep("_maternal$", rownames(merged_haps), value = TRUE)
          
          dpc_assigned <- merged_haps[identified_dpc_haplo, , drop = FALSE]
          maternal_assigned <- merged_haps[identified_maternal_haplo, , drop = FALSE]
          
          pat_results[[j]] <- as.data.frame(dpc_assigned)
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
  
  Route2_True_FLIP <- Route2_flipping(perfect_haplotypes = TRUE, Data_type = "True", pedigree = simulated_pedigree, method = "power_mean")
  Route2_GE_phased_FLIP <- Route2_flipping(perfect_haplotypes = FALSE, Data_type = "Phased_GE", pedigree = simulated_pedigree, method = "power_mean")
  Route2_nGE_phased_FLIP <- Route2_flipping(perfect_haplotypes = FALSE, Data_type = "Phased_nGE", pedigree = simulated_pedigree, method = "power_mean")
  
  Route2_Slov_FLIP <- Route2_flipping(perfect_haplotypes = FALSE, Data_type = "Slov_real_data", pedigree = Slov_pedigree, method = "power_mean")
 
}






#######################################################################################################################
#########   Step 7: Gametic Mendelian Sampling Value ################################################
#######################################################################################################################
# To assess the accuracy of haplotype parent-of-origin assignments, we calculated the 
# gametic Mendelian sampling values using both simulated and real data.
#
# According to Wright's pedigree model (Wright, 1921), the genetic value of an individual (g_i)
# is the sum of the average of their parental genetic values and a Mendelian sampling term (r_i):
#
#   g_i = 0.5 * g_m(i) + 0.5 * g_p(i) + r_i
#
# Here, g_m(i) and g_p(i) represent the genetic values of the mother and father, respectively.
# The term r_i captures the deviation due to Mendelian sampling (random inheritance).
#
# For our analysis, we extended this model to the haplotype level to evaluate deviations more precisely.
# We separated the maternal and paternal contributions as follows:
#
#   g_i_maternal = 0.5 * g_m(i)_1 + 0.5 * g_m(i)_2 + r_i_maternal
#   g_i_paternal = 0.5 * g_p(i)_1 + 0.5 * g_p(i)_2 + r_i_paternal
#
# This allows us to isolate the Mendelian sampling deviations (r_i_maternal and r_i_paternal)
# for each parent, using phased haplotype information.

#For haplotype p-o ROUTE1 (pedigree used)
{
  calculate_gametic_relatedness_R1 <- function(sorted_offspring_haplotypes, all_haplotypes, pedigree) {
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
  
  colnames(Route1_True_FLIP$real_results_flipped) <- colnames(true_haplotypes)
  Route1_Gametic_real <- calculate_gametic_relatedness_R1(sorted_offspring_haplotypes = Route1_True_FLIP$real_results_flipped , all_haplotypes = true_haplotypes, pedigree = simulated_pedigree)
  Route1_Gametic_real$Phasing <- "True"
  
  Route1_Gametic_nGE_phased <- calculate_gametic_relatedness_R1(sorted_offspring_haplotypes = Route1_nGE_phased_FLIP$results_flipped, all_haplotypes = nGE_phasedhaplotypes_phasedwithPed, pedigree = simulated_pedigree)
  Route1_Gametic_nGE_phased$Phasing <- "Phased_nGE"
  
  Route1_Gametic_GE_phased <- calculate_gametic_relatedness_R1(sorted_offspring_haplotypes = Route1_GE_phased_FLIP$results_flipped, all_haplotypes = GE_phasedhaplotypes_phasedwithPed, pedigree = simulated_pedigree)
  Route1_Gametic_GE_phased$Phasing <- "Phased_GE"
  
  Route1_Gametic_Slov <- calculate_gametic_relatedness_R1(sorted_offspring_haplotypes = Route1_Slov_FLIP$results_flipped, all_haplotypes = Slov_phasedhaplotypes_phasedwithPed, pedigree = Slov_pedigree_withDPC)
  Route1_Gametic_Slov$Phasing <- "Real"
  
  Gametic_df <- rbind(Route1_Gametic_real, Route1_Gametic_nGE_phased,  Route1_Gametic_GE_phased)
  
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
  Plotting_Gametic(df = Route1_Gametic_Slov, phased_type = c("Real"), plotting_styles = c("density"))
}

#For haplotype p-o ROUTE2 (mother info only)
{
  calculate_gametic_relatedness_R2 <- function(sorted_offspring_haplotypes, all_haplotypes, pedigree) {
    # Initialize an empty data frame to store results
    gametic_results <- data.frame(
      Offspring_ID = character(),
      Mother_ID = character(),
      Maternal_Mendelian = numeric(),
      stringsAsFactors = FALSE
    )
    
    cols <- colnames(sorted_offspring_haplotypes)
    all_haplotypes_cols <- all_haplotypes[, cols, drop = FALSE]
    
    # Loop through each offspring in the pedigree file
    for (i in 1:nrow(pedigree)) {
      offspring_id <- pedigree$id[i]
      mother_id <- pedigree$mother[i]
      
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
      
      # Calculate parental average and Mendelian sampling for each haplotype
      ma_hap <- 0.5 * (Maternal_Hap1 + Maternal_Hap2)
      ri_hap1 <- offspring_maternal - ma_hap
      ri_hap1_sum <- sum(ri_hap1)
      
      # Add the values to the results data frame
      gametic_results <- rbind(gametic_results, data.frame(
        Offspring_ID = offspring_id,
        Mother_ID = mother_id,
        Maternal_Mendelian = ri_hap1_sum,
        stringsAsFactors = FALSE
      ))
    }
    
    return(gametic_results)
  }
  
  
  colnames(Route2_True_FLIP$real_results_flipped) <- colnames(true_haplotypes)
  Route2_Gametic_TRUE <- calculate_gametic_relatedness_R2(sorted_offspring_haplotypes = Route2_True_FLIP$real_results_flipped , all_haplotypes = true_haplotypes, pedigree = simulated_pedigree)
  Route2_Gametic_TRUE$Phasing <- "True"
  
  Route2_Gametic_nGE_phased <- calculate_gametic_relatedness_R2(sorted_offspring_haplotypes = Route2_nGE_phased_FLIP$phased_results_flipped, all_haplotypes = nGE_phasedhaplotypes_phasedNoPed, pedigree = simulated_pedigree)
  Route2_Gametic_nGE_phased$Phasing <- "Phased_nGE"
  
  Route2_Gametic_GE_phased <- calculate_gametic_relatedness_R2(sorted_offspring_haplotypes = Route2_GE_phased_FLIP$phased_results_flipped, all_haplotypes = GE_phasedhaplotypes_phasedNoPed, pedigree = simulated_pedigree)
  Route2_Gametic_GE_phased$Phasing <- "Phased_GE"
  
  Route2_Gametic_Slov <- calculate_gametic_relatedness_R2(sorted_offspring_haplotypes = Route2_Slov_FLIP$phased_results_flipped, all_haplotypes = Slov_phasedhaplotypes_phasedNoPed, pedigree = Slov_pedigree)
  Route2_Gametic_Slov$Phasing <- "Real"
  
  Gametic_df_R2 <- rbind(Route2_Gametic_TRUE, Route2_Gametic_nGE_phased,  Route2_Gametic_GE_phased)
  
  Plotting_Gametic_R2 <- function(df, phased_type, plotting_styles) {
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
                                      values = c("Maternal" = "blue"),
                                      labels = c("Maternal"))
    fill_scale <- scale_fill_manual(name = "Assigned parent haplotype",
                                    values = c("Maternal" = "blue"),
                                    labels = c("Maternal"))
    
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
      # Create a plot to extract the legend from
      legend_plot <- ggplot(df[df$Phasing == "Phased_GE", ], aes(x = Maternal_Mendelian, fill = "Maternal")) +
        geom_histogram(aes(y = ..count..), binwidth = 1, alpha = 0.7, position = 'identity') +
        fill_scale +  # Apply fill scale
        theme_minimal() +
        theme(legend.position = "right",
              legend.title = element_text(size = 16),   # Increase legend title font size
              legend.text = element_text(size = 14))   # Increase legend text font size
      
      # Extract legend using get_legend from cowplot
      legend <- cowplot::get_legend(legend_plot)
      
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
  
  
  Plotting_Gametic_R2(df = Gametic_df_R2, phased_type = c("True", "Phased_nGE", "Phased_GE"), plotting_styles = c("density"))
  Plotting_Gametic_R2(df = Route2_Gametic_Slov, phased_type = c("Real"), plotting_styles = c("density"))
  
}







#######################################################################################################################
#########   Step 8: Patriline determination ################################################
#######################################################################################################################
# Estimates the number of patrilines by comparing the paternally assigned haplotypes workers

#Again I'm taking two routes here: 

# - ROUTE 1: Using actual drone information using the true haplotypes to test the accuracy of the threshold method
#            - checks the paternity accuracy really closely but is unrealistic unless you have father-drone information. 
#            - uses sister thresholds and father_thresholds to visualise the most accurate thresholds required

{

# ------------------------------------------------------------
# KEY STEPS IN `calc_nPaternity_Accuracy()`:
# ------------------------------------------------------------
# 1. Extract paternal haplotypes for all workers.
# 2. Loop through each queen (colony mother) to group her offspring.
# 3. For each queen:
#     - Compare worker haplotypes to true father haplotypes using a similarity threshold.
#     - Build a similarity matrix between all workers and fathers.
#     - Assign workers to father groups if similarity exceeds threshold.
#     - For unassigned workers, check if they resemble any assigned sisters.
#     - Still unmatched workers are compared with each other to form new groups.
#     - Remaining unmatched individuals are placed in their own group.
# 4. Count and store:
#     - Actual number of fathers
#     - Estimated number of fathers from grouping
#     - Number of correctly identified fathers
#
# ------------------------------------------------------------
# `run_paternity_tests()`:
# ------------------------------------------------------------
# - Accepts a matrix of paternal haplotypes (`results_arg`), a list of `sister_thresholds`, 
#   and a list of `father_test_thresholds`.
# - For each combination of thresholds, runs `calc_nPaternity_Accuracy()`.
# - Stores results in a named list (`results_list`) for downstream analysis.


#Load your Pop object that contains the fathers and the pedigree with the fathers id in there 
father_id_all <- simulated_pedigree$father
father_id_1 <- paste(father_id_all, "1", sep = "_")
all_perfect_haplotypes <- pullSnpHaplo(PopMerged, simParam = SP)

#since they're coded as identical diploids we can just use one of them
father_haplotypes_1 <- all_perfect_haplotypes[rownames(all_perfect_haplotypes) %in% father_id_1,]

#replace _1 with _paternal so that it can all be compared 
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
sister_thresholds <- c(1.0, 0.95, 0.90, 0.85, 0.80, 0.75)
father_test_thresholds <- c(0.95, 1)

#Run on the real/nGE/GE for haplotype p-o route 1 (don't think using either route will change the outcome)
{
  # Run the function and specify the results dataset to use
  #Route1 
  r1_Real_FatherTest_Results <- run_paternity_tests(
    results_arg = Route1_True_FLIP$real_results_flipped, 
    sister_thresholds = sister_thresholds,
    father_test_thresholds = father_test_thresholds,
    father_haplotypes = father_haplotype_pat
  )

}


plot_paternity_number_grid_DRONES <- function(Results) {
  
  # Combine the list of results into one data frame
  df <- dplyr::bind_rows(Results, .id = "threshold_combination")
  
  # Extract thresholds from the list names if not already in the data
  df <- df %>%
    tidyr::separate(threshold_combination, into = c("sis", "fat"), sep = "_") %>%
    dplyr::mutate(
      sister_threshold = as.numeric(gsub("sis", "", sis)),
      father_accuracy_threshold = as.numeric(gsub("fat", "", fat))
    ) %>%
    dplyr::select(-sis, -fat)
  
  # Pivot to long format, excluding num_workers
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = c(num_fathers_actual, num_fathers_estimated),
      names_to = "measure_type",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      measure_type = dplyr::recode(measure_type,
                                   num_fathers_actual = "Actual nPatrilines",
                                   num_fathers_estimated = "Estimated nPatrilines"),
      measure_type = factor(measure_type, levels = c("Actual nPatrilines", "Estimated nPatrilines"))
    )
  
  # Custom labeller function for facets
  custom_labeller <- function(variable, value) {
    if (variable == "sister_threshold") {
      return(paste0("sis_threshold: ", value))
    } else if (variable == "father_accuracy_threshold") {
      return(paste0("drone_threshold: ", value))
    } else {
      return(as.character(value))
    }
  }
  
  # Plot
  plot <- ggplot(df_long, aes(x = as.factor(queen_id), y = count, color = measure_type, shape = measure_type)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
    facet_grid(
      rows = vars(father_accuracy_threshold),
      cols = vars(sister_threshold),
      labeller = labeller(
        sister_threshold = function(x) paste0("sis_threshold: ", x),
        father_accuracy_threshold = function(x) paste0("drone_threshold: ", x)
      )
    ) +
    scale_color_manual(name = "Measure Type",
                       values = c("Actual nPatrilines" = "red", "Estimated nPatrilines" = "blue")) +
    scale_shape_manual(name = "Measure Type",
                       values = c("Actual nPatrilines" = 16, "Estimated nPatrilines" = 17)) +
    labs(x = "Colony ID", y = "Number of Patrilines", title = "Estimated vs Actual number of Patrilines") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.text.x = element_text(size = 12),
      strip.text.y = element_text(size = 12),
      strip.placement = "outside",
      strip.background = element_rect(fill = NA),
      plot.title = element_text(size = 18, face = "bold")
    ) +
    labs(
      subtitle = NULL,
      caption = NULL
    ) 
  
  return(plot)
}

plot_paternity_number_grid_DRONES(r1_Real_FatherTest_Results)

}

# - ROUTE 2: Use the paternally assigned haplotypes from this script and ONLY sister thresholds to determine patriline numbers
#           - again using multiple sister thresholds to determine which of the Slov real data results are the most reliable. 
{
{SisterONLY_calc_nPaternity_Accuracy <- function(results, sister_threshold, pedigree, simulated = NULL) {
  
  # Extract paternal haplotypes (still named '_paternal' in the results)
  results_paternal <- rownames(results)[grep(paste0("_paternal"), rownames(results))]
  results_paternal <- results[results_paternal, , drop = FALSE]
  
  # Create a loop grouping together sisters using pedigree 
  queen_ids <- unique(pedigree$mother)
  
  # Initialize a data frame to store the results
  if(simulated == TRUE){
    nPaternity <- data.frame(queen_id = queen_ids, 
                             num_workers = rep(0, length(queen_ids)),
                             num_sister_groups_estimated = rep(0, length(queen_ids)),
                             stringsAsFactors = FALSE,
                             sister_threshold = sister_threshold,
                             actual_number_fathers = rep(0, length(queen_ids)))
  }
  else if(simulated == FALSE){
    nPaternity <- data.frame(queen_id = queen_ids, 
                             num_workers = rep(0, length(queen_ids)),
                             num_sister_groups_estimated = rep(0, length(queen_ids)),
                             stringsAsFactors = FALSE,
                             sister_threshold = sister_threshold)
  }
  
  for (i in 1:length(queen_ids)) {
    print(i)
    queen_id <- queen_ids[i]
    sister_worker_ids <- t(pedigree[pedigree$mother == queen_id, "id"])
    
    if(simulated == TRUE){
      queen_group <- pedigree[pedigree$mother == queen_id,]
      father_ids <- t(unique(pedigree[pedigree$mother == queen_id, "father"]))
      nPaternity$actual_number_fathers[i] <- length(father_ids)
    }
    
    # Pull out the results_paternal with the worker IDs
    pattern <- paste0("^(", paste(sister_worker_ids, collapse = "|"), ")_paternal$")
    matching_indices <- grep(pattern, rownames(results_paternal))
    sister_paternal <- results_paternal[matching_indices, , drop = FALSE]
    
    # Store the number of workers for this queen
    nPaternity$num_workers[i] <- nrow(sister_paternal)
    
    # Ensure we have data to work with
    if (nrow(sister_paternal) == 0) {
      next
    }
    
    # Initialize a list to hold sister groups (clusters)
    sister_groups <- list()
    still_unmatched <- rownames(sister_paternal)
    
    # Define the similarity function to compare workers
    calculate_similarity <- function(hap1, hap2) {
      valid_indices <- !is.na(hap1) & !is.na(hap2)
      
      if (sum(valid_indices) == 0) {
        return(NA)  # Return NA if no valid comparisons can be made
      }
      
      sum(hap1[valid_indices] == hap2[valid_indices]) / length(hap1[valid_indices])
    }
    
    # Step 1: Group workers based on similarity (sister threshold)
    for (a in 1:(length(still_unmatched) - 1)) {
      for (b in (a + 1):length(still_unmatched)) {
        worker_a <- still_unmatched[a]
        worker_b <- still_unmatched[b]
        
        similarity_score <- calculate_similarity(sister_paternal[worker_a, ], sister_paternal[worker_b, ])
        
        # If they meet the sister threshold, group them together
        if (!is.na(similarity_score) && similarity_score >= sister_threshold) {
          group_found <- FALSE
          # Check if worker_a or worker_b already belongs to a group
          for (group_name in names(sister_groups)) {
            if (worker_a %in% sister_groups[[group_name]] || worker_b %in% sister_groups[[group_name]]) {
              sister_groups[[group_name]] <- unique(c(sister_groups[[group_name]], worker_a, worker_b))
              group_found <- TRUE
              break
            }
          }
          # If no group was found, create a new group for them
          if (!group_found) {
            new_group_name <- paste0("group_", length(sister_groups) + 1)
            sister_groups[[new_group_name]] <- c(worker_a, worker_b)
          }
        }
      }
    }
    
    # Step 2: Ensure all remaining workers (if any) are assigned to their own groups
    unmatched_workers <- setdiff(still_unmatched, unlist(sister_groups))
    for (worker in unmatched_workers) {
      new_group_name <- paste0("group_", length(sister_groups) + 1)
      sister_groups[[new_group_name]] <- c(worker)
    }
    
    # Store the number of estimated sister groups (clusters)
    nPaternity$num_sister_groups_estimated[i] <- length(sister_groups)
  }
  
  return(nPaternity)
}
  
  run_sister_clustering_tests <- function(results_arg, sister_thresholds, pedigree, simulated = NULL) {
    
    # Initialize a list to store the results
    results_list <- list()
    
    # Loop over each sister_threshold
    for (sister_threshold in sister_thresholds) {
      
      # Create a descriptive name for the result
      result_name <- paste0("sis", sister_threshold)
      
      # Print to ensure the threshold is being passed correctly
      print(paste("Running for sister_threshold:", sister_threshold))
      
      # Run the function for each sister_threshold
      result <- SisterONLY_calc_nPaternity_Accuracy(
        results = results_arg,      # Pass the dynamic results argument here
        sister_threshold = sister_threshold,  # sister threshold
        pedigree = pedigree,# Pass the pedigree information
        simulated = simulated
      )
      
      # Store the result with the unique name
      results_list[[result_name]] <- result
      
      # Print to confirm each result has been stored
      print(paste("Finished:", result_name))
    }
    
    # Return the list of results
    return(results_list)
  }
  
  {
    #results from haplotype p-o Route1 
    pr2_True_FatherTest_Results_hapr1 <- run_sister_clustering_tests(
      results_arg = Route1_True_FLIP$real_results_flipped, 
      sister_thresholds = sister_thresholds,
      pedigree = simulated_pedigree,
      simulated = TRUE
    )
    pr2_nGE_Father_Test_Results_hapr1 <- run_sister_clustering_tests(
      results_arg = Route1_nGE_phased_FLIP$results_flipped, 
      sister_thresholds = sister_thresholds,
      pedigree = simulated_pedigree,
      simulated = TRUE
    )
    pr2_GE_Father_Test_Results_hapr1 <- run_sister_clustering_tests(
      results_arg = Route1_GE_phased_FLIP$results_flipped,  
      sister_thresholds = sister_thresholds,
      pedigree = simulated_pedigree,
      simulated = TRUE
    )
    
    
    
    #results for haplotype p-o. Route2 
    
    pr2_True_FatherTest_Results_hapr2 <- run_sister_clustering_tests(
      results_arg = Route2_True_FLIP$real_results_flipped, 
      sister_thresholds = sister_thresholds,
      pedigree = simulated_pedigree,
      simulated = TRUE
    )
    
    pr2_nGE_Father_Test_Results_hapr2 <- run_sister_clustering_tests(
      results_arg = Route2_nGE_phased_FLIP$phased_results_flipped, 
      sister_thresholds = sister_thresholds,
      pedigree = simulated_pedigree,
      simulated = TRUE
    )
    
    pr2_GE_Father_Test_Results_hapr2 <- run_sister_clustering_tests(
      results_arg = Route2_GE_phased_FLIP$phased_results_flipped,  
      sister_thresholds = sister_thresholds,
      pedigree = simulated_pedigree,
      simulated = TRUE
    )
    
  }
  
  
  
  
  Slov_SisterONLY_calc_nPaternity_Accuracy <- function(results, sister_threshold, pedigree) {
    
    # Extract paternal haplotypes (still named '_paternal' in the results)
    results_paternal <- rownames(results)[grep(paste0("_paternal"), rownames(results))]
    results_paternal <- results[results_paternal, , drop = FALSE]
    
    # Create a loop grouping together sisters using pedigree 
    queen_ids <- unique(pedigree$mother)
    
    # Initialize a data frame to store the results
    
    nPaternity <- data.frame(queen_id = queen_ids, 
                             num_workers = rep(0, length(queen_ids)),
                             num_sister_groups_estimated = rep(0, length(queen_ids)),
                             stringsAsFactors = FALSE,
                             sister_threshold = sister_threshold)
    
    for (i in 1:length(queen_ids)) {
      print(i)
      queen_id <- queen_ids[i]
      sister_worker_ids <- t(pedigree[pedigree$mother == queen_id, "id"])
      
      # Pull out the results_paternal with the worker IDs
      pattern <- paste0("^(", paste(sister_worker_ids, collapse = "|"), ")_paternal$")
      matching_indices <- grep(pattern, rownames(results_paternal))
      sister_paternal <- results_paternal[matching_indices, , drop = FALSE]
      
      # Store the number of workers for this queen
      nPaternity$num_workers[i] <- nrow(sister_paternal)
      
      # Ensure we have data to work with
      if (nrow(sister_paternal) == 0) {
        next
      }
      
      # Initialize a list to hold sister groups (clusters)
      sister_groups <- list()
      still_unmatched <- rownames(sister_paternal)
      
      # Define the similarity function to compare workers
      calculate_similarity <- function(hap1, hap2) {
        valid_indices <- !is.na(hap1) & !is.na(hap2)
        
        if (sum(valid_indices) == 0) {
          return(NA)  # Return NA if no valid comparisons can be made
        }
        
        sum(hap1[valid_indices] == hap2[valid_indices]) / length(hap1[valid_indices])
      }
      
      # Step 1: Group workers based on similarity (sister threshold)
      for (a in 1:(length(still_unmatched) - 1)) {
        for (b in (a + 1):length(still_unmatched)) {
          worker_a <- still_unmatched[a]
          worker_b <- still_unmatched[b]
          
          similarity_score <- calculate_similarity(sister_paternal[worker_a, ], sister_paternal[worker_b, ])
          
          # If they meet the sister threshold, group them together
          if (!is.na(similarity_score) && similarity_score >= sister_threshold) {
            group_found <- FALSE
            # Check if worker_a or worker_b already belongs to a group
            for (group_name in names(sister_groups)) {
              if (worker_a %in% sister_groups[[group_name]] || worker_b %in% sister_groups[[group_name]]) {
                sister_groups[[group_name]] <- unique(c(sister_groups[[group_name]], worker_a, worker_b))
                group_found <- TRUE
                break
              }
            }
            # If no group was found, create a new group for them
            if (!group_found) {
              new_group_name <- paste0("group_", length(sister_groups) + 1)
              sister_groups[[new_group_name]] <- c(worker_a, worker_b)
            }
          }
        }
      }
      
      # Step 2: Ensure all remaining workers (if any) are assigned to their own groups
      unmatched_workers <- setdiff(still_unmatched, unlist(sister_groups))
      for (worker in unmatched_workers) {
        new_group_name <- paste0("group_", length(sister_groups) + 1)
        sister_groups[[new_group_name]] <- c(worker)
      }
      
      # Store the number of estimated sister groups (clusters)
      nPaternity$num_sister_groups_estimated[i] <- length(sister_groups)
    }
    
    return(nPaternity)
  }
  
  Slov_run_sister_clustering_tests <- function(results_arg, sister_thresholds, pedigree) {
    
    # Initialize a list to store the results
    results_list <- list()
    
    # Loop over each sister_threshold
    for (sister_threshold in sister_thresholds) {
      
      # Create a descriptive name for the result
      result_name <- paste0("sis", sister_threshold)
      
      # Print to ensure the threshold is being passed correctly
      print(paste("Running for sister_threshold:", sister_threshold))
      
      # Run the function for each sister_threshold
      result <- Slov_SisterONLY_calc_nPaternity_Accuracy(
        results = results_arg,      # Pass the dynamic results argument here
        sister_threshold = sister_threshold,  # sister threshold
        pedigree = pedigree# Pass the pedigree information
      )
      
      # Store the result with the unique name
      results_list[[result_name]] <- result
      
      # Print to confirm each result has been stored
      print(paste("Finished:", result_name))
    }
    
    # Return the list of results
    return(results_list)
  }
  
  pr2_Slov_FatherTest_results_hapr1 <- Slov_run_sister_clustering_tests(
    results_arg = Route1_Slov_FLIP$results_flipped,  
    sister_thresholds = sister_thresholds,
    pedigree = Slov_pedigree
  )
  
  pr2_Slov_FatherTest_results_hapr2 <- Slov_run_sister_clustering_tests(
    results_arg = Route2_Slov_FLIP$phased_results_flipped,  
    sister_thresholds = sister_thresholds,
    pedigree = Slov_pedigree
  )}

#MAY25th25 Plot: change the function names(results_sets) for each datasets
plot_paternity_number_grid_SISTERONLY <- function(...) {
  result_sets <- list(...)
  
  # Assign meaningful dataset labels
  names(result_sets) <- c("True_Route2", "nGE_Route2", "GE_Route2")
  
  # Combine all into one dataframe
  df <- purrr::imap_dfr(result_sets, function(results_list, dataset_name) {
    dplyr::bind_rows(results_list, .id = "threshold_combination") %>%
      mutate(dataset = dataset_name)
  })
  
  # Process threshold
  df <- df %>%
    tidyr::separate(threshold_combination, into = c("sis"), sep = "_", fill = "right", extra = "drop") %>%
    dplyr::mutate(
      sister_threshold = as.numeric(gsub("sis", "", sis))
    ) %>%
    dplyr::mutate(
      sister_threshold = factor(sister_threshold, levels = sort(unique(sister_threshold), decreasing = TRUE))
    ) %>%
    dplyr::select(-sis)
  
  # Long format
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = c(actual_number_fathers, num_sister_groups_estimated),
      names_to = "measure_type",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      measure_type = dplyr::recode(measure_type,
                                   actual_number_fathers = "Actual nPatrilines",
                                   num_sister_groups_estimated = "Estimated nPatrilines"),
      measure_type = factor(measure_type, levels = c("Actual nPatrilines", "Estimated nPatrilines"))
    )
  
  # Plot
  plot <- ggplot(df_long, aes(x = as.factor(queen_id), y = count, color = measure_type, shape = measure_type)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
    facet_grid(
      rows = vars(dataset),
      cols = vars(sister_threshold),
      labeller = labeller(
        sister_threshold = function(x) paste0("sis_threshold: ", x),
        dataset = function(x) paste0(x)
      )
    ) +
    scale_color_manual(name = "Measure Type",
                       values = c("Actual nPatrilines" = "red", "Estimated nPatrilines" = "blue")) +
    scale_shape_manual(name = "Measure Type",
                       values = c("Actual nPatrilines" = 16, "Estimated nPatrilines" = 17)) +
    labs(x = "Colony ID", y = "Number of Patrilines", title = "Estimated vs Actual number of Patrilines") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.text.x = element_text(size = 12),
      strip.text.y = element_text(size = 12),
      strip.placement = "outside",
      strip.background = element_rect(fill = NA),
      plot.title = element_text(size = 18, face = "bold")
    )
  
  return(plot)
}
plot_paternity_number_grid_SISTERONLY(
  Dataset1 = pr2_True_FatherTest_Results_hapr1,
  Dataset2 = pr2_nGE_Father_Test_Results_hapr1,
  Dataset3= pr2_GE_Father_Test_Results_hapr1)


plot_paternity_number_grid_SISTERONLY(
  Dataset1 = pr2_True_FatherTest_Results_hapr2,
  Dataset2 = pr2_nGE_Father_Test_Results_hapr2,
  Dataset3= pr2_GE_Father_Test_Results_hapr2)


Slov_plot_paternity_number_grid_SISTERONLY <- function(...) {
  result_sets <- list(...)
  
  # Assign meaningful dataset labels
  names(result_sets) <- c("Real_Route1", "Real_Route2")
  
  # Combine all into one dataframe
  df <- purrr::imap_dfr(result_sets, function(results_list, dataset_name) {
    dplyr::bind_rows(results_list, .id = "threshold_combination") %>%
      mutate(dataset = dataset_name)
  })
  
  # Process threshold
  df <- df %>%
    tidyr::separate(threshold_combination, into = c("sis"), sep = "_", fill = "right", extra = "drop") %>%
    dplyr::mutate(
      sister_threshold = as.numeric(gsub("sis", "", sis))
    ) %>%
    dplyr::mutate(
      sister_threshold = factor(sister_threshold, levels = sort(unique(sister_threshold), decreasing = TRUE))
    ) %>%
    dplyr::select(-sis)
  
  # Long format
  df_long <- df %>%
    dplyr::group_by(dataset, sister_threshold) %>%
    dplyr::mutate(queen_index = dplyr::row_number()) %>%  # <- add this line
    tidyr::pivot_longer(
      cols = c(num_workers, num_sister_groups_estimated),
      names_to = "measure_type",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      measure_type = dplyr::recode(measure_type,
                                   num_workers = "nWorkers",
                                   num_sister_groups_estimated = "Estimated nPatrilines"),
      measure_type = factor(measure_type, levels = c("nWorkers", "Estimated nPatrilines"))
    ) %>%
    dplyr::ungroup()
  
  
  # Plot
  plot <- ggplot(df_long, aes(x = queen_index, y = count, color = measure_type, shape = measure_type)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
    facet_grid(
      rows = vars(dataset),
      cols = vars(sister_threshold),
      labeller = labeller(
        sister_threshold = function(x) paste0("sis_threshold: ", x),
        dataset = function(x) paste0(x)
      )
    ) +
    scale_x_continuous(breaks = df_long$queen_index) +
    scale_color_manual(name = "Measure Type",
                       values = c("nWorkers" = "red", "Estimated nPatrilines" = "blue")) +
    scale_shape_manual(name = "Measure Type",
                       values = c("nWorkers" = 16, "Estimated nPatrilines" = 17)) +
    labs(x = "Colony ID", y = "Number of Patrilines", title = "Estimated number of Patrilines") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.text.x = element_text(size = 12),
      strip.text.y = element_text(size = 12),
      strip.placement = "outside",
      strip.background = element_rect(fill = NA),
      plot.title = element_text(size = 18, face = "bold")
    )
  
  return(plot)
}
Slov_plot_paternity_number_grid_SISTERONLY(
  Dataset1 = pr2_Slov_FatherTest_results_hapr1,
  Dataset2 = pr2_Slov_FatherTest_results_hapr2
)
}




