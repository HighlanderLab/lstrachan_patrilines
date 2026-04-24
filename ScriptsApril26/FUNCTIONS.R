#######################################################################################################################
#**                          FUNCTIONS                                       **
#######################################################################################################################
#Gathering all of the functions together and providing summaries as to what each function does

#######################################################################################################################
#**•• 1_Preparing_RealData_forPhasing Functions ••******************************
#######################################################################################################################
extract_parts <- function(x) {
  parts <- unlist(strsplit(x, "[._]"))
  chromosome <- paste(parts[2], paste(parts[3], parts[4], sep="."), sep="_")
  position <- as.numeric(parts[5])
  ref_allele <- parts[6]
  data.frame(Chromosome = chromosome, Position = position, RefAllele = ref_allele, stringsAsFactors = FALSE)
}
# Explanation of `extract_parts()`:
{
  # INPUT: 
  # - A character string (x) containing identifiers separated by "." or "_".
  #
  # WHAT'S HAPPENING:
  # - Splits the string into a vector of parts.
  # - Recombines specific segments to form a standardized Chromosome name.
  # - Extracts and converts the numerical Position and Reference Allele.
  #
  # OUTPUT:
  # - A data frame with three columns: Chromosome, Position, and RefAllele.
}

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
# Explanation of `AB_to_12()`:
{
  # INPUT: 
  # - ped_file: A PLINK-format pedigree file (A/B genotypes).
  # - map_ref_file: A reference file containing the "RefAllele" for each SNP.
  # - output_file_prefix: The string used to name the resulting .ped and .map files.
  #
  # WHAT'S HAPPENING:
  # 1. Data Cleaning: Loads files and standardizes alleles in the reference data 
  #    (converting "C" to "A" and "D" to "B") to ensure consistency.
  #
  # 2. Reference-Based Recoding: Iterates through SNP columns in pairs (alleles).
  #    - If RefAllele is "A": Recodes A -> 1 (reference) and B -> 2 (alternative).
  #    - If RefAllele is "B": Recodes B -> 1 (reference) and A -> 2 (alternative).
  #    - Purpose: Standardizes genotypes relative to a specific reference genome 
  #      rather than just arbitrary A/B labels.
  #
  # 3. Metadata Generation: Creates a new PLINK .map file using the chromosome 
  #    and position data from the reference file.
  #
  # OUTPUT:
  # - A recoded .ped file (genotypes as 1s and 2s).
  # - A corresponding .map file formatted for PLINK analysis.
}

filter_chromosome_nc <- function(df) {
  # Filter the dataframe to keep only rows where the Chromosome column starts with "NC"
  filtered_df <- df[grep("^NC", df$Chromosome), ]
  return(filtered_df)
}
# Explanation of `filter_chromosome_nc()`:
{
  # INPUT: 
  # - A data frame (df) containing a "Chromosome" column.
  #
  # WHAT'S HAPPENING:
  # - Only Chromosomes beginning with NC are the actual chromosomes
  # - Searches the "Chromosome" column for values starting with "NC".
  # - Subsets the data frame to retain only those matching rows.
  #
  # OUTPUT:
  # - A filtered data frame containing only "NC" (RefSeq) chromosome entries.
}

#######################################################################################################################
#**•• 2_InbreedBasePop Functions ••*********************************************
#######################################################################################################################
createArray= function(array_name, array_number, nChr, segSites, nSNPPerChr, pop){
snpArray = vector("list", nChr)
for (chr in 1:nChr) {
  #get the haplotypes for all sites within the pop
  x = pullSegSiteHaplo(pop, chr = chr)
  #pop allele frequency at all sites
  alleleFreq = apply(X = x, MARGIN = 2, FUN = mean)
  #Will pick SNP based on allele frequency distribution
  # tmp = runif_from_nonunif(x = data.frame(id = 1:segSites[[chr]], value = alleleFreq), n = nSNPPerChr[[chr]]) #run with runif_from_nonunif
  
  tmp_data <- data.frame(id = 1:length(alleleFreq), value = alleleFreq)
  #Select SNP at random
  tmp <- sample(nrow(tmp_data), size = nSNPPerChr[[chr]], replace = FALSE) # run without runif_from_nonunif
  #save selected SNP
  sel <- tmp_data[tmp, "id"] #without runif_from_nonunif
  #sel = tmp$id #with runif_from_nonunif
  #save SNP in array
  snpArray[[chr]] = sort(sel) # Must be sorted
}

snpArray = do.call("c", snpArray) # Collapse list to vector
snpArray = new(
  Class = "LociMap",
  nLoci = sum(as.integer(nSNPPerChr)),
  lociPerChr = as.integer(nSNPPerChr),
  lociLoc = snpArray,
  name = array_name
)
SP$snpChips[[array_number]] = snpArray

}
# Explanation of `createArray()`:
{
  # INPUT: 
  # - array_name: A string used to name the new SNP chip.
  # - array_number: The index/slot within the global simulation parameters (SP) 
  #     where the array will be stored.
  # - nChr: The total number of chromosomes in the population.
  # - segSites: A list containing the indices of all segregating sites for each chromosome.
  # - nSNPPerChr: A list or vector specifying how many SNPs should be selected 
  #     for the array from each chromosome.
  # - pop: The population object from which the genetic data is pulled.
  #
  # WHAT'S HAPPENING:
  # 1. Iterative SNP Selection: Loops through each chromosome to identify SNPs:
  #     - Haplotype Extraction: Retrieves raw haplotype data for all segregating 
  #       sites on the current chromosome.
  #     - Frequency Calculation: Calculates allele frequency for every site to 
  #       characterize available genetic variation.
  #     - Random Sampling: Randomly selects a specific number of SNPs (defined 
  #       by nSNPPerChr) without replacement.
  #     - Sorting: Numerically sorts selected SNP indices to meet downstream requirements.
  #
  # 2. Object Construction: Collapses the chromosome-specific lists into a single 
  #    vector and initializes a formal "LociMap" S4 object. This object stores 
  #    total SNP counts, distribution across chromosomes, and specific locations.
  #
  # 3. Global Integration: Assigns the newly created LociMap to the global 
  #    simulation object (SP$snpChips) at the specified array_number index.
  #
  # OUTPUT:
  # - This function does not return a value; it has the side effect of updating 
  #   the global SP object with a defined SNP chip for use in genotyping.
}

runif_from_nonunif <- function(x, n, n_bins = 100) {
  samples_min <- min(x$value)
  samples_max <- max(x$value)
  bin_size <- (samples_max - samples_min) / n_bins
  bin_seq <- seq(from = samples_min, to = samples_max, by = bin_size)
  x$bin <- cut(x = x$value, breaks = bin_seq)
  bin_freq <- as.data.frame(table(x$bin))
  colnames(bin_freq) <- c("bin", "freq")
  x <- merge(x = x, y = bin_freq)
  # Sample without replacement and up-weight low frequency values so that
  # once we sample these out, we can then move to more common & high frequency
  # values
  # TODO: maybe this should be done differently/better by sampling bins at
  # random and then randomly within a bin?
  sel <- sample.int(n = nrow(x), size = n, prob = 1 / x$freq, replace = FALSE)
  return(x[sel, c("id", "value")])
}
# Explanation of `runif_from_nonunif()`:
{
  # INPUT:
  # - x: A data frame containing 'id' and 'value' (e.g., allele frequencies).
  # - n: The number of items to sample.
  # - n_bins: The number of groups to divide the data into for frequency analysis.
  #
  # WHAT'S HAPPENING:
  # 1. Binning the Data: Calculates the range of values and divides them into 
  #    equally sized segments (bins). Each row in the input is then assigned 
  #    to a bin based on its "value".
  #
  # 2. Frequency Calculation: Determines how many data points fall into each 
  #    bin. This identifies which "values" are common and which are rare 
  #    within the dataset.
  #
  # 3. Probability Weighting: Assigns a sampling probability to each row that 
  #    is inversely proportional to its bin frequency (1 / freq). 
  #    - Purpose: This "up-weights" rare values and "down-weights" common ones.
  #
  # 4. Balanced Sampling: Uses the calculated probabilities to select 'n' 
  #    items. This forces the resulting sample to be more representative of 
  #    the entire range of values (uniform-like) rather than being biased 
  #    toward the most common values in the original distribution.
  #
  # OUTPUT:
  # - A data frame containing the 'id' and 'value' of the sampled items, 
  #   providing a subset that is distributed more evenly than the input.
}


#######################################################################################################################
#**•• 3_PreparingSimulatedfilesForAlphaAssign Functions ••*********************************************
#######################################################################################################################

generateGenoErr <- function(geno, error, error2, sampling_error, missing) {
  nLoci <- ncol(geno)
  nInd <- nrow(geno)
  #error for sampling error
  for (ind in 1:nInd){
    for (locus in 1:nLoci){
      if (geno[ind, locus] != 9 && rbinom(n = 1, size = 1, prob = sampling_error) == 1) {
        geno[ind, locus] = 9 }
    }}
  #error for genotyping error
  for (ind in 1:nInd){
    for (locus in 1:nLoci){
      if (geno[ind, locus] != 9 && rbinom(n = 1, size = 1, prob = error) == 1) {
        if (geno[ind, locus] == 0) {
          geno[ind, locus] <- sample(c(1, 2, 9), size = 1, prob = c(1 - error2 - missing, error2, missing))
        } else if (geno[ind, locus] == 1) {
          geno[ind, locus] <- sample(c(0, 2, 9), size = 1, prob = c(1 - error2 - missing, 1 - error2 - missing, missing))
        } else if (geno[ind, locus] == 2) {
          geno[ind, locus] <- sample(c(0, 1, 9), size = 1, prob = c(error2, 1 - error2 - missing, missing))
        }
      }
    }
  }
  return(geno)
}
# Explanation of `generateGenoErr()`:
{
  # INPUT:
  # - geno: A matrix of genotypes (typically 0, 1, 2 representing allele counts).
  # - error: The overall probability that a genotyping error occurs at a locus.
  # - error2: A specific probability parameter used to weight the direction of 
  #     the error (e.g., probability of a homozygous flip).
  # - sampling_error: The probability that a genotype is lost/missing (set to 9) 
  #     during an initial sampling pass.
  # - missing: The probability that an error results specifically in a missing 
  #     value (9) during the secondary genotyping pass.
  #
  # WHAT'S HAPPENING:
  # 1. Initial Sampling Error Pass: Iterates through every individual and locus. 
  #    Using a binomial distribution (rbinom), it identifies random sites to 
  #    be marked as missing (9) based on the sampling_error rate.
  #
  # 2. Genotyping Error Pass: Iterates again to simulate "technical" errors:
  #    - Check: Only applies to sites not already marked as missing.
  #    - Trigger: An error is triggered based on the 'error' probability.
  #    - Recoding Logic: If an error occurs, the original genotype is replaced 
  #      using a random sample of the other possible states:
  #        - From 0 (Homozygous Ref): Can become 1, 2, or 9.
  #        - From 1 (Heterozygous): Can become 0, 2, or 9.
  #        - From 2 (Homozygous Alt): Can become 0, 1, or 9.
  #    - Probability Weighting: Uses error2 and missing to determine how 
  #      likely the genotype is to "flip" to a specific alternative state.
  #
  # OUTPUT:
  # - A genotype matrix of the same dimensions as the input, but with 
  #   simulated noise, miscalls, and missing data points (9s).
}

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
# Explanation of `switchToNA()`:
{
  # INPUT:
  # - genoErr: A numeric genotype matrix.
  #
  # WHAT'S HAPPENING:
  # The function validates that the input is a matrix and then standardizes 
  # missing data representation. it identifies numeric placeholders commonly 
  # used in bioinformatics (9 and -9) and converts them into R's formal 
  # NA (Not Available) value to ensure compatibility with statistical functions.
  #
  # OUTPUT:
  # - A matrix where all instances of 9 and -9 are replaced by NA.
}

genotypes_to_haplotypes <- function(genotypes) {
  inferred_hap1 <- matrix(0, nrow = nrow(genotypes), ncol = ncol(genotypes))
  inferred_hap2 <- matrix(0, nrow = nrow(genotypes), ncol = ncol(genotypes))
  
  for (i in 1:nrow(genotypes)) {
    for (j in 1:ncol(genotypes)) {
      if (genotypes[i, j] == 9) {
        inferred_hap1[i, j] <- 9
        inferred_hap2[i, j] <- 9
      } else if (genotypes[i, j] == 2) {
        inferred_hap1[i, j] <- 1
        inferred_hap2[i, j] <- 1
      } else if (genotypes[i, j] == 0) {
        inferred_hap1[i, j] <- 0
        inferred_hap2[i, j] <- 0
      } else if (genotypes[i, j] == 1) {
        tmp <- sample(c(0,1), 1)
        inferred_hap1[i, j] <- tmp
        inferred_hap2[i, j] <- 1 - tmp
      }
    }
  }
  
  list(haplotype_1 = inferred_hap1, haplotype_2 = inferred_hap2)
}
# Explanation of `genotypes_to_haplotypes()`:
{
  # INPUT:
  # - genotypes: A matrix where rows represent individuals and columns represent 
  #     SNPs (coded as 0, 1, 2 for allele counts, or 9 for missing).
  #
  # WHAT'S HAPPENING:
  # 1. Initialization: Creates two empty matrices (`inferred_hap1` and `inferred_hap2`) 
  #    of the same size as the input to store the phased alleles.
  #
  # 2. Allele Phasing Logic: Iterates through every genotype and splits it:
  #     - Missing (9): Assigns 9 to both haplotypes.
  #     - Homozygous Alt (2): Assigns 1 (the alternative allele) to both haplotypes.
  #     - Homozygous Ref (0): Assigns 0 (the reference allele) to both haplotypes.
  #     - Heterozygous (1): Since it is unknown which allele came from which 
  #       parent, it randomly assigns a 0 to one haplotype and a 1 to the 
  #       other (using `sample`).
  #
  # 3. Structural Transformation: Converts the combined diploid counts back 
  #    into two distinct linear strands of genetic information.
  #
  # OUTPUT:
  # - A list containing two matrices: `haplotype_1` and `haplotype_2`, 
  #   representing the simulated phased haplotypes for all individuals.
}

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
# Explanation of `combine_haplotypes()`:
{
  # INPUT:
  # - haplotype_matrix: A matrix where individual haplotypes are stored as 
  #     separate rows (e.g., "ID1_hap1" and "ID1_hap2").
  #
  # WHAT'S HAPPENING:
  # 1. ID Extraction: Identifies unique individual IDs by stripping away the 
  #    haplotype suffixes (everything after the underscore) from the row names.
  #
  # 2. Header Preparation: Generates a new set of column names where every 
  #    original SNP column is split into two adjacent columns (suffixed _1 and _2). 
  #    This prepares the structure for a "long" genotype format.
  #
  # 3. Matrix Restructuring:
  #     - Initializes a new matrix with one row per unique individual and 
  #       twice as many columns as the input.
  #     - Iterates through each unique ID to find its corresponding haplotype rows.
  #     - Reorganizes the data so that the two alleles for each SNP are placed 
  #       side-by-side in the new row.
  #
  # OUTPUT:
  # - A consolidated matrix where each row represents a single individual, 
  #   and each pair of columns represents the two alleles of a specific locus.
}

#######################################################################################################################
#**•• 4_AlphaAssign Functions ••*********************************************
#######################################################################################################################

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
      nOffspring = nOffspring,
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
# Explanation of `Process_AlphaAssign_output()`:
{
  # INPUT:
  # - GE: Logical (TRUE/FALSE) indicating if Genotyping Error was included.
  # - True_pedigree: A reference data frame containing the actual known parentage.
  # - nOffspring: The total number of offspring used in the simulation.
  #
  # WHAT'S HAPPENING:
  # 1. Parameter Setup: Validates the 'GE' input and sets up file naming strings 
  #    ("WithGE" vs "NoGE") to ensure the function reads the correct external files.
  #
  # 2. Iterative File Processing: Loops through five different SNP group results:
  #     - File Reading: Loads the ".sires" output file generated by AlphaAssign.
  #     - Filtering: Extracts only the records where a sire was successfully 
  #       "chosen" (chosen == 1).
  #     - Pedigree Comparison: Merges the assigned sires with the "True_pedigree" 
  #       table based on the offspring ID.
  #     - Accuracy Calculation: Counts how many assigned sires match the known 
  #       biological sires.
  #
  # 3. Data Compilation: For each of the five groups, it stores the results 
  #    (counts of assignments and correct matches) into a structured data frame.
  #
  # 4. Result Aggregation: Combines the individual group data frames into one 
  #    comprehensive table for final analysis.
  #
  # OUTPUT:
  # - A single data frame summarizing the assignment accuracy and performance 
  #   of AlphaAssign across all tested SNP groups.
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
# Explanation of `ped_to_raw()`:
{
  # INPUT:
  # - ped_file: A PLINK-format .ped file (genotypes usually as two columns per SNP).
  # - map_file: A corresponding .map file containing SNP identifiers.
  # - output_file: The destination path for the processed data.
  #
  # WHAT'S HAPPENING:
  # 1. Data Loading: Reads the pedigree and map files into R, separating the 
  #    first six metadata columns (FID, IID, etc.) from the genotype columns.
  #
  # 2. Allele Pairing: Iterates through the genotype columns in pairs. It 
  #    concatenates the two alleles for each SNP into a single string (e.g., 
  #    "1" and "2" becomes "12").
  #
  # 3. Numeric Recoding: Uses a lookup map to convert the paired strings into 
  #    standard dosage format:
  #     - "11" (Homozygous Ref) -> 0
  #     - "12" or "21" (Heterozygous) -> 1
  #     - "22" (Homozygous Alt) -> 2
  #     - "00" (Missing) -> 9
  #
  # 4. Final Assembly: Renames the columns using the SNP IDs from the map file 
  #    and recombines the recoded genotypes with the original individual metadata.
  #
  # OUTPUT:
  # - A text file containing a consolidated matrix of metadata and numeric 
  #   genotype counts (0, 1, 2, 9).
}

#######################################################################################################################
#**•• 5_Prep_andRun_Beagle4 Functions ••*********************************************
#######################################################################################################################

convert_alleles <- function(allele) {
  if (allele == "1") {
    return("A")
  } else if (allele == "2") {
    return("C")
  } else if (allele == "9"){
    return("0")
    
  } else {
    return(allele)
  } }
# Explanation of `convert_alleles()`:
{
  # INPUT:
  # - allele: A single character or numeric value representing a genetic allele.
  #
  # WHAT'S HAPPENING:
  # The function performs a simple remapping of allele nomenclature. It checks 
  # the input against specific codes and translates them:
  #  - "1" is converted to "A".
  #  - "2" is converted to "C".
  #  - "9" (missing) is converted to "0".
  #  - Any other value is returned as-is without changes.
  #
  # OUTPUT:
  # - A single character value representing the converted allele.
}

convert_ped_genotypes <- function(ped_data, GE = NULL) {
  # Applying the conversion function to each allele in the genotype columns
  genotype_data <- ped_data[, 7:ncol(ped_data)]
  
  # Convert each allele using the sapply function
  corrected_genotype_data <- apply(genotype_data, 2, function(column) sapply(column, convert_alleles))
  
  
  # Replace the original genotype data with the converted data
  ped_data[, 7:ncol(ped_data)] <- corrected_genotype_data
  
  return(ped_data)
}
# Explanation of `convert_ped_genotypes()`:
{
  # INPUT:
  # - ped_data: A data frame in PLINK .ped format (where genotypes start from 
  #     column 7 onwards).
  # - GE: (Optional) A parameter for genotyping error, though not explicitly 
  #     used in the current logic of this function.
  #
  # WHAT'S HAPPENING:
  # 1. Data Subsetting: Isolates the genotype portion of the data frame by 
  #    selecting all columns from index 7 to the end, leaving the first 6 
  #    metadata columns untouched.
  #
  # 2. Element-wise Transformation: Uses `apply` and `sapply` to iterate 
  #    through every individual cell in the genotype matrix. It passes each 
  #    allele to the `convert_alleles` helper function to remap values 
  #    (e.g., changing "1" to "A" and "2" to "C").
  #
  # 3. Data Reintegration: Overwrites the original genotype columns in the 
  #    input data frame with the newly converted character-based alleles.
  #
  # OUTPUT:
  # - The modified `ped_data` data frame with standardized allele labels 
  #   (A, C, 0) in the genotype fields.
}

#######################################################################################################################
#**•• 6_Converting_PhasedVCF Functions ••*********************************************
#######################################################################################################################

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
# Explanation of `get_out_haplotypes()`:
{
  # INPUT:
  # - ped_matrix: A matrix where each cell contains a genotype string (e.g., "A C").
  # - ind_id_1: A vector of IDs for the first set of haplotypes.
  # - ind_id_2: A vector of IDs for the second set of haplotypes.
  # - realData: A logical flag determining which sorting function to use.
  #
  # WHAT'S HAPPENING:
  # 1. Matrix Initialization: Creates two empty matrices (`haplotype_matrix_0` 
  #    and `haplotype_matrix_1`) to store the separated alleles.
  #
  # 2. Haplotype Splitting: Iterates through every row and column of the input 
  #    matrix. For each genotype string, it splits the characters (e.g., 
  #    splitting "A C" into "A" and "C") and assigns the first allele to 
  #    matrix 0 and the second allele to matrix 1.
  #
  # 3. Labeling: Assigns the original SNP names to columns and the provided 
  #    individual IDs (`ind_id_1` and `ind_id_2`) to the rows of the 
  #    respective matrices.
  #
  # 4. Aggregation and Sorting: Vertically stacks the two matrices (using `rbind`) 
  #    to create a single long-format matrix. It then calls a helper function 
  #    (`order_by_prefix` or `order_by_prefix_realData`) to organize the rows 
  #    based on their ID prefixes.
  #
  # OUTPUT:
  # - A single combined matrix containing all separated haplotypes, sorted 
  #   sequentially by individual.
}

convert_genotypes <- function(genotypes) {
  #genotypes[genotypes == '1'] <- 0
  genotypes[genotypes == '0'] <- NA
  #genotypes[genotypes == '2'] <- 1
  genotypes[genotypes == 'A'] <- 0
  genotypes[genotypes == 'C'] <- 1
  
  return(as.numeric(genotypes))
}
# Explanation of `convert_genotypes()`:
{
  # INPUT:
  # - genotypes: A vector or matrix of genotype values (typically character strings).
  #
  # WHAT'S HAPPENING:
  # The function performs a specific remapping of genotype labels to numeric 
  # values. It scans the input and applies the following logic:
  #  - Converts "0" to NA (Not Available) to handle missing data.
  #  - Converts "A" to 0 (typically representing the reference allele).
  #  - Converts "C" to 1 (typically representing the alternative allele).
  #  - Note: Commented-out lines suggest it was previously designed to handle 
  #    numeric "1" and "2" as well.
  #
  # 2. Type Casting: Forces the final result into a numeric format using 
  #    `as.numeric()` to ensure it is ready for mathematical modeling.
  #
  # OUTPUT:
  # - A numeric vector or matrix where allele labels are replaced by 0, 1, or NA.
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
# Explanation of `convert_VCF()`:
{
  # INPUT:
  # - vcf_file: Path to a phased VCF file containing genetic variants.
  # - map_file: A reference file containing SNP identifiers/positions.
  #
  # WHAT'S HAPPENING:
  # 1. External Data Extraction: Uses system calls to `bcftools` to extract 
  #    essential information from the VCF:
  #     - "PhasedGT.txt": Contains the Reference allele, Alternative allele, 
  #       and the phased genotypes (e.g., 0|1).
  #     - "SampleIDs.txt": Contains the list of individual sample names.
  #
  # 2. Genotype Decoding: 
  #     - Loads the extracted data and transposes the genotype matrix so 
  #       individuals are rows and SNPs are columns.
  #     - Iterates through each SNP and replaces the numeric codes with actual 
  #       nucleotides: '0' becomes the REF allele and '1' becomes the ALT allele.
  #     - Replaces the phasing pipe symbol ("|") with a space to separate the 
  #       two alleles (e.g., "0|1" becomes "A T").
  #
  # 3. Labeling and Cleanup:
  #     - Cleans individual IDs by removing specific prefixes (e.g., "AMEL_").
  #     - Assigns the cleaned IDs to the rows and the SNP names from the 
  #       map file to the columns.
  #
  # OUTPUT:
  # - A character matrix of phased genotypes where cells contain space-separated 
  #   nucleotides (e.g., "A C") for each individual and locus.
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
# Explanation of `Haplotype_using_pedigree()`:
{
  # INPUT:
  # - GenErr: A logical or numeric indicator of genotyping error (validated but
  #     not directly used in this specific snippet).
  # - n: A parameter likely representing the iteration or group number.
  # - ped_recon: The reconstructed pedigree data object.
  # - pedigree_name: File path to the pedigree metadata (ID, Sire, Dam).
  # - haplo_name: The string name of the global object containing the 
  #     full haplotype data.
  #
  # WHAT'S HAPPENING:
  # 1. Input Validation & Loading: Ensures all required parameters are present, 
  #    loads the pedigree file, and retrieves the large haplotype matrix 
  #    from the global environment using `get()`.
  #
  # 2. Population Subsetting: Filters the master haplotype list into three 
  #    distinct groups based on the pedigree:
  #     - Mothers (Dams)
  #     - DPCs (Sires/Drones)
  #     - Workers (Offspring)
  #    It then recombines these into a single matrix containing only the 
  #    relevant individuals for the current analysis.
  #
  # 3. ID Generation: Creates unique identifiers for the two physical strands 
  #    of each individual by appending "_1" and "_2" to their base IDs.
  #
  # 4. Haplotype Processing:
  #     - Uses `get_out_haplotypes` to separate the diploid genotype strings 
  #       into individual alleles.
  #     - Applies `convert_genotypes` to recode character alleles (A, C) 
  #       into a numeric format (0, 1).
  #
  # OUTPUT:
  # - A numeric matrix of separated haplotypes where rows are identified by 
  #   individual ID and strand number (e.g., "ID1_1"), formatted for 
  #   downstream genetic analysis.
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
# Explanation of `Haplotype_using_pedigree_realData()`:
{
  # INPUT:
  # - pedigree_name: File path to the pedigree metadata containing 'id', 
  #     'dpc' (sire), and 'mother'.
  # - haplo_name: The string name of the global object containing the 
  #     full haplotype data.
  #
  # WHAT'S HAPPENING:
  # 1. Pedigree Alignment: Loads the pedigree file and filters it to include 
  #    only individuals that actually exist within the provided haplotype 
  #    dataset, ensuring the analysis doesn't crash on missing records.
  #
  # 2. Data Retrieval: Fetches the master haplotype matrix from the global 
  #    environment and cleans up row names by removing unnecessary prefixes.
  #
  # 3. Targeted Subsetting: Extracts the specific haplotypes for mothers, 
  #    dpc/sires, and offspring (workers) defined in the pedigree, then 
  #    combines them into a focused dataset.
  #
  # 4. Haplotype Extraction:
  #     - Calls `get_out_haplotypes` with the `realData = TRUE` flag, 
  #       which triggers a specific ordering logic optimized for real-world 
  #       genomic datasets.
  #     - Converts the resulting character alleles (e.g., A, C) into numeric 
  #       dosage format (0, 1) via `convert_genotypes`.
  #
  # OUTPUT:
  # - A numeric matrix of separated and phased haplotypes (two rows per 
  #   individual), sorted and labeled according to real-data prefixing 
  #   conventions.
}

order_by_prefix <- function(df) {
  # Extract the part before the underscore
  prefix <- as.numeric(sapply(rownames(df), function(x) strsplit(x, "_")[[1]][1]))
  
  # Order the data frame based on the extracted prefix
  df_ordered <- df[order(prefix, rownames(df)), , drop = FALSE]
  
  return(df_ordered)
}
# Explanation of `order_by_prefix()`:
{
  # INPUT:
  # - df: A data frame or matrix where row names contain numeric IDs followed 
  #     by an underscore (e.g., "10_1", "2_2").
  #
  # WHAT'S HAPPENING:
  # 1. Prefix Extraction: It splits each row name at the underscore and 
  #    extracts the first part. It converts these strings into actual 
  #    numbers so they can be sorted mathematically (ensuring "2" comes 
  #    before "10").
  #
  # 2. Multi-Level Sorting: It reorganizes the rows using two criteria:
  #     - Primary: Sorts by the numeric value of the prefix.
  #     - Secondary: Sorts by the full row name string (to keep "_1" and 
  #       "_2" in the correct order for the same ID).
  #
  # 3. Structure Preservation: Uses `drop = FALSE` to ensure the object 
  #    remains a data frame/matrix even if the result only contains one row.
  #
  # OUTPUT:
  # - An ordered version of the input data frame where individuals are 
  #   sequenced numerically by their ID.
}

order_by_prefix_realData <- function(df) {
  # Extract the part before the underscore
  prefix <- sapply(rownames(df), function(x) strsplit(x, "_")[[1]][1])
  
  # Order the data frame based on the extracted prefix
  df_ordered <- df[order(prefix, rownames(df)), , drop = FALSE]
  
  return(df_ordered)
}
# Explanation of `order_by_prefix_realData()`:
{
  # INPUT:
  # - df: A data frame or matrix where row names contain alphanumeric prefixes 
  #     followed by an underscore (e.g., "SampleA_1", "SampleB_2").
  #
  # WHAT'S HAPPENING:
  # 1. String Extraction: Splits each row name at the underscore and extracts 
  #    the first part as a character string (prefix). Unlike the standard 
  #    version, this does not convert to numeric, allowing it to handle 
  #    non-numeric sample names.
  #
  # 2. Alphabetical Sorting: Reorganizes the rows using a two-tier sort:
  #     - Primary: Sorts alphabetically by the extracted string prefix.
  #     - Secondary: Sorts by the full row name to ensure haplotype pairs 
  #       (_1 and _2) remain grouped and ordered correctly.
  #
  # 3. Structure Preservation: Uses `drop = FALSE` to prevent the matrix from 
  #    coercing into a vector if the resulting data has only one row.
  #
  # OUTPUT:
  # - A data frame or matrix ordered lexicographically by the sample 
  #   identifiers found in the row names.
}

#######################################################################################################################
#**•• 7_Haplotype_ParentageAssignments Functions ••*********************************************
#######################################################################################################################

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
# Explanation of `name_genotypes()`:
{
  # INPUT:
  # - map_file: A data frame representing a PLINK .map file, where column V2 
  #     contains the SNP/marker names.
  #
  # WHAT'S HAPPENING:
  # 1. Row Duplication: The function iterates through each marker in the map 
  #    file. For every single marker, it creates two identical copies.
  #
  # 2. Suffix Assignment: It modifies the marker name (column V2) for each 
  #    copy to distinguish between the two alleles of a diploid genotype:
  #     - The first copy gets a "_1" suffix (e.g., "SNP101_1").
  #     - The second copy gets a "_2" suffix (e.g., "SNP101_2").
  #
  # 3. List Expansion: These modified rows are stored in a list, effectively 
  #    doubling the total number of entries.
  #
  # 4. Vector Extraction: It recombines the list into a "long" data frame 
  #    and isolates the newly created marker names.
  #
  # OUTPUT:
  # - A character vector of marker names that is twice as long as the input, 
  #   formatted for phased haplotype or long-format genotype data.
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
# Explanation of `transform_genotype()`:
{
  # INPUT:
  # - genotype_matrix: A matrix where individual genotypes are currently 
  #     stored in a "wide" format with paired columns (e.g., "SNP1_1" and "SNP1_2").
  #
  # WHAT'S HAPPENING:
  # 1. Coordinate Identification: Identifies individual IDs from the row names 
  #    and unique marker names by stripping the "_1/_2" suffixes from 
  #    the column headers.
  #
  # 2. Data Deconstruction (Tidying):
  #     - Iterates through every individual and every marker.
  #     - Extracts the two separate alleles (Value 1 and Value 2) for each marker.
  #     - Creates new identifiers for the haplotypes by appending "_1" and 
  #       "_2" to the individual’s ID.
  #     - Stores these in a "long" list format (ID, Marker, Value).
  #
  # 3. Data Reconstruction (Reshaping):
  #     - Uses `reshape()` to pivot the data from a long list back into 
  #       a wide matrix.
  #     - The new structure places the separated haplotypes (e.g., "Ind1_1", 
  #       "Ind1_2") as unique rows and the markers as columns.
  #
  # 4. Cleanup: Standardizes column names by removing the "Value." prefix 
  #    generated during reshaping and sets the final row names.
  #
  # OUTPUT:
  # - A transformed matrix where each row represents a single haplotype 
  #   strand and each column represents a specific SNP marker.
}

get_genotypes <- function(ped, id) {
  row <- ped %>% filter(IID == id)
  genotypes <- row[,7:ncol(row)]
  return(genotypes)
}
# Explanation of `get_genotypes()`:
{
  # INPUT:
  # - ped: A data frame in PLINK .ped format (containing metadata in the first 
  #     6 columns and genotypes from column 7 onwards).
  # - id: The specific individual identifier (IID) to be extracted.
  #
  # WHAT'S HAPPENING:
  # 1. Filtering: Uses the `dplyr` pipe (`%>%`) and `filter` to search the 
  #    dataset for the row where the Individual ID (IID) matches the provided 'id'.
  #
  # 2. Column Selection: Once the correct row is found, it slices the data 
  #    to remove the first 6 metadata columns (Family ID, Sex, Phenotype, etc.).
  #
  # 3. Isolation: It focuses exclusively on the genetic markers (columns 7 
  #    to the end of the data frame).
  #
  # OUTPUT:
  # - A single-row data frame containing only the genotype data for the 
  #   specified individual.
}

Hap1 <- function(df, map) {
  # Get the column names that end in _2
  columns_ending_in_1 <- grep("_1$", colnames(df), value = TRUE)
  
  # Subset the dataframe with the selected columns
  new_df <- df[, columns_ending_in_1, drop = FALSE]
  
  colnames(new_df) <- map[,2]
  return(new_df)
}
# Explanation of `Hap1()`:
{
  # INPUT:
  # - df: A data frame or matrix of genotypes in "wide" format (where columns 
  #     are named with _1 and _2 suffixes, e.g., "SNP1_1", "SNP1_2").
  # - map: A reference map data frame where the second column (V2) contains 
  #     the original SNP names.
  #
  # WHAT'S HAPPENING:
  # 1. Pattern Matching: Uses `grep` to identify every column name that 
  #    ends specifically with the suffix "_1".
  #
  # 2. Haplotype Isolation: Subsets the original data frame to keep only 
  #    those columns, effectively extracting the first physical strand 
  #    (Haplotype 1) of the genetic data.
  #
  # 3. Label Resetting: Replaces the suffixed column names (e.g., "SNP1_1") 
  #    with the clean, original SNP names from the map file to make the 
  #    data readable.
  #
  # OUTPUT:
  # - A data frame containing only the first haplotype for all individuals, 
  #   with standardized SNP names as headers.
}

Hap2 <- function(df, map) {
  # Get the column names that end in _2
  columns_ending_in_2 <- grep("_2$", colnames(df), value = TRUE)
  
  # Subset the dataframe with the selected columns
  new_df <- df[, columns_ending_in_2, drop = FALSE]
  
  colnames(new_df) <- map[,2]
  return(new_df)
}
# Explanation of `Hap2()`:
{
  # INPUT:
  # - df: A data frame or matrix of genotypes in "wide" format (columns 
  #     ending in _1 and _2, representing the two paired alleles).
  # - map: A reference map data frame where the second column (V2) contains 
  #     the original SNP names.
  #
  # WHAT'S HAPPENING:
  # 1. Pattern Matching: Uses `grep` to search the column headers for every 
  #    name that ends specifically with the suffix "_2".
  #
  # 2. Haplotype Isolation: Filters the data frame to retain only those 
  #    columns, effectively extracting the second physical strand 
  #    (Haplotype 2) of the genetic data.
  #
  # 3. Label Resetting: Renames the subsetted columns (e.g., changing 
  #    "SNP101_2" to "SNP101") using the clean SNP IDs from the map file.
  #
  # OUTPUT:
  # - A data frame containing only the second haplotype for all individuals, 
  #   with original SNP names as column headers.
}

calcGeno <- function(Hap1, Hap2, map){
  Hap1 <- t(Hap1)
  Hap2 <- t(Hap2)
  Geno <- matrix(data = Hap1 + Hap2, nrow = 1)
  colnames(Geno) <- map[,2]
  
  return(Geno)
}
# Explanation of `calcGeno()`:
{
  # INPUT:
  # - Hap1: A vector or single-row data frame representing the first haplotype 
  #     (alleles 0 or 1).
  # - Hap2: A vector or single-row data frame representing the second haplotype.
  # - map: A reference map data frame where column 2 contains SNP names.
  #
  # WHAT'S HAPPENING:
  # 1. Transposition: Converts both input haplotypes into a vertical format (t) 
  #    to ensure they are treated as vectors for addition.
  #
  # 2. Additive Genotype Calculation: Adds the two haplotypes together. Since 
  #    alleles are coded as 0 and 1, this creates standard "dosage" values:
  #     - 0 + 0 = 0 (Homozygous Reference)
  #     - 0 + 1 or 1 + 0 = 1 (Heterozygous)
  #     - 1 + 1 = 2 (Homozygous Alternative)
  #
  # 3. Formatting: Wraps the resulting values into a single-row matrix and 
  #    assigns the correct SNP identifiers from the map file to the columns.
  #
  # OUTPUT:
  # - A single-row matrix representing the full diploid genotype (0, 1, 2) 
  #   for an individual across all markers.
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
# Explanation of `assign_parent_haplo()`:
{
  # INPUT:
  # - df: A data frame or matrix representing a single haplotype strand.
  # - maternal: A logical or object flag indicating the source is the mother.
  # - paternal: A logical or object flag indicating the source is the father.
  # - offspring_id: The ID of the child/offspring being processed.
  #
  # WHAT'S HAPPENING:
  # 1. Labeling Logic: The function determines which parent the genetic material 
  #    came from to create a specific row name:
  #     - If 'maternal' is provided, it creates a label like "ChildID_maternal".
  #     - If 'paternal' is provided, it creates a label like "ChildID_paternal".
  #    *Note: If both are provided, the paternal label will overwrite the 
  #    maternal one due to the sequential "if" structure.*
  #
  # 2. Row Name Assignment: It applies this new descriptive label to the row 
  #    names of the input data frame 'df'.
  #
  # OUTPUT:
  # - The modified data frame 'df' with updated row names, clarifying the 
  #   parental origin of that specific haplotype.
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
# Explanation of `get_phased_haplotypes()`:
{
  # INPUT:
  # - Geno_Error: A character string, either "GE" (Genotyping Error) 
  #     or "nGE" (No Genotyping Error).
  #
  # WHAT'S HAPPENING:
  # 1. Environment Setup: Changes the working directory to locate the 
  #    required pedigree and genomic files.
  #
  # 2. Conditional File Loading:
  #     - If "GE": Loads files containing simulated genotyping errors.
  #     - If "nGE": Loads the clean version of the phased dataset.
  #     - Throws an error (stop) if neither flag is provided correctly.
  #
  # 3. Numeric Conversion:
  #     - Strips the metadata (columns 1-6).
  #     - Uses the previously defined `convert_genotypes` function to 
  #       change allele characters (A, C, etc.) into 0 and 1.
  #
  # 4. Labeling & Formatting:
  #     - Uses `name_genotypes` to create paired marker headers (e.g., SNP_1, SNP_2).
  #     - Sets the row names to the individual IDs (V2 from the .ped file).
  #
  # 5. Data Restructuring: 
  #    - Calls `transform_genotype` to pivot the data from a wide format 
  #      (one row per person) to a long/phased format (two rows per person, 
  #      representing the two individual haplotype strands).
  #
  # OUTPUT:
  # - A structured matrix/data frame of phased haplotypes ready for 
  #   downstream analysis.
}

order_by_prefix <- function(df) {
  # Extract the part before the underscore
  prefix <- as.numeric(sapply(rownames(df), function(x) strsplit(x, "_")[[1]][1]))
  
  # Order the data frame based on the extracted prefix
  df_ordered <- df[order(prefix, rownames(df)), , drop = FALSE]
  
  return(df_ordered)
}
# Explanation of `order_by_prefix()`:
{
  # INPUT:
  # - df: A data frame or matrix where row names contain numeric IDs followed 
  #     by an underscore and a haplotype index (e.g., "10_1", "10_2", "2_1").
  #
  # WHAT'S HAPPENING:
  # 1. Prefix Extraction: It uses `strsplit` to isolate the text before the 
  #    underscore. Crucially, it wraps this in `as.numeric()`. This ensures 
  #    that "10" is treated as a larger number than "2" (standard string 
  #    sorting would put "10" before "2").
  #
  # 2. Sequential Ordering: The `order()` function is called with two arguments:
  #     - Primary: The numeric `prefix` (sorts by Individual ID).
  #     - Secondary: The original `rownames(df)` (sorts by the suffix _1 or _2).
  #
  # 3. Matrix Maintenance: Uses `drop = FALSE` to prevent R from automatically 
  #    converting a single-row result into a vector, preserving the data 
  #    frame structure.
  #
  # OUTPUT:
  # - A data frame sorted logically by ID: 
  #   (e.g., 1_1, 1_2, 2_1, 2_2, 10_1, 10_2).
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
# Explanation of `check_haplotype()`:
{
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
}

weighing_haplotypes <- function(off_hap, par_geno, method = NULL, result = NULL, j){
  

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
# Explanation of `weighing_haplotypes()`:
{
  # This function acts as a genetic scoring engine that determines how well a potential parent matches an offspring.
  # INPUT:
  # - off_hap: A numeric vector (0/1) representing the offspring's haplotype.
  # - par_geno: A numeric vector (0/1/2) representing the potential parent's genotype.
  # - method: A string defining the mathematical scoring model (e.g., "quadratic").
  # - j: An identifier for the current chromosome being analyzed.
  #
  # WHAT'S HAPPENING:
  # 1. Conflict Detection:
  #    It creates a 'logical_vec' to find impossible genetic combinations.
  #    If (Offspring Haplotype - Parent Genotype) is -2 or 1, it implies a 
  #    mismatch (e.g., offspring has allele '1' but parent is '0/0').
  #
  # 2. Run Length Encoding (RLE):
  #    The function uses `rle()` to group consecutive matches (FALSE) and 
  #    consecutive mismatches (TRUE). In genetics, long "runs" of matches 
  #    suggest a true biological relationship, while frequent breaks or 
  #    long runs of mismatches suggest the individual is not the parent.
  #
  # 3. Statistical Scoring:
  #    Based on the chosen 'method', it calculates scores for:
  #     - false_score: Penalty for mismatches (higher = less likely to be parent).
  #     - true_score: Reward for match runs (higher = more likely to be parent).
  #     - overall_score: Usually (False - True). A lower (or more negative) 
  #       score indicates a better parental candidate.
  #
  # 4. Math Models:
  #    - method options :  c("mean","sum","quadratic","logarithmic" ,"geometric","composite", "exponential", "harmonic","power_mean") 
  #    - "mean/sum": Standard linear averages.
  #    - "quadratic/exponential": Penalizes long mismatch segments heavily.
  #    - "logarithmic/harmonic": Less sensitive to extreme outliers in run lengths.
  #
  # OUTPUT:
  # - A data frame containing the mismatch penalty, the match reward, the 
  #   final consolidated score, the method used, and the chromosome ID.
}

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
# Explanation of `check_dist()`:
{
  # SUMMARY:
  # This function acts as a diagnostic tool to map the distribution of 
  # genetic matches and mismatches across a chromosome.
  #
  # INPUT:
  # - off_hap: The offspring's haplotype (0/1).
  # - par_geno: The potential parent's genotype (0/1/2).
  # - j: The identifier for the current chromosome.
  #
  # WHAT'S HAPPENING:
  # 1. Conflict Identification:
  #    It creates a 'logical_vec' where TRUE represents a biological mismatch 
  #    (impossible inheritance) and FALSE represents a valid match.
  #
  # 2. Segment Analysis (RLE):
  #    Uses Run Length Encoding to group consecutive matches and mismatches. 
  #    It identifies how many markers in a row stay the same before a change 
  #    occurs (e.g., a "run" of 50 matches followed by a "run" of 2 mismatches).
  #
  # 3. Data Structuring:
  #    - Extracts lengths of "FALSE" runs (matches) into a data frame.
  #    - Extracts lengths of "TRUE" runs (mismatches) into a data frame.
  #    - Includes safety checks to assign 'NA' if no runs of a certain type exist.
  #
  # 4. Consolidation:
  #    Combines these into a single table, labeled with the chromosome ID (j).
  #
  # OUTPUT:
  # - A data frame listing every continuous genetic segment, its length, 
  #   and whether it was a Match (FALSE) or Mismatch (TRUE).
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
# Explanation of `Route1_flipping()`:
{
  # SUMMARY:
  # This function identifies whether an offspring's haplotypes came 
  # from the mother or the father and labels them accordingly.
  #
  # INPUT:
  # - Data_type: The specific dataset being processed.
  # - pedigree: A table linking offspring to their parents.
  # - perfect_haplotypes: Boolean; toggles between simulated "true" data 
  #     and experimental phased data.
  # - method: The mathematical scoring model used to test parental matches.
  #
  # WHAT'S HAPPENING:
  # 1. Data Setup: For every family in the pedigree, it gathers the 
  #    haplotypes for the offspring and both parents.
  #
  # 2. Comparison: It iterates through the chromosomes, comparing each of 
  #    the offspring's two strands against the mother's and father's genotypes.
  #
  # 3. Decision Logic: It analyzes the match/mismatch scores to determine 
  #    parental origin. For example, if Strand A matches the father much 
  #    better than the mother, it is labeled "paternal."
  #
  # 4. Error Tracking: It runs `check_dist` during the process to keep 
  #    track of the lengths of matching and mismatching segments.
  #
  # 5. Organization: It collects all the labeled strands and statistical 
  #    results into a final structured list.
  #
  # OUTPUT:
  # - A list containing the final labeled haplotypes ("_maternal" or 
  #   "_paternal") and the scoring data used to make those decisions.
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
      offspring_id <- as.character(pedigree$id[i])
      dam_id <- as.character(pedigree$mother[i])
      
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
      dam_id <- pedigree$mother[i]
      
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
      
      for (j in 1){#unique(map[, 1])) {
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
# Explanation of `Route2_flipping()`:
{
  # SUMMARY:
  # This function determines which offspring haplotype is maternal by 
  # comparing them against only the mother's genotype; it then assumes 
  # the remaining strand is paternal.
  #
  # INPUT:
  # - Data_type: The dataset being used (e.g., "Real_Slov_data").
  # - pedigree: A table linking offspring to their mothers.
  # - perfect_haplotypes: Boolean; toggles between "true" simulated data 
  #     and experimental phased data.
  # - method: The mathematical scoring model used to test for matches.
  #
  # WHAT'S HAPPENING:
  # 1. Targeted Comparison: Unlike Route 1, this function only pulls the 
  #    maternal genotype for comparison. It calculates two scores: 
  #    (Offspring Hap1 vs Mom) and (Offspring Hap2 vs Mom).
  #
  # 2. Origin by Exclusion:
  #    - The strand with the significantly better match to the mother is 
  #      formally labeled "_maternal".
  #    - The other strand is automatically labeled "_paternal" by process 
  #      of elimination.
  #
  # 3. Decision Logic: It checks the "max difference" between the two 
  #    maternal scores. If the scores are tied, the function stops with 
  #    an error, as it cannot distinguish the source.
  #
  # 4. Data Collection: Similar to Route 1, it gathers diagnostic run-length 
  #    data (`check_dist`) and scoring statistics across all chromosomes.
  #
  # OUTPUT:
  # - A list containing the labeled maternal and paternal haplotypes and 
  #   the statistical data used to justify the maternal assignment.
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
# Explanation of `check_haplotype_postFlip()`:
{
  # SUMMARY:
  # This function acts as a validation tool that calculates the accuracy 
  # of the flipping process by comparing the assigned haplotypes against 
  # a "ground truth" dataset.
  #
  # INPUT:
  # - complete_haplotypes: The "truth" dataset (correctly labeled haplotypes).
  # - results: The output from Route 1 or Route 2 flipping functions.
  # - pedigree: A table of offspring IDs to be checked.
  #
  # WHAT'S HAPPENING:
  # 1. Alignment: It filters the "truth" data to match the column names 
  #    (markers) present in your results.
  #
  # 2. Iterative Comparison: For each offspring, it isolates the 
  #    maternal and paternal strands from both the "truth" and the "results."
  #
  # 3. Accuracy Calculation:
  #    - It performs a marker-by-marker comparison (results == truth).
  #    - It calculates a percentage score for both strands (e.g., "98% match").
  #
  # 4. Statistical Summary: After checking everyone, it calculates and 
  #    prints the mean accuracy and the range (min/max) for both 
  #    maternal and paternal assignments.
  #
  # OUTPUT:
  # - A data frame containing the percentage accuracy scores for 
  #   every individual in the pedigree.
}


#######################################################################################################################
#**•• 8_Checking_Gametic_MendelianSamplingValues Functions ••*********************************************
#**- NOT YET ALL TEST AND MAY NEED UPDATED -**
#######################################################################################################################

Plotting_Gametic_R1 <- function(df, phased_type, plotting_styles) {
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
# Explanation of `Plotting_Gametic_R1()`:
{
  # SUMMARY:
  # This function generates visual reports to compare Mendelian sampling values 
  # between maternal and paternal haplotypes across different phasing scenarios.
  #
  # INPUT:
  # - df: The data frame containing Mendelian sampling values and phasing labels.
  # - phased_type: The categories to plot (e.g., "Phased_GE", "True", "Real").
  # - plotting_styles: The type of charts to create ("histogram", "density").
  #
  # WHAT'S HAPPENING:
  # 1. Validation: It checks that your phasing types and requested styles 
  #    are valid before attempting to draw anything.
  #
  # 2. Color Standardization: It defines a strict color code—Blue for Maternal 
  #    and Red for Paternal—to ensure all charts are easy to read consistently.
  #
  # 3. Plot Generation Loop:
  #    - It iterates through each Phasing category.
  #    - For each category, it builds the requested chart (Histogram or Density).
  #    - It overlays the Maternal and Paternal data on the same axes to 
  #      highlight differences or overlaps.
  #
  # 4. Layout & cowplot: It uses the `cowplot` package to stitch individual 
  #    charts together into a single, clean multi-panel figure.
  #
  # 5. Legend Management: It extracts the legend from a single plot and 
  #    attaches it to the side of the final combined image to avoid cluttering 
  #    individual panels.
  #
  # OUTPUT:
  # - A single combined graphic displaying the distribution of Mendelian 
  #   sampling across the specified datasets and styles.
}

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
# Explanation of `Plotting_Gametic_R2()`:
{
  # SUMMARY:
  # This function generates visual reports specifically for maternal Mendelian 
  # sampling values across different phasing categories.
  #
  # INPUT:
  # - df: The data frame containing Mendelian sampling values and phasing labels.
  # - phased_type: The categories to plot (e.g., "Phased_GE", "True", "Real").
  # - plotting_styles: The type of charts to create ("histogram", "density").
  #
  # WHAT'S HAPPENING:
  # 1. Validation: It confirms that the requested phasing categories and 
  #    plotting styles are valid before execution.
  #
  # 2. Targeted Focus: Unlike R1, this function focuses exclusively on 
  #    the Maternal data, applying a consistent blue color scale.
  #
  # 3. Chart Generation:
  #    - It loops through each phasing type and style.
  #    - It builds either histograms (for counts) or density plots (for 
  #      distribution shape) of maternal Mendelian sampling.
  #
  # 4. Multi-Panel Layout: It uses `cowplot` to arrange the different 
  #    charts into a single row per phasing type, then combines those rows 
  #    into one final dashboard.
  #
  # 5. Clean Legend: It extracts a single legend and places it on the 
  #    right side of the final grid to keep the individual charts uncluttered.
  #
  # OUTPUT:
  # - A single combined graphic displaying maternal Mendelian sampling 
  #   distributions for the specified datasets and styles.
}

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
# Explanation of `calculate_gametic_relatedness_R1()`:
{
  # SUMMARY:
  # This function calculates "Mendelian Sampling" values for both maternal 
  # and paternal haplotypes to see how much an offspring deviates from 
  # the average of its parents' DNA.
  #
  # INPUT:
  # - sorted_offspring_haplotypes: The offspring strands already labeled 
  #     as "_maternal" and "_paternal".
  # - all_haplotypes: The master list of all individual haplotypes.
  # - pedigree: A table linking offspring to their specific mothers and fathers.
  #
  # WHAT'S HAPPENING:
  # 1. Family Retrieval: For every offspring in the pedigree, it pulls the 
  #    offspring's sorted strands and both parent's original haplotypes.
  #
  # 2. Parental Averaging: 
  #    It calculates the "Parental Average" (mid-parent value) for each parent 
  #    by averaging their two haplotypes. This represents the expected 
  #    value if inheritance were perfectly uniform.
  #
  # 3. Mendelian Sampling Calculation:
  #    - It subtracts the Mother's average from the offspring's maternal strand.
  #    - It subtracts the Father's average from the offspring's paternal strand.
  #    - Formula: $ri = OffspringHap - 0.5(ParentHap1 + ParentHap2)$
  #
  # 4. Summation: It sums these differences across all markers to get a 
  #    single "Mendelian sampling" value for each parent-offspring side.
  #
  # OUTPUT:
  # - A data frame containing the Offspring ID and the calculated Mendelian 
  #   sampling deviations for both the Maternal and Paternal sides.
}

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
# Explanation of `calculate_gametic_relatedness_R2()`:
{
  # SUMMARY:
  # This function calculates the Mendelian Sampling value for the maternal 
  # haplotype only, identifying how much the inherited maternal strand 
  # deviates from the mother's average genotype.
  #
  # INPUT:
  # - sorted_offspring_haplotypes: Offspring strands labeled by the 
  #     Route 2 "exclusion" method.
  # - all_haplotypes: The master list of all individual haplotypes.
  # - pedigree: A table linking offspring to their mothers.
  #
  # WHAT'S HAPPENING:
  # 1. Targeted Retrieval: For each offspring, it pulls the identified 
  #    maternal strand and the mother's two original haplotypes.
  #
  # 2. Maternal Averaging: 
  #    It calculates the expected mid-parent value for the mother by 
  #    averaging her two strands ($0.5 \times (Hap1 + Hap2)$).
  #
  # 3. Deviation Calculation:
  #    It subtracts the mother's average from the offspring's maternal 
  #    strand. This numerical difference represents the "sampling" 
  #    effect—essentially which alleles the offspring actually received 
  #    versus the 50/50 statistical expectation.
  #
  # 4. Results Consolidation: 
  #    It sums these differences across all markers and saves the 
  #    total for each offspring-mother pair.
  #
  # OUTPUT:
  # - A data frame containing the Offspring ID, Mother ID, and the 
  #   single calculated Maternal Mendelian sampling value.
}

#######################################################################################################################
#**•• 9_DeterminingPatrilines Functions ••*********************************************
#**- NOT YET ALL TEST AND MAY NEED UPDATED -**
#######################################################################################################################

#Route 1 - comparing to drones with sister and father thresholds 
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
# Explanation of `calc_nPaternity_Accuracy()`:
{
  # SUMMARY:
  # This function evaluates how accurately the system can identify the number 
  # of unique fathers (drones) for a set of sisters (workers) based on 
  # their paternal haplotypes.
  #
  # INPUT:
  # - results: The paternal haplotypes identified for each worker.
  # - sister_threshold: How similar two workers' paternal DNA must be 
  #     to consider them "full sisters" (from the same father).
  # - father_test_threshold: How similar a worker's DNA must be to a known 
  #     father's DNA to confirm a direct match.
  # - father_haplotypes: The master list of "true" father DNA.
  #
  # WHAT'S HAPPENING:
  # 1. Queen-Based Grouping: It processes the population one queen at a 
  #    time, pulling all her offspring (workers) and the known fathers 
  #    she mated with.
  #
  # 2. Match Phase (Worker vs. Known Father):
  #    It compares every worker's paternal strand against the actual 
  #    fathers. If the similarity exceeds the `father_test_threshold`, 
  #    the worker is assigned to that father’s group.
  #
  # 3. Clustering Phase (Sister vs. Sister):
  #    For workers who didn't clearly match a known father, the function 
  #    compares them to each other. If they are very similar 
  #    (`sister_threshold`), it clusters them together into a 
  #    "predicted" father group.
  #
  # 4. Accuracy Check:
  #    It compares the "Estimated" number of fathers (based on the DNA clusters) 
  #    against the "Actual" number of fathers from the pedigree.
  #
  # OUTPUT:
  # - A table summarizing for each queen: how many workers she had, how many 
  #   fathers were expected, how many were estimated, and how many of those 
  #   estimates were perfectly correct.
}

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
# Explanation of `run_paternity_tests()`:
{
  # SUMMARY:
  # This function is a "batch runner." It automates the paternity analysis 
  # by testing multiple different sensitivity settings (thresholds) to 
  # find which ones yield the most accurate results.
  #
  # INPUT:
  # - results_arg: The paternal haplotypes identified in previous steps.
  # - sister_thresholds: A list of values to test for sibling similarity.
  # - father_test_thresholds: A list of values to test for father similarity.
  # - father_haplotypes: The master list of "true" father DNA.
  #
  # WHAT'S HAPPENING:
  # 1. Nested Looping: The function uses a "grid search" approach. It 
  #    pairs every possible sister threshold with every possible father 
  #    threshold (e.g., trying 90% sister match with 95% father match, 
  #    then 90% with 98%, etc.).
  #
  # 2. Dynamic Naming: For each run, it creates a unique ID based on the 
  #    settings used (e.g., "sis0.9_fat0.95"). This prevents data from 
  #    being overwritten.
  #
  # 3. Execution: It calls `calc_nPaternity_Accuracy` for every combination, 
  #    passing the specific thresholds and the DNA data into the engine.
  #
  # 4. Storage: It collects the results of every single combination and 
  #    saves them into one large master list.
  #
  # OUTPUT:
  # - A comprehensive list of results. Each entry in the list shows how 
  #   well the paternity estimation performed under that specific 
  #   combination of settings.
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
# Explanation of `plot_paternity_number_grid_DRONES()`:
{
  # SUMMARY:
  # This function visualizes the results of your paternity tests by creating 
  # a grid of charts. It compares the "Actual" number of fathers to the 
  # "Estimated" number across every threshold setting you tested.
  #
  # INPUT:
  # - Results: The master list generated by `run_paternity_tests()`.
  #
  # WHAT'S HAPPENING:
  # 1. Data Merging: It takes the list of results and flattens them into 
  #    one large table. It extracts the threshold values (sister and drone) 
  #    directly from the list names.
  #
  # 2. Reshaping: It pivots the data into a "long format" so that "Actual" 
  #    and "Estimated" counts can be easily mapped to different colors 
  #    and shapes on the same plot.
  #
  # 3. Creating the Grid (Faceting): 
  #    This is the core of the function. It uses `facet_grid` to create 
  #    a matrix of plots:
  #    - Columns represent different Sister Thresholds.
  #    - Rows represent different Drone (Father) Thresholds.
  #
  # 4. Visualization:
  #    - Red Circles: Represent the ground truth (Actual nPatrilines).
  #    - Blue Triangles: Represent what the algorithm predicted (Estimated).
  #    - The points are slightly "dodged" (shifted) so they don't hide 
  #      each other if the estimate is perfectly accurate.
  #
  # OUTPUT:
  # - A comprehensive ggplot object. This visual diagnostic helps you 
  #   instantly see which combination of thresholds makes the blue 
  #   triangles land closest to the red circles.
}

#Route 2 - using only sister thresholds to determine partilines 
SisterONLY_calc_nPaternity_Accuracy <- function(results, sister_threshold, pedigree, simulated = NULL) {
  
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
# Explanation of `SisterONLY_calc_nPaternity_Accuracy()`:
{
  # SUMMARY:
  # This function estimates the number of fathers (patrilines) by clustering 
  # offspring together based solely on how similar their paternal DNA is 
  # to one another, without comparing them to known father DNA.
  #
  # INPUT:
  # - results: The paternal haplotypes identified for the offspring.
  # - sister_threshold: The minimum similarity required to group two 
  #     workers as "full sisters" (sharing the same father).
  # - pedigree: The family tree used to group offspring by their mother.
  # - simulated: Boolean; if TRUE, it compares the estimate against 
  #     known "actual" father counts from a simulation.
  #
  # WHAT'S HAPPENING:
  # 1. Family Sorting: Like the previous versions, it groups offspring 
  #    by their mother (Queen ID) to analyze one colony at a time.
  #
  # 2. Pairwise Comparison: It compares every worker in the colony against 
  #    every other worker. It calculates a similarity score based on 
  #    matching alleles.
  #
  # 3. Agglomerative Clustering:
  #    - If Worker A and Worker B are similar enough (>= sister_threshold), 
  #      they are put into a "Sister Group."
  #    - If Worker C matches either A or B, Worker C is added to that same group.
  #    - Workers who don't match anyone are eventually placed in their own 
  #      individual groups (representing a father who only appeared once).
  #
  # 4. Counting: The total number of unique clusters formed represents the 
  #    estimated number of fathers who contributed to that colony.
  #
  # OUTPUT:
  # - A data frame showing the Queen ID, the number of workers analyzed, 
  #   and the final count of estimated sister groups (patrilines).
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
# Explanation of `run_sister_clustering_tests()`:
{
  # SUMMARY:
  # This function is the "automated supervisor" for sister-only clustering. 
  # It iterates through various similarity settings to help you find the 
  # "sweet spot" for grouping sisters without using father DNA.
  #
  # INPUT:
  # - results_arg: The paternal haplotypes for the population.
  # - sister_thresholds: A vector of similarity values to test (e.g., 0.75, 0.90).
  # - pedigree: The family tree data.
  # - simulated: Boolean; tells the engine whether to check estimates 
  #     against a known "truth" dataset.
  #
  # WHAT'S HAPPENING:
  # 1. Parameter Sweeping: It loops through every value provided in 
  #    `sister_thresholds`. This allows you to see how your father count 
  #    estimates change if you make the "sisterhood" rules stricter or looser.
  #
  # 2. Execution: For each threshold, it calls the `SisterONLY` engine. 
  #    It names each run (e.g., "sis0.85") so you can compare them later.
  #
  # 3. Progress Tracking: It prints status updates to the console so 
  #    you know which threshold is currently being calculated and when 
  #    the entire batch is finished.
  #
  # 4. Result Bundling: It collects all the different data frames 
  #    (one for each threshold) into a single master list.
  #
  # OUTPUT:
  # - A list of data frames. Each entry contains the estimated number 
  #   of patrilines for that specific sister-similarity threshold.
}

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
# Explanation of `plot_paternity_number_grid_SISTERONLY()`:
{
  # SUMMARY:
  # This function creates a visual "stress test" grid. It compares actual 
  # vs. estimated father counts across three different datasets and 
  # multiple sister-similarity thresholds simultaneously.
  #
  # INPUT:
  # - ... : Accepts multiple result lists (specifically True, nGE, and GE 
  #     datasets) generated by the sister clustering tests.
  #
  # WHAT'S HAPPENING:
  # 1. Multi-Dataset Consolidation: It takes three distinct result sets 
  #    and merges them into a single master data frame, tagging each row 
  #    with its source (e.g., "True_Route2", "nGE_Route2", "GE_Route2").
  #
  # 2. Threshold Extraction: It cleans up the "sis0.XX" labels to convert 
  #    them into numeric values for proper sorting and plotting.
  #
  # 3. Comparative Faceting (The Grid):
  #    It builds a powerful visual matrix using `facet_grid`:
  #    - ROWS: Represent the different Datasets (how does the model perform 
  #      on "True" data vs. "GE" data?).
  #    - COLUMNS: Represent the Sister Thresholds (how does changing the 
  #      similarity requirement affect the accuracy?).
  #
  # 4. Visual Comparison:
  #    - Red Circles: The actual number of fathers.
  #    - Blue Triangles: The estimated number of sister groups.
  #    - Like the drone plot, it "dodges" the points so you can see exactly 
  #      where the estimate and the reality diverge.
  #
  # OUTPUT:
  # - A large, multi-paneled ggplot. This is the "final verdict" graphic 
  #   that tells you which dataset and which threshold produce the most 
  #   reliable paternity counts.
}














