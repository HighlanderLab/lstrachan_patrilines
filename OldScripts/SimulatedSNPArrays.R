# Clean workspace
rm(list = ls())

#library
library(AlphaSimR)
library(SIMplyBee)
library(readr)
library(genio)

# Make 4 SNP chip sizes with an SNP array #######

#Modified functions from Audrey's nesting functions 
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




# Create/download the founder genomes from SIMplyBee 
# founderGenomes <- simulateHoneyBeeGenomes(nCar = 10,
#                                           nChr = 16,
#                                           nSegSites = 1000)

load("~/Desktop/Slovenia data/Attempt2/founderGenomes_10ind.Rdata")

#Set up the SP #
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
write_csv(col, file = "worker_pedigree.csv", col_names = TRUE)


#Save the general files at this point (just incase we need to come back to them)
save(queen_colonies, file = "queen_colonies.Rdata")
save(DPQs, file = "DPQs.Rdata")
save(matingStation_drones, file = "matingStation_drones.Rdata")
save(fathers, file = "real_fathers.Rdata")
save(PopMerged, file = "Pop_withFathers.Rdata")
save(PopMerged_noFathers, file = "Pop_withNoFathers.Rdata")
save(SP, file = "SP_object.Rdata")

#Pull out the whole GenMap for this simulation (all of the sexes)
GenMap_full <- getGenMap(SP, sex = "A")
GenMap_femaleOnly <- getGenMap(SP, sex = "F")
write.table(GenMap_full, file = "GenMap_full.txt", sep = " ", quote = F, col.names = T, row.names = F)
write.table(GenMap_femaleOnly, file = "GenMap_femaleOnly.txt", sep = " ", quote = F, col.names = T, row.names = F)


#We are showing genotyping errors and non-genotyping error version 
#If we are adding in genotyping errors - these are the parameters 
error  = 0.0001
error2 = 0.00001
missing = 0.005 #using this as it matches the missingness in the real Slovenian Data 
sampling_error = 0.005

# Genotyping error functions
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


#Modifying the ped files #####
#For Beagle the genotypes need to be in the format AGCT so lets change them 

# Function to convert allele codes
convert_alleles <- function(allele) {
  if (allele == "2") {
    return("C")
  } else if (allele == "1") {
    return("A")
  } else if (allele == "0"){
    return("0")
  
  } else {
    return(allele)
  }
}

convert_ped_genotypes <- function(ped_data) {
  # Applying the conversion function to each allele in the genotype columns
  genotype_data <- ped_data[, 7:ncol(ped_data)]
  
  # Convert each allele using the sapply function
  corrected_genotype_data <- apply(genotype_data, 2, function(column) sapply(column, convert_alleles))
  
  # Replace the original genotype data with the converted data
  ped_data[, 7:ncol(ped_data)] <- corrected_genotype_data
  
  return(ped_data)
}

GE_ped_12 <- convert_ped_genotypes(ped_data = GE_ped_AC)

# Convert alleles
Slov_ped_filtered2 <- convert_ped_genotypes(ped_data = Slov_ped_filtered)
GE_ped_AC <- convert_ped_genotypes(ped_data = GE_ped)
nGE_ped_AC <- convert_ped_genotypes(ped_data = nGE_ped)





#################################################################################

# Use PLINK for quality control and conversion to a vcf ready for phasing

#################################################################################

# take the corrected ped file and the map file: 
#./plink --file  --make-bed  --geno 0.1 --mind 0.1 --maf 0.01 --out  #quality control 
#./plink --bfile  --recode vcf --out  #make a vcf of the QC for phasing 
#./plink --vcf  "vcf file name" --recode --out  #get ped files back for pedigree reconstruction

