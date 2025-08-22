# Clean workspace
rm(list = ls())

#library
library(AlphaSimR)
library(SIMplyBee)
library(readr)
library(genio)

# Create/download the founder genomes from SIMplyBee 
# founderGenomes <- simulateHoneyBeeGenomes(nCar = 10,
#                                           nChr = 16,
#                                           nSegSites = 1000)

load("~/Desktop/Slovenia data/Attempt2/founderGenomes_10ind.Rdata")



#Set up the SP ####
SP <- SimParamBee$new(founderGenomes, csdChr = ifelse(16 >= 3, 3, 1), nCsdAlleles = 128)
SP$nWorkers <- 30
# Track the pedigree
SP$setTrackPed(TRUE)
# Track the recombination
SP$setTrackRec(TRUE)

SP$addSnpChip(nSnpPerChr = 10) # CHANGE THE SNP SIZE HERE!!!!! (ADD A FEW EXTRA FOR QC THEN SAMPLE THE CORRECT SIZES)# define csd chromomsome
csdChr <- SP$csdChr
# Add traits - taken from the QuantGen vignettte
mean <- c(20, 0)
varA <- c(1, 1 / SP$nWorkers)
corA <- matrix(data = c( 1.0, -0.5,
                         -0.5,  1.0), nrow = 2, byrow = TRUE)
SP$addTraitA(nQtlPerChr = 100, mean = mean, var = varA, corA = corA,
             name = c("queenTrait", "workersTrait"))

varE <- c(3, 3 / SP$nWorkers)
# TODO: what is a reasonable environmental correlation between queen and worker effects?
corE <- matrix(data = c(1.0, 0.3,
                        0.3, 1.0), nrow = 2, byrow = TRUE)
SP$setVarE(varE = varE, corE = corE)

#Create base population
basePop <- createVirginQueens(founderGenomes)
basePop1 <- randCross(basePop)
basePop2 <- randCross(basePop1)
#Create drones to create dpc (THAT ARE RELATED BECAUSE MATING STATION!!!)
pFathers <- nFathersPoisson
nDronesPerQueen = 50
drones <- createDrones(basePop[1], nInd = nDronesPerQueen)
fathers <- pullDroneGroupsFromDCA(drones, n = 1, nDrones = pFathers)

#Create a colony and cross it
colony1 <- createColony(x = basePop[2])
colony1 <- cross(colony1, drones = fathers[[1]])


#create Mating station Drones where the dpcs are
DPQs <- createVirginQueens(colony1, nInd = 4)
matingStation_drones <- createDrones(DPQs, nInd = nDronesPerQueen)

# Introduce 8 visting virgin queens
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


#If we are adding in genotyping errors - these are the parameters 
error  = 0.0001
error2 = 0.00001
missing = 0.005 #using this as it matches the missingness in the real Slovenian Data 
sampling_error = 0.005

#Genotyping error functions
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



#Using PopMerged_noFathers since the genotypic info of the drones is not available in the real Slov data 
geno_noError = pullSnpGeno(PopMerged_noFathers, chr = NULL, asRaw = FALSE, simParam=SP)
haplo_noError = pullSnpHaplo(PopMerged_noFathers, chr = NULL, asRaw = FALSE, simParam=SP)
geno_WithError = generateGenoErr(geno = geno_noError, error = error, error2 = error2, sampling_error = sampling_error, missing = missing)
SNPMap = getSnpMap(snpChip = 1, simParam = SP)

write.table(SNPMap[c("chr", "id", "site")], paste0("SNPMapFile_2.txt"), sep = " ", na = "0", quote = F, row.names = TRUE, col.names = TRUE)
write.table(geno_noError, paste0("SNP2_geno_noError.txt"), sep = " ", na = "0", quote = F, row.names = TRUE, col.names = TRUE)
write.table(haplo_noError, paste0("SNP2_haplo_noError.txt"), sep = " ", na = "0", quote = F, row.names = TRUE, col.names = TRUE)
write.table(geno_WithError, paste0("SNP2_geno_WithError.txt"), sep = " ", na = "0", quote = F, row.names = TRUE, col.names = TRUE)

#write plink files (ped and map) for each SNP array size 
#we use the population without fathers to mimic the real data since we would not have the drone data 
writePlink(PopMerged_noFathers, baseName = paste0('SNP2_NoGenoError'), simParam = SP)

#Next we need to phase the genotyped data 
#We are going to take the largest array size, phase that then take out the phased smaller arrays from that (finger crossed)

#Load the ped file of the largest array size 
Bigarrayped <- read.table("SNP2_NoGenoError.ped")
#The genotypes here have no genotyping error so they need to be replaced with the ones that do 
BigarrayGenoError <- read.table("SNP2_geno_WithError.txt")


transform_genotypes <- function(genotype_df) {
  # Define the transformation mapping
  mapping <- c("0" = "1 1", "1" = "1 2", "2" = "2 2", "9" = "0 0")
  
  # Initialize an empty list to store the transformed columns
  transformed_list <- list()
  
  # Loop through each column in the dataframe
  for (col_name in colnames(genotype_df)) {
    # Apply the transformation to each element in the column
    transformed_col <- sapply(genotype_df[[col_name]], function(x) mapping[as.character(x)])
    
    # Split the transformed values into two separate columns
    split_col <- do.call(rbind, strsplit(transformed_col, " "))
    
    # Create new column names
    new_col_names <- paste(col_name, c("1", "2"), sep = "_")
    
    # Add the new columns to the list
    transformed_list[[new_col_names[1]]] <- as.numeric(split_col[, 1])
    transformed_list[[new_col_names[2]]] <- as.numeric(split_col[, 2])
  }
  
  # Convert the list to a dataframe
  transformed_df <- as.data.frame(transformed_list)
  
  # Preserve the row names (individual IDs)
  rownames(transformed_df) <- rownames(genotype_df)
  
  return(transformed_df)
}

transformed_df <- transform_genotypes(BigarrayGenoError)

Bigarray_ped_withGenoError <- cbind(Bigarrayped[,c(1:6)], transformed_df)
#lets save the 1 2 0 version
write.table(Bigarray_ped_withGenoError, "SNP2_WithGenoError_12format.ped", sep = " ", quote = F, col.names = T, row.names = F)

#For Beagle the genotypes need to be in the format AGCT so lets change them 

# Function to convert allele codes
convert_alleles <- function(allele) {
  if (allele == "1") {
    return("A")
  } else if (allele == "2") {
    return("C")
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

# Convert alleles
converted_ped_data <- convert_ped_genotypes(ped_data = Bigarray_ped_withGenoError )


# Write the corrected PED file
write.table(converted_ped_data, file = "SNP2_withGenoError_ACformat.ped", quote = FALSE, sep = " ", row.names = FALSE, col.names = F)



