# Clean workspace
rm(list = ls())

#library
library(AlphaSimR)
library(SIMplyBee)
library(readr)
library(genio)
library(forecast)

args = commandArgs(trailingOnly=TRUE)
Rep = args[1]
workingDir = args[2]
softwareDir = args[3]

repDir = paste0(workingDir, "/SimRep", Rep, "/")


# Set working directory
setwd(repDir)


nSNP_array <- rbind(c(3125, 5L), c(108, 4L), c(50, 3L), c(10, 2L), c(1,1L))
# Make 4 SNP chip sizes with an SNP array #######

# Create/download the founder genomes from SIMplyBee
# founderGenomes <- simulateHoneyBeeGenomes(nCar = 10,
#                                           nChr = 16,
#                                           nSegSites = 1000)

#setwd("~/github/lstrachan_patrilines/")

pathToPlink=softwareDir
#workingDir = "/home/jana/github/lstrachan_patrilines/"
#setwd(workingDir)
load("Pipeline/0_FounderPop_inbred.RData")

matingStation_drones <- createDrones(DPQs, nInd = nDronesPerQueen)

# Introduce 8 visiting virgin queens
virgin_queens <- basePop[3:10]

# cross mating station DCA  with virgin queens
nCastePoisson <- SIMplyBee::nDronesPoisson
queens <- SIMplyBee::cross(x = virgin_queens, drones = pullDroneGroupsFromDCA(DCA = matingStation_drones, n = 8, nDrones = pFathers))
queen_colonies <- createMultiColony(queens)
#create 30 workers per queen
queen_colonies <- buildUp(queen_colonies, nWorkers = 50, exact = TRUE)

for (colony in 1:nColonies(queen_colonies)) {
  queen_colonies@colonies[[colony]]@workers <- selectInd(queen_colonies@colonies[[colony]]@workers, nInd = 30, use = "rand")
}

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

write_csv(col, file = "Data/worker_pedigree.csv", col_names = TRUE)


#Save the general files at this point (just incase we need to come back to them)
save(queen_colonies, file = "Data/queen_colonies.Rdata")
save(DPQs, file = "Data/DPQs.Rdata")
save(matingStation_drones, file = "Data/matingStation_drones.Rdata")
save(fathers, file = "Data/real_fathers.Rdata")
save(PopMerged, file = "Data/Pop_withFathers.Rdata")
save(PopMerged_noFathers, file = "Data/Pop_withNoFathers.Rdata")
save(SP, file = "Data/SP_object.Rdata")

#Pull out the whole GenMap for this simulation (all of the sexes)
GenMap_full <- getGenMap(SP, sex = "A")
GenMap_femaleOnly <- getGenMap(SP, sex = "F")
write.table(GenMap_full, file = "Data/GenMap_full.txt", sep = " ", quote = F, col.names = T, row.names = F)
write.table(GenMap_femaleOnly, file = "Data/GenMap_femaleOnly.txt", sep = " ", quote = F, col.names = T, row.names = F)


#We are showing genotyping errors and non-genotyping error version
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

#nest arrayssss
for (n in 2:6){ #we create the first array (larger one) outside of the loop and subset it to create the others
  tmp1 = vector()
  for (chr in 1:16){
    tmp2 = sample(colnames(pullSnpGeno(pop = founderGenomes, snpChip = (n-1), chr = chr, simParam = SP)), size = nSNP_array[(n-1),1], replace = FALSE)
    
    tmp1 = c(tmp1, tmp2)
  }
  SP$addSnpChipByName(tmp1, name=paste0('SNP_', nSNP_array[n-1,2]))
  write.table(tmp1, paste0("Data/SNP_", nSNP_array[(n-1),2], "_array.txt"), sep = " ", na = "NA", quote = F, row.names = FALSE, col.names = FALSE)
}

#this function will go through the genotypes with errors and use the true haplotypes as a reference for the genotype with value 1.
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




dir.create("Data/Sim_NoGE")
dir.create("Data/Sim_WithGE")


#Now we make a loop for all of the different arrays, going through and saving the needed info for each one.
#We are saving the genotype matrix (both full version and one with genotypic errors), along with the map files of the SNP chips and the haplotypes for the full genotypic data
for (n in 1:nrow(nSNP_array)){
  
  print(paste0('SNP_', nSNP_array[n,2], '_array.txt'))
  SNP_array = as.vector(read_table(paste0('Data/SNP_', nSNP_array[n,2], '_array.txt'), show_col_types = FALSE, col_names = FALSE))
  SP$addSnpChipByName(SNP_array$X1, name=paste0('array_', nSNP_array[n,2]))
  
  #Using PopMerged_noFathers since the genotypic info of the drones is not available in the real Slov data
  print("Adding genotyping errors")
  geno_noError = pullSnpGeno(PopMerged_noFathers, snpChip= paste0('array_', nSNP_array[n,2]), chr = NULL, asRaw = FALSE, simParam=SP)
  haplo_noError = pullSnpHaplo(PopMerged_noFathers, snpChip= paste0('array_', nSNP_array[n,2]), chr = NULL, asRaw = FALSE, simParam=SP)
  SNPmap_file = getSnpMap(snpChip= paste0('array_', nSNP_array[n,2]), simParam=SP)
  geno_WithError = generateGenoErr(geno = geno_noError, error = error, error2 = error2, sampling_error = sampling_error, missing = missing)
  
  print("Writing files")
  write.table(SNPmap_file[c("chr", "id", "site")], paste0("Data/Sim_NoGE/SNPMapFile_", nSNP_array[n,2] , ".txt"), sep = " ", na = "0", quote = F, row.names = TRUE, col.names = TRUE)
  write.table(SNPmap_file[c("chr", "id", "site")], paste0("Data/Sim_WithGE/SNPMapFile_", nSNP_array[n,2] , ".txt"), sep = " ", na = "0", quote = F, row.names = TRUE, col.names = TRUE)
  
  write.table(geno_noError, paste0("Data/Sim_NoGE/SNP", nSNP_array[n,2], "_geno_noError.txt"), sep = " ", na = "0", quote = F, row.names = TRUE, col.names = TRUE)
  write.table(haplo_noError, paste0("Data/Sim_NoGE/SNP", nSNP_array[n,2], "_haplo_noError.txt"), sep = " ", na = "0", quote = F, row.names = TRUE, col.names = TRUE)
  write.table(geno_WithError, paste0("Data/Sim_WithGE/SNP", nSNP_array[n,2], "_geno_WithError.txt"), sep = " ", na = "0", quote = F, row.names = TRUE, col.names = TRUE)
  
  #write plink files (ped and map) for each SNP array size
  #we use the population without fathers to mimic the real data since we would not have the drone data
  print("Write Plink files")
  writePlink(PopMerged_noFathers, baseName = paste0('Data/Sim_NoGE/SNP_', nSNP_array[n,2], '_NoGE'), snpChip = paste0('array_', nSNP_array[n,2]), simParam = SP)
  
  #we have the genotyping errors, but we need to make them into a plink ped file and that is done below:
  print("Preparing plink files + vcfs")
  #1. first we take the true haplotypes so that we can use them for column names
  true_haplotype1 <- haplo_noError[grep("_1$", rownames(haplo_noError)),]
  true_haplotype2 <- haplo_noError[grep("_2$", rownames(haplo_noError)),]
  
  #2. Because we pull genotypes from SimplyBee, with pullSnpGeno(0,1,2), they need to be converted back to haploid form 
  genoError_haplotypes <- genotypes_to_haplotypes(genotypes = geno_WithError)
  genoError_hap1 <- genoError_haplotypes[["haplotype_1"]]
  rownames(genoError_hap1) <- rownames(true_haplotype1)
  genoError_hap2 <- genoError_haplotypes[["haplotype_2"]]
  rownames(genoError_hap2) <- rownames(true_haplotype2)
  
  #3. This joins _1 alleles and then _2 alleles 
  combined_matrix <- rbind(genoError_hap1, genoError_hap2)
  colnames(combined_matrix) <- colnames(haplo_noError)
  
  #4. Now we need to order it by individual
  #Extract the numeric part and whether it's _1 or _2
  numeric_part <- as.numeric(gsub("_.*", "", rownames(combined_matrix)))
  suffix_part <- gsub(".*_", "", rownames(combined_matrix))
  # Create a data frame to order the rows
  order_df <- data.frame(numeric_part, suffix_part)
  # Order first by the numeric part, then by the suffix (_1 before _2)
  ordered_indices <- order(order_df$numeric_part, order_df$suffix_part)
  # Reorder the combined matrix
  ordered_matrix <- combined_matrix[ordered_indices, ]
  
  #5. This create a plink ped file with 2 columns for 1 SNP 
  GE_ped_format <- combine_haplotypes(ordered_matrix)
  GE_ped_format[GE_ped_format == 1] <- 2
  GE_ped_format[GE_ped_format == 0] <- 1
  GE_ped_format[GE_ped_format == 9] <- 0
  
  #6. Read plink file in to get the first 6 columns to make plink ped for one with genotyping errors
  first_six_columns = data.frame(FamID = "AMEL", ID = rownames(GE_ped_format), FID = 0, MID = 0, Sex = 0, Phenotype = -9)
  GE_ped <- cbind(first_six_columns, GE_ped_format)
  write.table(GE_ped, file= paste0("Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE.ped"), quote=F, row.names=F, col.names=F, sep = " ")
  
  #Need to modify the ped files for BEAGLE (genotypes need to being the 0ACGT format not 0,1,2,9) 
  #print("Converting to AC format for Beagle")
  #NoGE_ped_AC <- convert_ped_genotypes(ped_data = NoGE_ped, GE = FALSE)
  #write.table(NoGE_ped_AC, file = paste0("Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE_ACformat.ped"), quote = FALSE, sep = " ", row.names = FALSE, col.names = F)
  #system(paste0("cp Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE.map Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE_ACformat.map")) 

  #WithGE_ped_AC <- convert_ped_genotypes(ped_data = GE_ped, GE = TRUE)
  #write.table(WithGE_ped_AC, file = paste0("Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE_ACformat.ped"), quote = FALSE, sep = " ", row.names = FALSE, col.names = F)
  #system(paste0("cp Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE.map Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE_ACformat.map"))
  system(paste0("cp Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE.map Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE.map"))

  
  
  #Now you can put things into plink for quality control 
  print("Running PLINK QC")
  # take the corrected ped file and the map file:
  #system(paste0(pathToPlink,"/plink --file Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE_ACformat --make-bed  --geno 0.1 --mind 0.1 --maf 0.01 --out Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE_ACformat_QC")) #quality control
  system(paste0(pathToPlink,"/plink --file Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE --make-bed  --geno 0.1 --mind 0.1 --maf 0.01 --missing --out Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE_QC")) #quality control
  #system(paste0(pathToPlink,"/plink --file Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE_ACformat --make-bed  --geno 0.1 --mind 0.1 --maf 0.01 --out Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE_ACformat_QC")) #quality control
  system(paste0(pathToPlink,"/plink --file Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE --make-bed  --geno 0.1 --mind 0.1 --maf 0.01 --missing --out Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE_QC")) #quality control


  #system(paste0(pathToPlink,"/plink --bfile Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE_ACformat_QC  --recode --out Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE_ACformat_QC "))  #get the ped and map files back
  system(paste0(pathToPlink,"/plink --bfile Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE_QC  --recode --out Data/Sim_NoGE/SNP_", nSNP_array[n,2], "_NoGE_QC "))  #get the ped and map files back
  #system(paste0(pathToPlink,"/plink --bfile Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE_ACformat_QC --recode --out Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE_ACformat_QC"))   #get the ped and map files back
  system(paste0(pathToPlink,"/plink --bfile Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE_QC --recode --out Data/Sim_WithGE/SNP_", nSNP_array[n,2], "_WithGE_QC"))   #get the ped and map files back
}
  

save.image(file = paste0(repDir, "/Pipeline/1_SimulatedData_prepared.Rdata"))
