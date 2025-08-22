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

#I THINK MY GENOTYPE-to-HAPLOTYPE function has corrected this problem 
#You cant have A 9 in the ped file for the vcf it needs to be  0 0 so lets change it 
# Function to check and correct half-missing genotypes
# correct_half_missing <- function(genotype1, genotype2) {
#   if (genotype1 == 9) genotype1 <- 0
#   if (genotype2 == 9) genotype2 <- 0
#   if ((genotype1 == 0 && genotype2 != 0) || (genotype1 != 0 && genotype2 == 0)) {
#     return(c(0, 0)) # Replace half-missing with full-missing
#   } else {
#     return(c(genotype1, genotype2))
#   }
# }
# 
# # Apply the function to each genotype pair in the dataframe to preQC dataframe 
# for (i in seq(7, ncol(GE_nQC_ped), by = 2)) { # Start from the 7th column (first genotype) and step by 2
#   GE_nQC_ped[, c(i, i+1)] <- t(apply(GE_nQC_ped[, c(i, i+1)], 1, function(x) correct_half_missing(x[1], x[2])))
# }
# Write the corrected PED file
#write.table(converted_ped_data, file = "SNP4_ped_withGenoError_ACformat.ped", quote = FALSE, sep = " ", row.names = FALSE, col.names = F)


#Now you can put things into plink 
# take the corrected ped file and the map file: 
./plink --file  --make-bed  --geno 0.1 --mind 0.1 --maf 0.01 --out  #quality control 
./plink --bfile  --recode vcf --out  #make a vcf of the QC for phasing 
./plink --bfile  --recode --out  #get the ped and map files back 


#convert the vcf back into a useable ped file
./plink --vcf phased_test.vcf --recode --out phased_SNP5_nested



#THIS FOR THE SOFTWARE INPUTS ##############################

#Setting up the input files #####
pedigree_file <- read.csv("worker_pedigree.csv")

# Sequoia ####

#Life History - need ID, Sex, BirthYear in a csv file 
Workers <- data.frame(ID = pedigree_file$id,
                      Sex = rep(1, length(pedigree_file$id)),
                      BirthYear = rep(2024, length(pedigree_file$id)))
#make the parents birthyear the year after offspring 
Mothers <- data.frame(ID = unique(pedigree_file$mother),
                      Sex = rep(1, length(unique(pedigree_file$mother))),
                      BirthYear = rep(2023, length(unique(pedigree_file$mother))))

#make the dpc sex 2 = male to make sequoia think they're fathers 
Dpc <- data.frame(ID = unique(pedigree_file$dpc),
                      Sex = rep(2, length(unique(pedigree_file$dpc))),
                      BirthYear = rep(2023, length(unique(pedigree_file$dpc))))

LifeHistory <- rbind(Mothers, Dpc, Workers)
write.csv(LifeHistory, file = "LifeHistory.csv", sep = ",", quote = F, col.names = T, row.names = F)

#Known Dpcs
Worker_known <- data.frame(id = pedigree_file$id,
                           dam = pedigree_file$mother,
                           sire = pedigree_file$dpc)
Mothers_known <- data.frame(id = unique(pedigree_file$mother),
                            dam = rep(NA, length(unique(pedigree_file$mother))),
                            sire = rep(NA, length(unique(pedigree_file$mother))))

Dpc_known <- data.frame(id = unique(pedigree_file$dpc),
                        dam = rep(NA, length(unique(pedigree_file$dpc))),
                        sire = rep(NA, length(unique(pedigree_file$dpc))))

Known_Dpc <- rbind(Mothers_known, Dpc_known, Worker_known)
write.csv(Known_Dpc, file = "Known_Dpc.csv", sep = ",", quote = F, col.names = T, row.names = F)
#Unknown Dpcs 
Worker_Unknown <- data.frame(id = pedigree_file$id,
                           dam = pedigree_file$mother,
                           sire = rep(NA, length(pedigree_file$id)))

Mothers_Unknown <- data.frame(id = unique(pedigree_file$mother),
                            dam = rep(NA, length(unique(pedigree_file$mother))),
                            sire = rep(NA, length(unique(pedigree_file$mother))))

Dpc_Unknown <- data.frame(id = unique(pedigree_file$dpc),
                        dam = rep(NA, length(unique(pedigree_file$dpc))),
                        sire = rep(NA, length(unique(pedigree_file$dpc))))

Unknown_Dpc <- rbind(Mothers_Unknown, Dpc_Unknown, Worker_Unknown)
write.csv(Unknown_Dpc, file = "Unknown_Dpc.csv", sep = ",", quote = F, col.names = T, row.names = F)


library(sequoia)
#convert ped file to sequoia format 
map <- read.table("phased_SNP2_ped_withGenoError_ACformat.map", header = F)
ped <- read.table("phased_SNP2_ped_withGenoError_ACformat.ped", header = F)

SNP_names <- map[,2]
SNP_names_new <- unlist(lapply(SNP_names, function(name) c(paste(name, "1", sep="_"), paste(name, "2", sep="_"))))
colnames(ped)[7:ncol(ped)] <- SNP_names_new


#add Output file if you are running it on Eddie 
Sequoia_ped <- GenoConvert(InData = ped, InFormat = 'ped', Missing = '0', IDcol = 2, header = F)
#Sequoia_ped <- GenoConvert(InData = ped, InFormat = 'ped', Missing = '0', IDcol = 2, header = F, OutFile = "Sequoia_ped.txt")

print("Run sequoia")
SequoiaOutPut <- sequoia(GenoM = Sequoia_ped2, LifeHistData = LifeHistory, Module = "ped", Plot = TRUE)
save(SequoiaOutPut, file = "SequoiaOutPut_NestSNP2_phasedGE.Rdata")

rm(... = map, ped, SNP_names, SNP_names_new, Sequoia_ped, SequoiaOutPut)

#Count opposite homozygous (OH) loci between parent-offspring pairs and Mendelian errors (ME) between parent-parent-offspring trios, and calculate the parental log-likelihood ratios (LLR)
CalcOHLLR_unknown <- CalcOHLLR(Pedigree = Unknown_Dpc, GenoM = Sequoia_ped, CalcLLR = TRUE, LifeHistData = LifeHistory, SeqList = SequoiaOutPut)
CalcOHLLR_known <- CalcOHLLR(Pedigree = Known_Dpc, GenoM = Sequoia_ped, CalcLLR = TRUE, LifeHistData = LifeHistory, SeqList = SequoiaOutPut)

write.table(CalcOHLLR_unknown, file = "CalcOHLLR_unknown.txt", quote = F, sep = " ", col.names = T, row.names = F, na = "NA")
write.table(CalcOHLLR_known, file = "CalcOHLLR_known.txt", quote = F, sep = " ", col.names = T, row.names = F, na = "NA")

#Calculate Pedigree relatedness 
PedRel_Pedpar <- CalcRped(Pedigree = SequoiaOutPut$PedigreePar, OUT = 'DF')


#Organise the Sequoia summary file 
#Compare with known parents
Known_Dpc <- read_csv("Known_Dpc.csv")
Known_Dpc <- as.data.frame(Known_Dpc)
load("SequoiaOutPut_NestSNP4_nGE.Rdata")

PC_par <- PedCompare(Ped1 = Known_Dpc[, c("id", "dam", "sire")],
                     Ped2 = SequoiaOutPut$PedigreePar)

nSires_assigned <- sum(!is.na(SequoiaOutPut$PedigreePar$sire))
nCorrect_sires <- sum(PC_par[["MergedPed"]][["sire.class"]] == "Match")

Sequoia_file <- data.frame(Data_Group = "Nested",
                           Test = "Non_GenoErr",
                           nOffspring = 240,
                           SNP_group = 4,
                           nSires_assigned = nSires_assigned,
                           nCorrect_sires = nCorrect_sires,
                           Software = "Sequoia")

write.table(Sequoia_file, file = "SNP4_nGE_Seq.txt", sep = " ", quote = F, col.names = T, row.names = F)


tmp1 <- read.table("SNP4_nGE_Seq.txt", header = T)
tmp2 <- read.table("SNP3_nGE_Seq.txt", header = T)
tmp3 <- read.table("SNP2_nGE_Seq.txt", header = T)
tmp4 <- read.table("SNP1_nGE_Seq.txt", header = T)


SNP1_summary <- rbind(tmp1, tmp2, tmp3, tmp4)
write.table(SNP1_summary, file = "nGE_summary.txt", sep = " ", quote = F, col.names = T, row.names = F)

SNP4_summary <- read.table("phasedGE_summary.txt", header = T)
SNP3_summary <- read.table("nGE_summary.txt", header = T)
SNP2_summary <- read.table("GE_summary.txt", header = T)

Sequoia_summary <- rbind(SNP2_summary, SNP3_summary, SNP4_summary)

write.table(Sequoia_summary, file = "Sequoia_summary.txt", sep = " ", quote = F, col.names = T, row.names = F)


#On the slov data
Sequoia_PedPar <- SequoiaOutPut$PedigreePar


Sequoia_file <- data.frame(Data_Group = "Slov",
                           Test = "GenoErr",
                           nOffspring = 235,
                           SNP_group = NA,
                           nSires_assigned = sum(!is.na(Sequoia_PedPar$sire)),
                           nCorrect_sires = NA,
                           Software = "Sequoia")

write.table(Sequoia_file, file = "Sequoia_summary.txt", sep = " ", quote = F, col.names = T, row.names = F)


#Create summary of them all (all in their own folders)
Slov_Sequoia <- read.table("Sequoia_summary.txt", header = T)
Non_Nested_Sequoia <- read.table("Sequoia_summary.txt", header = T)
Nested_Sequoia <- read.table("Sequoia_summary.txt", header = T)

Sequoia_summary <- rbind(Slov_Sequoia, Non_Nested_Sequoia, Nested_Sequoia)

write.table(Sequoia_summary, file = "Sequoia_summary.txt", sep = " ", quote = F, col.names = T, row.names = F)

# AlphaAssign ####
pedigree_file <- read.csv("worker_pedigree.csv")

Alpha_pedigree <- data.frame(id = pedigree_file$id,
                             sire = rep(0, length(pedigree_file$id)),
                             dam = pedigree_file$mother)
write.table(Alpha_pedigree, file = "Pedigree.txt", sep = " ", quote = F, col.names = F, row.names = F)

Potential_fathers <- data.frame(id = pedigree_file$id,
                                Dpc1 = rep(unique(pedigree_file$dpc)[1], length(pedigree_file$dpc)),
                                Dpc2 = rep(unique(pedigree_file$dpc)[2], length(pedigree_file$dpc)),
                                Dpc3 = rep(unique(pedigree_file$dpc)[3], length(pedigree_file$dpc)),
                                Dpc4 = rep(unique(pedigree_file$dpc)[4], length(pedigree_file$dpc)))
write.table(Potential_fathers, file = "PotentialFathers.list", sep = " ",  quote = F, col.names = F, row.names = F)


#format the ped/map files into AlphaAssign format -using recode A 
#If there are genotyping errors (NA) you need to change them to 9 
rm(list = ls())
AlphaPed <- read.table("SNP2_phasedGE_recodeA.raw", header = T)
AlphaGeno <- AlphaPed[,7:ncol(AlphaPed)]

AlphaGeno[is.na(AlphaGeno)] <- 9
AlphaGeno_id <- cbind(AlphaPed$IID, AlphaGeno)
write.table(AlphaGeno_id, file = "AlphaGenoSNP2_phasedGE.txt", sep = " ", quote = F, col.names = F, row.names = F)


#Create AlphaAssign summary 
Known_Dpc <- read.csv("Known_Dpc.csv")

Alpha_file <- read.table("Alpha_SNP4_phasedGE_NonNested.sires", header = T)

Sires_assigned <- Alpha_file[Alpha_file$chosen == 1, ]
nSires_assigned <- nrow(Sires_assigned)

Pairwise <- Sires_assigned[, c(1,2)]
colnames(Pairwise) <- c("id", "sire")
merged_df2 <- merge(Pairwise, Known_Dpc, by="id", suffixes=c("_pairwise", "_known"))

# Count the number of matches and mismatches
nCorrect_sires <- sum(merged_df2$sire_pairwise == merged_df2$sire_known)


SNP4_phasedGE <- data.frame(Data_Group = "Non_Nested",
                           Test = "Phased_GenoErr",
                           nOffspring = 240,
                           SNP_group = 4,
                           nSires_assigned = nSires_assigned,
                           nCorrect_sires = nCorrect_sires,
                           Software = "AlphaAssign")

Alpha_summary <- rbind(SNP1_nGE, SNP1_GE,
                       SNP2_nGE, SNP2_GE, SNP2_phasedGE,
                       SNP3_nGE, SNP3_GE, SNP3_phasedGE,
                       SNP4_nGE, SNP4_GE, SNP4_phasedGE)

write.table(Alpha_summary, file = "Alpha_summary.txt", sep = " ", quote = F, col.names = T, row.names = F)


#Colony ####
#Set up the colony files to be put into the input file 
pedigree_file <- read.csv("worker_pedigree.csv")
#Known Mothers 
Known_Mothers <- data.frame(worker_id = pedigree_file$id,
                            mother_id = pedigree_file$mother)

write.table(Known_Mothers, file = "Known_Mothers.txt", sep = " ", quote = F, col.names = F, row.names = F)

#Excluded mothers 
create_excluded_mothers <- function(Known_Mothers) {
  
  # Step 1: Extract unique worker_ids and mother_ids
  unique_worker_ids <- unique(Known_Mothers$worker_id)
  unique_mother_ids <- unique(Known_Mothers$mother_id)
  
  # Initialize an empty list to store the result
  result_list <- list()
  
  # Step 2: For each worker_id, find the non-mothers and store in result_list
  for (worker in unique_worker_ids) {
    non_mothers <- unique_mother_ids[!(unique_mother_ids %in% Known_Mothers$mother_id[Known_Mothers$worker_id == worker])]
    total_excluded <- length(non_mothers)
    result_list[[as.character(worker)]] <- c(worker, total_excluded, non_mothers)
  }
  
  # Convert the result_list to a data frame
  result_df <- do.call(rbind, result_list)
  result_df <- as.data.frame(result_df)
  
  # Set column names
  colnames(result_df) <- c("worker_id", "total_excluded", paste0("non_mother", 1:(ncol(result_df)-2)))
  
  return(result_df)
}


excluded_mother <- create_excluded_mothers(Known_Mothers)

write.table(excluded_mother, file = "Excluded_mothers.txt", sep = " ", quote = F, col.names = F, row.names = F)


#Excluded siblings 
# Load necessary libraries

# Define the function for excluded siblings
create_excluded_siblings <- function(Known_Mothers) {
  
  # Step 1: Extract unique worker_ids
  unique_worker_ids <- unique(Known_Mothers$worker_id)
  
  # Initialize an empty list to store the result
  result_list <- list()
  
  # Step 2: For each worker_id, find the non-siblings and store in result_list
  for (worker in unique_worker_ids) {
    # Find the mother_id of the current worker
    worker_mother_id <- Known_Mothers$mother_id[Known_Mothers$worker_id == worker]
    
    # Find all workers who do not share the same mother_id
    non_siblings <- Known_Mothers$worker_id[Known_Mothers$mother_id != worker_mother_id]
    
    # Calculate the total number of non-siblings
    total_siblings <- length(non_siblings)
    
    # Store the result in the list
    result_list[[as.character(worker)]] <- c(worker, total_siblings, non_siblings)
  }
  
  # Convert the result_list to a data frame
  result_df <- do.call(rbind, result_list)
  result_df <- as.data.frame(result_df)
  
  # Set column names
  colnames(result_df) <- c("worker_id", "total_siblings", paste0("non_sibling", 1:(ncol(result_df)-2)))
  
  return(result_df)
}

excluded_siblings <- create_excluded_siblings(Known_Mothers)

write.table(excluded_siblings, file = "Excluded_siblings.txt", sep = " ", quote = F, col.names = F, row.names = F)


#Worker genotypes
rm(list = ls())
Ped <- read.table("SNP2_phasedGE_recode12.ped")
ColonyGeno <- Ped[,-c(1,3,4,5,6)] 

WorkerGeno <- ColonyGeno[-c(1:12), ]
write.table(WorkerGeno, file = "Col_workerGeno_SNP2_phasedGE.txt", sep = " ", quote = F, col.names = F, row.names = F)

#Mother genotypes 
MotherGeno <- ColonyGeno[c(1:8), ]
write.table(MotherGeno, file = "Col_motherGeno_SNP2_phasedGE.txt", sep = " ", quote = F, col.names = F, row.names = F)


#Dpc genotypes 
DpcGeno <- ColonyGeno[c(9:12), ]
write.table(DpcGeno, file = "Col_dpcGeno_SNP2_phasedGE.txt", sep = " ", quote = F, col.names = F, row.names = F)



#Create a summary file of the Colony outputs for the Pairwise Paternity file
setwd("~/Desktop/Slovenia data/Attempt2/Nested/Colony/PairwisePaternity_output")
Nest_SNP1_nGE <- read.csv("NestedSNP1_nGE.PairwisePaternity")
Nest_SNP1_GE <- read.csv("NestedSNP1_GE.PairwisePaternity")
Nest_SNP2_nGE <- read.csv("NestedSNP2_nGE.PairwisePaternity")
Nest_SNP2_GE <- read.csv("NestSNP2_GE.PairwisePaternity")
Nest_SNP2_phasedGE <- read.csv("NestedSNP2_phasedGE.PairwisePaternity")
Nest_SNP3_nGE <- read.csv("NestedSNP3_nGE.PairwisePaternity")
Nest_SNP3_GE <- read.csv("NestedSNP3_GE.PairwisePaternity")
Nest_SNP3_phasedGE <- read.csv("NestedSNP3_phasedGE.PairwisePaternity")
Nest_SNP4_nGE <- read.csv("NestedSNP4_nGE.PairwisePaternity")
Nest_SNP4_GE <- read.csv("NestedSNP4_GE.PairwisePaternity")
Nest_SNP4_phasedGE <- read.csv("NestedSNP4_phasedGE.PairwisePaternity")

Nest_SNP1_nGE$Test <- rep("No_GenoErr")
Nest_SNP2_nGE$Test <- rep("No_GenoErr")
Nest_SNP3_nGE$Test <- rep("No_GenoErr")
Nest_SNP4_nGE$Test <- rep("No_GenoErr")

Nest_SNP1_GE$Test <- rep("GenoErr")
Nest_SNP2_GE$Test <- rep("GenoErr")
Nest_SNP3_GE$Test <- rep("GenoErr")
Nest_SNP4_GE$Test <- rep("GenoErr")

Nest_SNP2_phasedGE$Test <- rep("Phased_GenoErr")
Nest_SNP3_phasedGE$Test <- rep("Phased_GenoErr")
Nest_SNP4_phasedGE$Test <- rep("Phased_GenoErr")

Nest_SNP1_GE$Data_Group <- rep(1)
Nest_SNP1_nGE$Data_Group <-rep(1)
Nest_SNP2_nGE$Data_Group <- rep(2)
Nest_SNP2_GE$Data_Group <- rep(2)
Nest_SNP2_phasedGE$Data_Group <- rep(2)
Nest_SNP3_nGE$Data_Group <- rep(3)
Nest_SNP3_GE$Data_Group <- rep(3)
Nest_SNP3_phasedGE$Data_Group <- rep(3)
Nest_SNP4_nGE$Data_Group <- rep(4)
Nest_SNP4_GE$Data_Group <- rep(4)
Nest_SNP4_phasedGE$Data_Group <- rep(4)


SNP1 <- rbind(Nest_SNP1_nGE, Nest_SNP1_GE)
SNP2 <- rbind(Nest_SNP2_GE, Nest_SNP2_nGE, Nest_SNP2_phasedGE)
SNP3 <- rbind( Nest_SNP3_nGE, Nest_SNP3_GE, Nest_SNP3_phasedGE)
SNP4 <- rbind(Nest_SNP4_nGE, Nest_SNP4_GE, Nest_SNP4_phasedGE)


Known_Dpc <- read_csv("Known_Dpc.csv")
colnames(Known_Dpc) <- c("OffspringID", "Mother", "CandidateID")
merged_df2 <- merge(SNP1, Known_Dpc, by="OffspringID", suffixes=c("_sim", "_known"))
merged_df2$CorrectSires <- ifelse(merged_df2$CandidateID_sim == merged_df2$CandidateID_known, yes = TRUE, no = FALSE)

SNP1 <- merged_df2


ColonySummary <- rbind(SNP1, SNP2, SNP3, SNP4)

write.table(ColonySummary, file = "Colony_summary_NESTED_confidence.txt", quote = F, sep = " ", col.names = T, row.names = F )


Sires_assigned <- Nest_SNP1_phasedGE$CandidateID
nSires_assigned <- length(Sires_assigned)

Pairwise <- Nest_SNP1_phasedGE[, c(1,2)]
colnames(Pairwise) <- c("id", "sire")
merged_df2 <- merge(Pairwise, Known_Dpc, by="id", suffixes=c("_pairwise", "_known"))

# Count the number of matches and mismatches
nCorrect_sires <- sum(merged_df2$sire_pairwise == merged_df2$sire_known)


SNP1_phasedGE <- data.frame(Data_Group = "Non_Nested",
                      Test = "Phased_GenoErr",
                      SNP_group = 1,
                      nSires_assigned = nSires_assigned,
                      nCorrect_sires = nCorrect_sires,
                      Software = "Colony")

Colony_group_summary <- rbind(SNP1_nGE, SNP1_GE,
                              SNP2_nGE, SNP2_GE, SNP2_phasedGE,
                              SNP3_phasedGE, SNP3_GE, SNP3_nGE,
                              SNP4_nGE, SNP4_GE, SNP4_phasedGE)

write.table(Colony_group_summary, file = "Colony_nSires_summary_NON.txt", quote = F, sep = " ", col.names = T, row.names = F)









#KING #####
#Create a summary file for KING 
#Create a summary file of the Colony outputs for the Pairwise Paternity file
setwd("~/Desktop/Slovenia data/Attempt2/Non-Nested/KINGsummary")
Nest_SNP1_nGE <- read.table("KING_SNP1_nGE.kin0")
Nest_SNP1_GE <- read.table("KING_SNP1_GE.kin0")
Nest_SNP2_nGE <- read.table("KING_SNP2_nGE.kin0")
Nest_SNP2_GE <- read.table("KING_SNP2_GE.kin0")
Nest_SNP2_phasedGE <- read.table("KING_SNP2_phasedGE.kin0")
Nest_SNP3_nGE <- read.table("KING_SNP3_nGE.kin0")
Nest_SNP3_GE <- read.table("KING_SNP3_GE.kin0")
Nest_SNP3_phasedGE <- read.table("KING_SNP3_phasedGE.kin0")
Nest_SNP4_nGE <- read.table("KING_SNP4_nGE.kin0")
Nest_SNP4_GE <- read.table("KING_SNP4_GE.kin0")
Nest_SNP4_phasedGE <- read.table("KING_SNP4_phasedGE.kin0")

Nest_SNP1_nGE$Test <- rep("No_GenoErr")
Nest_SNP2_nGE$Test <- rep("No_GenoErr")
Nest_SNP3_nGE$Test <- rep("No_GenoErr")
Nest_SNP4_nGE$Test <- rep("No_GenoErr")

Nest_SNP1_GE$Test <- rep("GenoErr")
Nest_SNP2_GE$Test <- rep("GenoErr")
Nest_SNP3_GE$Test <- rep("GenoErr")
Nest_SNP4_GE$Test <- rep("GenoErr")

Nest_SNP2_phasedGE$Test <- rep("Phased_GenoErr")
Nest_SNP3_phasedGE$Test <- rep("Phased_GenoErr")
Nest_SNP4_phasedGE$Test <- rep("Phased_GenoErr")

Nest_SNP1_GE$SNP_group <- rep(1)
Nest_SNP1_nGE$SNP_group <-rep(1)
Nest_SNP2_nGE$SNP_group <- rep(2)
Nest_SNP2_GE$SNP_group <- rep(2)
Nest_SNP2_phasedGE$SNP_group <- rep(2)
Nest_SNP3_nGE$SNP_group <- rep(3)
Nest_SNP3_GE$SNP_group <- rep(3)
Nest_SNP3_phasedGE$SNP_group <- rep(3)
Nest_SNP4_nGE$SNP_group <- rep(4)
Nest_SNP4_GE$SNP_group <- rep(4)
Nest_SNP4_phasedGE$SNP_group <- rep(4)


SNP1 <- rbind(Nest_SNP1_nGE, Nest_SNP1_GE)
SNP2 <- rbind(Nest_SNP2_GE, Nest_SNP2_nGE, Nest_SNP2_phasedGE)
SNP3 <- rbind( Nest_SNP3_nGE, Nest_SNP3_GE, Nest_SNP3_phasedGE)
SNP4 <- rbind(Nest_SNP4_nGE, Nest_SNP4_GE, Nest_SNP4_phasedGE)

KING_summary <- rbind(SNP1, SNP2, SNP3, SNP4)
colnames(KING_summary) <- c("FID1","IID1",	"FID2",	"IID2",	"NSNP",	"HETHET",	"IBS0",	"KINSHIP", "Test", "SNP_group")
KING_summary$Data_Group <- rep("Non_Nested")

write.table(KING_summary, file = "KING_summary.txt", quote = F, sep = " ", col.names = T, row.names = F)


#Now lets have a look to see if we can see any fathers 
# IBS0: The value in this column indicates the fraction of genetic markers where the two individuals share no alleles. 
# This can be used to infer the degree of relatedness between individuals.
# For instance, a low IBS0 value typically indicates a close relationship (such as parent-child or full siblings),
# while a high IBS0 value suggests a more distant relationship or even unrelated individuals.

#The KINSHIP coefficient is a measure of the genetic relatedness or the proportion of alleles shared identical-by-descent (IBD) between two individuals.
# It estimates the fraction of the genome where two individuals are expected to share alleles inherited from a common ancestor.
# KINSHIP values range from 0 (unrelated individuals) to 0.5 (identical twins). For example:
# Parent-Offspring: KINSHIP ≈ 0.25
# Full Siblings: KINSHIP ≈ 0.25
# Half Siblings: KINSHIP ≈ 0.125
# Unrelated Individuals: KINSHIP ≈ 0

#So we're going to plot IBS0 against KINSHIP to see if we can spot any father clumping 

library(ggplot2)

plot_ibs0_kinship <- function(data) {
  # Check if the necessary columns are present in the dataframe
  if (!all(c("IBS0", "KINSHIP") %in% colnames(data))) {
    stop("Dataframe must contain 'IBS0' and 'KINSHIP' columns.")
  }
  
  # Create the scatter plot using ggplot2
  p <- ggplot(data, aes(x = KINSHIP, y = IBS0)) +
    geom_point() +
    labs(title = "Scatter Plot of IBS0 vs KINSHIP",
         x = "KINSHIP",
         y = "IBS0") +
    theme_minimal() +
    scale_y_continuous(expand = c(0, 0)) # Ensure the y-axis starts at 0
  
  # Print the plot
  print(p)
}

plot_ibs0_kinship(KING_summary)

# Select those with a Kinship >0.2 to more accuratetly find the fathers 
KING_0.2 <- KING_summary[KING_summary$KINSHIP >= 0.2,]
#Also take only the ones where dpc has been assigned 
KING_sires <- KING_0.2[KING_0.2$IID2 %in% c("61","62","63","64"), ]
KING_sires <- KING_sires[!KING_sires$IID1 %in% c("61","62","63","64", "3", "4", "5", "6","7","8","9","10"), ]

plot_ibs0_kinship(KING_sires)

write.table(KING_sires, file = "Filtered_KING_sires_NON.txt", quote = F, sep = " ", col.names = T, row.names = F)


#Check if the sires are correct
Known_Dpc <- read.csv("Known_Dpc.csv")

SNP <- Nest_KING_summary[Nest_KING_summary$SNP_group == 3 & Nest_KING_summary$Test == "Non_GenoErr", ]


nSires_assigned <- nrow(SNP)

Pairwise <- SNP[, c(2,4)]
colnames(Pairwise) <- c("id", "sire")
merged_df2 <- merge(Pairwise, Known_Dpc, by="id", suffixes=c("_pairwise", "_known"))

SNP3_nGE_kinship <- cbind(merged_df2[,c(1,2,4)], SNP[,c(7:11)])
SNP3_nGE_kinship$Match <- ifelse(SNP3_nGE_kinship$sire_pairwise == SNP3_nGE_kinship$sire_known, yes = TRUE, no = FALSE)


# Count the number of matches and mismatches
nCorrect_sires <- sum(merged_df2$sire_pairwise == merged_df2$sire_known)


SNP4_nGE <- data.frame(Data_Group = "Non_Nested",
                            Test = "No_GenoErr",
                            nOffspring = 240,
                            SNP_group = 4,
                            nSires_assigned = nSires_assigned,
                            nCorrect_sires = nCorrect_sires,
                            Software = "KING")

rm(SNP, merged_df2, Pairwise, nSires_assigned, nCorrect_sires)



KING_summary <- rbind(SNP1_nGE_kinship, SNP1_GE_kinship,
                       SNP2_nGE_kinship, SNP2_GE_kinship, SNP2_phasedGE_kinship,
                       SNP3_nGE_kinship, SNP3_GE_kinship, SNP3_phasedGE_kinship,
                       SNP4_nGE_kinship, SNP4_GE_kinship, SNP4_phasedGE_kinship)

write.table(KING_summary, file = "KING_kinship_NEST_summary.txt", sep = " ", quote = F, col.names = T, row.names = F)




#when the slov fathers get assigned via a software assignment, update the pedigree
# Assuming your dataframes are named AlphaSires and Slov_pedigree
# Load necessary library
library(dplyr)
library(tidyr)

# Merge the dataframes on the ID column
merged_data <- Slov_pedigree %>%
  left_join(AlphaSires, by = "ID", suffix = c(".old", ".new"))

merged_data <- merged_data[,-c(2)]
merged_data <- merged_data[,c(1,3,2)]
merged_data[is.na(merged_data)] <- 0
merged_data$FamilyID <- rep("APIS")
merged_data <- merged_data[,c(4,1,2,3)]
colnames(merged_data) <- c("FamilyID", "ID", "SIRE", "DAM")
# View the updated Slov_pedigree
print(Slov_pedigree)




#####genotyping errors need to be added to the ped file accurately so that the mendelian haplotype stuff works correctly 

# Load necessary library
library(data.table)

# Parameters for genotyping errors
error  = 0.0001
error2 = 0.00001
missing = 0.005
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




filter_king <- function(KING_file){
  #filter Kinship >= 0.2 and IBS < 0.005
  KING_file <- KING_file[KING_file$KINSHIP >= 0.2 & KING_file$IBS0 < 0.005, ]
  return(KING_file)
}





#Check the kinship coefficients

check_king <- function(KING_file, pedigree){

  check_king <- list()
  
  for (i in 1:nrow(pedigree)){
  
  dpc_id <- pedigree$dpc[i]
  worker_id <- pedigree$id[i]
  
  Known_pair <- KING_file[KING_file[,2] == worker_id & KING_file[,4] == dpc_id,]
  
  check_king[[i]] <- Known_pair
  }

  check_king <- do.call(rbind,check_king)
  tmp <- nrow(check_king)
return(tmp)
}


check_king_wrong <- function(KING_file, pedigree){
  
  check_king <- list()
  
  for (i in 1:nrow(pedigree)){
    
    dpc_id <- pedigree$dpc[i]
    worker_id <- pedigree$id[i]
    
    Known_pair <- KING_file[KING_file[,2] == worker_id & KING_file[,4] != dpc_id,]
    
    check_king[[i]] <- Known_pair
  }
  
  check_king <- do.call(rbind,check_king)
  tmp <- nrow(check_king)
  return(tmp)
}

SNP3_GE_check_wrong <- check_king(KING_file = KING_SNP3_GE, pedigree = worker_pedigree)



