#################################################################################

# Prepare file for Colony pedigree reconstruction software

#################################################################################
rm(list = ls())


args = commandArgs(trailingOnly=TRUE)
Rep = args[1]
workingDir = args[2]
softwareDir = args[3]
pathToPlink <- softwareDir
pathToBeagle <- softwareDir

repDir = paste0(workingDir, "/SimRep", Rep, "/")

# Set working directory
setwd(repDir)

library(Eagle)
library(tidyr)
library(dplyr)

#Excluded siblings - For each worker_id, find the non-siblings and store in result_list
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

#Excluded mothers - For each worker_id, find the non-mothers and store in result_list
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

dir.create("Data/Colony", showWarnings = F)
#****** 1. SIMULATED *******************************************************

#Setting up the input files *****************
pedigree_file <- read.csv("Data/worker_pedigree.csv")


#Genotypes 
# •• NO GE ••

worker_ids <- unique(pedigree_file$id)
queen_ids <- unique(pedigree_file$mother)
dpc_ids <- unique(pedigree_file$dpc)


if (method == "NoGE") {
  Genotypes <- read.table(paste0("Data/Sim_NoGE/SNP_",n,"_NoGE_QC.ped"))
  # If a queen has been removed in the QC, remove the workers as well
  removedQueen = unique(pedigree_file$mother[!pedigree_file$mother %in% Genotypes$V2])
  offspring_to_remove = pedigree_file$id[pedigree_file$mother %in% removedQueen]
  Genotypes <- Genotypes[!Genotypes$V2 %in% offspring_to_remove,]

  Worker_genotypes <- Genotypes[Genotypes$V2 %in% worker_ids,]
  Worker_genotypes <- Worker_genotypes[, c(2, 7:ncol(Worker_genotypes))]
  Worker_genotypes[Worker_genotypes == 0] <- -9

  Queen_genotypes <- Genotypes[Genotypes$V2 %in% queen_ids,]
  Queen_genotypes <- Queen_genotypes[, c(2, 7:ncol(Queen_genotypes))]
  Queen_genotypes[Queen_genotypes == 0] <- -9


  Dpc_genotypes <- Genotypes[Genotypes$V2 %in% dpc_ids,]
  Dpc_genotypes <- Dpc_genotypes[, c(2, 7:ncol(Dpc_genotypes))]
  Dpc_genotypes[Dpc_genotypes == 0] <- -9

}

#•• WITH GE ••
if (method == "WithGE")  {
  Genotypes <- read.table(paste0("Data/Sim_WithGE/SNP_",n,"_WithGE_QC.ped"))
  # If a queen has been removed in the QC, remove the workers as well
  removedQueen = unique(pedigree_file$mother[!pedigree_file$mother %in% Genotypes$V2])
  offspring_to_remove = pedigree_file$id[pedigree_file$mother %in% removedQueen]
  Genotypes <- Genotypes[!Genotypes$V2 %in% offspring_to_remove,]

  Worker_genotypes <- Genotypes[Genotypes$V2 %in% worker_ids,]
  Worker_genotypes <- Worker_genotypes[, c(2, 7:ncol(Worker_genotypes))]


  Queen_genotypes <- Genotypes[Genotypes$V2 %in% queen_ids,]
  Queen_genotypes <- Queen_genotypes[, c(2, 7:ncol(Queen_genotypes))]


  Dpc_genotypes <- Genotypes[Genotypes$V2 %in% dpc_ids,]
  Dpc_genotypes <- Dpc_genotypes[, c(2, 7:ncol(Dpc_genotypes))]

}

#Known Mothers 
Known_Mothers <- data.frame(worker_id = pedigree_file$id[pedigree_file$id %in% Genotypes$V2],
                            mother_id = pedigree_file$mother[pedigree_file$mother %in% Genotypes$V2])
excluded_mother <- create_excluded_mothers(Known_Mothers)
excluded_siblings <- create_excluded_siblings(Known_Mothers)




nLoci = (ncol(Worker_genotypes)-1) / 2



sink(paste0("Data/Colony/colony2_SNP", n, "_", method, ".dat"))
cat(paste0("'", method, "_SNP", n, "'   !Dataset name"))
cat("\n")
cat(paste0("'", method, "_SNP", n, "'   !Output file name"))
cat("\n")
cat(paste0(nrow(Worker_genotypes), "         ! Number of offspring in the sample"))
cat("\n")
cat(paste0(nLoci, "        ! Number of loci / same as genotype markers? (make sure same order is everywhere)"))

cat("
123         ! Seed for random number generator (default used)
1           ! 0/1=Not updating/updating allele frequency
2           ! 2/1=Dioecious/Monoecious species
0           ! 0/1=No inbreeding/inbreeding
0           ! 0/1=Diploid species/HaploDiploid species
0  0        ! 0/1=Polygamy/Monogamy for males & females
0           ! 0/1=Clone inference =No/Yes
1           ! 0/1=Full sibship size scaling =No/Yes
1 1.0 1.0   ! 0,1,2,3=No,weak,medium,strong sibship size prior; mean paternal; materal sibship size (defaults used)
0           ! 0/1=Unknown/Known population allele frequency
1           ! Number of runs
2           ! Length of run
0           ! 0/1=Monitor method by Iterate#/Time in second (default)
100000      ! Monitor interval in Iterate# / in seconds (default)
0           ! non-Windows version
1           ! Analysis method = Full-likelihood
1           ! 0/1/2/3=low/medium/high/v high Precision for Full-likelihood
MK@          !Marker names 
0@           !Marker types, 0/1 = codominant/dominant
0.1@      !Allelic dropout rate (default)
0.1@      !false allele rate (default)")
cat("\n")
cat("\n")
write.table(Worker_genotypes,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

cat("
1.0 1.0     !prob. of dad/mum included in the candidates (we definitely know mothers but not all fathers)")
cat("\n")
cat("\n")
cat(paste0(
nrow(Dpc_genotypes), "   ", nrow(Queen_genotypes), "       !numbers of candidate males & females (this assumes founders which are unrelated and are either parents or unrelated to offspring"))
cat("\n")
write.table(Dpc_genotypes,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
cat("\n")

write.table(Queen_genotypes,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
cat("\n")
cat(" ")
cat("\n")
cat("
0   0        !number of known father-offspring; paternity exclusion threshold
")
cat("\n")
cat(
  paste0(nrow(Worker_genotypes), " 0       !number of known mother-offspring, maternity exclusion threshold (If known maternity is larger than zero, the on each row list a single relationship containing the offpsring ID followed by it's known mother ID)
"))
cat("\n")
write.table(Known_Mothers,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
cat("\n")
cat("
0            !#known paternal sibship with unknown fathers (we don't know sibling type if fathers are same) 

0            !#known maternal sibship with unknown mothers (no mothers are unknown)

0            !#known paternity exclusions 

")
cat("\n")
cat(paste0(
  nrow(Worker_genotypes), "          !#known maternity exclusions (some candiate females are definitely impossible to be the mother of  a specific offspring) (Offspring ID followed by known NOT mothers' id)
"))
write.table(excluded_mother,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
cat("\n")
cat("
0           !#known paternal sibship exclusions

")
cat("\n")
cat(paste0(
  nrow(Worker_genotypes), "           !#known maternal sibship exclusions (provide the number of offspring that are known to each has at least one excluded offspring as its maternal sibling) (Exculded this for now becuase it would be suuuupppper long (SEE IF  THIS DATA file works first))
"))
write.table(excluded_siblings,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
sink()


cmd = paste0(softwareDir, "/colony/colony2s.ifort.out IFN:Data/Colony/colony2_SNP", n, "_", method, ".dat")
system(cmd)


  