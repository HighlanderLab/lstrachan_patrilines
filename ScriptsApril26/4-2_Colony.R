#################################################################################

# Prepare file for Colony pedigree reconstruction software

#################################################################################
rm(list = ls())

library(Eagle)
library(tidyr)
library(dplyr)

#****** 1. SIMULATED *******************************************************

pathToPlink <- "~/Desktop/PLINK/./"
workingDir = "~/Desktop/lstrachan_patrilines"
setwd(workingDir)

#Setting up the input files *****************

pedigree_file <- read.csv("Data/worker_pedigree.csv")
#Known Mothers 
Known_Mothers <- data.frame(worker_id = pedigree_file$id,
                            mother_id = pedigree_file$mother)
write.table(Known_Mothers, file = "Known_Mothers.txt", sep = " ", quote = F, col.names = F, row.names = F)

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

excluded_mother <- create_excluded_mothers(Known_Mothers)

write.table(excluded_mother, file = "Excluded_mothers.txt", sep = " ", quote = F, col.names = F, row.names = F)

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

excluded_siblings <- create_excluded_siblings(Known_Mothers)

write.table(excluded_siblings, file = "Excluded_siblings.txt", sep = " ", quote = F, col.names = F, row.names = F)

#Genotypes 
# •• NO GE ••
for (n in 1:5){
Ped <- read.table("SNP",n,"_Genotypes.ped") # TODO : Add in correct file name
Genotypes <- read.table("Slov_fM_AC_QC.ped")
worker_ids <- unique(pedigree_file$id)
queen_ids <- unique(pedigree_file$mother)
dpc_ids <- unique(pedigree_file$dpc)

Worker_genotypes <- Genotypes[Genotypes$V2 %in% worker_ids,]
write.table(WorkerGeno, file = "Col_workerGeno_SNP",n,"_NoGE.txt", sep = " ", quote = F, col.names = F, row.names = F)

Queen_genotypes <- Genotypes[Genotypes$V2 %in% queen_ids,]
write.table(Queen_genotypes, file = "Col_queenGeno_SNP",n,"_NoGE.txt", sep = " ", quote = F, col.names = F, row.names = F)

Dpc_genotypes <- Genotypes[Genotypes$V2 %in% dpc_ids,]
write.table(Queen_genotypes, file = "Col_dpcGeno_SNP",n,"_NoGE.txt", sep = " ", quote = F, col.names = F, row.names = F)

}
#•• WITH GE ••
for (n in 1:5){
  Genotypes <- read.table("ADD PED FILE NAME HERE")
  worker_ids <- unique(pedigree_file$id)
  queen_ids <- unique(pedigree_file$mother)
  dpc_ids <- unique(pedigree_file$dpc)
  
  Worker_genotypes <- Genotypes[Genotypes$V2 %in% worker_ids,]
  write.table(WorkerGeno, file = "Col_workerGeno_SNP",n,"_WithGE.txt", sep = " ", quote = F, col.names = F, row.names = F)
  
  Queen_genotypes <- Genotypes[Genotypes$V2 %in% queen_ids,]
  write.table(Queen_genotypes, file = "Col_queenGeno_SNP",n,"_WithGE.txt", sep = " ", quote = F, col.names = F, row.names = F)
  
  Dpc_genotypes <- Genotypes[Genotypes$V2 %in% dpc_ids,]
  write.table(Queen_genotypes, file = "Col_dpcGeno_SNP",n,"_WithGE.txt", sep = " ", quote = F, col.names = F, row.names = F)
  
}



#****** 2. REAL DATA *******************************************************
pathToPlink <- "~/Desktop/PLINK/./"
workingDir = "~/Desktop/lstrachan_patrilines"
setwd(workingDir)

pedigree_file <- read.csv("Data/Real_data/Real_Data_pedigree.csv")
colnames(pedigree_file) <- c("id", "sire", "dam")

#Known Mothers 
Known_Mothers <- data.frame(worker_id = pedigree_file$id,
                            mother_id = pedigree_file$dam)

write.table(Known_Mothers, file = "Known_Mothers.txt", sep = " ", quote = F, col.names = F, row.names = F)

excluded_mother <- create_excluded_mothers(Known_Mothers)
write.table(excluded_mother, file = "Excluded_mothers.txt", sep = " ", quote = F, col.names = F, row.names = F)


excluded_siblings <- create_excluded_siblings(Known_Mothers)
write.table(excluded_siblings, file = "Excluded_siblings.txt", sep = " ", quote = F, col.names = F, row.names = F)

queen_ids <- unique(pedigree_file$dam)
all_ids <- unique(pedigree_file$id)
workers_ids <- setdiff(all_ids, queen_ids)

Ped <- read.table("ADD PED NAME HERE")

Worker_genotypes <- Genotypes[Genotypes$V2 %in% worker_ids,]
write.table(WorkerGeno, file = "Col_workerGeno_SNP",n,"_WithGE.txt", sep = " ", quote = F, col.names = F, row.names = F)

Queen_genotypes <- Genotypes[Genotypes$V2 %in% queen_ids,]
write.table(Queen_genotypes, file = "Col_queenGeno_SNP",n,"_WithGE.txt", sep = " ", quote = F, col.names = F, row.names = F)

