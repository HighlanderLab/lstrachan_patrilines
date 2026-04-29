


#******************************************************************************
# Prepare files for AlphaAssign pedigree reconstruction software
#******************************************************************************

#                     1. SIMULATED DATA 
#******************************************************************************

##############.  Prepping the input files ##################
rm(list = ls()) #Lets clear workspace just incase something sneaks in

args = commandArgs(trailingOnly=TRUE)
workingDir = args[1]
softwareDir = args[2]
pathToPlink <- softwareDir

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



#******************************************************************************
#                     2. REAL DATA
#******************************************************************************

##############.  Prepping the input files ##################
#rm(list = ls()) #Lets clear workspace just incase something sneaks in 

setwd(paste0(workingDir, "/Real_data/"))
dir.create("Data/AlphaAssign", showWarnings = FALSE)

pedigree_file_real <- read.csv("Data/Real_Data_pedigree.csv")
colnames(pedigree_file_real) <- c("id", "sire", "dam")

removed_queen_id <- "B0003710" #Id of the queen removed during QC 
pedigree_file_real <- pedigree_file_real[pedigree_file_real$id != removed_queen_id,]
pedigree_file_real <- pedigree_file_real[pedigree_file_real$dam != removed_queen_id,] #remove all of her offspring from the pedigree too


Alpha_pedigree_real <- data.frame(id = pedigree_file_real$id,
                             sire = pedigree_file_real$sire,
                             dam = pedigree_file_real$dam)


write.table(Alpha_pedigree_real, file = "Data/AlphaAssign/Pedigree.txt", sep = " ", quote = F, col.names = F, row.names = F)

SNP_samples <- read.csv("Data/SNP_samples_2022.csv", header= T)
workers_IDs_beforeQC <- SNP_samples$snp_id[SNP_samples$biotype == "worker"]
SNP_samples <- SNP_samples[SNP_samples$biotype == "dpc", ]
SNP_samples <- SNP_samples$snp_id

n <- nrow(pedigree_file_real)

Potential_fathers <- data.frame(
  id = pedigree_file_real$id,
  Dpc1 = rep(SNP_samples[1], n),
  Dpc2 = rep(SNP_samples[2], n),
  Dpc3 = rep(SNP_samples[3], n),
  Dpc4 = rep(SNP_samples[4], n),
  Dpc5 = rep(SNP_samples[5], n)
)
write.table(Potential_fathers, file = "Data/AlphaAssign/PotentialFathers.list", sep = " ",  quote = F, col.names = F, row.names = F)


#Recode quality controlled files for AlphaAssign
ped_to_raw(ped_file = "Data/Slov_fM_QC.ped", map_file = "Data/Slov_fM_QC.map", output_file = "Data/Slov_fM_QC.raw")

#Process the genotypes
AlphaPed <- read.table("Data/Slov_fM_QC.raw", header=TRUE)
AlphaGeno <- AlphaPed[,7:ncol(AlphaPed)]; AlphaGeno[is.na(AlphaGeno)] <- 9
AlphaGeno_id <- cbind(AlphaPed$IID, AlphaGeno)

write.table(AlphaGeno_id, file=paste0("Data/AlphaAssign/AlphaGeno_RealData.txt"), sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE)

############### Running AlphAssign in terminal################

#create a for loop here for all of the SNP array sizes <-------------- 
dir.create("Outputs/AlphaAssign", showWarnings = FALSE)
system(paste0("bash ", workingDir, "/ScriptsApril26/RunAlphaAssign_RealData.sh AlphaGeno_RealData ", workingDir, "/Real_data/Data/AlphaAssign ", workingDir, "/Real_data/Outputs/AlphaAssign"))

#Summarise the dataset     
Alpha_output <- read.table("Outputs/AlphaAssign/AlphaGeno_RealData.sires", header = TRUE)

Sires_assigned <- Alpha_output[Alpha_output$chosen == 1, ]
nSires_assigned <- nrow(Sires_assigned)

Offspring_and_candidateParent <- Sires_assigned[, c(1,2)]
colnames(Offspring_and_candidateParent) <- c("id", "sire")

    
Sires_assigned <- Alpha_output[Alpha_output$chosen == 1, ]
nSires_assigned <- nrow(Sires_assigned)

Offspring_and_candidateParent <- Sires_assigned[, c(1,2)]
colnames(Offspring_and_candidateParent) <- c("id", "sire")

Worker_known <- sum(AlphaPed$IID %in% workers_IDs_beforeQC)
Real_Alpha_df <- data.frame(
  Test = "Real",
  nOffspring = Worker_known,
  SNP_group = "1.7k",
  nSires_assigned = nSires_assigned,
  nCorrect_sires = NA,
  Software = "AlphaAssign")


write.table(Real_Alpha_df, file = paste0(workingDir, "/Real_data/Outputs/AlphaAssign/Alpha_summary.txt"), sep = " ", quote = F, col.names = T, row.names = F)

##################################################
# Update the pedigree
# Match ids to get corresponding sires from WithGE_50k_pedigree
idx <- match(Alpha_pedigree_real$id, Offspring_and_candidateParent$id)
# Replace only where sire == 0 and a match exists
update_sires <- Alpha_pedigree_real$sire == 0 & !is.na(idx)
Alpha_pedigree_real$sire[update_sires] <- Offspring_and_candidateParent$sire[idx[update_sires]]
write.table(Alpha_pedigree_real, file = "Outputs/AlphaAssign/Alpha_pedigree_Real.txt", sep = " ", quote = F, col.names = F, row.names = F)

save.image(file = paste0(workingDir, "/Real_data/Pipeline/2_AlphaAssign.Rdata"))
