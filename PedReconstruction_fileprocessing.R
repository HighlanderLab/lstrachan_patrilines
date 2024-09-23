#################################################################################

# Prepare and run Sequoia pedigree reconstruction software

#################################################################################

#Setting up the input files #####
pedigree_file <- read.csv("worker_pedigree.csv")

#Life History - need ID, Sex, BirthYear in a csv file 
Workers <- data.frame(ID = pedigree_file$id,
                      Sex = rep(1, length(pedigree_file$id)),
                      BirthYear = rep(2024, length(pedigree_file$id)))

#make the parents birth year the year after offspring 
Mothers <- data.frame(ID = unique(pedigree_file$mother),
                      Sex = rep(1, length(unique(pedigree_file$mother))),
                      BirthYear = rep(2023, length(unique(pedigree_file$mother))))

#make the dpc sex 2 = male to make sequoia think they're fathers 
Dpc <- data.frame(ID = unique(pedigree_file$dpc),
                  Sex = rep(2, length(unique(pedigree_file$dpc))),
                  BirthYear = rep(2023, length(unique(pedigree_file$dpc))))

LifeHistory <- rbind(Mothers, Dpc, Workers)
write.csv(LifeHistory, file = "LifeHistory.csv", sep = ",", quote = F, col.names = T, row.names = F)

#################################################################################
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

#################################################################################
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


##################################################################################
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


#On the Real data
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


#################################################################################

# Prepare files for AlphaAssign pedigree reconstruction software

#################################################################################

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

#################################################################################

# Prepare file for Colony pedigree reconstruction software

#################################################################################

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





#################################################################################

# Summarise the files from the KING pedigree reconstruction software

#################################################################################

#Create a summary file for KING 
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



