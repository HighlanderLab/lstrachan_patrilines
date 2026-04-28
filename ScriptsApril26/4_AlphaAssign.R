#******************************************************************************
# Prepare files for AlphaAssign pedigree reconstruction software
#******************************************************************************

#                     1. SIMULATED DATA 
#******************************************************************************

##############.  Prepping the input files ##################
rm(list = ls()) #Lets clear workspace just incase something sneaks in
pathToPlink <- "/home/jana/bin/"
workingDir = "/home/jana/github/lstrachan_patrilines/"

setwd(workingDir)
source("ScriptsApril26/Ped_to_raw.R")

pedigree_file <- read.csv("Data/worker_pedigree.csv")

Alpha_pedigree_sim <- data.frame(id = pedigree_file$id,
                             sire = rep(0, length(pedigree_file$id)),
                             dam = pedigree_file$mother)
f = data.frame(id = unique(Alpha_pedigree_sim$sire), sire = 0, dam = 0)
m = data.frame(id = unique(Alpha_pedigree_sim$dam), sire = 0, dam = 0)
Alpha_pedigree_sim_full <- rbind(f, m, Alpha_pedigree_sim)
write.table(Alpha_pedigree_sim_full, file = "Data/SimData_Pedigree_Full_Maternal.txt", quote = F, row.names = F, col.names = F)

dir.create("Data/AlphaAssign", showWarnings = FALSE)
write.table(Alpha_pedigree_sim, file = "Data/AlphaAssign/SimData_Pedigree.txt", sep = " ", quote = F, col.names = F, row.names = F)

Potential_fathers <- data.frame(id = pedigree_file$id,
                                Dpc1 = rep(unique(pedigree_file$dpc)[1], length(pedigree_file$dpc)),
                                Dpc2 = rep(unique(pedigree_file$dpc)[2], length(pedigree_file$dpc)),
                                Dpc3 = rep(unique(pedigree_file$dpc)[3], length(pedigree_file$dpc)),
                                Dpc4 = rep(unique(pedigree_file$dpc)[4], length(pedigree_file$dpc)))

write.table(Potential_fathers, file = "Data/AlphaAssign/SimData_PotentialFathers.list", sep = " ",  quote = F, col.names = F, row.names = F)

#Get known DPQs from pedigree 
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

True_pedigree <- rbind(Mothers_known, Dpc_known, Worker_known)
# write.table(True_pedigree, file = "Data/AlphaAssign/SimData_TruePedigree.csv", sep = ",", quote = F, col.names = T, row.names = F)

#With genotyping errors: 
setwd("Data/Sim_WithGE")
for (n in 1:5){ 
  print(n)
  ped_to_raw(ped_file = paste0("SNP_",n,"_WithGE_QC.ped"), map_file = paste0("SNP_",n,"_WithGE_QC.map"), output_file = paste0("SNP_",n,"_WithGE_QC_RecodeA.raw"))

  AlphaPed <- read.table(paste0("SNP_",n,"_WithGE_QC_RecodeA.raw"), header=TRUE)
  AlphaGeno <- AlphaPed[,7:ncol(AlphaPed)]; AlphaGeno[is.na(AlphaGeno)] <- 9
  AlphaGeno_id <- cbind(AlphaPed$IID, AlphaGeno)
  
  write.table(AlphaGeno_id, file=paste0(workingDir, "Data/AlphaAssign/AlphaGeno_SNP",n,"_WithGE.txt"), sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE) 
}
setwd(workingDir)
  
#No genotyping errors:   
setwd("Data/Sim_NoGE")
for (n in 1:5){ 
  ped_to_raw(ped_file = paste0("SNP_",n,"_NoGE_QC.ped"), map_file = paste0("SNP_",n,"_NoGE_QC.map"), output_file = paste0("SNP_",n,"_NoGE_QC_RecodeA.raw"))
  
  AlphaPed <- read.table(paste0("SNP_",n,"_NoGE_QC_RecodeA.raw"), header=TRUE)
  AlphaGeno <- AlphaPed[,7:ncol(AlphaPed)]; AlphaGeno[is.na(AlphaGeno)] <- 9
  AlphaGeno_id <- cbind(AlphaPed$IID, AlphaGeno)
  
  write.table(AlphaGeno_id, file=paste0(workingDir, "Data/AlphaAssign/AlphaGeno_SNP",n,"_NoGE.txt"), sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE) 
}
setwd(workingDir)


############### Running AlphAssign in terminal################
  #create a for loop here for all of the SNP array sizes <-------------- 
dir.create("Outputs/", showWarnings = FALSE)
dir.create("Outputs/AlphaAssign", showWarnings = FALSE)

#With genotyping errors:
for (n in 1:5){
  print(n)
  system(paste0("bash ScriptsApril26/RunAlphaAssign_SimData.sh AlphaGeno_SNP", n, "_WithGE ", workingDir, "Data/AlphaAssign/ ", workingDir, "Outputs/AlphaAssign/"))
}

#No genotyping errors:
for (n in 1:5){
  print(n)
  system(paste0("bash ScriptsApril26/RunAlphaAssign_SimData.sh AlphaGeno_SNP", n, "_NoGE ", workingDir, "Data/AlphaAssign/ ", workingDir, "Outputs/AlphaAssign/"))
}


############### Processing AlphaAssign output #################

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

setwd(paste0(workingDir, "Outputs/AlphaAssign/"))
NoGE_Alpha_output <- Process_AlphaAssign_output(GE = FALSE, True_pedigree = True_pedigree, nOffspring = nrow(Worker_known))
WithGE_Alpha_output <- Process_AlphaAssign_output(GE = TRUE, True_pedigree = True_pedigree, nOffspring = nrow(Worker_known))



# Updating the pedigree 
#Original pedigree = Alpha_pedigree
# NoGE updated pedigree ----------------------------------------------------
#Using the 1.7k SNP and the 50k SNP
NoGE_2k_pedigree <- read.table("AlphaGeno_SNP4_NoGE.sires", header=T)
NoGE_2k_pedigree <- NoGE_2k_pedigree[NoGE_2k_pedigree$chosen == 1, ]
NoGE_2k_pedigree <- NoGE_2k_pedigree[, c(1,2)]
colnames(NoGE_2k_pedigree) <- c("id", "sire")

# Match ids to get corresponding sires from NoGE_2k_pedigree
idx <- match(Alpha_pedigree_sim$id, NoGE_2k_pedigree$id)
# Replace only where sire == 0 and a match exists
update_sires <- Alpha_pedigree_sim$sire == 0 & !is.na(idx)
Alpha_pedigree_sim_2K <- Alpha_pedigree_sim
Alpha_pedigree_sim_2K$sire[update_sires] <- NoGE_2k_pedigree$sire[idx[update_sires]]
write.table(Alpha_pedigree_sim_2K, file = paste0(workingDir, "Outputs/AlphaAssign/Alpha_pedigree_2k_NoGE.txt"), sep = " ", quote = F, col.names = F, row.names = F)


NoGE_50k_pedigree <- read.table("AlphaGeno_SNP5_NoGE.sires", header=T)
NoGE_50k_pedigree <- NoGE_50k_pedigree[NoGE_50k_pedigree$chosen == 1, ]
NoGE_50k_pedigree <- NoGE_50k_pedigree[, c(1,2)]
colnames(NoGE_50k_pedigree) <- c("id", "sire")

# Match ids to get corresponding sires from NoGE_50k_pedigree
idx <- match(Alpha_pedigree_sim$id, NoGE_50k_pedigree$id)
# Replace only where sire == 0 and a match exists
update_sires <- Alpha_pedigree_sim$sire == 0 & !is.na(idx)
Alpha_pedigree_sim_50K <- Alpha_pedigree_sim
Alpha_pedigree_sim_50K$sire[update_sires] <- NoGE_50k_pedigree$sire[idx[update_sires]]
write.table(Alpha_pedigree_sim_50K, file = paste0(workingDir, "Outputs/AlphaAssign/Alpha_pedigree_50k_NoGE.txt"), sep = " ", quote = F, col.names = F, row.names = F)


# WithGE updated pedigree ----------------------------------------------------
#Using the 1.7k SNP and the 50k SNP
WithGE_2k_pedigree <- read.table("AlphaGeno_SNP4_WithGE.sires", header=T)
WithGE_2k_pedigree <- WithGE_2k_pedigree[WithGE_2k_pedigree$chosen == 1, ]
WithGE_2k_pedigree <- WithGE_2k_pedigree[, c(1,2)]
colnames(WithGE_2k_pedigree) <- c("id", "sire")

# Match ids to get corresponding sires from WithGE_2k_pedigree
idx <- match(Alpha_pedigree_sim$id, WithGE_2k_pedigree$id)
# Replace only where sire == 0 and a match exists
Alpha_pedigree_sim_2K_withGE <- Alpha_pedigree_sim
update_sires <- Alpha_pedigree_sim$sire == 0 & !is.na(idx)
Alpha_pedigree_sim_2K_withGE$sire[update_sires] <- WithGE_2k_pedigree$sire[idx[update_sires]]
write.table(Alpha_pedigree_sim_2K_withGE, file = paste0(workingDir, "Outputs/AlphaAssign/Alpha_pedigree_2k_WithGE.txt"), sep = " ", quote = F, col.names = F, row.names = F)


WithGE_50k_pedigree <- read.table("AlphaGeno_SNP5_WithGE.sires", header=T)
WithGE_50k_pedigree <- WithGE_50k_pedigree[WithGE_50k_pedigree$chosen == 1, ]
WithGE_50k_pedigree <- WithGE_50k_pedigree[, c(1,2)]
colnames(WithGE_50k_pedigree) <- c("id", "sire")

# Match ids to get corresponding sires from WithGE_50k_pedigree
idx <- match(Alpha_pedigree_sim$id, WithGE_50k_pedigree$id)
# Replace only where sire == 0 and a match exists
update_sires <- Alpha_pedigree_sim$sire == 0 & !is.na(idx)
Alpha_pedigree_sim_50K_withGE <- Alpha_pedigree_sim
Alpha_pedigree_sim_50K_withGE$sire[update_sires] <- WithGE_50k_pedigree$sire[idx[update_sires]]
write.table(Alpha_pedigree_sim_50K_withGE, file = paste0(workingDir, "Outputs/AlphaAssign/Alpha_pedigree_50k_WithGE.txt"), sep = " ", quote = F, col.names = F, row.names = F)


#******************************************************************************
#                     2. REAL DATA
#******************************************************************************

##############.  Prepping the input files ##################
#rm(list = ls()) #Lets clear workspace just incase something sneaks in 
pathToPlink <- "~/Desktop/PLINK/./"
workingDir = "~/Desktop/lstrachan_patrilines"

pathToPlink <- "/home/jana/bin/"
workingDir <- "/home/jana/github/lstrachan_patrilines/"

setwd(workingDir)
dir.create("Data/Real_data/AlphaAssign", showWarnings = FALSE)

pedigree_file_real <- read.csv("Data/Real_data/Real_Data_pedigree.csv")
colnames(pedigree_file_real) <- c("id", "sire", "dam")

removed_queen_id <- "B0003710" #Id of the queen removed during QC 
pedigree_file_real <- pedigree_file_real[pedigree_file_real$id != removed_queen_id,]
pedigree_file_real <- pedigree_file_real[pedigree_file_real$dam != removed_queen_id,] #remove all of her offspring from the pedigree too


Alpha_pedigree_real <- data.frame(id = pedigree_file_real$id,
                             sire = pedigree_file_real$sire,
                             dam = pedigree_file_real$dam)


write.table(Alpha_pedigree_real, file = "Data/Real_data/AlphaAssign/Pedigree.txt", sep = " ", quote = F, col.names = F, row.names = F)

SNP_samples <- read.csv("Data/Real_data/SNP_samples_2022.csv", header= T)
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
write.table(Potential_fathers, file = "Data/Real_data/AlphaAssign/PotentialFathers.list", sep = " ",  quote = F, col.names = F, row.names = F)


#Recode quality controlled files for AlphaAssign
ped_to_raw(ped_file = "Data/Real_data/Slov_fM_QC.ped", map_file = "Data/Real_data/Slov_fM_QC.map", output_file = "Data/Real_data/Slov_fM_QC.raw")

#Process the genotypes
AlphaPed <- read.table("Data/Real_data/Slov_fM_QC.raw", header=TRUE)
AlphaGeno <- AlphaPed[,7:ncol(AlphaPed)]; AlphaGeno[is.na(AlphaGeno)] <- 9
AlphaGeno_id <- cbind(AlphaPed$IID, AlphaGeno)

write.table(AlphaGeno_id, file=paste0("Data/Real_data/AlphaGeno_RealData.txt"), sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE) 

############### Running AlphAssign in terminal################
setwd(workingDir)
#create a for loop here for all of the SNP array sizes <-------------- 
system(paste0("bash ScriptsApril26/RunAlphaAssign_RealData.sh AlphaGeno_RealData ", workingDir, "/Data/Real_data/AlphaAssign ", workingDir, "/Outputs/AlphaAssign"))

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


AlphaAssign_Summary <- rbind(NoGE_Alpha_output, WithGE_Alpha_output, Real_Alpha_df)
write.table(AlphaAssign_Summary, file = paste0(workingDir, "Outputs/AlphaAssign/Alpha_summary.txt"), sep = " ", quote = F, col.names = T, row.names = F)

##################################################
# Update the pedigree
# Match ids to get corresponding sires from WithGE_50k_pedigree
idx <- match(Alpha_pedigree_real$id, Offspring_and_candidateParent$id)
# Replace only where sire == 0 and a match exists
update_sires <- Alpha_pedigree_real$sire == 0 & !is.na(idx)
Alpha_pedigree_real$sire[update_sires] <- Offspring_and_candidateParent$sire[idx[update_sires]]
write.table(Alpha_pedigree_real, file = "Outputs/AlphaAssign/Alpha_pedigree_Real.txt", sep = " ", quote = F, col.names = F, row.names = F)

save.image(file = paste0(workingDir, "/Data/Pipeline/4_AlphaAssign.Rdata"))