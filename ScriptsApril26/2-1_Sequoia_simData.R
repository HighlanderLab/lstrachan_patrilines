#################################################################################

# Prepare and run Sequoia pedigree reconstruction software
rm(list = ls())

library(Eagle)
library(tidyr)
library(dplyr)

#################################################################################

#****** 1. SIMULATED *******************************************************


args = commandArgs(trailingOnly=TRUE)
Rep = args[1]
workingDir = args[2]
softwareDir = args[3]
pathToPlink <- softwareDir
pathToBeagle <- softwareDir

repDir = paste0(workingDir, "/SimRep", Rep, "/")

# Set working directory
setwd(repDir)

#Setting up the input files 
# *************** LIFE HISTORY ********************************
pedigree_file <- read.csv("Data/worker_pedigree.csv")

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

dir.create("Data/Sequoia")
setwd("Data/Sequoia")
write.table(LifeHistory, file = "LifeHistory_simulated.csv", sep = ",", quote = F, col.names = T, row.names = F)


# *********** Known DPCs **************************************
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
write.table(Known_Dpc, file = "Known_Dpc_simulated.csv", sep = ",", quote = F, col.names = T, row.names = F)

# •••••••••••• RUNNING SEQUOIA ••••••••••••••••••••••••••••••••
library(sequoia)

setwd(repDir)
dir.create("Outputs/Sequoia", showWarnings = FALSE)
for (n in 1:5){ 
  print(n)

  print("NoGE")
  # NoGE
  AlphaPed <- read.table(paste0(repDir, "/Data/Sim_NoGE/SNP_",n,"_NoGE_QC_RecodeA.raw"), header=TRUE)
  AlphaGeno <- AlphaPed[,7:ncol(AlphaPed)]
  AlphaGeno[AlphaGeno == 9] <- NA
  AlphaGeno <- as.matrix(AlphaGeno)
  rownames(AlphaGeno) <- AlphaPed$IID


  print("Run sequoia")
  SequoiaOutPut <- sequoia(GenoM = AlphaGeno, LifeHistData = LifeHistory, Module = "ped", Plot = TRUE)
  save(SequoiaOutPut, file = paste0("Outputs/Sequoia/SNP", n, "_NoGE.Rdata"))

  print("WithGE")
  # WithGE
  AlphaPed <- read.table(paste0(repDir, "/Data/Sim_WithGE/SNP_",n,"_WithGE_QC_RecodeA.raw"), header=TRUE)
  AlphaGeno <- AlphaPed[,7:ncol(AlphaPed)]
  AlphaGeno[AlphaGeno == 9] <- NA
  AlphaGeno <- as.matrix(AlphaGeno)
  rownames(AlphaGeno) <- AlphaPed$IID


  print("Run sequoia")
  SequoiaOutPut <- sequoia(GenoM = AlphaGeno, LifeHistData = LifeHistory, Module = "ped", Plot = TRUE)
  save(SequoiaOutPut, file = paste0("Outputs/Sequoia/SNP", n, "_WithGE.Rdata"))

}

#****** 3. SUMMARISE OUTPUTS *******************************************************
setwd(repDir)
#Load Simulated files: 
#•• NO GE ••
results_list <- list()
for (n in 1:5){
  load(paste0("Outputs/Sequoia/SNP",n,"_NoGE.Rdata"))
  Known_Dpc <- read.csv("Data/Sequoia/Known_Dpc_simulated.csv")
  Known_Dpc <- as.data.frame(Known_Dpc)
  
  PC_par <- PedCompare(Ped1 = Known_Dpc[, c("id", "dam", "sire")],
                       Ped2 = SequoiaOutPut$PedigreePar)
  
  nSires_assigned <- sum(!is.na(SequoiaOutPut$PedigreePar$sire))
  nCorrect_sires <- sum(PC_par[["MergedPed"]][["sire.class"]] == "Match")
  nCorrect_dams <- sum(PC_par[["MergedPed"]][["dam.class"]] == "Match")
  nMismatch <- sum((mergedPed$sire.class == "Match" & mergedPed$dam.class != "Match" ))
  
  Sequoia_file <- data.frame(Test = "NoGE",
                             nOffspring = 240,
                             SNP_group = n,
                             nSires_assigned = nSires_assigned,
                             nCorrect_sires = nCorrect_sires,
                             nCorrect_dam = nCorrect_dams,
                             nDamSireMismatch = nMismatch,
                             Software = "Sequoia")
  results_list[[n]] <- Sequoia_file
}
NoGE_SequoiaTable <- do.call(rbind, results_list)


#•• WITH GE ••
results_list <- list()
for (n in 1:5){
  load(paste0("Outputs/Sequoia/SNP",n,"_WithGE.Rdata"))
  Known_Dpc <- read.csv("Data/Sequoia/Known_Dpc_simulated.csv")
  Known_Dpc <- as.data.frame(Known_Dpc)
  
  PC_par <- PedCompare(Ped1 = Known_Dpc[, c("id", "dam", "sire")],
                       Ped2 = SequoiaOutPut$PedigreePar)
  
  nSires_assigned <- sum(!is.na(SequoiaOutPut$PedigreePar$sire))
  mergedPed <- PC_par[["MergedPed"]]
  nCorrect_sires <- sum(mergedPed$sire.class == "Match")
  nCorrect_dams <- sum(mergedPed$dam.class == "Match")
  nMismatch <- sum((mergedPed$sire.class == "Match" & mergedPed$dam.class != "Match" ))
  
  Sequoia_file <- data.frame(Test = "WithGE",
                             nOffspring = 240,
                             SNP_group = n,
                             nSires_assigned = nSires_assigned,
                             nCorrect_sires = nCorrect_sires,
                             nCorrect_dam = nCorrect_dams,
                             nDamSireMismatch = nMismatch,
                             Software = "Sequoia")
  results_list[[n]] <- Sequoia_file
}
WithGE_SequoiaTable <- do.call(rbind, results_list)


write.csv(NoGE_SequoiaTable, "Outputs/Sequoia/NoGE_SequoiaTable.csv", quote=F, row.names=F)
write.csv(WithGE_SequoiaTable, "Outputs/Sequoia/WithGE_SequoiaTable.csv", quote=F, row.names=F)
