#################################################################################

# Prepare and run Sequoia pedigree reconstruction software
rm(list = ls())

library(Eagle)
library(tidyr)
library(dplyr)

#################################################################################

#****** 1. SIMULATED *******************************************************

pathToPlink <- "~/Desktop/PLINK/./"
workingDir = "~/Desktop/lstrachan_patrilines"
setwd(workingDir)

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
write.csv(LifeHistory, file = "LifeHistory_simulated.csv", sep = ",", quote = F, col.names = T, row.names = F)


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
write.csv(Known_Dpc, file = "Known_Dpc_simulated.csv", sep = ",", quote = F, col.names = T, row.names = F)

# •••••••••••• RUNNING SEQUOIA ••••••••••••••••••••••••••••••••
library(sequoia)

#1. No Genotyping Errors
#TODO : Make a loop here for each SNP size <----------------------
#convert ped file to sequoia format 
map <- read.table("SNPnames_NOGE_ACformat.map", header = F) #TODO : Need to fix name 
ped <- read.table("SNPnames_NOGE_ACformat.ped", header = F)

SNP_names <- map[,2]
SNP_names_new <- unlist(lapply(SNP_names, function(name) c(paste(name, "1", sep="_"), paste(name, "2", sep="_"))))
colnames(ped)[7:ncol(ped)] <- SNP_names_new

#add Output file if you are running it on Eddie 
Sequoia_ped <- GenoConvert(InData = ped, InFormat = 'ped', Missing = '0', IDcol = 2, header = F)
#Sequoia_ped <- GenoConvert(InData = ped, InFormat = 'ped', Missing = '0', IDcol = 2, header = F, OutFile = "Sequoia_ped.txt")

print("Run sequoia")
SequoiaOutPut <- sequoia(GenoM = Sequoia_ped2, LifeHistData = LifeHistory, Module = "ped", Plot = TRUE)
save(SequoiaOutPut, file = "SequoiaOutPut_SNPn_NoGE.Rdata")

rm(... = map, ped, SNP_names, SNP_names_new, Sequoia_ped, SequoiaOutPut)

#2. With Genotyping Errors
#TODO : Make a loop here for each SNP size <----------------------
#convert ped file to sequoia format 
map <- read.table("SNPnames_WithGE_ACformat.map", header = F) #TODO : Need to fix name 
ped <- read.table("SNPnames_WithGE_ACformat.ped", header = F)

SNP_names <- map[,2]
SNP_names_new <- unlist(lapply(SNP_names, function(name) c(paste(name, "1", sep="_"), paste(name, "2", sep="_"))))
colnames(ped)[7:ncol(ped)] <- SNP_names_new

#add Output file if you are running it on Eddie 
Sequoia_ped <- GenoConvert(InData = ped, InFormat = 'ped', Missing = '0', IDcol = 2, header = F)
#Sequoia_ped <- GenoConvert(InData = ped, InFormat = 'ped', Missing = '0', IDcol = 2, header = F, OutFile = "Sequoia_ped.txt")

print("Run sequoia")
SequoiaOutPut <- sequoia(GenoM = Sequoia_ped, LifeHistData = LifeHistory, Module = "ped", Plot = TRUE)
save(SequoiaOutPut, file = "SequoiaOutPut_SNPn_NoGE.Rdata")

rm(... = map, ped, SNP_names, SNP_names_new, Sequoia_ped, SequoiaOutPut)




#****** 2. REAL DATA *******************************************************
rm(list = ls())
pathToPlink <- "~/Desktop/PLINK/./"
workingDir = "~/Desktop/lstrachan_patrilines"
setwd(workingDir)

# ******* LIFE HISTORY *****************
pedigree_file <- read.csv("~/Desktop/lstrachan_patrilines/Data/Real_data/Real_Data_pedigree.csv")
colnames(pedigree_file) <- c("id", "sire", "dam")

#Life History - need ID, Sex, BirthYear in a csv file 
Workers <- data.frame(ID = pedigree_file$id,
                      Sex = rep(1, length(pedigree_file$id)),
                      BirthYear = rep(2024, length(pedigree_file$id)))

#make the parents birth year the year after offspring 
Mothers <- data.frame(ID = unique(pedigree_file$dam),
                      Sex = rep(1, length(unique(pedigree_file$dam))),
                      BirthYear = rep(2023, length(unique(pedigree_file$dam))))


#Do we still need to add DPCs to the LifeHistory or will it be confused? 
Dpc <- data.frame(ID = unique(pedigree_file$sire),
                  Sex = rep(2, length(unique(pedigree_file$sire))),
                  BirthYear = rep(2023, length(unique(pedigree_file$sire))))

LifeHistory <- rbind(Mothers, Dpc, Workers)
write.csv(LifeHistory, file = "LifeHistory_Real.csv", sep = ",", quote = F, col.names = T, row.names = F)

# No known dpqs 

# •••••••• RUN SEQUOIA ••••••••••••••••••••••••
library(sequoia)

#convert ped file to sequoia format 
map <- read.table("processedRealfile.map", header = F) # TODO : Add in the right names for the files
ped <- read.table("processedRealfile.ped", header = F)

SNP_names <- map[,2]
SNP_names_new <- unlist(lapply(SNP_names, function(name) c(paste(name, "1", sep="_"), paste(name, "2", sep="_"))))
colnames(ped)[7:ncol(ped)] <- SNP_names_new

#add Output file if you are running it on Eddie 
Sequoia_ped <- GenoConvert(InData = ped, InFormat = 'ped', Missing = '0', IDcol = 2, header = F)
#Sequoia_ped <- GenoConvert(InData = ped, InFormat = 'ped', Missing = '0', IDcol = 2, header = F, OutFile = "Sequoia_ped.txt")

print("Run sequoia")
SequoiaOutPut <- sequoia(GenoM = Sequoia_ped, LifeHistData = LifeHistory, Module = "ped", Plot = TRUE)
save(SequoiaOutPut, file = "SequoiaOutPut_Real.Rdata")

rm(... = map, ped, SNP_names, SNP_names_new, Sequoia_ped, SequoiaOutPut)



#****** 3. SUMMARISE OUTPUTS *******************************************************
setwd(workingDir)
#Load Simulated files: 
#•• NO GE ••
results_list <- list()
for (n in 1:5){
  SequoiaOutPut <- load("Data/Sequoia/SequoiaOutPut_SNP",n,"_NoGE.Rdata")
  Known_Dpc <- read.csv("Data/Sequoia/Known_Dpc.csv")
  Known_Dpc <- as.data.frame(Known_Dpc)
  
  PC_par <- PedCompare(Ped1 = Known_Dpc[, c("id", "dam", "sire")],
                       Ped2 = SequoiaOutPut$PedigreePar)
  
  nSires_assigned <- sum(!is.na(SequoiaOutPut$PedigreePar$sire))
  nCorrect_sires <- sum(PC_par[["MergedPed"]][["sire.class"]] == "Match")
  
  Sequoia_file <- data.frame(Test = "NoGE",
                             nOffspring = 240,
                             SNP_group = n,
                             nSires_assigned = nSires_assigned,
                             nCorrect_sires = nCorrect_sires,
                             Software = "Sequoia")
  results_list[[n]] <- Sequoia_file
}
NoGE_SequoiaTable <- do.call(rbind, results_list)


#•• WITH GE ••
results_list <- list()
for (n in 1:5){
  SequoiaOutPut <- load("Data/Sequoia/SequoiaOutPut_SNP",n,"_WithGE.Rdata")
  Known_Dpc <- read.csv("Data/Sequoia/Known_Dpc.csv")
  Known_Dpc <- as.data.frame(Known_Dpc)
  
  PC_par <- PedCompare(Ped1 = Known_Dpc[, c("id", "dam", "sire")],
                       Ped2 = SequoiaOutPut$PedigreePar)
  
  nSires_assigned <- sum(!is.na(SequoiaOutPut$PedigreePar$sire))
  nCorrect_sires <- sum(PC_par[["MergedPed"]][["sire.class"]] == "Match")
  
  Sequoia_file <- data.frame(Test = "WithGE",
                             nOffspring = 240,
                             SNP_group = n,
                             nSires_assigned = nSires_assigned,
                             nCorrect_sires = nCorrect_sires,
                             Software = "Sequoia")
  results_list[[n]] <- Sequoia_file
}
WithGE_SequoiaTable <- do.call(rbind, results_list)


#•• REAL DATA •• 
SequoiaOutPut <- load("Data/Sequoia/SequoiaOutPut_Real.Rdata")
nSires_assigned <- sum(!is.na(SequoiaOutPut$PedigreePar$sire))

Real_SequoiaTable <- data.frame(Test = "Real",
                           nOffspring = 235,
                           SNP_group = 4,
                           nSires_assigned = nSires_assigned,
                           nCorrect_sires = NA,
                           Software = "Sequoia")

Sequoia_Summary <- rbind(Real_SequoiaTable, NoGE_SequoiaTable, WithGE_SequoiaTable)
Sequoia_Summary

