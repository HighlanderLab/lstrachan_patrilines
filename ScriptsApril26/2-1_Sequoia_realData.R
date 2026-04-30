#################################################################################

# Prepare and run Sequoia pedigree reconstruction software
rm(list = ls())

library(Eagle)
library(tidyr)
library(dplyr)
library(sequoia)


args = commandArgs(trailingOnly=TRUE)
workingDir = args[1]
softwareDir = args[2]
pathToPlink <- softwareDir
pathToBeagle <- softwareDir


# Set working directory
setwd(paste0(workingDir, "/Real_data/"))
#################################################################################


# ******* LIFE HISTORY *****************
pedigree_file <- read.csv("Data/Real_Data_pedigree.csv")
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

dir.create("Data/Sequoia", showWarnings = F)
write.table(LifeHistory, file = "Data/Sequoia/LifeHistory_Real.csv", sep = ",", quote = F, col.names = T, row.names = F)

# No known dpqs 

# •••••••• RUN SEQUOIA ••••••••••••••••••••••••

AlphaPed <- read.table("Data/Slov_fM_QC.raw", header=TRUE)
AlphaGeno <- AlphaPed[,7:ncol(AlphaPed)]
AlphaGeno[AlphaGeno == 9] <- NA
AlphaGeno <- as.matrix(AlphaGeno)
rownames(AlphaGeno) <- AlphaPed$IID

print("Run sequoia")
SequoiaOutPut <- sequoia(GenoM = AlphaGeno, LifeHistData = LifeHistory, Module = "ped", Plot = TRUE)

dir.create("Outputs/Sequoia", showWarnings = F)
save(SequoiaOutPut, file = "Outputs/Sequoia/Real.Rdata")



#•• REAL DATA •• 
nSires_assigned <- sum(!is.na(SequoiaOutPut$PedigreePar$sire))

Real_SequoiaTable <- data.frame(Test = "Real",
                           nOffspring = 235,
                           SNP_group = 4,
                           nSires_assigned = nSires_assigned,
                           nCorrect_sires = NA,
                           Software = "Sequoia")


write.csv(Real_SequoiaTable, "Outputs/Sequoia/Real_SequoiaTable.csv", quote=F, row.names=F)

