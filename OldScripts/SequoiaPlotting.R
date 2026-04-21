
rm(list = ls())

#Packages
library(readr)
library(Eagle)
library(sequoia)


setwd("~/Desktop/Slovenia data/Jana's data /PostQCwithMaf0.01")
load("GMatrix_queenlessWorkersRemoved.Rdata")
#Convert it back to sequoia format
GenoMatrix[is.na(GenoMatrix)] <- -9

SNP_samplesQC_maf0_01 <- read_csv("~/Desktop/Slovenia data/Jana's data /PostQCwithMaf0.01/SNP_samplesQC_maf0.01_queenlessRemoved.csv")

LF <- SNP_samplesQC_maf0_01[, c(1,10)]
tmp <- LF[LF$biotype == c("queen"),]
tmp$BirthYear <- 2023 #example for now for pedigree reconstruction
tmp$Sex <- 1

tmp2 <- LF[LF$biotype == "worker", ]
tmp2$BirthYear <- 2024
tmp2$Sex <- 1


tmp3 <- LF[LF$biotype == "dpc", ]
tmp3$BirthYear <- 2023
tmp3$Sex <- 2

LifeHistory <- rbind(tmp, tmp2, tmp3)
LifeHistory$ID <- LifeHistory$snp_id
LifeHistory <- LifeHistory[, -c(1,2)]
LifeHistory <- LifeHistory[, c(3,2,1)]
LifeHistory$Sex <- as.integer(LifeHistory$Sex)
LifeHistory$BirthYear <- as.integer(LifeHistory$BirthYear)

#Have to run this for some reason to make LH work
LifeHistory$Sex <- as.integer(LifeHistory$Sex)
LifeHistory$BirthYear <- as.integer(LifeHistory$BirthYear)
LifeHistory <- as.data.frame(LifeHistory)

write_csv(LifeHistory, file = "LifeHistoryQueenlessremoved.csv")


#duplicate check and parentage assignment
ParOUT <- sequoia(GenoM = GenoMatrix, LifeHistData = LifeHistory, Module = "par", quiet = FALSE, Plot = TRUE)

#Compare with known parents
KnownParents <- read_csv("KnownParents_dpc.csv")
KnownParents <- as.data.frame(KnownParents)

PC_par <- PedCompare(Ped1 = KnownParents[, c("id", "dam", "sire")],
                     Ped2 = SequoiaOutPut$PedigreePar)


#G = genotyped
#D = dummy
#T = totals


#Since there are no sires in known pedigree - have a look at the dams only
PC_par$Counts[,,"dam"]

#Lets looks at the mismatches (dam1 is Ped1 (known parents))
MisMatches <- PC_par$MergedPed[which(PC_par$MergedPed$dam.class == "Mismatch"), c("id", "dam.1", "dam.2")]

#Can look at the stats for each queen
MergedPed <-  PC_par$MergedPed[, c("id", "dam.1", "dam.2")]
MisMatchedQueen1 <- MergedPed[which(MergedPed$dam.1 == "B0003389"),] #does having a different queen suggest the queens are perhaps related?


load("~/Desktop/Slovenia data/Sequoia/PostQC_maf0.01/SequoiaOutput.RData")

SummarySeq(SequoiaOutPut$PedigreePar)
SummarySeq(SequoiaOutPut$Pedigree)


PedPar <- SequoiaOutPut$PedigreePar
Number_Dams <- sum(!is.na(PedPar$dam))

SeqPed <- SequoiaOutPut$Pedigree
B_ids <- SeqPed$dam[grep("^B", SeqPed$dam, ignore.case = FALSE)]                 #Find out number of observed dam-offspring
DummyDam_ids <- SeqPed$dam[grep("^F", SeqPed$dam, ignore.case = FALSE)]          #find out number of dummy dam-offspring
DummySire_ids <- SeqPed$sire[grep("^M", SeqPed$sire, ignore.case = FALSE)]       #find out number of dummy sire-offspring










