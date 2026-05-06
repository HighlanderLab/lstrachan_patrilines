#################################################################################

# Prepare and run Sequoia pedigree reconstruction software
rm(list = ls())

library(Eagle)
library(tidyr)
library(dplyr)
library(sequoia)
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


# •••••••••••• RUNNING SEQUOIA ••••••••••••••••••••••••••••••••

setwd(repDir)
dir.create("Data/Sequoia", showWarnings = FALSE)
dir.create("Outputs/Sequoia", showWarnings = FALSE)
timing = data.frame()

for (n in 1:5){ 
  for (method in c("NoGE", "WithGE")) {
    setwd(repDir)
    print(n)
    print(method)

    # Read in genotypes
    AlphaPed <- read.table(paste0(repDir, "/Data/Sim_", method, "/SNP_",n,"_", method, "_QC_RecodeA.raw"), header=TRUE)


    #Setting up the input files 
    # *************** LIFE HISTORY ********************************
    pedigree_file <- read.csv("Data/worker_pedigree.csv")
    pedigree_file <- pedigree_file[pedigree_file$id %in% AlphaPed$IID,]
    pedigree_file <- pedigree_file[pedigree_file$mother %in% AlphaPed$IID,]
    pedigree_file <- pedigree_file[pedigree_file$dpc %in% AlphaPed$IID,]

    AlphaPed <- AlphaPed[AlphaPed$IID %in% c(pedigree_file$id, pedigree_file$dpc, pedigree_file$mother),]
    AlphaGeno <- AlphaPed[,7:ncol(AlphaPed)]
    AlphaGeno[AlphaGeno == 9] <- NA
    AlphaGeno <- as.matrix(AlphaGeno)
    rownames(AlphaGeno) <- AlphaPed$IID


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


    setwd("Data/Sequoia")
    #write.table(LifeHistory, file = paste0("LifeHistory_simulated_", method, "_SNP", n, ".csv"), sep = ",", quote = F, col.names = T, row.names = F)


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

    write.table(Known_Dpc, file = paste0("Known_Dpc_simulated_", method, "_SNP", n, ".csv"), sep = ",", quote = F, col.names = T, row.names = F)

    print("Run sequoia")
    setwd(repDir)
    start = Sys.time()
    SequoiaOutPut <- sequoia(GenoM = AlphaGeno, LifeHistData = LifeHistory, Module = "ped", Plot = TRUE)
    end = Sys.time()
    save(SequoiaOutPut, file = paste0("Outputs/Sequoia/SNP", n, "_", method, ".Rdata"))

    timing <- rbind(timing, c(method = method, n = n, Software = "Sequoia", Rep = Rep, time = as.numeric(as.difftime(end-start, units = "s"))))


  }
}

write.csv(timing, "Outputs/Sequoia/Timing.csv", quote=F, row.names=F)

#****** 3. SUMMARISE OUTPUTS *******************************************************
setwd(repDir)
#Load Simulated files: 
#•• NO GE ••
results_list <- list()
x = 1
for (n in 1:5){
  for (method in c("NoGE", "WithGE")) {
    print(n)
    print(method)

    load(paste0("Outputs/Sequoia/SNP",n,"_", method, ".Rdata"))
    Known_Dpc <- read.csv(paste0("Data/Sequoia/Known_Dpc_simulated_", method, "_SNP", n, ".csv"))
    Known_Dpc <- as.data.frame(Known_Dpc)
    
    PC_par <- PedCompare(Ped1 = Known_Dpc[, c("id", "dam", "sire")],
                        Ped2 = SequoiaOutPut$PedigreePar)
    
    nSires_assigned <- sum(!is.na(SequoiaOutPut$PedigreePar$sire))
    mergedPed <- PC_par[["MergedPed"]]
    nCorrect_sires <- sum(PC_par[["MergedPed"]][["sire.class"]] == "Match")
    nCorrect_dams <- sum(PC_par[["MergedPed"]][["dam.class"]] == "Match")
    nMismatch <- sum((mergedPed$sire.class == "Match" & mergedPed$dam.class != "Match" ))
    
    Sequoia_file <- data.frame(Test = method,
                              nOffspring = nrow(Known_Dpc),
                              SNP_group = n,
                              nSires_assigned = nSires_assigned,
                              nCorrect_sires = nCorrect_sires,
                              nCorrect_dam = nCorrect_dams,
                              nDamSireMismatch = nMismatch,
                              Software = "Sequoia")
    results_list[[x]] <- Sequoia_file
    x = x+1
  }
}
  
SequoiaTable <- do.call(rbind, results_list)




write.csv(SequoiaTable, "Outputs/Sequoia/SequoiaTable.csv", quote=F, row.names=F)
  