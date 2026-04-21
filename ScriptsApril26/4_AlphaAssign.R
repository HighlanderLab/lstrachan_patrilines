#******************************************************************************
# Prepare files for AlphaAssign pedigree reconstruction software
#******************************************************************************

#                     1. SIMULATED DATA 
#******************************************************************************

##############.  Prepping the input files ##################
pathToPlink <- "/home/jana/bin/"
workingDir = "/home/jana/github/lstrachan_patrilines/"
setwd(workingDir)
pedigree_file <- read.csv("Data/worker_pedigree.csv")

Alpha_pedigree <- data.frame(id = pedigree_file$id,
                             sire = rep(0, length(pedigree_file$id)),
                             dam = pedigree_file$mother)

dir.create("Data/AlphaAssign", showWarnings = FALSE)
write.table(Alpha_pedigree, file = "Data/AlphaAssign/Pedigree.txt", sep = " ", quote = F, col.names = F, row.names = F)

Potential_fathers <- data.frame(id = pedigree_file$id,
                                Dpc1 = rep(unique(pedigree_file$dpc)[1], length(pedigree_file$dpc)),
                                Dpc2 = rep(unique(pedigree_file$dpc)[2], length(pedigree_file$dpc)),
                                Dpc3 = rep(unique(pedigree_file$dpc)[3], length(pedigree_file$dpc)),
                                Dpc4 = rep(unique(pedigree_file$dpc)[4], length(pedigree_file$dpc)))

write.table(Potential_fathers, file = "Data/AlphaAssign/PotentialFathers.list", sep = " ",  quote = F, col.names = F, row.names = F)

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

Known_Dpc <- rbind(Mothers_known, Dpc_known, Worker_known)
write.table(Known_Dpc, file = "Data/AlphaAssign/Known_Dpc.csv", sep = ",", quote = F, col.names = T, row.names = F)

#With genotyping errors: 
setwd("Data/Sim_WithGE")
for (n in 1:5){ 
  print(n)
  system(paste0(pathToPlink, "plink --file SNP_",n,"_WithGE_ACformat_QC --recode A --out SNP_",n,"_WithGE_QC_RecodeA"))
  
  AlphaPed <- read.table(paste0("SNP_",n,"_WithGE_QC_RecodeA.raw"), header=TRUE)
  AlphaGeno <- AlphaPed[,7:ncol(AlphaPed)]; AlphaGeno[is.na(AlphaGeno)] <- 9
  AlphaGeno_id <- cbind(AlphaPed$IID, AlphaGeno)
  
  write.table(AlphaGeno_id, file=paste0(workingDir, "Data/AlphaAssign/AlphaGeno_SNP",n,"_WithGE.txt"), sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE) 
}
setwd(workingDir)
  
#No genotyping errors:   
setwd("Data/Sim_NoGE")
for (n in 1:5){ 
  system(paste0(pathToPlink, "plink --file SNP_",n,"_NoGE_ACformat_QC --recode A --out SNP_",n,"_NoGE_QC_RecodeA"))
  
  AlphaPed <- read.table(paste0("SNP_",n,"_NoGE_QC_RecodeA.raw"), header=TRUE)
  AlphaGeno <- AlphaPed[,7:ncol(AlphaPed)]; AlphaGeno[is.na(AlphaGeno)] <- 9
  AlphaGeno_id <- cbind(AlphaPed$IID, AlphaGeno)
  
  write.table(AlphaGeno_id, file=paste0(workingDir, "Data/AlphaAssign/AlphaGeno_SNP",n,"_NoGE.txt"), sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE) 
}
setwd(workingDir)


############### Running AlphAssign in terminal################
  #create a for loop here for all of the SNP array sizes <-------------- 
system(paste0("bash ~/Desktop/Slovenia_data/April26/RunAlphaAssign.sh AlphaGeno_SNP4_WithGE"))







############### Processing AlphaAssign output #################

Process_AlphaAssign_output <- function(GE = NULL){
  
  if (isTRUE(GE)) {
    setwd("~/Desktop/Slovenia data/April26/Simulated/Data/Sim_WithGE")
    test <- "WithGE"
    filetype <- "_WithGE"
  } else if (isFALSE(GE)) {
    setwd("~/Desktop/Slovenia data/April26/Simulated/Data/Sim_NoGE")
    test <- "NoGE"
    filetype <- "_NoGE"
  } else {
    stop("GE must be TRUE or FALSE")
  }
  
  results_list <- list()
  
  for (n in 1:5){
    
    # FIX: proper file name construction
    file_name <- paste0("AlphaGeno_SNP", n, filetype)
    
    Alpha_output <- read.table(file_name, header = TRUE)
    
    Sires_assigned <- Alpha_output[Alpha_output$chosen == 1, ]
    nSires_assigned <- nrow(Sires_assigned)
    
    Offspring_and_candidateParent <- Sires_assigned[, c(1,2)]
    colnames(Offspring_and_candidateParent) <- c("id", "sire")
    
    compare_sires <- merge(
      Offspring_and_candidateParent,
      Known_Dpc,
      by = "id",
      suffixes = c("_assigned", "_known")
    )
    
    nCorrect_sires <- sum(
      compare_sires$sire_assigned == compare_sires$sire_known,
      na.rm = TRUE
    )
    
    df <- data.frame(
      Test = test,
      nOffspring = nrow(Worker_known),
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

NoGE_Alpha_output <- Process_AlphaAssign_output(GE = FALSE)
WithGE_Alpha_output <- Process_AlphaAssign_output(GE = TRUE)

AlphaAssign_Summary <- rbind(NoGE_Alpha_output, WithGE_Alpha_output)
write.table(AlphaSummary, file = "Alpha_summary.txt", sep = " ", quote = F, col.names = T, row.names = F)
