
args = commandArgs(trailingOnly=TRUE)
Rep = args[1]
workingDir = args[2]
softwareDir = args[3]
pathToPlink <- softwareDir
pathToBeagle <- softwareDir

repDir = paste0(workingDir, "/SimRep", Rep, "/")

# Set working directory
setwd(repDir)

colony_summary <- data.frame()

for (n in 1:5) {
  for (method in c("NoGE", "WithGE")) {
    if (file.exists(paste0(repDir, "/Data/Colony/", method, "_SNP", n, ".PairwisePaternity"))) {
      paternity <- read.csv(paste0(repDir, "/Data/Colony/", method, "_SNP", n, ".PairwisePaternity"))
      maternity <- read.csv(paste0(repDir, "/Data/Colony/", method, "_SNP", n, ".PairwiseMaternity"))

      Known_Dpc <- read.csv(paste0("Data/Sequoia/Known_Dpc_simulated_", method, "_SNP", n, ".csv"))
      colnames(Known_Dpc) <- c("OffspringID", "Mother", "CandidateID")

      merged_df <- merge(paternity, Known_Dpc, by="OffspringID", suffixes=c("_sim", "_known"))
      nCorrectSires <- merged_df$CandidateID_sim == merged_df$CandidateID_known

      colony_summary <- rbind(colony_summary, 
        data.frame(Rep = Rep, 
          Method = method, 
          SNP = n, 
          nOffspring = nrow(maternity),
          nAssigned_Sires = nrow(paternity),
          nCorrect_Sires = sum(nCorrectSires)))
        }
    }
}

write.csv(colony_summary, file = paste0(repDir, "/Outputs/Colony/Colony_summary.csv"), quote=FALSE, row.names = FALSE)