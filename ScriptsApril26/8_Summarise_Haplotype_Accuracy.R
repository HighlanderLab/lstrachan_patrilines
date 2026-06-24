args = commandArgs(trailingOnly=TRUE)

# Set working directory
setwd(workingDir)

phasing_accuracy <- data.frame()

nRep = 10
for (rep in 1:nRep) {
  print(rep)
  load(paste0(workingDir, "/SimRep", rep, "/Pipeline/5_Haplotype_ParentAssignments.RData"))

  for (snp_array in c("2k", "50k")) {
    print(snp_array)
    for (route in c("recPed", "matPed")) {
      true_vs_NoGEphasedSNP <- get(paste0("true_vs_NoGEphasedSNP", snp_array, route))
      tmp <- true_vs_NoGEphasedSNP$summary
      tmp$SNPArray <- snp_array
      tmp$Route <- route
      phasing_accuracy <- rbind(phasing_accuracy, tmp)
    }
  }
}