#******************************************************************************
# Running and analysing KING pedigree reconstruction software - SIMULATED DATA
#******************************************************************************

pathToPlink2 <- "~/Desktop/PLINK2"
workingDir = "~/Desktop/lstrachan_patrilines/SimRep2/Data/"


worker_pedigree <- read_csv("~/Desktop/lstrachan_patrilines/SimRep2/Data/worker_pedigree.csv")
dpq_ids <- unique(worker_pedigree$dpc)
worker_ids <- unique(worker_pedigree$id)


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

# •• No GE ••
setwd(paste0(workingDir, "Sim_NoGE/"))
for (n in 1:5){ 
  print("NoGE_SNP", n)
  system(paste0(pathToPlink2,"/plink2 --bfile SNP_",n,"_NoGE_QC --make-king-table --out NoGE_SNP",n,"_KINGoutput"))
}

# •• With GE ••
setwd(paste0(workingDir, "Sim_WithGE/"))
for (n in 1:5){ 
  print("SNP", n)
  system(paste0(pathToPlink2,"/plink2 --bfile SNP_",n,"_WithGE_QC --make-king-table --out WithGE_SNP",n,"_KINGoutput"))
}

tests <- c("NoGE", "WithGE")
snp_range <- 1:5 

results_list <- list()

for (t in tests) {
  
  # Construct the specific test directory path
  # This targets "~/Desktop/lstrachan_patrilines/SimRep2/Data/Sim_NoGE/" etc.
  testDir <- paste0(workingDir, "Sim_", t, "/")
  
  # Change the working directory to the current test folder
  if (dir.exists(testDir)) {
    setwd(testDir)
    message(paste("Processing directory:", getwd()))
  } else {
    warning(paste("Directory not found, skipping:", testDir))
    next 
  }
  
  for (s in snp_range) {
    
    file_name <- paste0(t, "_SNP", s, "_KINGoutput.kin0")
    
    if (file.exists(file_name)) {
      # Read the data
      temp_data <- read.table(file_name, header = TRUE)
      
      # Select first 8 columns and rename
      temp_data <- temp_data[, 1:8]
      colnames(temp_data) <- c("FID1","IID1", "FID2", "IID2", "NSNP", "HETHET", "IBS0", "KINSHIP")
      
      # Filter for Worker (IID1) and Queen/Drone (IID2) relationships
      filtered_WD <- temp_data[temp_data$KINSHIP >= 0.2 & 
                                 temp_data$IBS0 < 0.005 & 
                                 temp_data$IID1 %in% worker_ids & 
                                 temp_data$IID2 %in% dpq_ids, ]
      
      # Run matching function
      matches <- check_king(filtered_WD, worker_pedigree)
      
      # Create the data frame for this iteration
      res_df <- data.frame(
        Test = t,
        SNP_group = s,
        nWorkers = length(worker_ids),
        nDPQs_assigned = nrow(filtered_WD),
        Multi_perWorker = length(unique(filtered_WD$IID1)) != nrow(filtered_WD),
        nDPQs_correct = matches,
        pDPQs_correct = ifelse(nrow(filtered_WD) > 0, (matches / nrow(filtered_WD)) * 100, 0),
        Software = "KING",
        stringsAsFactors = FALSE
      )
      
      # Store in list
      results_list[[paste0(t, "_SNP", s)]] <- res_df
      
    } else {
      message(paste("  File not found:", file_name))
    }
  }
}

# Combine all results into one table
final_results_table <- do.call(rbind, results_list)
rownames(final_results_table) <- NULL

write.table(final_results_table, file = "../KING/KING_summary_Simulated.txt", col.names = T, row.names = F)
save.image(file = "../KING/KING_simulated.Rdata")
