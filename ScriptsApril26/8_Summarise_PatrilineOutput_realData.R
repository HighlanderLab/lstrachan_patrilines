#######################################################################################################################
#************************************* Summarising Patriline determination *******************************
#######################################################################################################################
#Creating workable tables to view and summarise results of the patriline determination 

rm(list = ls())
# --- Libraries ---

args = commandArgs(trailingOnly=TRUE)
workingDir = args[1]
softwareDir = args[2]
pathToPlink <- softwareDir
pathToBeagle <- softwareDir


# Set working directory
setwd(paste0(workingDir, "/Real_data/"))

load("Pipeline/7_DeterminingPatrilines.Rdata")


#**•• REAL ••**

summarize_sisterPaternity_results_Real <- function(df_list) {
  
  # 1. Initialize a list to collect each summary row
  summary_collection <- list()
  
  for (i in 1:length(df_list)) {
    
    # Extract the element
    tmp <- df_list[[i]]
    
    # 2. Safety check: only process if the loop hit a successful dataframe
    if (is.data.frame(tmp)) {
      
      # YOUR EXACT DATAFRAME STRUCTURE
      summary <- data.frame( 
        num_queens = length(tmp$queen_id), 
        avg_dpqs_recon = mean(tmp$num_dpqs_recon_ped, na.rm = TRUE), 
        avg_nworkers_actual = mean(tmp$num_workers_actual, na.rm = TRUE), 
        avg_nworkers_recon = mean(tmp$num_workers_recon, na.rm = TRUE), 
        avg_sistergroups_est = mean(tmp$num_sister_groups_estimated, na.rm = TRUE), 
        sister_thresholds = unique(tmp$sister_threshold), 
        data_type = unique(tmp$data_type), 
        haplo_assignment_type = unique(tmp$haplo_assignment_type) 
      )
      
      # 3. Store the summary row
      summary_collection[[i]] <- summary
      
    } else {
      # If the element was a "FAILED" message from tryCatch, we skip it
      message(paste("Skipping index", i, "as it is not a dataframe."))
    }
  }
  
  # 4. Rbind everything together at the end
  # Filter out any NULL entries in the list before rbinding
  summary_collection <- summary_collection[!sapply(summary_collection, is.null)]
  final_summary_df <- do.call(rbind, summary_collection)
  
  return(final_summary_df)
}

real_df1 <-  summarize_sisterPaternity_results_Real(PatR2_Real_Haplo1)
real_df2 <- summarize_sisterPaternity_results_Real(PatR2_Real_Haplo2)

PatR2_real_summary <-  rbind(real_df1, real_df2)

dir.create("Outputs/PatrilineDetermination")
write.csv(PatR2_real_summary, "Outputs/PatrilineDetermination/Patriline_summary.csv", quote=F, row.names=F)
save.image(file = "Pipeline/8_SummarisingData.RData")



