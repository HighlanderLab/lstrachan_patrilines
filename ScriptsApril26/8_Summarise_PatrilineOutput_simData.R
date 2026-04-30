#######################################################################################################################
#************************************* Summarising Patriline determination *******************************
#######################################################################################################################
#Creating workable tables to view and summarise results of the patriline determination 


load("Data/Pipeline/7_Haplotype_ParentAssignments.Rdata")

#** ••• SIMULATED ••• **
#Father + sister thresholds 

summarize_paternity_results <- function(df_list) {
  
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
        avg_dpqs_actual = mean(tmp$num_dpqs_actual, na.rm = TRUE), 
        avg_dpqs_recon = mean(tmp$num_dpqs_recon_ped, na.rm = TRUE), 
        avg_nfathers_actual = mean(tmp$num_fathers_actual, na.rm = TRUE), 
        avg_nworkers_actual = mean(tmp$num_workers_actual, na.rm = TRUE), 
        avg_nworkers_recon = mean(tmp$num_workers_recon, na.rm = TRUE), 
        avg_nfathers_est = mean(tmp$num_fathers_estimated, na.rm = TRUE), 
        avg_nfathers_correct = mean(tmp$num_fathers_correct, na.rm = TRUE), 
        father_thresholds = unique(tmp$father_accuracy_threshold), 
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

df1 <- summarize_paternity_results(PatR1_SimTrue_Haplo1)
df2 <- summarize_paternity_results(PatR1_NoGE_SNP2k_Haplo1)
df3 <- summarize_paternity_results(PatR1_NoGE_SNP50k_Haplo1)
df4 <- summarize_paternity_results(PatR1_WithGE_SNP2k_Haplo1)
df5 <- summarize_paternity_results(PatR1_WithGE_SNP50k_Haplo1)
PatR1_Haplo1_summary <- rbind(df1,df2,df3,df4,df5)


df6 <- summarize_paternity_results(PatR1_SimTrue_Haplo2)
df7 <- summarize_paternity_results(PatR1_NoGE_SNP2k_Haplo2)
df8 <- summarize_paternity_results(PatR1_NoGE_SNP50k_Haplo2)
df9 <- summarize_paternity_results(PatR1_WithGE_SNP2k_Haplo2)
df10 <- summarize_paternity_results(PatR1_WithGE_SNP50k_Haplo2)
PatR1_Haplo2_summary <- rbind(df6,df7,df8,df9,df10)

PatR1_summary <- rbind(PatR1_Haplo1_summary, PatR1_Haplo2_summary)



#Only Sister thresholds 
summarize_sisterPaternity_results <- function(df_list) {
  
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
        avg_dpqs_actual = mean(tmp$num_dpqs_actual, na.rm = TRUE), 
        avg_dpqs_recon = mean(tmp$num_dpqs_recon_ped, na.rm = TRUE), 
        avg_nworkers_actual = mean(tmp$num_workers_actual, na.rm = TRUE), 
        avg_nworkers_recon = mean(tmp$num_workers_recon, na.rm = TRUE), 
        avg_sistergroups_est = mean(tmp$num_sister_groups_estimated, na.rm = TRUE), 
        avg_nfathers_actual = mean(tmp$actual_number_fathers, na.rm = TRUE), 
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
df11 <- summarize_sisterPaternity_results(PatR2_SimTrue_Haplo1)
df12 <- summarize_sisterPaternity_results(PatR2_NoGE_SNP2k_Haplo1)
df13 <- summarize_sisterPaternity_results(PatR2_NoGE_SNP50k_Haplo1)
df14 <- summarize_sisterPaternity_results(PatR2_NoGE_SNP2k_Haplo1)
df15 <- summarize_sisterPaternity_results(PatR2_WithGE_SNP2k_Haplo1)
PatR2_Haplo1_summary <- rbind(df11,df12,df13,df14,df15)



df16 <- summarize_sisterPaternity_results(PatR2_SimTrue_Haplo2)
df17 <- summarize_sisterPaternity_results(PatR2_NoGE_SNP2k_Haplo2)
df18 <- summarize_sisterPaternity_results(PatR2_NoGE_SNP50k_Haplo2)
df19 <- summarize_sisterPaternity_results(PatR2_WithGE_SNP2k_Haplo2)
df20 <- summarize_sisterPaternity_results(PatR2_WithGE_SNP50k_Haplo2)
PatR2_Haplo2_summary <- rbind(df16,df17,df18,df19,df20)


PatR2_summary <- rbind(PatR2_Haplo1_summary, PatR2_Haplo2_summary)


save.image(file = "Pipeline/8_SummarisingData.RData")
