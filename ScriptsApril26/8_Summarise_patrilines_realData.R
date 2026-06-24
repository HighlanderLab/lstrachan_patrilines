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
    
    # 2. Safety check: only process if the loop hit a successful dataframe
    if (is.data.frame(df_list[[i]][[1]])) {
      tmp <- df_list[[i]]$nPaternity
      tmp$patrilines_per_worker = tmp$num_sister_groups_estimated / tmp$num_workers_actual 
      tmp$patrilines_times_dpq = tmp$num_sister_groups_estimated / tmp$num_dpqs_recon_ped

      # YOUR EXACT DATAFRAME STRUCTURE
      summary <- data.frame( 
        queen = tmp$queen_id, 
        number_patrilines = tmp$num_sister_groups_estimated,
        avg_patrilines_per_worker = tmp$patrilines_per_worker, 
        avg_patrilines_times_dpq = tmp$patrilines_times_dpq,
        sister_thresholds = tmp$sister_threshold, 
        data_type = tmp$data_type, 
        haplo_assignment_type = tmp$haplo_assignment_type 
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

library(tidyr)
library(dplyr)
library(ggplot2)
real_df1 <-  summarize_sisterPaternity_results_Real(PatSis_Real_Route1)
real_df2 <- summarize_sisterPaternity_results_Real(PatSis_Real_Route2)

PatSis_real_summary <-  rbind(real_df1, real_df2)

PatSis_real_summary <- PatSis_real_summary %>% 
  mutate(haplo_assignment_type = 
    recode(haplo_assignment_type, "Used mat pedigree" = "matPed", "Used recon pedigree" = "recPed"))

PatSis_real_summary_L <- PatSis_real_summary |> 
  pivot_longer(cols = c(number_patrilines, avg_patrilines_per_worker, avg_patrilines_times_dpq),
  values_to = "Mean", names_to = "Measure")

PatSis_real_summary_L <- PatSis_real_summary_L |> mutate(Measure_plot = recode(Measure, 
  "number_patrilines" = "nPatrilines",
   "avg_patrilines_times_dpq" = "Patriline/\nDPQ",
   "father_per_sister_group" = "Father/\nSister Group"
))

plot_dat <- PatSis_real_summary_L %>%
  filter(Measure %in% c("number_patrilines", "avg_patrilines_times_dpq")) %>%
  filter(
    sister_thresholds %in% c(0.75, 0.85, 0.9, 0.95, 1.0)
  ) %>%
    mutate(Sister_threshold = sister_thresholds)

palette <- c(
  "#4E79A7",  # muted blue
  "#F28E2B",  # muted orange
  "#59A14F",  # muted green
  "#B07AA1",  # muted purple
  "#E15759"   # muted red
)

plot_dat$facet_cols = paste0(plot_dat$haplo_assignment_type)
plot_dat$hline <- case_when(
  plot_dat$Measure %in% c(
    "correct_n_patrilines",
    "correct_patrilines"
  ) ~ 100,
  TRUE ~ 1
)

plot_dat$Sister_threshold <- as.factor(plot_dat$Sister_threshold)
plot_dat$Colony <- as.factor(as.numeric(as.factor(plot_dat$queen)))

fatsis_plot_real <- ggplot(plot_dat,
  aes(x = Colony, y = Mean, 
  fill = Sister_threshold, colour = Sister_threshold)) + 
  #geom_point(size = 4)+
  geom_col(position = "dodge") + 
  facet_grid(cols = vars(haplo_assignment_type), rows = vars(Measure_plot), 
scales = "free", space = "free_x") + 
  theme_bw(base_size = 24) + 
  scale_fill_manual(values = palette) + 
  scale_colour_manual(values = palette)  + 
  geom_hline(aes(yintercept = hline), linetype = "dashed", colour = "#464a4a", linewidth = 1) +
  xlab("Colony") + ylab("") + 
  theme(legend.position = "bottom", legend.box="vertical") 

png(paste0(workingDir, "/Patriline_determination_PaSis_real.png"), width = 1200, height = 800)
print(fatsis_plot_real)
dev.off()

PatSis_real_summary |> group_by(sister_thresholds) |> summarise(mean_patriline = mean(number_patrilines),
                                                                min_patriline = min(number_patrilines), 
                                                                max_patriline = max(number_patrilines))

dir.create("Outputs/PatrilineDetermination", showWarnings = FALSE)
write.csv(PatSis_real_summary, "Outputs/PatrilineDetermination/Patriline_summary.csv", quote=F, row.names=F)
save.image(file = "Pipeline/8_SummarisingData.RData")



