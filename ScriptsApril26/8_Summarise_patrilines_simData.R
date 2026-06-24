###########################################################################
###########################################################################
# Second, summarise the patriline results
###########################################################################
###########################################################################
# Father and sister threshold
summarize_paternity_results <- function(df_list) {
  
  # 1. Initialize a list to collect each summary row
  summary_collection <- list()
  
  for (i in 1:length(df_list)) {
    print(paste("Processing index:", names(df_list)[i]))
   
    
    # 2. Safety check: only process if the loop hit a successful dataframe
    if (is.data.frame(df_list[[i]][[1]])) {
      # This is numbers of assigned
      tmp <- df_list[[i]]$nPaternity
      tmp$per_dpq = tmp$num_dpqs_recon_ped / tmp$num_dpqs_actual * 100
      tmp$per_patrilines = tmp$num_fathers_estimated / tmp$num_fathers_actual * 100
      tmp$patrilines_per_worker = tmp$num_fathers_estimated / tmp$num_workers_actual
      tmp$patrilines_times_dpq = tmp$num_fathers_estimated / tmp$num_dpqs_recon_ped

      # This is checking how many are correct
      tmp1 <- df_list[[i]]$worker_assignments
      tmp1_summary <- tmp1 |> group_by(queen_id) |> 
        summarise(correct_father = sum(actual_father == assigned_father) / length(actual_father) * 100,
                  correct_dpq = sum(actual_dpq == assigned_dpq) / length(actual_dpq) * 100,
                .groups = "drop")
      
      
      # YOUR EXACT DATAFRAME STRUCTURE
      summary <- data.frame( 
        num_queens = length(unique(tmp$queen_id)), 
        avg_correct_n_dpq = mean(tmp$per_dpq, na.rm = TRUE), 
        avg_correct_n_patrilines = mean(tmp$per_patrilines, na.rm = TRUE), 
        avg_patrilines_per_worker = mean(tmp$patrilines_per_worker, na.rm = TRUE), 
        avg_patrilines_times_dpq = mean(tmp$patrilines_times_dpq, na.rm = TRUE),
        avg_correct_patrilines = mean(tmp1_summary$correct_father, na.rm = TRUE),
        avg_correct_dpq = mean(tmp1_summary$correct_dpq, na.rm = TRUE),
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

#Only Sister thresholds 
summarize_sisterPaternity_results <- function(df_list) {
  
  # 1. Initialize a list to collect each summary row
  summary_collection <- list()
  
  for (i in 1:length(df_list)) {
    
    
    # 2. Safety check: only process if the loop hit a successful dataframe
    if (is.data.frame(df_list[[i]][[1]])) {
      tmp <- df_list[[i]]$nPaternity
      tmp$per_patrilines = tmp$num_sister_groups_estimated / tmp$actual_number_fathers * 100
      tmp$patrilines_per_worker = tmp$num_sister_groups_estimated / tmp$num_workers_actual 
      tmp$patrilines_times_dpq = tmp$num_sister_groups_estimated / tmp$num_dpqs_recon_ped

      # This is checking how many are correct
      tmp1 <- df_list[[i]]$worker_assignments
      tmp1_summary <- tmp1|> group_by(queen_id, sister_group) |> 
        summarise(nFather = length(unique(actual_father)), .groups = "drop") 
      
      # YOUR EXACT DATAFRAME STRUCTURE
      summary <- data.frame( 
        num_queens = length(unique(tmp$queen_id)), 
        avg_correct_n_patrilines = mean(tmp$per_patrilines, na.rm = TRUE), 
        avg_patrilines_per_worker = mean(tmp$patrilines_per_worker, na.rm = TRUE), 
        avg_patrilines_times_dpq = mean(tmp$patrilines_times_dpq, na.rm = TRUE),
        avg_nFather_per_sisterGroup = mean(tmp1_summary$nFather, na.rm = TRUE),
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

library(dplyr)
library(ggplot2)

workingDir = getwd()
nRep = 10
PatFatSis_summary = data.frame()
PatSis_summary = data.frame()

for (Rep in 1:nRep) {
  print(paste0("Rep: ", Rep))
  repDir = paste0(workingDir, "/SimRep", Rep, "/")
  setwd(repDir)
  # List all files that start with Pipeline/7_Pat

  files = Sys.glob("Pipeline/7_Pat*")
  i = 1
  for (file in files) {
    #print(abs(difftime(file.info(file)$mtime, Sys.time(), units = "hours")))
      name = basename(file)
      name = gsub(".Rdata", "", name)
      name = gsub("7_", "", name)
      print(name)

      print(i)
      load(file)
      assign(name, result)
    
    i = i + 1
  }


  #** ••• SIMULATED ••• **
  #Father + sister thresholds 
  df1 <- summarize_paternity_results(PatFatSis_SimTrue_Route1)
  df2 <- summarize_paternity_results(PatFatSis_NoGE_SNP2k_Route1)
  df3 <- summarize_paternity_results(PatFatSis_NoGE_SNP50k_Route1)
  df4 <- summarize_paternity_results(PatFatSis_WithGE_SNP2k_Route1)
  df5 <- summarize_paternity_results(PatFatSis_WithGE_SNP50k_Route1)
  PatFatSis_Route1_summary <- rbind(df1,df2,df3,df4,df5)
  
  

  df6 <- summarize_paternity_results(PatFatSis_SimTrue_Route2)
  df7 <- summarize_paternity_results(PatFatSis_NoGE_SNP2k_Route2)
  df8 <- summarize_paternity_results(PatFatSis_NoGE_SNP50k_Route2)
  df9 <- summarize_paternity_results(PatFatSis_WithGE_SNP2k_Route2)
  df10 <- summarize_paternity_results(PatFatSis_WithGE_SNP50k_Route2)
  PatFatSis_Route2_summary <- rbind(df6,df7,df8,df9,df10)
  

  PatFatSis_summary_rep <- rbind(PatFatSis_Route1_summary, PatFatSis_Route2_summary)
  PatFatSis_summary_rep$Rep <- Rep
  PatFatSis_summary <- rbind(PatFatSis_summary, PatFatSis_summary_rep)

  df11 <- summarize_sisterPaternity_results(PatSis_SimTrue_Route1)
  df12 <- summarize_sisterPaternity_results(PatSis_NoGE_SNP2k_Route1)
  df13 <- summarize_sisterPaternity_results(PatSis_NoGE_SNP50k_Route1)
  df14 <- summarize_sisterPaternity_results(PatSis_WithGE_SNP2k_Route1)
  df15 <- summarize_sisterPaternity_results(PatSis_WithGE_SNP50k_Route1)
  PatSis_Route1_summary <- rbind(df11,df12,df13,df14,df15)

  df16 <- summarize_sisterPaternity_results(PatSis_SimTrue_Route2)
  df17 <- summarize_sisterPaternity_results(PatSis_NoGE_SNP2k_Route2)
  df18 <- summarize_sisterPaternity_results(PatSis_NoGE_SNP50k_Route2)
  df19 <- summarize_sisterPaternity_results(PatSis_WithGE_SNP2k_Route2)
  df20 <- summarize_sisterPaternity_results(PatSis_WithGE_SNP50k_Route2)
  PatSis_Route2_summary <- rbind(df16,df17,df18,df19,df20)

  PatSis_summary_rep <- rbind(PatSis_Route1_summary, PatSis_Route2_summary)
  PatSis_summary_rep$Rep <- Rep
  PatSis_summary <- rbind(PatSis_summary, PatSis_summary_rep)

  dir.create("Outputs/PatrilineDetermination", showWarnings = FALSE)
  write.csv(PatSis_summary_rep, "Outputs/PatrilineDetermination/Patrilines_PatSis_summary.csv", quote=F, row.names=F)
  write.csv(PatFatSis_summary_rep, "Outputs/PatrilineDetermination/Patrilines_PatFatSis_summary.csv", quote=F, row.names=F)
  #save.image(file = "Pipeline/8_SummarisingData.RData")
 }

write.csv(PatSis_summary, paste0(workingDir, "/PatSis_summary.csv"), quote=F, row.names=F)
write.csv(PatFatSis_summary, paste0(workingDir, "/PatFatSis_summary.csv"), quote=F, row.names=F)

PatSis_summary <- read.csv(paste0(workingDir, "/PatSis_summary.csv"))
PatFatSis_summary <- read.csv(paste0(workingDir, "/PatFatSis_summary.csv"))

PatSis_summaryA <- PatSis_summary %>%
  group_by(sister_thresholds, data_type, haplo_assignment_type) %>%
  summarise(
    mean_correct_n_patrilines = mean(avg_correct_n_patrilines, na.rm = TRUE),
    sd_correct_n_patrilines = sd(avg_correct_n_patrilines, na.rm = TRUE),
    mean_patrilines_times_dpq = mean(avg_patrilines_times_dpq, na.rm = TRUE),
    sd_patrilines_times_dpq = sd(avg_patrilines_times_dpq, na.rm = TRUE),
    mean_father_per_sister_group = mean(avg_nFather_per_sisterGroup, na.rm = TRUE),
    sd_father_per_sister_group = sd(avg_nFather_per_sisterGroup, na.rm = TRUE),
    .groups = "drop"
  )


PatFatSis_summaryA <- PatFatSis_summary %>%
  group_by(father_thresholds, sister_thresholds, data_type, haplo_assignment_type) %>%
  summarise(
    mean_correct_n_patrilines = mean(avg_correct_n_patrilines, na.rm = TRUE),
    sd_correct_n_patrilines = sd(avg_correct_n_patrilines, na.rm = TRUE),
    mean_patrilines_times_dpq = mean(avg_patrilines_times_dpq, na.rm = TRUE),
    sd_patrilines_times_dpq = sd(avg_patrilines_times_dpq, na.rm = TRUE),
    mean_correct_patrilines = mean(avg_correct_patrilines, na.rm = TRUE),
    sd_correct_patrilines = sd(avg_correct_patrilines, na.rm = TRUE),
    .groups = "drop"
  )

PatSis_summaryA$GEType = sapply(strsplit(PatSis_summaryA$data_type, "_"), '[', 1)
PatSis_summaryA$SNPArray = sapply(strsplit(PatSis_summaryA$data_type, "_"), '[', 2)
PatFatSis_summaryA$GEType = sapply(strsplit(PatFatSis_summaryA$data_type, "_"), '[', 1)
PatFatSis_summaryA$SNPArray = sapply(strsplit(PatFatSis_summaryA$data_type, "_"), '[', 2)
PatSis_summaryA <- PatSis_summaryA %>% mutate(haplo_assignment_type = recode(haplo_assignment_type, "Used mat pedigree" = "matPed", "Used recon pedigree" = "recPed"))
PatFatSis_summaryA <- PatFatSis_summaryA %>% mutate(haplo_assignment_type = recode(haplo_assignment_type, "Used mat pedigree" = "matPed", "Used recon pedigree" = "recPed"))



#################### --- PATSIS --- ################################3
library(tidyr)
PatSis_summaryA_meanL <- PatSis_summaryA |> 
  select(-c(sd_correct_n_patrilines, sd_patrilines_times_dpq, 
    sd_father_per_sister_group)) |> 
  pivot_longer(cols = 
  c(mean_correct_n_patrilines, mean_patrilines_times_dpq, 
    mean_father_per_sister_group), values_to = "Mean", names_to = "Measure")
PatSis_summaryA_meanL$Measure <- gsub("mean_", "", PatSis_summaryA_meanL$Measure)

PatSis_summaryA_sdL <- PatSis_summaryA |> 
  select(-c(mean_correct_n_patrilines, mean_patrilines_times_dpq, mean_father_per_sister_group)) |> 
  pivot_longer(cols = 
  c(sd_correct_n_patrilines, sd_patrilines_times_dpq, sd_father_per_sister_group), values_to = "Sd", names_to = "Measure")
PatSis_summaryA_sdL$Measure <- gsub("sd_", "", PatSis_summaryA_sdL$Measure)

table(PatSis_summaryA_meanL$Measure)
PatSis_summaryA_L <- merge(PatSis_summaryA_meanL, PatSis_summaryA_sdL, 
  by = c("sister_thresholds", "data_type", "haplo_assignment_type", "GEType", "SNPArray", "Measure"))
PatSis_summaryA_L$sister_thresholds <- as.factor(PatSis_summaryA_L$sister_thresholds)

PatSis_summaryA_L$SNPArray[PatSis_summaryA_L$data_type == "True"] <- "True"
PatSis_summaryA_L$SNPArray <- factor(PatSis_summaryA_L$SNPArray, levels = c("True", "2k", "50k"))
PatSis_summaryA_L$GEType <- factor(PatSis_summaryA_L$GEType, levels = c("True", "NoGE", "WithGE"))

PatSis_summaryA_L <- PatSis_summaryA_L |> mutate(Measure_plot = recode(Measure, 
  "correct_n_patrilines" = "%\nnPatrilines",
   "patrilines_times_dpq" = "Patriline/\nDPQ",
   "father_per_sister_group" = "Father/\nSister Group"
))

plot_dat <- PatSis_summaryA_L %>%
  filter(
    sister_thresholds %in% c(0.75, 0.85, 0.9, 0.95, 1.0)
  ) %>%
  filter(
    (GEType == "True" & SNPArray == "True") |
    (GEType %in% c("NoGE", "WithGE") & SNPArray %in% c("2k", "50k"))
  ) |> 
  mutate(Sister_threshold = sister_thresholds)

palette <- c(
  "#4E79A7",  # muted blue
  "#F28E2B",  # muted orange
  "#59A14F",  # muted green
  "#B07AA1",  # muted purple
  "#E15759"   # muted red
)
                   
plot_dat$facet_cols = paste0(plot_dat$haplo_assignment_type, "\n", plot_dat$GEType)
plot_dat$hline <- case_when(
  plot_dat$Measure %in% c(
    "father_per_sister_group",
    "patrilines_times_dpq"
  ) ~ 1,
  TRUE ~ 100
)
plot_dat <- plot_dat |> mutate(GETypePlot = recode(GEType, "NoGE" = "Sim_NoGE", "WithGE" = "Sim_WithGE", "True" = "Sim_True"))
patsis_plot <- ggplot(plot_dat,
  aes(x = SNPArray, y = Mean, 
  fill = Sister_threshold, colour = Sister_threshold)) + 
  #geom_point(size = 4)+
  geom_col(position = "dodge") + 
  facet_grid(cols = vars(haplo_assignment_type, GETypePlot), rows = vars(Measure_plot), 
scales = "free", space = "free_x") + 
  theme_bw(base_size = 24) + 
  scale_fill_manual(values = palette) + 
  scale_colour_manual(values = palette) + 
  geom_hline(aes(yintercept = hline), linetype = "dashed", colour = "#464a4a", linewidth = 1) +
  xlab("Data") + ylab("") + 
  theme(legend.position = "bottom")

setwd(workingDir)
png("Patriline_determination_PatSis.png", width = 1800, height = 800)
patsis_plot
dev.off()


############# --- PATFATSIS --- ###################################
library(tidyr)
PatFatSis_summaryA_meanL <- PatFatSis_summaryA |> 
  select(-c(sd_correct_n_patrilines, sd_patrilines_times_dpq, sd_correct_patrilines)) |> 
  pivot_longer(cols = 
  c(mean_correct_n_patrilines, mean_patrilines_times_dpq, mean_correct_patrilines), values_to = "Mean", names_to = "Measure")
PatFatSis_summaryA_meanL$Measure <- gsub("mean_", "", PatFatSis_summaryA_meanL$Measure)

PatFatSis_summaryA_sdL <- PatFatSis_summaryA |> 
  select(-c(mean_correct_n_patrilines, mean_correct_patrilines, mean_patrilines_times_dpq)) |> 
  pivot_longer(cols = 
  c(sd_correct_n_patrilines, sd_correct_patrilines, sd_patrilines_times_dpq), values_to = "Sd", names_to = "Measure")
PatFatSis_summaryA_sdL$Measure <- gsub("sd_", "", PatFatSis_summaryA_sdL$Measure)

table(PatFatSis_summaryA_meanL$Measure)
PatFatSis_summaryA_L <- merge(PatFatSis_summaryA_meanL, PatFatSis_summaryA_sdL, 
  by = c("sister_thresholds", "father_thresholds", "data_type", "haplo_assignment_type", "GEType", "SNPArray", "Measure"))
PatFatSis_summaryA_L$sister_thresholds <- as.factor(PatFatSis_summaryA_L$sister_thresholds)
PatFatSis_summaryA_L$father_thresholds <- as.factor(PatFatSis_summaryA_L$father_thresholds)

PatFatSis_summaryA_L$SNPArray[PatFatSis_summaryA_L$data_type == "True"] <- "True"
PatFatSis_summaryA_L$SNPArray <- factor(PatFatSis_summaryA_L$SNPArray, levels = c("True", "2k", "50k"))
PatFatSis_summaryA_L$GEType <- factor(PatFatSis_summaryA_L$GEType, levels = c("True", "NoGE", "WithGE"))

PatFatSis_summaryA_L <- PatFatSis_summaryA_L |> mutate(Measure_plot = recode(Measure, 
  "correct_n_patrilines" = "%\nnPatrilines",
   "patrilines_times_dpq" = "Patriline/\nDPQ",
  "correct_patrilines" = "Accuracy\nPatriline"))

plot_dat <- PatFatSis_summaryA_L %>%
  filter(
    sister_thresholds %in% c(0.75, 0.85, 0.9, 0.95, 1.0)
  ) %>%
  filter(
    (GEType == "True" & SNPArray == "True") |
    (GEType %in% c("NoGE", "WithGE") & SNPArray %in% c("2k", "50k"))
  ) |> 
  mutate(Sister_threshold = sister_thresholds) |> 
  mutate(Father_threshold = father_thresholds)

palette <- c(
  "#4E79A7",  # muted blue
  "#F28E2B",  # muted orange
  "#59A14F",  # muted green
  "#B07AA1",  # muted purple
  "#E15759"   # muted red
)

plot_dat$facet_cols = paste0(plot_dat$haplo_assignment_type, "\n", plot_dat$GEType)
plot_dat$hline <- case_when(
  plot_dat$Measure %in% c(
    "correct_n_patrilines",
    "correct_patrilines"
  ) ~ 100,
  TRUE ~ 1
)
plot_dat <- plot_dat |> mutate(GETypePlot = recode(GEType, "NoGE" = "Sim_NoGE", "WithGE" = "Sim_WithGE", "True" = "Sim_True"))
patfatsis_plot <- ggplot(plot_dat,
  aes(x = SNPArray, y = Mean, 
  fill = Sister_threshold, colour = Sister_threshold, alpha = Father_threshold)) + 
  #geom_point(size = 4)+
  geom_col(position = "dodge") + 
  facet_grid(cols = vars(haplo_assignment_type, GETypePlot), rows = vars(Measure_plot), 
scales = "free", space = "free_x") + 
  theme_bw(base_size = 24) + 
  scale_fill_manual(values = palette) + 
  scale_colour_manual(values = palette)  + 
  geom_hline(aes(yintercept = hline), linetype = "dashed", colour = "#464a4a", linewidth = 1) +
  xlab("Data") + ylab("") + 
  theme(legend.position = "bottom", legend.box="vertical") 


setwd(workingDir)
png("Patriline_determination_PatFatSis.png", width = 1800, height = 800)
patfatsis_plot
dev.off()
