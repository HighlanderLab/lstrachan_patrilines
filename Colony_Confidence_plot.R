
#COLONY ####

# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(dplyr)

# Function to create a density plot with frequency on y-axis
plot_density <- function(data, data_group) {
  # Filter the data based on Data_Group
  filtered_data <- data %>% filter(Data_Group == data_group)
  
  # Convert Confidence from character to numeric
  filtered_data$Confidence <- as.numeric(gsub("%", "", filtered_data$Confidence))
  
  # Create the histogram plot
  ggplot(filtered_data, aes(x = factor(Confidence), fill = as.factor(SNP_group))) +
    geom_bar() +
    labs(x = "Confidence (%)", y = "Frequency", fill = "SNP Group") +
    scale_fill_brewer(palette = "Set2") +
    facet_wrap(~ SNP_group, scales = "free_x", nrow = 2, ncol = 2) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top",
          strip.text = element_blank())+
    ggtitle(paste("Data Group:", data_group))
}

# Example of using the function
 plot_density(Colony_confidence, "Nested")

 
 KnownParents_dpc <- read_csv("~/Desktop/Slovenia data/Simulated /nested/data/KnownParents_dpc.csv")
 KnownParents_dpc <- KnownParents_dpc[,-2]
 colnames(KnownParents_dpc) <- c("id", "KnownDPC")
 
 # Merge the data frames
 CheckMatch <- merge(Colony_confidence,Known_Dpc, by.x = "OffspringID", by.y = "id", all.x = TRUE)
 CheckMatch$nCorrect <- ifelse(CheckMatch$CandidateID_sim == CheckMatch$KnownDPC, "match", "mismatch")
 
 
 
 # Count the number of matches and mismatches
 num_matches <- sum(UnCol_160$CandidateID == UnCol_160$KnownDPC)
 num_mismatches <- sum(UnCol_160$CandidateID != UnCol_160$KnownDPC)
 

 
 plot_density_matches <- function(data, data_group, SNP_group) {
   
   # Check if only one SNP_group is chosen
   if (length(SNP_group) == 1) {
     # Filter the data based on Data_Group and SNP_group
     filtered_data <- data[data$SNP_group %in% SNP_group & data$Data_Group %in% data_group, ]
     
     # Convert Confidence from character to numeric
     filtered_data$Confidence <- as.numeric(gsub("%", "", filtered_data$Confidence))
     
     # Add the nCorrect column based on the condition
     filtered_data$nCorrect <- ifelse(filtered_data$CandidateID_sim == filtered_data$CandidateID_known, "match", "mismatch")
     
     # Create the histogram plot
     ggplot(filtered_data, aes(x = factor(Confidence), fill = as.factor(nCorrect))) +
       geom_bar() +
       labs(x = "Confidence (%)", y = "Frequency", fill = "nCorrect") +
       scale_fill_brewer(palette = "Set2") +
       facet_grid(SNP_group ~ nCorrect, scales = "free_x", space = "free_x", labeller = labeller(nCorrect = c("match" = "Match", "mismatch" = "Mismatch"))) +
       theme_minimal(base_size = 14) +
       theme(legend.position = "none",
             strip.text = element_text(size = 12, face = "bold"),
             strip.background = element_blank(),
             strip.placement = "outside",
             strip.text.x = element_text(margin = margin(b = 10)),
             strip.text.y = element_blank(),
             panel.border = element_rect(color = "black", fill = NA, size = 1)) +
       ggtitle(paste("Data Group:", data_group, "with SNP Size:", SNP_group))
   } else {
     # Filter the data based on Data_Group and SNP_group
     filtered_data <- data[data$SNP_group %in% SNP_group & data$Data_Group %in% data_group, ]
     
     # Convert Confidence from character to numeric
     filtered_data$Confidence <- as.numeric(gsub("%", "", filtered_data$Confidence))
     
     # Add the nCorrect column based on the condition
     filtered_data$nCorrect <- ifelse(filtered_data$CandidateID_sim == filtered_data$CandidateID_known, "match", "mismatch")
     
     
     # Use the original function provided
     # Create the histogram plot
     ggplot(filtered_data, aes(x = factor(Confidence), fill = as.factor(SNP_group))) +
       geom_bar() +
       labs(x = "Confidence (%)", y = "Frequency", fill = "SNP Size") +
       scale_fill_brewer(palette = "Set2") +
       facet_grid(SNP_group ~ nCorrect, scales = "free_x", space = "free_x", labeller = labeller(nCorrect = c("match" = "Match", "mismatch" = "Mismatch"))) +
       theme_minimal(base_size = 14) +
       theme(legend.position = "top",
             strip.text = element_text(size = 12, face = "bold"),
             strip.background = element_blank(),
             strip.placement = "outside",
             strip.text.x = element_text(margin = margin(b = 10)),
             strip.text.y = element_blank(),
             panel.border = element_rect(color = "black", fill = NA, size = 1)) +
       ggtitle(paste("Data Group:", data_group))
   }
   
 }
 
 plot_density_matches(CheckMatch, data_group = c("Nested","Non-Nested"), SNP_group = c(1))
 

 
 library(ggplot2)
 library(dplyr)
 library(RColorBrewer)
 
 plot_density_matches <- function(data, data_group, SNP_group) {
   
   # Filter the data based on Data_Group and SNP_group
   filtered_data <- data[data$SNP_group %in% SNP_group & data$Data_Group %in% data_group, ]
   
   # Convert Confidence from character to numeric
   filtered_data$Confidence <- as.numeric(gsub("%", "", filtered_data$Confidence))
   
   # Add the nCorrect column based on the condition
   filtered_data$nCorrect <- ifelse(filtered_data$CandidateID_sim == filtered_data$CandidateID_known, "match", "mismatch")
   
   # Check if only one SNP_group is chosen and multiple data_groups are chosen
   if (length(SNP_group) == 1 && length(data_group) > 1) {
     ggplot(filtered_data, aes(x = factor(Confidence), fill = as.factor(nCorrect))) +
       geom_bar() +
       labs(x = "Confidence (%)", y = "Frequency", fill = "nCorrect") +
       scale_fill_brewer(palette = "Set2") +
       facet_grid(Data_Group ~ nCorrect, scales = "free_x", space = "free_x", labeller = labeller(nCorrect = c("match" = "Match", "mismatch" = "Mismatch"))) +
       theme_minimal(base_size = 14) +
       theme(legend.position = "none",
             strip.text = element_text(size = 12, face = "bold"),
             strip.background = element_blank(),
             strip.placement = "outside",
             strip.text.x = element_text(margin = margin(b = 10)),
             strip.text.y = element_text(margin = margin(r = 10)),
             panel.border = element_rect(color = "black", fill = NA, size = 1))
   } else {
     # Check if only one SNP_group is chosen
     if (length(SNP_group) == 1) {
       ggplot(filtered_data, aes(x = factor(Confidence), fill = as.factor(nCorrect))) +
         geom_bar() +
         labs(x = "Confidence (%)", y = "Frequency", fill = "nCorrect") +
         scale_fill_brewer(palette = "Set2") +
         facet_grid(SNP_group ~ nCorrect, scales = "free_x", space = "free_x", labeller = labeller(nCorrect = c("match" = "Match", "mismatch" = "Mismatch"))) +
         theme_minimal(base_size = 14) +
         theme(legend.position = "none",
               strip.text = element_text(size = 12, face = "bold"),
               strip.background = element_blank(),
               strip.placement = "outside",
               strip.text.x = element_text(margin = margin(b = 10)),
               strip.text.y = element_blank(),
               panel.border = element_rect(color = "black", fill = NA, size = 1)) +
         ggtitle(paste("Data Group:", data_group, "with SNP Size:", SNP_group))
     } else {
       # For multiple SNP_groups
       ggplot(filtered_data, aes(x = factor(Confidence), fill = as.factor(SNP_group))) +
         geom_bar() +
         labs(x = "Confidence (%)", y = "Frequency", fill = "SNP Size") +
         scale_fill_brewer(palette = "Set2") +
         facet_grid(SNP_group ~ nCorrect, scales = "free_x", space = "free_x", labeller = labeller(nCorrect = c("match" = "Match", "mismatch" = "Mismatch"))) +
         theme_minimal(base_size = 14) +
         theme(legend.position = "top",
               strip.text = element_text(size = 12, face = "bold"),
               strip.background = element_blank(),
               strip.placement = "outside",
               strip.text.x = element_text(margin = margin(b = 10)),
               strip.text.y = element_blank(),
               panel.border = element_rect(color = "black", fill = NA, size = 1)) +
         ggtitle(paste("Data Group:", data_group))
     }
   }
 }
 
 # Example usage
 # plot_density_matches(data, c("Group1", "Group2"), 1)
 