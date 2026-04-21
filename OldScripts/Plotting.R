#download libraries 
library(readr)
library(ggplot2)
library(viridisLite)
library(patchwork)
library(dplyr)
library(gridExtra)


#load all the summary data 
setwd("~/Desktop/Slovenia data/Attempt2/Summary ")

#### Plotting all together ####

Software_summary <- bind_rows(Sequoia, Colony_nSire, AlphaAssign, KING_nSire)
save(Software_summary, file = "Software_summary_Simulated.Rdata")

plot_ALL_percentage <- function(data, Data_Group = NULL, Software = NULL, Test = NULL) {
  # Check if required columns exist in the data
  required_columns <- c("Data_Group", "Test", "SNP_group", "nCorrect_sires", "nSires_assigned", "Software")
  missing_columns <- setdiff(required_columns, colnames(data))
  
  if (length(missing_columns) > 0) {
    stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
  }
  
  # Convert columns to appropriate types if necessary
  data <- data %>%
    mutate(across(c(Data_Group, Software, Test, SNP_group), as.factor)) %>%
    mutate(across(c(nSires_assigned, nCorrect_sires), as.numeric))
  
  # Filter data based on provided arguments
  if (!is.null(Data_Group)) {
    if (!all(Data_Group %in% levels(data$Data_Group))) {
      stop("Invalid Data_Group values provided.")
    }
    data <- data %>% filter(Data_Group %in% Data_Group)
  }
  
  if (!is.null(Software)) {
    if (!all(Software %in% levels(data$Software))) {
      stop("Invalid Software values provided.")
    }
    data <- data %>% filter(Software %in% Software)
  }
  
  if (!is.null(Test)) {
    if (!all(Test %in% levels(data$Test))) {
      stop("Invalid Test values provided.")
    }
    data <- data %>% filter(Test %in% Test)
  }
  
  if (nrow(data) == 0) {
    stop("No data left after filtering. Check your filter criteria.")
  }
  
  # Calculate percentage of correct sire assignments
  data <- data %>%
    mutate(percent_correct_sires = ifelse(nSires_assigned > 0, (nCorrect_sires / nSires_assigned) * 100, NA))
  
  # Split data based on Test
  plot_data_genoerr <- data %>% filter(Test == "GenoErr")
  plot_data_no_genoerr <- data %>% filter(Test == "No_GenoErr")
  
  if (nrow(plot_data_genoerr) == 0) {
    stop("No data for 'GenoErr'.")
  }
  if (nrow(plot_data_no_genoerr) == 0) {
    stop("No data for 'No_GenoErr'.")
  }
  
  # Create plot for nSires_assigned
  p_genoerr_sires <- ggplot(plot_data_genoerr, aes(x = factor(SNP_group), y = nSires_assigned, colour = Software, shape = Software)) +
    geom_point(size = 5, stroke = 2) +
    facet_grid(Data_Group ~ ., scales = "free_y") +
    labs(title = "With Genotyping Errors", x = element_blank(), y = "Number Sires Assigned") +
    scale_y_continuous(limits = c(-10, 300), expand = c(0, 0), breaks = seq(0, 300, by = 20)) +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "lightgrey"),
      panel.grid.minor.y = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title.x = element_text(color = "black", size = 14),
      axis.title.y = element_text(color = "black", size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right",
      text = element_text(size = 14),
      strip.text = element_blank()
    ) +
    scale_shape_manual(values = c("AlphaAssign" = 15, "Colony" = 16, "Sequoia" = 17, "KING" = 18)) +
    scale_color_manual(values = c("AlphaAssign" = "blue", "Colony" = "red", "Sequoia" = "lightgreen", "KING" = "orange"))
  
  p_no_genoerr_sires <- ggplot(plot_data_no_genoerr, aes(x = factor(SNP_group), y = nSires_assigned, colour = Software, shape = Software)) +
    geom_point(size = 5, stroke = 2) +
    facet_grid(Data_Group ~ ., scales = "free_y") +
    labs(title = "Without Genotyping Errors", x = element_blank(), y = "Number Sires Assigned") +
    scale_y_continuous(limits = c(-10, 300), expand = c(0, 0), breaks = seq(0, 300, by = 20)) +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "lightgrey"),
      panel.grid.minor.y = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title.x = element_text(color = "black", size = 14),
      axis.title.y = element_text(color = "black", size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right",
      text = element_text(size = 14),
      strip.text = element_blank()
    ) +
    scale_shape_manual(values = c("AlphaAssign" = 15, "Colony" = 16, "Sequoia" = 17, "KING" = 18)) +
    scale_color_manual(values = c("AlphaAssign" = "blue", "Colony" = "red", "Sequoia" = "lightgreen", "KING" = "orange"))
  
  # Create plot for percent_correct_sires
  p_genoerr_correct <- ggplot(plot_data_genoerr %>% filter(Data_Group %in% c("Simulated")), aes(x = factor(SNP_group), y = percent_correct_sires, colour = Software, shape = Software)) +
    geom_point(size = 5, stroke = 2) +
    facet_grid(Data_Group ~ ., scales = "free_y") +
    labs(x = "SNP group", y = "Accuracy (%)") +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0), breaks = seq(0, 105, by = 10)) +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "lightgrey"),
      panel.grid.minor.y = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title.x = element_text(color = "black", size = 14),
      axis.title.y = element_text(color = "black", size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right",
      text = element_text(size = 14),
      strip.text = element_blank()
    ) +
    scale_shape_manual(values = c("AlphaAssign" = 15, "Colony" = 16, "Sequoia" = 17, "KING" = 18)) +
    scale_color_manual(values = c("AlphaAssign" = "blue", "Colony" = "red", "Sequoia" = "lightgreen", "KING" = "orange"))
  
  p_no_genoerr_correct <- ggplot(plot_data_no_genoerr %>% filter(Data_Group %in% c("Simulated")), aes(x = factor(SNP_group), y = percent_correct_sires, colour = Software, shape = Software)) +
    geom_point(size = 5, stroke = 2) +
    facet_grid(Data_Group ~ ., scales = "free_y") +
    labs(x = "SNP group", y = "Accuracy (%)") +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0), breaks = seq(0, 105, by = 10)) +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "lightgrey"),
      panel.grid.minor.y = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title.x = element_text(color = "black", size = 14),
      axis.title.y = element_text(color = "black", size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right",
      text = element_text(size = 14),
      strip.text = element_blank()
    ) +
    scale_shape_manual(values = c("AlphaAssign" = 15, "Colony" = 16, "Sequoia" = 17, "KING" = 18)) +
    scale_color_manual(values = c("AlphaAssign" = "blue", "Colony" = "red", "Sequoia" = "lightgreen", "KING" = "orange"))
  
  # Combine plots into a single layout with Test as columns
  combined_plot <- (p_no_genoerr_sires | p_genoerr_sires) / (p_no_genoerr_correct | p_genoerr_correct) +
    plot_layout(guides = "collect")
  
  return(combined_plot)
}

plot_ALL_Real <- function(data, Data_Group = NULL, Software = NULL, Test = NULL) {
  # Check if required columns exist in the data
  required_columns <- c("Data_Group", "Test", "SNP_group", "nCorrect_sires", "nSires_assigned", "Software")
  missing_columns <- setdiff(required_columns, colnames(data))
  
  if (length(missing_columns) > 0) {
    stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
  }
  
  # Convert columns to appropriate types if necessary
  data <- data %>%
    mutate(across(c(Data_Group, Software, Test, SNP_group), as.factor)) %>%
    mutate(across(c(nSires_assigned, nCorrect_sires), as.numeric))
  
  # Filter data based on provided arguments
  if (!is.null(Data_Group)) {
    if (!all(Data_Group %in% levels(data$Data_Group))) {
      stop("Invalid Data_Group values provided.")
    }
    data <- data %>% filter(Data_Group %in% Data_Group)
  }
  
  if (!is.null(Software)) {
    if (!all(Software %in% levels(data$Software))) {
      stop("Invalid Software values provided.")
    }
    data <- data %>% filter(Software %in% Software)
  }
  
  if (!is.null(Test)) {
    if (!all(Test %in% levels(data$Test))) {
      stop("Invalid Test values provided.")
    }
    data <- data %>% filter(Test %in% Test)
  }
  
  if (nrow(data) == 0) {
    stop("No data left after filtering. Check your filter criteria.")
  }
  
  # Split data based on Test
  plot_data_genoerr <- data %>% filter(Test == "GenoErr")
  
  # Create bar plot for nSires_assigned
  p_genoerr_sires <- ggplot(plot_data_genoerr, aes(x = Software, y = nSires_assigned, fill = Software)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    facet_grid(Data_Group ~ ., scales = "free_y") +
    labs(title = element_blank(), x = "Software", y = "Number Sires Assigned") +
    scale_y_continuous(limits = c(0, 130), expand = c(0, 0), breaks = seq(0, 300, by = 10)) +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "lightgrey"),
      panel.grid.minor.y = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title.x = element_text(color = "black", size = 16),
      axis.title.y = element_text(color = "black", size = 16),
      axis.text = element_text(size = 16),
      legend.position = "none", # Remove the legend
      text = element_text(size = 16),
      strip.text = element_blank()
    ) +
    scale_fill_manual(values = c("AlphaAssign" = "blue", "Colony" = "red", "Sequoia" = "lightgreen", "KING" = "orange"))
  
  return(p_genoerr_sires)
}

plot_ALL_cluster <- function(data, Data_Group = NULL, Software = NULL, Test = NULL) {
  # Check if required columns exist in the data
  required_columns <- c("Data_Group", "Test", "SNP_group", "nCorrect_sires", "nSires_assigned", "Software")
  missing_columns <- setdiff(required_columns, colnames(data))
  
  if (length(missing_columns) > 0) {
    stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
  }
  
  # Convert columns to appropriate types if necessary
  data <- data %>%
    mutate(across(c(Data_Group, Software, Test, SNP_group), as.factor)) %>%
    mutate(across(c(nSires_assigned, nCorrect_sires), as.numeric))
  
  # Filter data based on provided arguments
  if (!is.null(Data_Group)) {
    if (!all(Data_Group %in% levels(data$Data_Group))) {
      stop("Invalid Data_Group values provided.")
    }
    data <- data %>% filter(Data_Group %in% Data_Group)
  }
  
  if (!is.null(Software)) {
    if (!all(Software %in% levels(data$Software))) {
      stop("Invalid Software values provided.")
    }
    data <- data %>% filter(Software %in% Software)
  }
  
  if (!is.null(Test)) {
    if (!all(Test %in% levels(data$Test))) {
      stop("Invalid Test values provided.")
    }
    data <- data %>% filter(Test %in% Test)
  }
  
  if (nrow(data) == 0) {
    stop("No data left after filtering. Check your filter criteria.")
  }
  
  # Calculate percentage of correct sire assignments
  data <- data %>%
    mutate(percent_correct_sires = ifelse(nSires_assigned > 0, (nCorrect_sires / nSires_assigned) * 100, NA))
  
  # Split data based on Test
  plot_data_genoerr <- data %>% filter(Test == "GenoErr")
  plot_data_no_genoerr <- data %>% filter(Test == "No_GenoErr")
  
  # Create clustered bar plot for nSires_assigned with genotyping errors
  p_genoerr_sires <- ggplot(plot_data_genoerr, aes(x = factor(SNP_group), y = nSires_assigned, fill = Software)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
    facet_grid(Data_Group ~ ., scales = "free_y") +
    labs(title = "With Genotyping Errors", x = "SNP group", y = "Number Sires Assigned") +
    scale_y_continuous(limits = c(-10, 300), expand = c(0, 0), breaks = seq(0, 300, by = 40)) +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "lightgrey"),
      panel.grid.minor.y = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title.x = element_text(color = "black", size = 14),
      axis.title.y = element_text(color = "black", size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right",
      text = element_text(size = 14),
      strip.text = element_blank()
    ) +
    scale_fill_manual(values = c("AlphaAssign" = "blue", "Colony" = "red", "Sequoia" = "lightgreen", "KING" = "orange"))
  
  # Create clustered bar plot for nSires_assigned without genotyping errors
  p_no_genoerr_sires <- ggplot(plot_data_no_genoerr, aes(x = factor(SNP_group), y = nSires_assigned, fill = Software)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
    facet_grid(Data_Group ~ ., scales = "free_y") +
    labs(title = "Without Genotyping Errors", x = "SNP group", y = "Number Sires Assigned") +
    scale_y_continuous(limits = c(-10, 300), expand = c(0, 0), breaks = seq(0, 300, by = 40)) +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "lightgrey"),
      panel.grid.minor.y = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title.x = element_text(color = "black", size = 14),
      axis.title.y = element_text(color = "black", size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right",
      text = element_text(size = 14),
      strip.text = element_blank()
    ) +
    scale_fill_manual(values = c("AlphaAssign" = "blue", "Colony" = "red", "Sequoia" = "lightgreen", "KING" = "orange"))
  
  # Create clustered bar plot for percent_correct_sires with genotyping errors
  p_genoerr_correct <- ggplot(plot_data_genoerr, aes(x = factor(SNP_group), y = percent_correct_sires, fill = Software)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
    facet_grid(Data_Group ~ ., scales = "free_y") +
    labs(x = "SNP group", y = "Accuracy (%)") +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0), breaks = seq(0, 105, by = 10)) +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "lightgrey"),
      panel.grid.minor.y = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title.x = element_text(color = "black", size = 14),
      axis.title.y = element_text(color = "black", size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right",
      text = element_text(size = 14),
      strip.text = element_blank()
    ) +
    scale_fill_manual(values = c("AlphaAssign" = "blue", "Colony" = "red", "Sequoia" = "lightgreen", "KING" = "orange"))
  
  # Create clustered bar plot for percent_correct_sires without genotyping errors
  p_no_genoerr_correct <- ggplot(plot_data_no_genoerr, aes(x = factor(SNP_group), y = percent_correct_sires, fill = Software)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
    facet_grid(Data_Group ~ ., scales = "free_y") +
    labs(x = "SNP group", y = "Accuracy (%)") +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0), breaks = seq(0, 105, by = 10)) +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "lightgrey"),
      panel.grid.minor.y = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title.x = element_text(color = "black", size = 14),
      axis.title.y = element_text(color = "black", size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right",
      text = element_text(size = 14),
      strip.text = element_blank()
    ) +
    scale_fill_manual(values = c("AlphaAssign" = "blue", "Colony" = "red", "Sequoia" = "lightgreen", "KING" = "orange"))
  
  # Combine plots into a single layout with Test as columns
  combined_plot <- (p_no_genoerr_sires | p_genoerr_sires) / (p_no_genoerr_correct | p_genoerr_correct) +
    plot_layout(guides = "collect")
  
  return(combined_plot)
}

plot_ALL_cluster(data = Software_nested, Data_Group = c("Simulated"), Test = c("GenoErr", "No_GenoErr"), Software = c("Sequoia","Colony", "AlphaAssign", "KING"))


plot_ALL_Real(data = Software_Slov, Data_Group = c("Real"), Test = c("GenoErr"), Software = c("Sequoia","Colony", "AlphaAssign", "KING"))

##### Plot Colony the timings ######

plot_colony_time <- function(data) {
  
  # Filter out rows with NA in Total_time
  data <- data[!is.na(data$Total_time), ]
  
  # Arrange SNP_group in ascending order
  data$SNP_group <- factor(data$SNP_group, levels = unique(data$SNP_group)[order(as.numeric(unique(data$SNP_group)))])
  
  # Reorder Data_Group so that "Simulated" is before "Real"
  data$Data_Group <- factor(data$Data_Group, levels = c("Simulated", "Real"))
  
  # Define color and shape mappings for Data_Group
  color_map <- c( 
    "Real" = "red",
    "Simulated" = "#56B4E9" # Light Blue
  )
  shape_map <- c("Real" = 18, "Simulated" = 16)  # You can adjust shapes as needed
  
  # Define color and shape scales
  color_scale <- scale_color_manual(values = color_map)
  shape_scale <- scale_shape_manual(values = shape_map)
  
  # Create the base plot with separate layers for each Data_Group
  base_plot <- ggplot(data, aes(x = SNP_group, y = Total_time)) +
    geom_point(data = subset(data, Data_Group == "Simulated"), aes(color = Data_Group, shape = Data_Group), size = 8) +
    geom_point(data = subset(data, Data_Group == "Real"), aes(color = Data_Group, shape = Data_Group), size = 8) +
    labs(x = "SNP Group", y = "Total time (hrs)", color = "Data Group:", shape = "Data Group:") +
    color_scale + shape_scale +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = "lightblue"),  # Light blue strip background
      axis.line = element_line(color = "black"),
      strip.text.x = element_text(color = "black"),  # Color of facet labels
      panel.grid = element_line(color = "lightgrey"),
      panel.grid.minor.y = element_blank(),
      axis.title.x = element_text(color = "black", size = 16),  # Adjust size as needed
      axis.title.y = element_text(color = "black", size = 16),  # Adjust size as needed
      axis.text = element_text(size = 16),  # Adjust size as needed
      axis.title = element_text(size = 16),  # Adjust size as needed
      legend.position = "right",  # Show legend on the right
      text = element_text(size = 16),  # Adjust size as needed
      strip.text = element_text(size = 16)
    )
  
  # Custom labels for facets
  custom_labeller <- labeller(Test = c(
    "nGE" = "Without Genotyping Errors",
    "GE" = "With Genotyping Errors"
  ))
  
  # Add facet_wrap to the base plot
  final_plot <- base_plot + facet_wrap(~ Test, scales = "free_x", strip.position = "top", labeller = custom_labeller)
  
  # Print the final plot
  print(final_plot)
}

plot_colony_time(Timing)


create_confidence_plot <- function(data) {
  # Transform the Confidence column to a factor for ordered plotting
  data$Confidence <- as.numeric(gsub("%", "", data$Confidence))
  data$CorrectSires <- factor(data$CorrectSires, levels = c(TRUE, FALSE))
  
  # Create the plot
  plot <- ggplot(data, aes(x = factor(Confidence), fill = CorrectSires)) +
    geom_bar(stat = "count", position = "dodge") +
    facet_wrap(~ CorrectSires, scales = "free_y", labeller = labeller(CorrectSires = c("TRUE" = "Correct Assignments", "FALSE" = "Incorrect Assignments"))) +
    labs(
         x = "Confidence (%)",
         y = "Frequency") +
    scale_fill_manual(values = c("TRUE" = "seagreen3", "FALSE" = "coral"), 
                      labels = c("Correct Assignments", "Incorrect Assignments")) +
    theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightblue"),  # Light blue strip background
    axis.line = element_line(color = "black"),
    strip.text.x = element_text(color = "black"),  # Color of facet labels
    panel.grid = element_line(color = "lightgrey"),
    panel.grid.minor.y = element_blank(),
    axis.title.x = element_text(color = "black", size = 16),  # Adjust size as needed
    axis.title.y = element_text(color = "black", size = 16),  # Adjust size as needed
    axis.text = element_text(size = 16),  # Adjust size as needed
    axis.title = element_text(size = 16),  # Adjust size as needed
    legend.position = "none",  # Show legend on the right
    text = element_text(size = 16),  # Adjust size as needed
    strip.text = element_text(size = 16)
  )
  
  # Print the plot
  print(plot)
  
  
}

create_confidence_plot(Colony_confidence_SNP1)


#### Plot number of patrilines determined #########

# Data for Table A and Table B
data <- data.frame(
  Known_Patrilines = c(13, 14, 12, 17, 15, 10, 13, 16),
  True_1 = c(13, 14, 12, 17, 15, 10, 13, 16),
  True_095 = c(13, 14, 12, 17, 15, 10, 13, 16),
  True_085 = c(7, 8, 8, 10, 9, 7, 10, 12),
  Phased_nGE_1 = c(30, 28, 30, 29, 28, 30, 30, 30),
  Phased_nGE_095 = c(13, 14, 12, 17, 15, 10, 13, 16),
  Phased_nGE_085 = c(9, 10, 8, 11, 10, 8, 11, 13),
  Phased_GE_1 = c(30, 30, 29, 30, 30, 30, 30, 30),
  Phased_GE_095 = c(13, 14, 12, 17, 15, 10, 13, 16),
  Phased_GE_085 = c(8, 10, 9, 10, 12, 7, 11, 11)
)

breaks_y <- seq(0, 20, by = 10)  # Adjust the range and step size as needed

# Boxplot with specific y-axis ticks and updated legend title
ggplot(data_long, aes(x=Threshold, y=Difference, fill=Group)) +
  geom_boxplot(outlier.size = 2, outlier.colour = "red") +
  labs(title=element_blank(),
       x="Similarity Threshold",
       y="Difference from Known Patrilines",
       fill="Dataset") +  # Set legend title to "Dataset"
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title.x = element_text(size = 18, margin = margin(t = 15)),
    axis.title.y = element_text(size = 18, margin = margin(r = 15)),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 18, hjust = 0.5, margin = margin(b = 20)),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  coord_cartesian(ylim = range_diff) +  # Adjust y-axis limits based on calculated range
  scale_y_continuous(breaks = breaks_y)  # Set specific y-axis ticks
