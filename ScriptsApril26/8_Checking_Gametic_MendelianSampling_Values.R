#######################################################################################################################
#************ Calculating Gametic Mendelian Sampling Values *************
#######################################################################################################################

# To assess the accuracy of haplotype parent-of-origin assignments, we calculated the 
# gametic Mendelian sampling values using both simulated and real data.
#
# According to Wright's pedigree model (Wright, 1921), the genetic value of an individual (g_i)
# is the sum of the average of their parental genetic values and a Mendelian sampling term (r_i):
#
#   g_i = 0.5 * g_m(i) + 0.5 * g_p(i) + r_i
#
# Here, g_m(i) and g_p(i) represent the genetic values of the mother and father, respectively.
# The term r_i captures the deviation due to Mendelian sampling (random inheritance).
#
# For our analysis, we extended this model to the haplotype level to evaluate deviations more precisely.
# We separated the maternal and paternal contributions as follows:
#
#   g_i_maternal = 0.5 * g_m(i)_1 + 0.5 * g_m(i)_2 + r_i_maternal
#   g_i_paternal = 0.5 * g_p(i)_1 + 0.5 * g_p(i)_2 + r_i_paternal
#
# This allows us to isolate the Mendelian sampling deviations (r_i_maternal and r_i_paternal)
# for each parent, using phased haplotype information.

############################################################################################################
#**••• Functions and library  •••**
############################################################################################################

library(ggplot2)
library(gridExtra)
library(cowplot) 

Plotting_Gametic_R1 <- function(df, phased_type, plotting_styles) {
  # Ensure the phased_type is valid
  if (!all(phased_type %in% c("True", "Phased_NoGE_2k", "Phased_WithGE_2k", "Phased_NoGE_50k", "Phased_WithGE_50k", "Real"))) {
    stop("Invalid phased_type. Choose either True, Phased_NoGE_2k, Phased_WithGE_2k, Phased_NoGE_50k, Phased_WithGE_50k or 'Real'.")
  }
  
  # Ensure the plotting_styles are valid
  valid_styles <- c("histogram", "density", "scatter")
  if (!all(plotting_styles %in% valid_styles)) {
    stop("Invalid plotting_styles. Choose either 'histogram', 'density', or 'scatter'.")
  }
  
  # Define color and fill scales for consistency
  color_scale <- scale_color_manual(name = "Assigned parent haplotype",
                                    values = c("Maternal" = "blue", "Paternal" = "red"),
                                    labels = c("Maternal", "Paternal"))
  fill_scale <- scale_fill_manual(name = "Assigned parent haplotype",
                                  values = c("Maternal" = "blue", "Paternal" = "red"),
                                  labels = c("Maternal", "Paternal"))
  
  # Initialize an empty list to store plots
  plots <- list()
  
  # Generate plots for each phased_type
  for (phase in phased_type) {
    # Initialize list to store plots for current phased_type
    phase_plots <- list()
    # Filter data for this specific phase loop
    curr_df <- df[df$Phasing == phase, ]
    # Calculate means for the vertical lines
    m_mean <- mean(curr_df$Maternal_Mendelian, na.rm = TRUE)
    p_mean <- mean(curr_df$Dpc_Mendelian, na.rm = TRUE)
    # Generate plots for each plotting_style
    for (style in plotting_styles) {
      if (style == "histogram") {
        p <- ggplot(df[df$Phasing == phase, ], aes(x = Maternal_Mendelian, fill = "Maternal")) +
          geom_histogram(aes(y = ..count..), binwidth = 1, alpha = 0.7, position = 'identity') +
          geom_histogram(aes(x = Dpc_Mendelian, y = ..count.., fill = "Paternal"), 
                         binwidth = 1, alpha = 0.7, position = 'identity') +
          labs(title = paste(phase),
               x = "Mendelian Sampling", y = "Count") +
          theme_minimal() +
          theme(strip.text = element_blank(),
                plot.title = element_text(size = 20),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text = element_text(size = 14),
                legend.position = "none") +  # Remove legend from individual plots
          fill_scale +  # Apply fill scale
          facet_wrap(~ Phasing, ncol = 1)  # Separate plots by Phasing type
        
      } else if (style == "density") {
        p <- ggplot(df[df$Phasing == phase, ], aes(x = Maternal_Mendelian, color = "Maternal")) +
          geom_density(aes(y = ..scaled..), alpha = 0.7) +
          geom_density(aes(x = Dpc_Mendelian, y = ..scaled.., color = "Paternal"), alpha = 0.7) +
          geom_vline(xintercept = m_mean, color = "blue", linetype = "dashed") +
          geom_vline(xintercept = p_mean, color = "red", linetype = "dashed") +
          labs(title = paste(phase),
               x = "Mendelian Sampling", y = "Density") +
          theme_minimal() +
          theme(strip.text = element_blank(),
                plot.title = element_text(size = 20),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text = element_text(size = 14),
                legend.position = "none") +  # Remove legend from individual plots
          color_scale +  # Apply color scale
          facet_wrap(~ Phasing, ncol = 1)  # Separate plots by Phasing type
        
      }
      # Add plot to current phased_type plots list
      phase_plots[[style]] <- p
    }
    
    # Combine plots for the current phased_type into a single row
    phase_plots_combined <- do.call(cowplot::plot_grid, c(phase_plots, list(nrow = 1)))
    
    # Add combined plots to the main plots list
    plots[[phase]] <- phase_plots_combined
  }
  
  # Combine all phased_type plots into a single row
  combined_plots <- do.call(cowplot::plot_grid, c(plots, list(nrow = 1)))
  
  # Add legend only if "Phased_GE" is in the phased_type
  if ("Phased_GE" %in% phased_type) {
    # Extract legend from one of the "Phased_GE" plots
    example_plot <- ggplot(df[df$Phasing == "Phased_GE", ], aes(x = Maternal_Mendelian, fill = "Maternal")) +
      geom_histogram(aes(y = ..count..), binwidth = 1, alpha = 0.7, position = 'identity') +
      geom_histogram(aes(x = Dpc_Mendelian, y = ..count.., fill = "Paternal"), 
                     binwidth = 1, alpha = 0.7, position = 'identity') +
      fill_scale +  # Apply fill scale
      theme_minimal()
    
    # Extract legend using get_legend from cowplot
    legend <- get_legend(
      example_plot + 
        theme(legend.position = "right",
              legend.title = element_text(size = 16),   # Increase legend title font size
              legend.text = element_text(size = 14))   # Increase legend text font size
    )
    
    # Combine plots and legend
    combined_plots_with_legend <- plot_grid(
      combined_plots,
      legend,
      ncol = 2,
      rel_widths = c(4, 1)  # Adjust widths as needed
    )
  } else {
    # Print combined plots without a legend
    combined_plots_with_legend <- combined_plots
  }
  
  # Print the combined plots with or without the legend
  print(combined_plots_with_legend)
}

Plotting_Gametic_Density_R1 <- function(df, phased_types) {
  
  # 1. Define the color scale
  color_scale <- scale_color_manual(
    name = "Assigned parent haplotype",
    values = c("Maternal" = "blue", "Paternal" = "red")
  )
  
  # 2. Initialize a list to store the individual plots
  plot_list <- list()
  
  # 3. Loop through each phased_type to create a separate plot
  for (phase in phased_types) {
    
    # Filter data for this specific phase
    phase_df <- df[df$Phasing == phase, ]
    
    # Skip if no data for this phase
    if (nrow(phase_df) == 0) next
    
    # Calculate means for the vertical lines
    mean_mat <- mean(phase_df$Maternal_Mendelian, na.rm = TRUE)
    mean_pat <- mean(phase_df$Dpc_Mendelian, na.rm = TRUE)
    
    # Create the plot
    p <- ggplot(phase_df) +
      # Maternal Density (adjust = 0.5 helps with the 'shifted peak' issue)
      geom_density(aes(x = Maternal_Mendelian, color = "Maternal", y = ..scaled..), 
                   size = 1, adjust = 0.5) +
      geom_vline(xintercept = mean_mat, color = "blue", linetype = "dashed", alpha = 0.5) +
      
      # Paternal Density
      geom_density(aes(x = Dpc_Mendelian, color = "Paternal", y = ..scaled..), 
                   size = 1, adjust = 0.5) +
      geom_vline(xintercept = mean_pat, color = "red", linetype = "dashed", alpha = 0.5) +
      
      labs(title = phase, x = "Mendelian Sampling", y = "Scaled Density") +
      theme_minimal() +
      theme(
        legend.position = "none", # Hide individual legends to save space
        plot.title = element_text(size = 14, face = "bold")
      ) +
      color_scale
    
    # Store plot in the list
    plot_list[[phase]] <- p
  }
  
  # 4. Combine all plots into a grid
  # We can set ncol to 2 or 3 depending on how many plots you have
  combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
  
  # 5. Extract a legend from one plot to show at the bottom
  # Create a dummy plot just for the legend
  legend_plot <- ggplot(phase_df) +
    geom_density(aes(x = Maternal_Mendelian, color = "Maternal")) +
    geom_density(aes(x = Dpc_Mendelian, color = "Paternal")) +
    color_scale + theme(legend.position = "bottom")
  
  legend <- cowplot::get_legend(legend_plot)
  
  # Return the grid with the legend at the bottom
  return(cowplot::plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1)))
}

Plotting_Gametic_R2 <- function(df, phased_type, plotting_styles) {
  # Ensure the phased_type is valid
  if (!all(phased_type %in% c("Phased_GE", "Phased_nGE",  "Pre_Phased_nGE", "Pre_Phased_GE","True", "Real"))) {
    stop("Invalid phased_type. Choose either 'Phased_GE', 'Phased_nGE',  'Pre_Phased_nGE', 'Pre_Phased_GE', 'True' or 'Real'.")
  }
  
  # Ensure the plotting_styles are valid
  valid_styles <- c("histogram", "density", "scatter")
  if (!all(plotting_styles %in% valid_styles)) {
    stop("Invalid plotting_styles. Choose either 'histogram', 'density', or 'scatter'.")
  }
  
  # Define color and fill scales for consistency
  color_scale <- scale_color_manual(name = "Assigned parent haplotype",
                                    values = c("Maternal" = "blue"),
                                    labels = c("Maternal"))
  fill_scale <- scale_fill_manual(name = "Assigned parent haplotype",
                                  values = c("Maternal" = "blue"),
                                  labels = c("Maternal"))
  
  # Initialize an empty list to store plots
  plots <- list()
  
  # Generate plots for each phased_type
  for (phase in phased_type) {
    # Initialize list to store plots for current phased_type
    phase_plots <- list()
    
    # Generate plots for each plotting_style
    for (style in plotting_styles) {
      if (style == "histogram") {
        p <- ggplot(df[df$Phasing == phase, ], aes(x = Maternal_Mendelian, fill = "Maternal")) +
          geom_histogram(aes(y = ..count..), binwidth = 1, alpha = 0.7, position = 'identity') +
          labs(title = paste(phase),
               x = "Mendelian Sampling", y = "Count") +
          theme_minimal() +
          theme(strip.text = element_blank(),
                plot.title = element_text(size = 20),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text = element_text(size = 14),
                legend.position = "none") +  # Remove legend from individual plots
          fill_scale +  # Apply fill scale
          facet_wrap(~ Phasing, ncol = 1)  # Separate plots by Phasing type
        
      } else if (style == "density") {
        p <- ggplot(df[df$Phasing == phase, ], aes(x = Maternal_Mendelian, color = "Maternal")) +
          geom_density(aes(y = ..scaled..), alpha = 0.7) +
          labs(title = paste(phase),
               x = "Mendelian Sampling", y = "Density") +
          theme_minimal() +
          theme(strip.text = element_blank(),
                plot.title = element_text(size = 20),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text = element_text(size = 14),
                legend.position = "none") +  # Remove legend from individual plots
          color_scale +  # Apply color scale
          facet_wrap(~ Phasing, ncol = 1)  # Separate plots by Phasing type
        
      }
      # Add plot to current phased_type plots list
      phase_plots[[style]] <- p
    }
    
    # Combine plots for the current phased_type into a single row
    phase_plots_combined <- do.call(cowplot::plot_grid, c(phase_plots, list(nrow = 1)))
    
    # Add combined plots to the main plots list
    plots[[phase]] <- phase_plots_combined
  }
  
  # Combine all phased_type plots into a single row
  combined_plots <- do.call(cowplot::plot_grid, c(plots, list(nrow = 1)))
  
  # Add legend only if "Phased_GE" is in the phased_type
  if ("Phased_GE" %in% phased_type) {
    # Create a plot to extract the legend from
    legend_plot <- ggplot(df[df$Phasing == "Phased_GE", ], aes(x = Maternal_Mendelian, fill = "Maternal")) +
      geom_histogram(aes(y = ..count..), binwidth = 1, alpha = 0.7, position = 'identity') +
      fill_scale +  # Apply fill scale
      theme_minimal() +
      theme(legend.position = "right",
            legend.title = element_text(size = 16),   # Increase legend title font size
            legend.text = element_text(size = 14))   # Increase legend text font size
    
    # Extract legend using get_legend from cowplot
    legend <- cowplot::get_legend(legend_plot)
    
    # Combine plots and legend
    combined_plots_with_legend <- plot_grid(
      combined_plots,
      legend,
      ncol = 2,
      rel_widths = c(4, 1)  # Adjust widths as needed
    )
  } else {
    # Print combined plots without a legend
    combined_plots_with_legend <- combined_plots
  }
  
  # Print the combined plots with or without the legend
  print(combined_plots_with_legend)
}


calculate_gametic_relatedness_R1 <- function(sorted_offspring_haplotypes, all_haplotypes, pedigree) {
  # Initialize an empty data frame to store results
  gametic_results <- data.frame(
    Offspring_ID = character(),
    Mother_ID = character(),
    DPC_ID = character(),
    Maternal_Mendelian = numeric(),
    Dpc_Mendelian = numeric(),
    stringsAsFactors = FALSE
  )
  
  cols <- colnames(sorted_offspring_haplotypes)
  all_haplotypes_cols <- all_haplotypes[, cols, drop = FALSE]
  
  # Loop through each offspring in the pedigree file
  for (i in 1:nrow(pedigree)) {
    print(i)
    offspring_id <- pedigree$id[i]
    mother_id <- pedigree$mother[i]
    dpc_id <- pedigree$dpc[i]
    
    
    #pull out the offspring maternal and paternal genotypes 
    offspring_maternal <- rownames(sorted_offspring_haplotypes)[grep(paste0("^", offspring_id, "_maternal"), rownames(sorted_offspring_haplotypes))]
    offspring_maternal <- sorted_offspring_haplotypes[offspring_maternal, , drop = FALSE]
    
    offspring_paternal <- rownames(sorted_offspring_haplotypes)[grep(paste0("^", offspring_id, "_paternal"), rownames(sorted_offspring_haplotypes))]
    offspring_paternal <- sorted_offspring_haplotypes[offspring_paternal, , drop = FALSE]
    
    # Separate parents into Haplotype 1 and Haplotype 2
    mother_row <- rownames(all_haplotypes_cols)[sapply(rownames(all_haplotypes_cols), function(x) strsplit(x,'_')[[1]][1]) %in% mother_id]
    mother_haplotypes_i <- all_haplotypes_cols[mother_row,]
    Maternal_Hap1 <- as.numeric(mother_haplotypes_i[1,])
    Maternal_Hap2 <- as.numeric(mother_haplotypes_i[2,])
    
    dpc_row <- rownames(all_haplotypes_cols)[sapply(rownames(all_haplotypes_cols), function(x) strsplit(x,'_')[[1]][1]) %in% dpc_id]
    dpc_haplotypes_i <- all_haplotypes_cols[dpc_row,]
    Dpc_Hap1 <- as.numeric(dpc_haplotypes_i[1,])
    Dpc_Hap2 <- as.numeric(dpc_haplotypes_i[2,])
    
    # Calculate parental average and Mendelian sampling for each haplotype
    ma_hap <- 0.5 * (Maternal_Hap1 + Maternal_Hap2)
    ri_hap1 <- offspring_maternal - ma_hap
    ri_hap1_sum <- sum(ri_hap1)
    
    pa_hap <- 0.5 * (Dpc_Hap1 + Dpc_Hap2)
    ri_hap2 <- offspring_paternal - pa_hap
    ri_hap2_sum <- sum(ri_hap2)
    
    # Add the values to the results data frame
    gametic_results <- rbind(gametic_results, data.frame(
      Offspring_ID = offspring_id,
      Mother_ID = mother_id,
      DPC_ID = dpc_id,
      Maternal_Mendelian = ri_hap1_sum,
      Dpc_Mendelian = ri_hap2_sum,
      stringsAsFactors = FALSE
    ))
  }
  
  return(gametic_results)
}

calculate_gametic_relatedness_R2 <- function(sorted_offspring_haplotypes, all_haplotypes, pedigree) {
  # Initialize an empty data frame to store results
  gametic_results <- data.frame(
    Offspring_ID = character(),
    Mother_ID = character(),
    Maternal_Mendelian = numeric(),
    stringsAsFactors = FALSE
  )
  
  cols <- colnames(sorted_offspring_haplotypes)
  all_haplotypes_cols <- all_haplotypes[, cols, drop = FALSE]
  
  # Loop through each offspring in the pedigree file
  for (i in 1:nrow(pedigree)) {
    print(i)
    offspring_id <- pedigree$id[i]
    mother_id <- pedigree$mother[i]
    
    #pull out the offspring maternal and paternal genotypes 
    offspring_maternal <- rownames(sorted_offspring_haplotypes)[grep(paste0("^", offspring_id, "_maternal"), rownames(sorted_offspring_haplotypes))]
    offspring_maternal <- sorted_offspring_haplotypes[offspring_maternal, , drop = FALSE]
    
    offspring_paternal <- rownames(sorted_offspring_haplotypes)[grep(paste0("^", offspring_id, "_paternal"), rownames(sorted_offspring_haplotypes))]
    offspring_paternal <- sorted_offspring_haplotypes[offspring_paternal, , drop = FALSE]
    
    # Separate parents into Haplotype 1 and Haplotype 2
    mother_row <- rownames(all_haplotypes_cols)[sapply(rownames(all_haplotypes_cols), function(x) strsplit(x,'_')[[1]][1]) %in% mother_id]
    mother_haplotypes_i <- all_haplotypes_cols[mother_row,]
    Maternal_Hap1 <- as.numeric(mother_haplotypes_i[1,])
    Maternal_Hap2 <- as.numeric(mother_haplotypes_i[2,])
    
    # Calculate parental average and Mendelian sampling for each haplotype
    ma_hap <- 0.5 * (Maternal_Hap1 + Maternal_Hap2)
    ri_hap1 <- offspring_maternal - ma_hap
    ri_hap1_sum <- sum(ri_hap1)
    
    # Add the values to the results data frame
    gametic_results <- rbind(gametic_results, data.frame(
      Offspring_ID = offspring_id,
      Mother_ID = mother_id,
      Maternal_Mendelian = ri_hap1_sum,
      stringsAsFactors = FALSE
    ))
  }
  
  return(gametic_results)
}


workingDir = "~/Desktop/lstrachan_patrilines"
workingDir = "/home/jana/github/lstrachan_patrilines"
setwd(workingDir)

############################################################################################################
#**••• Get the file we need to run •••**
############################################################################################################

#LOAD Rdata from script 7_Haplotype_ParentAssignments
load("Location where the Rdata is stored")
#Needs to contain the mat and rec pedigrees (including the filtered ones), true haplotypes, map files, and the assigned haplotypes for Route1 and Route2. 

#--- Real Pedigrees --- 

# Pedigree prior to ped reconstruction
colnames(Slov_pedigree_mat) <- c("id","dpc","mother")

#Pedigree After reconstruction
colnames(Slov_pedigree_rec) <- c("id","dpc","mother")


# --- Real assigned haplotypes from 7_Haplotype_ParentAssignments script 
Haplo_R1_Real <- Route1_Real
Haplo_R2_Real <- Route2_Real


# --- Simulated Pedigrees ---
#Pedigree prior to reconstruction 
Sim_pedigree_mat<- read.csv("/Data/worker_pedigree.csv")

#Pedigree post reconstruction (with AlphaAssign)
Alpha_pedigree_2k_NoGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_2k_NoGE.txt")
Alpha_pedigree_50k_NoGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_50k_NoGE.txt")

Alpha_pedigree_2k_WithGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_2k_WithGE.txt")
Alpha_pedigree_50k_WithGE <- read.table("Outputs/AlphaAssign/Alpha_pedigree_50k_WithGE.txt")


# --- Simulated assigned haplotypes from 7_Haplotype_ParentAssignments script 
Haplo_R1_SimTrue <-  Route1_SimTrue
Haplo_R1_NoGE_SNP2k <-  Route1_NoGE_SNP2k
Haplo_R1_WithGE_SNP2k <-  Route1_WithGE_SNP2k
Haplo_R1_NoGE_SNP50k <-  Route1_NoGE_SNP50k
Haplo_R1_WithGE_SNP50k <-  Route1_WithGE_SNP50k


Haplo_R2_SimTrue <-  Route2_SimTrue
Haplo_R2_NoGE_SNP2k <-  Route2_NoGE_SNP2k
Haplo_R2_WithGE_SNP2k <-  Route2_WithGE_SNP2k
Haplo_R2_NoGE_SNP50k <-  Route2_NoGE_SNP50k
Haplo_R2_WithGE_SNP50k <-  Route2_WithGE_SNP50k

############################################################################################################
#**••• Route 1 - with reconstructed pedigree and both parents info •••**
############################################################################################################

#** SIMULATED **|
colnames(Haplo_R1_SimTrue$real_results_flipped) <- colnames(true_haplotypes)
Route1_Gametic_True <- calculate_gametic_relatedness_R1(sorted_offspring_haplotypes = Haplo_R1_SimTrue$real_results_flipped , all_haplotypes = true_haplotypes, pedigree = Sim_pedigree_mat)
Route1_Gametic_True$Phasing <- "True"

Route1_Gametic_NoGE_phased_SNP2k <- calculate_gametic_relatedness_R1(sorted_offspring_haplotypes = Haplo_R1_NoGE_SNP2k$results_flipped, all_haplotypes = NoGE_SNP2k_PhasedHaplotypes_recPed, pedigree = Alpha_pedigree_2k_NoGE)
Route1_Gametic_NoGE_phased_SNP2k$Phasing <- "Phased_NoGE_2k"

Route1_Gametic_WithGE_phased_SNP2k <- calculate_gametic_relatedness_R1(sorted_offspring_haplotypes = Haplo_R1_WithGE_SNP2k$results_flipped, all_haplotypes = WithGE_SNP2k_PhasedHaplotypes_recPed, pedigree = Alpha_pedigree_2k_WithGE)
Route1_Gametic_WithGE_phased_SNP2k$Phasing <- "Phased_WithGE_2k"


Route1_Gametic_NoGE_phased_SNP50k <- calculate_gametic_relatedness_R1(sorted_offspring_haplotypes = Haplo_R1_NoGE_SNP50k$results_flipped, all_haplotypes = NoGE_SNP50k_PhasedHaplotypes_recPed, pedigree = Alpha_pedigree_50k_NoGE)
Route1_Gametic_NoGE_phased_SNP50k$Phasing <- "Phased_NoGE_50k"

Route1_Gametic_WithGE_phased_SNP50k <- calculate_gametic_relatedness_R1(sorted_offspring_haplotypes = Haplo_R1_WithGE_SNP50k$results_flipped, all_haplotypes = WithGE_SNP50k_PhasedHaplotypes_recPed, pedigree = Alpha_pedigree_50k_WithGE)
Route1_Gametic_WithGE_phased_SNP50k$Phasing <- "Phased_WithGE_50k"


Gametic_df <- rbind(Route1_Gametic_True, Route1_Gametic_NoGE_phased_SNP2k, Route1_Gametic_WithGE_phased_SNP2k, Route1_Gametic_NoGE_phased_SNP50k, Route1_Gametic_WithGE_phased_SNP50k)


#** REAL DATA **|
Slov_pedigree_rec_filtered <- Slov_pedigree_rec[Slov_pedigree_rec$dpc != 0 & Slov_pedigree_rec$mother != 0,]
Route1_Gametic_Slov <- calculate_gametic_relatedness_R1(sorted_offspring_haplotypes = Haplo_R1_Real$results_flipped, all_haplotypes = Slov_PhasedHaplotypes_recPed, pedigree = Slov_pedigree_rec_filtered)
Route1_Gametic_Slov$Phasing <- "Real"


#•• Plotting ••
Plotting_Gametic_R1(df = Gametic_df, phased_type = c("True", "Phased_NoGE_2k", "Phased_WithGE_2k", "Phased_NoGE_50k", "Phased_WithGE_50k"), plotting_styles = c("density"))
Plotting_Gametic_R1(df = Route1_Gametic_Slov, phased_type = c("Real"), plotting_styles = "density")


############################################################################################################
#**••• Route 2 - with maternal pedigree and maternal info only •••**
############################################################################################################

#** SIMULATED **
colnames(Route2_True_FLIP$real_results_flipped) <- colnames(true_haplotypes)
Route2_Gametic_TRUE <- calculate_gametic_relatedness_R2(sorted_offspring_haplotypes = Route2_True_FLIP$real_results_flipped , all_haplotypes = true_haplotypes, pedigree = Sim_pedigree_mat)
Route2_Gametic_TRUE$Phasing <- "True"

Route2_Gametic_NoGE_phased_2k <- calculate_gametic_relatedness_R2(sorted_offspring_haplotypes = Haplo_R2_NoGE_SNP2k$phased_results_flipped, all_haplotypes = NoGE_SNP2k_PhasedHaplotypes_matPed, pedigree = Sim_pedigree_mat)
Route2_Gametic_NoGE_phased_2k$Phasing <- "Phased_NoGE_2k"

Route2_Gametic_WithGE_phased_2k <- calculate_gametic_relatedness_R2(sorted_offspring_haplotypes = Haplo_R2_WithGE_SNP2k$phased_results_flipped, all_haplotypes = WithGE_SNP2k_PhasedHaplotypes_matPed, pedigree = Sim_pedigree_mat)
Route2_Gametic_WithGE_phased_2k$Phasing <- "Phased_WithGE_2k"

Route2_Gametic_NoGE_phased_50k <- calculate_gametic_relatedness_R2(sorted_offspring_haplotypes = Haplo_R2_NoGE_SNP50k$phased_results_flipped, all_haplotypes = NoGE_SNP50k_PhasedHaplotypes_matPed, pedigree = Sim_pedigree_mat)
Route2_Gametic_NoGE_phased_50k$Phasing <- "Phased_NoGE_50k"

Route2_Gametic_WithGE_phased_50k <- calculate_gametic_relatedness_R2(sorted_offspring_haplotypes = Haplo_R2_WithGE_SNP50k$phased_results_flipped, all_haplotypes = WithGE_SNP50k_PhasedHaplotypes_matPed, pedigree = Sim_pedigree_mat)
Route2_Gametic_WithGE_phased_50k$Phasing <- "Phased_WithGE_50k"


Gametic_df_R2 <- rbind(Route2_Gametic_True, Route2_Gametic_NoGE_phased_SNP2k, Route2_Gametic_WithGE_phased_SNP2k, Route2_Gametic_NoGE_phased_SNP50k, Route2_Gametic_WithGE_phased_SNP50k)


#** REAL DATA **
Slov_pedigree_mat_filtered <- Slov_pedigree_mat[Slov_pedigree_mat$mother != 0,] # Remove rows with unknown mothers 
tmp <-sub("_.*", "", rownames(Slov_PhasedHaplotypes_matPed))
tmp <- unique(tmp)
Slov_pedigree_mat_filtered <- Slov_pedigree_mat_filtered[Slov_pedigree_mat_filtered$id %in% tmp,]
Route2_Gametic_Slov <- calculate_gametic_relatedness_R2(sorted_offspring_haplotypes = Haplo_R2_Real$phased_results_flipped, all_haplotypes = Slov_PhasedHaplotypes_matPed, pedigree = Slov_pedigree_mat_filtered)
Route2_Gametic_Slov$Phasing <- "Real"

#•• Plotting ••
Plotting_Gametic_R2(df = Gametic_df_R2, phased_type = c("True", "Phased_NoGE_2k", "Phased_WithGE_2k", "Phased_NoGE_50k", "Phased_WithGE_50k"), plotting_styles = c("density"))
Plotting_Gametic_R2(df = Route2_Gametic_Slov, phased_type = c("Real"), plotting_styles = c("density"))
