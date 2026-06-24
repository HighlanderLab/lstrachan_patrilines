
Plotting_Gametic_R1 <- function(df, phased_type, plotting_styles) {
  # 1. Validation
  valid_phases <- c("True", "Phased_NoGE_2k", "Phased_WithGE_2k", 
                    "Phased_NoGE_50k", "Phased_WithGE_50k", "Real")
  if (!all(phased_type %in% valid_phases)) stop("Invalid phased_type.")
  
  # 2. Define scales with text wrapping (\n) to prevent cutoff
  color_scale <- scale_color_manual(name = "Assigned parent\nhaplotype",
                                    values = c("Maternal" = "blue", "Paternal" = "red"))
  
  fill_scale <- scale_fill_manual(name = "Assigned parent\nhaplotype",
                                  values = c("Maternal" = "blue", "Paternal" = "red"))
  
  # 3. Generate individual plots
  plots <- list()
  for (phase in phased_type) {
    phase_plots <- list()
    curr_df <- df[df$Phasing == phase, ]
    m_mean <- mean(curr_df$Maternal_Mendelian_mean, na.rm = TRUE)
    s_df <- sd(curr_df$Maternal_Mendelian_mean, na.rm = TRUE)
    p_mean <- mean(curr_df$Dpc_Mendelian_mean, na.rm = TRUE)
    p_sd <- sd(curr_df$Dpc_Mendelian_mean, na.rm = TRUE)
    
    for (style in plotting_styles) {
      if (style == "histogram") {
        p <- ggplot(curr_df, aes(x = Maternal_Mendelian, fill = "Maternal")) +
          geom_histogram(aes(y = after_stat(count)), binwidth = 1, alpha = 0.7, position = 'identity') +
          geom_histogram(aes(x = Dpc_Mendelian, y = after_stat(count), fill = "Paternal"), 
                         binwidth = 1, alpha = 0.7, position = 'identity') +
          geom_vline(xintercept = 0, color = "black", linetype = "solid", linewidth = 0.8) +
          labs(title = paste(phase), x = "Mendelian Sampling", y = "Count") +
          theme_minimal() +
          theme(plot.title = element_text(size = 20),
                axis.title = element_text(size = 16),
                legend.position = "none") +
          fill_scale
      } else if (style == "density") {
        p <- ggplot(curr_df, aes(x = Maternal_Mendelian_mean, color = "Maternal")) +
          geom_density(aes(y = after_stat(scaled)), alpha = 0.7) +
          geom_density(aes(x = Dpc_Mendelian_mean, y = after_stat(scaled), color = "Paternal"), alpha = 0.7) +
          geom_vline(xintercept = m_mean, color = "blue", linetype = "dashed") +
          geom_rect(aes(xmin = m_mean - s_df, xmax = m_mean + s_df), alpha = 0.3) +
          geom_vline(xintercept = p_mean, color = "red", linetype = "dashed") +
          geom_vline(xintercept = 0, color = "black", linetype = "solid", linewidth = 0.8) +
          labs(title = paste(phase), x = "Mendelian Sampling", y = "Density") +
          theme_minimal() +
          theme(plot.title = element_text(size = 20),
                axis.title = element_text(size = 16),
                legend.position = "none") +
          color_scale
      }
      phase_plots[[style]] <- p
    }
    plots[[phase]] <- do.call(cowplot::plot_grid, c(phase_plots, list(nrow = 1)))
  }
  
  # 4. Grid Logic
  num_rows <- if (length(phased_type) > 3) 2 else 1
  combined_plots <- do.call(cowplot::plot_grid, c(plots, list(nrow = num_rows)))
  
  # 5. Extract Legend with extra padding
  example_plot <- ggplot(df[df$Phasing == phased_type[1], ], aes(x = Maternal_Mendelian, fill = "Maternal")) +
    geom_histogram(alpha = 0.7) +
    geom_histogram(aes(x = Dpc_Mendelian, fill = "Paternal"), alpha = 0.7) +
    fill_scale + 
    theme_minimal() +
    theme(legend.position = "right",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          # Add margin to the right of the legend so it doesn't touch the edge
          legend.margin = margin(0, 20, 0, 0)) 
  
  legend <- cowplot::get_legend(example_plot)
  
  # Combine with more room for the legend (0.85 / 0.15)
  final_output <- cowplot::plot_grid(
    combined_plots,
    legend,
    ncol = 2,
    rel_widths = c(0.85, 0.15)
  )
  
  return(final_output)
}

Plotting_Gametic_R2 <- function(df, phased_type, plotting_styles) {
  # 1. Validation
  valid_phases <- c("True", "Phased_NoGE_2k", "Phased_WithGE_2k", 
                    "Phased_NoGE_50k", "Phased_WithGE_50k", "Real")
  if (!all(phased_type %in% valid_phases)) stop("Invalid phased_type.")
  
  # 2. Define scales with text wrapping (\n) to prevent cutoff
  color_scale <- scale_color_manual(name = "Assigned parent\nhaplotype",
                                    values = c("Maternal" = "blue"))
  
  fill_scale <- scale_fill_manual(name = "Assigned parent\nhaplotype",
                                  values = c("Maternal" = "blue"))
  
  # 3. Generate individual plots
  plots <- list()
  for (phase in phased_type) {
    phase_plots <- list()
    curr_df <- df[df$Phasing == phase, ]
    m_mean <- mean(curr_df$Maternal_Mendelian, na.rm = TRUE)
    
    for (style in plotting_styles) {
      if (style == "histogram") {
        p <- ggplot(curr_df, aes(x = Maternal_Mendelian, fill = "Maternal")) +
          geom_histogram(aes(y = after_stat(count)), binwidth = 1, alpha = 0.7, position = 'identity') +
          labs(title = paste(phase), x = "Mendelian Sampling", y = "Count") +
          theme_minimal() +
          theme(plot.title = element_text(size = 20),
                axis.title = element_text(size = 16),
                legend.position = "none") +
          fill_scale
      } else if (style == "density") {
        p <- ggplot(curr_df, aes(x = Maternal_Mendelian, color = "Maternal")) +
          geom_density(aes(y = after_stat(scaled)), alpha = 0.7) +
          geom_vline(xintercept = m_mean, color = "blue", linetype = "dashed") +
          geom_vline(xintercept = 0, color = "black", linetype = "solid", linewidth = 0.8) +
          labs(title = paste(phase), x = "Mendelian Sampling", y = "Density") +
          theme_minimal() +
          theme(plot.title = element_text(size = 20),
                axis.title = element_text(size = 16),
                legend.position = "none") +
          color_scale
      }
      phase_plots[[style]] <- p
    }
    plots[[phase]] <- do.call(cowplot::plot_grid, c(phase_plots, list(nrow = 1)))
  }
  
  # 4. Grid Logic
  num_rows <- if (length(phased_type) > 3) 2 else 1
  combined_plots <- do.call(cowplot::plot_grid, c(plots, list(nrow = num_rows)))
  
  # 5. Extract Legend with extra padding
  example_plot <- ggplot(df[df$Phasing == phased_type[1], ], aes(x = Maternal_Mendelian, fill = "Maternal")) +
    geom_histogram(alpha = 0.7) +
    fill_scale + 
    theme_minimal() +
    theme(legend.position = "right",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          # Add margin to the right of the legend so it doesn't touch the edge
          legend.margin = margin(0, 20, 0, 0)) 
  
  legend <- cowplot::get_legend(example_plot)
  
  # Combine with more room for the legend (0.85 / 0.15)
  final_output <- cowplot::plot_grid(
    combined_plots,
    legend,
    ncol = 2,
    rel_widths = c(0.85, 0.15)
  )
  
  return(final_output)
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
    ri_hap1_sum <- sum(ri_hap1) / ncol(ri_hap1)
    
    pa_hap <- 0.5 * (Dpc_Hap1 + Dpc_Hap2)
    ri_hap2 <- offspring_paternal - pa_hap
    ri_hap2_sum <- sum(ri_hap2) / ncol(ri_hap2)
    
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
    ri_hap1_sum <- sum(ri_hap1) / ncol(ri_hap1)
    
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