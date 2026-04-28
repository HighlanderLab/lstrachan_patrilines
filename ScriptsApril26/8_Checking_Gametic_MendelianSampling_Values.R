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
    m_mean <- mean(curr_df$Maternal_Mendelian, na.rm = TRUE)
    p_mean <- mean(curr_df$Dpc_Mendelian, na.rm = TRUE)
    
    for (style in plotting_styles) {
      if (style == "histogram") {
        p <- ggplot(curr_df, aes(x = Maternal_Mendelian, fill = "Maternal")) +
          geom_histogram(aes(y = after_stat(count)), binwidth = 1, alpha = 0.7, position = 'identity') +
          geom_histogram(aes(x = Dpc_Mendelian, y = after_stat(count), fill = "Paternal"), 
                         binwidth = 1, alpha = 0.7, position = 'identity') +
          labs(title = paste(phase), x = "Mendelian Sampling", y = "Count") +
          theme_minimal() +
          theme(plot.title = element_text(size = 20),
                axis.title = element_text(size = 16),
                legend.position = "none") +
          fill_scale
      } else if (style == "density") {
        p <- ggplot(curr_df, aes(x = Maternal_Mendelian, color = "Maternal")) +
          geom_density(aes(y = after_stat(scaled)), alpha = 0.7) +
          geom_density(aes(x = Dpc_Mendelian, y = after_stat(scaled), color = "Paternal"), alpha = 0.7) +
          geom_vline(xintercept = m_mean, color = "blue", linetype = "dashed") +
          geom_vline(xintercept = p_mean, color = "red", linetype = "dashed") +
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
,  # Initialize an empty data frame to store results
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

# Extract IDs from the target matrix (the smaller one)
target_ids <- sub("_.*", "", rownames(Haplo_R1_SimTrue$real_results_flipped))
# Get unique IDs to make the filtering more efficient
unique_target_ids <- unique(target_ids)

# Extract IDs from the large matrix
true_haplo_ids <- sub("_.*", "", rownames(true_haplotypes))

# Filter true_haplotypes
true_haplotypes_filtered <- true_haplotypes[true_haplo_ids %in% unique_target_ids, ]

# Verify the results
dim(true_haplotypes_filtered)



#** SIMULATED **|
colnames(Haplo_R1_SimTrue$real_results_flipped) <- colnames(true_haplotypes)[1:3200]
Route1_Gametic_True <- calculate_gametic_relatedness_R1(sorted_offspring_haplotypes = Haplo_R1_SimTrue$real_results_flipped , all_haplotypes = true_haplotypes[,1:3200], pedigree = Worker_pedigree)
Route1_Gametic_True$Phasing <- "True"

Route1_Gametic_NoGE_phased_SNP2k <- calculate_gametic_relatedness_R1(sorted_offspring_haplotypes = Haplo_R1_NoGE_SNP2k$results_flipped, all_haplotypes = NoGE_SNP2k_PhasedHaplotypes_recPed, pedigree = Alpha_pedigree_2k_NoGE)
Route1_Gametic_NoGE_phased_SNP2k$Phasing <- "Phased_NoGE_2k"

Route1_Gametic_WithGE_phased_SNP2k <- calculate_gametic_relatedness_R1(sorted_offspring_haplotypes = Haplo_R1_WithGE_SNP2k$results_flipped, all_haplotypes = WithGE_SNP2k_PhasedHaplotypes_recPed, pedigree = Alpha_pedigree_2k_WithGE)
Route1_Gametic_WithGE_phased_SNP2k$Phasing <- "Phased_WithGE_2k"


Route1_Gametic_NoGE_phased_SNP50k <- calculate_gametic_relatedness_R1(sorted_offspring_haplotypes = Haplo_R1_NoGE_SNP50k$results_flipped, all_haplotypes = NoGE_SNP50k_PhasedHaplotypes_recPed, pedigree = Alpha_pedigree_50k_NoGE)
Route1_Gametic_NoGE_phased_SNP50k$Phasing <- "Phased_NoGE_50k"

Route1_Gametic_WithGE_phased_SNP50k <- calculate_gametic_relatedness_R1(sorted_offspring_haplotypes = Haplo_R1_WithGE_SNP50k$results_flipped, all_haplotypes = WithGE_SNP50k_PhasedHaplotypes_recPed, pedigree = Alpha_pedigree_50k_WithGE)
Route1_Gametic_WithGE_phased_SNP50k$Phasing <- "Phased_WithGE_50k"


Route1_Mendelian_sim_df <- rbind(Route1_Gametic_True, Route1_Gametic_NoGE_phased_SNP2k, Route1_Gametic_WithGE_phased_SNP2k, Route1_Gametic_NoGE_phased_SNP50k, Route1_Gametic_WithGE_phased_SNP50k)
save(Route1_Mendelian_sim_df, file = "Outputs/MendelianSampling/Route1_Mend_Sim.Rdata")

#** REAL DATA **|
Slov_pedigree_rec_filtered <- Slov_pedigree_rec[Slov_pedigree_rec$dpc != 0 & Slov_pedigree_rec$mother != 0,]
Route1_Gametic_Slov <- calculate_gametic_relatedness_R1(sorted_offspring_haplotypes = Haplo_R1_Real$results_flipped, all_haplotypes = Slov_PhasedHaplotypes_recPed, pedigree = Slov_pedigree_rec_filtered)
Route1_Gametic_Slov$Phasing <- "Real"


#•• Plotting ••
Plotting_Gametic_R1(df = Route1_Mendelian_sim_df, phased_type = c("True", "Phased_NoGE_2k", "Phased_WithGE_2k", "Phased_NoGE_50k", "Phased_WithGE_50k"), plotting_styles = c("density"))
Plotting_Gametic_R1(df = Route1_Gametic_Slov, phased_type = c("Real"), plotting_styles = "density")
Route1_plot <- Plotting_Gametic_R1(df = rbind(Route1_Mendelian_sim_df, Route1_Gametic_Slov), phased_type = c("True", "Phased_NoGE_2k", "Phased_WithGE_2k", "Phased_NoGE_50k", "Phased_WithGE_50k", "Real"), plotting_styles = "density")

setwd(workingDir)
dir.create("Outputs/MendelianSampling")
ggsave(plot = Route1_plot, filename = "Outputs/MendelianSampling/Route1_test.png", width = 25, height = 10) #Something weird here - need to fix 

############################################################################################################
#**••• Route 2 - with maternal pedigree and maternal info only •••**
############################################################################################################

#** SIMULATED **
colnames(Haplo_R2_SimTrue$real_results_flipped) <- colnames(true_haplotypes)[1:3200]
Route2_Gametic_TRUE <- calculate_gametic_relatedness_R2(sorted_offspring_haplotypes = Haplo_R2_SimTrue$real_results_flipped , all_haplotypes = true_haplotypes[,1:3200], pedigree = Worker_pedigree)
Route2_Gametic_TRUE$Phasing <- "True"

Route2_Gametic_NoGE_phased_2k <- calculate_gametic_relatedness_R2(sorted_offspring_haplotypes = Haplo_R2_NoGE_SNP2k$phased_results_flipped, all_haplotypes = NoGE_SNP2k_PhasedHaplotypes_matPed, pedigree = Worker_pedigree)
Route2_Gametic_NoGE_phased_2k$Phasing <- "Phased_NoGE_2k"

Route2_Gametic_WithGE_phased_2k <- calculate_gametic_relatedness_R2(sorted_offspring_haplotypes = Haplo_R2_WithGE_SNP2k$phased_results_flipped, all_haplotypes = WithGE_SNP2k_PhasedHaplotypes_matPed, pedigree = Worker_pedigree)
Route2_Gametic_WithGE_phased_2k$Phasing <- "Phased_WithGE_2k"

Route2_Gametic_NoGE_phased_50k <- calculate_gametic_relatedness_R2(sorted_offspring_haplotypes = Haplo_R2_NoGE_SNP50k$phased_results_flipped, all_haplotypes = NoGE_SNP50k_PhasedHaplotypes_matPed, pedigree = Worker_pedigree)
Route2_Gametic_NoGE_phased_50k$Phasing <- "Phased_NoGE_50k"

Route2_Gametic_WithGE_phased_50k <- calculate_gametic_relatedness_R2(sorted_offspring_haplotypes = Haplo_R2_WithGE_SNP50k$phased_results_flipped, all_haplotypes = WithGE_SNP50k_PhasedHaplotypes_matPed, pedigree = Worker_pedigree)
Route2_Gametic_WithGE_phased_50k$Phasing <- "Phased_WithGE_50k"


Route2_Mendelian_sim_df <- rbind(Route2_Gametic_TRUE, Route2_Gametic_NoGE_phased_2k, Route2_Gametic_WithGE_phased_2k, Route2_Gametic_NoGE_phased_50k, Route2_Gametic_WithGE_phased_50k)


#** REAL DATA **
Slov_pedigree_mat_filtered <- Slov_pedigree_mat[Slov_pedigree_mat$dam != 0,] # Remove rows with unknown mothers 
tmp <-sub("_.*", "", rownames(Slov_PhasedHaplotypes_matPed))
tmp <- unique(tmp)
Slov_pedigree_mat_filtered <- Slov_pedigree_mat_filtered[Slov_pedigree_mat_filtered$id %in% tmp,]
colnames(Slov_pedigree_mat_filtered) <- c("id", "dpc", "mother")
Route2_Gametic_Slov <- calculate_gametic_relatedness_R2(sorted_offspring_haplotypes = Haplo_R2_Real$phased_results_flipped, all_haplotypes = Slov_PhasedHaplotypes_matPed, pedigree = Slov_pedigree_mat_filtered)
Route2_Gametic_Slov$Phasing <- "Real"

#•• Plotting ••
Plotting_Gametic_R2(df = Route2_Mendelian_sim_df, phased_type = c("Phased_NoGE_2k", "Phased_WithGE_2k", "Phased_NoGE_50k", "Phased_WithGE_50k"), plotting_styles = c("density"))
Plotting_Gametic_R2(df = Route2_Gametic_Slov, phased_type = c("Real"), plotting_styles = c("density"))
Route2_plot <- Plotting_Gametic_R2(df = rbind(Route2_Mendelian_sim_df, Route2_Gametic_Slov), phased_type = c("True","Phased_NoGE_2k", "Phased_WithGE_2k", "Phased_NoGE_50k", "Phased_WithGE_50k", "Real"), plotting_styles = c("density"))

save(Route2_Mendelian_sim_df, file = "Outputs/MendelianSampling/Route2_Mend_Sim.Rdata")


save.image(file = paste0(workingDir, "/Data/Pipeline/8_Checking_Gametic_MendelianSampling.Rdata"))