#######################################################################################################################
#************************************* Summarising Patriline determination *******************************
#######################################################################################################################
#Creating workable tables to view and summarise results of the patriline determination 

args = commandArgs(trailingOnly=TRUE)
Rep = args[1]
workingDir = args[2]
softwareDir = args[3]

library(ggplot2)
library(gridExtra)
library(cowplot) 
library(dplyr)
library(tidyr)

source(paste0(workingDir, "/ScriptsApril26/6_Checking_Gametic_MendelianSampling_Values_functions.R"))

# Set working directory
setwd(workingDir)


###########################################################################
# First summarise the MS results
###########################################################################

Route1_Mendelian_sim_df_all <- data.frame()
Route2_Mendelian_sim_df_all <- data.frame()

nRep = 10
for (rep in 1:nRep) {
  print(rep)
  load(paste0(workingDir, "/SimRep", rep, "/Outputs/MendelianSampling/Route1_Mend_Sim.Rdata"))
  load(paste0(workingDir, "/SimRep", rep, "/Outputs/MendelianSampling/Route2_Mend_Sim.Rdata"))
  Route1_Mendelian_sim_df$Rep <- rep
  Route1_Mendelian_sim_df_all <- rbind(Route1_Mendelian_sim_df_all, rbind(Route1_Mendelian_sim_df))
  Route2_Mendelian_sim_df$Rep <- rep
  Route2_Mendelian_sim_df_all <- rbind(Route2_Mendelian_sim_df_all, rbind(Route2_Mendelian_sim_df))
}


######### Route 1 ##################
# Convert to long format
df_long_R1 <- Route1_Mendelian_sim_df_all %>%
  pivot_longer(
    cols = c(Maternal_Mendelian, Dpc_Mendelian),
    names_to = "Path",
    values_to = "MS"
  )

grid <- seq(
  min(df_long_R1$MS, na.rm = TRUE),
  max(df_long_R1$MS, na.rm = TRUE),
  length.out = 512
)


# Compute density per replicate per variable
densities_R1 <- df_long_R1 %>%
  group_by(Phasing, Path, Rep) %>%
  do({

    d <- density(.$MS, from = min(grid), to = max(grid), n = length(grid))

    tibble(
      x = grid,   # <-- FORCE SAME X FOR ALL REPLICATES
      y = approx(d$x, d$y, xout = grid)$y
    )

  }) %>%
  ungroup()

# Summarize densities across replicates
summary_density_R1 <- densities_R1 %>%
  group_by(Phasing, Path, x) %>%
  summarise(
    mean_y = mean(y),
    sd = sd(y),
    sem = sd(y) / sqrt(n()),
    .groups = "drop"
  )

mean_lines_R1 <- df_long_R1 %>%
  group_by(Phasing, Path) %>%
  summarise(
    mean_value = mean(MS, na.rm = TRUE),
    sd_value = sd(MS, na.rm = TRUE),
    .groups = "drop"
  )


################### Route 2 #################3
# Convert to long format
df_long_R2 <- Route2_Mendelian_sim_df_all %>%
  pivot_longer(
    cols = c(Maternal_Mendelian),
    names_to = "Path",
    values_to = "MS"
  )

# Compute density per replicate per variable
densities_R2 <- df_long_R2 %>%
  group_by(Phasing, Path, Rep) %>%
  do({

    d <- density(.$MS, from = min(grid), to = max(grid), n = length(grid))

    tibble(
      x = grid,   # <-- FORCE SAME X FOR ALL REPLICATES
      y = approx(d$x, d$y, xout = grid)$y
    )

  }) %>%
  ungroup()

# Summarize densities across replicates
summary_density_R2 <- densities_R2 %>%
  group_by(Phasing, Path, x) %>%
  summarise(
    mean_y = mean(y),
    sd = sd(y),
    sem = sd(y) / sqrt(n()),
    .groups = "drop"
  )

mean_lines_R2 <- df_long_R2 %>%
  group_by(Phasing, Path) %>%
  summarise(
    mean_value = mean(MS, na.rm = TRUE),
    sd_value = sd(MS, na.rm = TRUE),
    .groups = "drop"
  )

############### Real data ####################
load(paste0(workingDir, "/Real_data/Outputs/MendelianSampling/Route1_Mend.Rdata")) #loads Route1_Gametic_Slov
load(paste0(workingDir, "/Real_data/Outputs/MendelianSampling/Route2_Mend.Rdata")) #loads Route2_Gametic_Slov

Route1_Gametic_Slov$Rep <- 1
df_long_R1_real <- Route1_Gametic_Slov %>%
  pivot_longer(
    cols = c(Maternal_Mendelian, Dpc_Mendelian),
    names_to = "Path",
    values_to = "MS"
  )

densities_R1_real <- df_long_R1_real %>%
  group_by(Phasing, Path, Rep) %>%
  do({
    d <- density(.$MS, from = min(grid), to = max(grid), n = length(grid))
    tibble(
      x = grid,   # <-- FORCE SAME X FOR ALL REPLICATES
      y = approx(d$x, d$y, xout = grid)$y
    )
  }) %>% ungroup()

# Summarize densities  - just to match sim datas
summary_density_R1_real <- densities_R1_real %>%
  group_by(Phasing, Path, x) %>%
  summarise(
    mean_y = mean(y),
    sd = sd(y),
    sem = sd(y) / sqrt(n()),
    .groups = "drop"
  )

mean_lines_R1_real <- df_long_R1_real %>%
  group_by(Phasing, Path) %>%
  summarise(
    mean_value = mean(MS, na.rm = TRUE),
    sd_value = sd(MS, na.rm = TRUE),
    .groups = "drop"
  )

###### Route 2
Route2_Gametic_Slov$Rep <- 1
df_long_R2_real <- Route2_Gametic_Slov %>%
  pivot_longer(
    cols = c(Maternal_Mendelian),
    names_to = "Path",
    values_to = "MS"
  )

densities_R2_real <- df_long_R2_real %>%
  group_by(Phasing, Path, Rep) %>%
  do({
    d <- density(.$MS, from = min(grid), to = max(grid), n = length(grid))
    tibble(
      x = grid,   # <-- FORCE SAME X FOR ALL REPLICATES
      y = approx(d$x, d$y, xout = grid)$y
    )
  }) %>%  ungroup()

# Summarize to match sim data
summary_density_R2_real <- densities_R2_real %>%
  group_by(Phasing, Path, x) %>%
  summarise(
    mean_y = mean(y),
    sd = sd(y),
    sem = sd(y) / sqrt(n()),
    .groups = "drop"
  )

mean_lines_R2_real <- df_long_R2_real %>%
  group_by(Phasing, Path) %>%
  summarise(
    mean_value = mean(MS, na.rm = TRUE),
    sd_value = sd(MS, na.rm = TRUE),
    .groups = "drop"
  )



####################### Combine for plotting ########################
densities_R1$Route <- "Route_recPed"
densities_R2$Route <- "Route_matPed"
densities_R1_real$Route <- "Route_recPed"
densities_R2_real$Route <- "Route_matPed"
densities <- do.call(rbind, list(densities_R1, densities_R2, densities_R1_real, densities_R2_real))
densities$Phasing <- gsub("Phased_", "", densities$Phasing)
table(densities$Phasing)
densities$Phasing <- factor(densities$Phasing,
                             levels = c("True", "NoGE_2k", "WithGE_2k", "NoGE_50k", "WithGE_50k", "Real"))

densities$Route <- factor(densities$Route, levels = c("Route_recPed", "Route_matPed"))

summary_density_R1$Route <- "Route_recPed"
summary_density_R2$Route <- "Route_matPed"
summary_density_R1_real$Route <- "Route_recPed"
summary_density_R2_real$Route <- "Route_matPed"
summary_density <- do.call(rbind, list(summary_density_R1, summary_density_R2, summary_density_R1_real, summary_density_R2_real))
summary_density$Phasing <- gsub("Phased_", "", summary_density$Phasing)
table(summary_density$Phasing)
summary_density$Phasing <- factor(summary_density$Phasing,
                             levels = c("True", "NoGE_2k", "WithGE_2k", "NoGE_50k", "WithGE_50k", "Real"))

mean_lines_R1$Route <- "Route_recPed"
mean_lines_R2$Route <- "Route_matPed"
mean_lines_R1_real$Route <- "Route_recPed"
mean_lines_R2_real$Route <- "Route_matPed"
mean_lines <- do.call(rbind, list(mean_lines_R1, mean_lines_R2, mean_lines_R1_real, mean_lines_R2_real))
mean_lines$Phasing <- gsub("Phased_", "", mean_lines$Phasing)  
mean_lines$Phasing <- factor(mean_lines$Phasing,
                             levels = c("True", "NoGE_2k", "WithGE_2k", "NoGE_50k", "WithGE_50k", "Real"))

mean_lines <- mean_lines |> mutate(PhasingPlot = recode(Phasing, 
                                        "NoGE_2k" = "Sim_NoGE_4", "NoGE_50k" = "Sim_NoGE_5", 
                                        "WithGE_2k" = "Sim_WithGE_4", "WithGE_50k" = "Sim_WithGE_5", 
                                        "True" = "Sim_True", "Real" = "Real"))
densities <- densities |> mutate(PhasingPlot = recode(Phasing, 
                                        "NoGE_2k" = "Sim_NoGE_4", "NoGE_50k" = "Sim_NoGE_5", 
                                        "WithGE_2k" = "Sim_WithGE_4", "WithGE_50k" = "Sim_WithGE_5", 
                                        "True" = "Sim_True", "Real" = "Real"))
summary_density <- summary_density |> mutate(PhasingPlot = recode(Phasing, 
                                        "NoGE_2k" = "Sim_NoGE_4", "NoGE_50k" = "Sim_NoGE_5", 
                                        "WithGE_2k" = "Sim_WithGE_4", "WithGE_50k" = "Sim_WithGE_5", 
                                        "True" = "Sim_True", "Real" = "Real"))

MS_plot <- ggplot() +
  geom_line(
    data = densities,
    aes(x = x, y = y,
        group = interaction(Rep, Path),
        color = Path),
    alpha = 0.3
  ) +
  facet_grid(rows = vars(Route), cols = vars(PhasingPlot)) +
    geom_vline(
    data = mean_lines,
    aes(xintercept = mean_value, color = Path),
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_ribbon(
    data = summary_density,
    aes(x = x,
        ymin = mean_y - sem,
        ymax = mean_y + sem,
        fill = Path),
    alpha = 0.25
  ) +
  geom_line(
    data = summary_density,
    aes(x = x, y = mean_y, color = Path),
    linewidth = 1.2
  ) + 
  scale_colour_manual("", values = c("#4E79A7", "#E15759"), labels = c("Paternal", "Maternal")) +
  scale_fill_manual("", values = c("#4E79A7", "#E15759"), labels = c("Paternal", "Maternal")) +
  #geom_vline(aes(xintercept = 0), colour = "black", linetype  = "solid") +
  labs(x = "Mendelian Sampling", y = "Density") +
  theme_bw(base_size = 18)

png("MS_plot.png", width = 1400, height = 400)
print(MS_plot)
dev.off()

