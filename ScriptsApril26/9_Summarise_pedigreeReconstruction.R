library(dplyr)

args = commandArgs(trailingOnly=TRUE)
workingDir = args[1]

setwd(workingDir)
nRep = 10

##########################################################################
#AlphaAssign
##########################################################################
alphaAssign_summary <- data.frame()
for (rep in 1:nRep) {
  print(rep)
  tmp <- read.table(paste0(workingDir, "/SimRep", rep, "/Outputs/AlphaAssign/Alpha_summary.txt"), header=T)
  tmp$Rep <- rep
  tmp$DataType = "Sim"
  alphaAssign_summary <- rbind(alphaAssign_summary, tmp)
}
nrow(alphaAssign_summary)

alphaAssign_summary$Percentage_assigned <- alphaAssign_summary$nSires_assigned / alphaAssign_summary$nOffspring * 100
alphaAssign_summary$Percentage_correct <- alphaAssign_summary$nCorrect_sires / alphaAssign_summary$nSires_assigned * 100

alphaAssign_summaryA <- 
alphaAssign_summary %>% filter(DataType == "Sim") %>% group_by(DataType, Test, SNP_group, Software) %>%
  summarise(Percentage_assigned_mean = mean(Percentage_assigned),
            Percentage_assigned_sd = sd(Percentage_assigned, na.rm=T),
            Percentage_correct_mean = mean(Percentage_correct),
            Percentage_correct_sd = sd(Percentage_correct, na.rm=T))

alphaAssign_real <- read.table("Real_data/Outputs/AlphaAssign/Alpha_summary.txt", header=T)

alphaAssign_results <- rbind(alphaAssign_summaryA, 
                      data.frame(DataType = "Real", Test = "Real", SNP_group = NA, Software = "AlphaAssign", 
                       Percentage_assigned_mean = alphaAssign_real$nSires_assigned / alphaAssign_real$nOffspring * 100, Percentage_assigned_sd = NA, 
                       Percentage_correct_mean = NA, Percentage_correct_sd = NA))


##########################################################################
#KING
##########################################################################
king_summary <- data.frame()
for (rep in 1:nRep) {
  print(rep)
  tmp <- read.table(paste0(workingDir, "/SimRep", rep, "/Outputs/KING/KING_summary.txt"), header=T)
  tmp$Rep <- rep
  tmp$DataType <- "Sim"
  king_summary <- rbind(king_summary, tmp)
}
nrow(king_summary)

king_summary$Percentage_assigned <- king_summary$nDPQs_assigned / king_summary$nWorkers * 100
king_summary$Percentage_correct <- king_summary$nDPQs_correct / king_summary$nDPQs_assigned * 100

king_summaryA <- king_summary %>% filter(DataType == "Sim") %>% group_by(DataType, Test, SNP_group, Software) %>%
  summarise(Percentage_assigned_mean = mean(Percentage_assigned),
            Percentage_assigned_sd = sd(Percentage_assigned, na.rm=T),
            Percentage_correct_mean = mean(Percentage_correct),
            Percentage_correct_sd = sd(Percentage_correct, na.rm=T))

king_real <- read.csv(paste0(workingDir, "/Real_data/Outputs/KING/KING_summary.txt"))

king_results <- rbind(king_summaryA, 
                      data.frame(DataType = "Real", Test = "Real", SNP_group = NA, Software = "KING", 
                       Percentage_assigned_mean = king_real$nDPQs_assigned / king_real$nOffspring * 100, 
                       Percentage_assigned_sd = NA, 
                       Percentage_correct_mean = NA, Percentage_correct_sd = NA))


##########################################################################
#Sequoia
##########################################################################
sequoia_summary <- data.frame()
for (rep in 1:nRep) {
  file = paste0(workingDir, "/SimRep", rep, "/Outputs/Sequoia/SequoiaTable.csv")
  if (file.exists(file)) {
    print(rep)
    tmp <- read.csv(file)
    tmp$Rep <- rep
    tmp$DataType = "Sim"
    sequoia_summary <- rbind(sequoia_summary, tmp)
  }
}
nrow(sequoia_summary)

sequoia_summary$Percentage_assigned <- sequoia_summary$nSires_assigned / sequoia_summary$nOffspring * 100
sequoia_summary$Percentage_correct <- sequoia_summary$nCorrect_sires / sequoia_summary$nSires_assigned * 100

sequoia_summaryA <- 
sequoia_summary %>% filter(DataType == "Sim") %>% group_by(DataType, Test, SNP_group, Software) %>%
  summarise(Percentage_assigned_mean = mean(Percentage_assigned),
            Percentage_assigned_sd = sd(Percentage_assigned, na.rm=T),
            Percentage_correct_mean = mean(Percentage_correct),
            Percentage_correct_sd = sd(Percentage_correct, na.rm=T))

sequoia_real <- read.csv(paste0(workingDir, "/Real_data/Outputs/Sequoia/Real_SequoiaTable.csv"), header=T)

sequoia_results <- rbind(sequoia_summaryA, data.frame(DataType = "Real", Test = "Real", SNP_group = NA, Software = "Sequoia", 
                                                     Percentage_assigned_mean = sequoia_real$nSires_assigned / sequoia_real$nOffspring * 100, 
                                                     Percentage_assigned_sd = NA, 
                                                     Percentage_correct_mean = NA, Percentage_correct_sd = NA))


##########################################################################
#Colony
##########################################################################
colony_summary <- data.frame()
for (rep in 1:nRep) {
  file = paste0(workingDir, "/SimRep", rep, "/Outputs/Colony/Colony_summary.csv")
  if (file.exists(file)) {
    print(rep)
    tmp <- read.csv(file)
    tmp$Rep <- rep
    tmp$DataType = "Sim"
    tmp$Software <- "Colony"
    colony_summary <- rbind(colony_summary, tmp)
  }
}
nrow(colony_summary)
colony_summary <- colony_summary |> rename(Test = Method) |> rename(SNP_group = SNP)

colony_summary$Percentage_assigned <- colony_summary$nAssigned_Sires / colony_summary$nOffspring * 100
colony_summary$Percentage_correct <- colony_summary$nCorrect_Sires / colony_summary$nAssigned_Sires * 100

colony_summaryA <- 
colony_summary %>% filter(DataType == "Sim") %>% group_by(DataType, Test, SNP_group, Software) %>%
  summarise(Percentage_assigned_mean = mean(Percentage_assigned),
            Percentage_assigned_sd = sd(Percentage_assigned),
            Percentage_correct_mean = mean(Percentage_correct),
            Percentage_correct_sd = sd(Percentage_correct))

#sequoia_real <- read.table("Real_data/Outputs/Sequoia/Sequoia_summary.txt", header=T)

colony_results <- rbind(colony_summaryA,
data.frame(DataType = "Real", Test = "Real", SNP_group = NA, Software = "Colony", 
                                                     Percentage_assigned_mean = NA, 
                                                     Percentage_assigned_sd = NA, 
                                                     Percentage_correct_mean = NA, 
                                                     Percentage_correct_sd = NA))
colony_results |> filter(DataType == "Sim" & is.na(SNP_group))



#####################################################################################3
#Put everything together
library(ggplot2)
library(viridis)
library(tidyr)
library(patchwork)

results <- rbind(alphaAssign_results, king_results, sequoia_results, colony_results)
resultsLongMean <- results %>% select(-c(Percentage_assigned_sd, Percentage_correct_sd)) |> 
  pivot_longer(cols = c(Percentage_assigned_mean, Percentage_correct_mean),names_to = "Measure", values_to = "Percentage_mean")
resultsLongMean$Measure <- gsub("_mean", "", resultsLongMean$Measure)

resultsLongSd <- results %>% select(-c(Percentage_assigned_mean, Percentage_correct_mean)) |> 
  pivot_longer(cols = c(Percentage_assigned_sd, Percentage_correct_sd),names_to = "Measure", values_to = "Percentage_sd")
resultsLongSd$Measure <- gsub("_sd", "", resultsLongSd$Measure)

resultsLong <- merge(resultsLongMean, resultsLongSd, by = c("DataType", "Test", "SNP_group", "Software", "Measure"))
resultsLong <- resultsLong |> mutate(
    Measure = recode(Measure,
      "Percentage_assigned" = "Assigned",
      "Percentage_correct" = "Correct"
    ))

resultsLong$SNP_group <- as.factor(resultsLong$SNP_group)
resultsLong <- resultsLong |> mutate(TestPlot = recode(Test, "NoGE" = "Sim_NoGE", "WithGE" = "Sim_WithGE", "Real" = "Real"))

palette <- c(
  "#4E79A7",  # muted blue
  "#F28E2B",  # muted orange
  "#59A14F",  # muted green
  "#B07AA1",  # muted purple
  "#E15759"   # muted red
)
table(resultsLong$Test)
pd <- position_dodge(width = 0.9)

p1 <- resultsLong %>% 
  filter(DataType == "Sim") %>% 
  ggplot(aes(x = SNP_group, y = Percentage_mean, fill = Software)) + 
  geom_col(position = pd) + 
  geom_errorbar(
    aes(
      ymin = Percentage_mean - Percentage_sd,
      ymax = Percentage_mean + Percentage_sd
    ),
    position = pd,
    width = 0.2
  ) +
  facet_grid(rows = vars(Measure), cols = vars(TestPlot)) + 
  scale_fill_manual(values = palette) + 
  theme_bw(base_size = 18) + 
  geom_hline(yintercept = 100,
             linetype = "dashed",
             colour = "#464a4a") +
  labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold"), legend.position = "none") + 
  xlab("SNP Array")

p2 <- resultsLong %>% 
  filter(DataType == "Real", Measure == "Assigned") %>%
  ggplot(aes(x = Software, y = Percentage_mean, fill = Software)) + 
  geom_col(position = "dodge") + 
  scale_fill_manual(values = palette) + 
  theme_bw(base_size = 18) + 
  geom_hline(yintercept = 100, linetype = "dashed", colour = "#464a4a") +
  labs(tag = "B") +
  facet_grid(cols = vars(Test)) +
  ylab("Percentage\nassigned") +
  theme(plot.tag = element_text(face = "bold"), legend.position = "none")


bottom_row <- (plot_spacer() | p2 | plot_spacer()) +
  plot_layout(widths = c(0.15, 0.7, 0.15))

png(paste0(workingDir, "/PedigreeReconstruction_results.png"), width = 1200, height = 800)
(p1 / bottom_row) +
  plot_layout(heights = c(3, 1), guides = "collect") &
  theme(legend.position = "top")
dev.off()
#################################################
#Timing
##################################################

nRep = 10
timing <- data.frame()
for (rep in 1:nRep) {
  print(rep)
  for (software in c("AlphaAssign", "KING", "Sequoia", "Colony")) {
    print(paste0(software, " ", rep))
    tmp <- read.csv(paste0(workingDir, "/SimRep", rep, "/Outputs/", software, "/Timing.csv"))
    timing <- rbind(timing, tmp)
  }
}
table(timing$Rep, timing$Software)

library(ggplot2)
timingA <- timing |> 
  group_by(Method, SNP_group, Software) |> 
  summarise(mean_time = mean(Time, na.rm=T), 
            sd_time = sd(Time, na.rm=T))

timingA$SNP_group <- as.factor(timingA$SNP_group)
timing$SNP_group <- as.factor(timing$SNP_group)
timing |> ggplot(aes(x = Software, y = Time, colour = SNP_group)) + 
  geom_boxplot() +
  facet_wrap(~ Method) +
  theme_minimal()


timingA |> filter(Software == "Colony") |> mutate(mean_time = round(mean_time, 2),
sd_time = round(sd_time, 2))

timingReal <- data.frame()
for (software in c("AlphaAssign", "KING", "Sequoia", "Colony")) {
  print(paste0(software))
  if (file.exists(paste0(workingDir, "/Real_data/Outputs/", software, "/Timing.csv"))) {
    tmp <- read.csv(paste0(workingDir, "/Real_data/Outputs/", software, "/Timing.csv"))
    print(head(tmp))
    timingReal <- rbind(timingReal, tmp)
  }
}
timingReal
