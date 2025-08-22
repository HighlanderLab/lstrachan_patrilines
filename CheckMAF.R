# Initialize a list to store allele frequencies
allele_frequencies <- list()

# Loop through each SNP
for (i in seq(1, ncol(GE_geno), by = 2)) {
  # Extract the alleles for the current SNP
  alleles <- as.vector(unlist(GE_geno[, i:(i+1)]))
  
  # Count the occurrences of each allele
  allele_counts <- table(alleles)
  
  # Calculate the total number of alleles (2 times the number of individuals)
  total_alleles <- sum(allele_counts)
  
  # Calculate allele frequencies
  freq <- allele_counts / total_alleles
  
  # Store the frequencies in the list
  allele_frequencies[[i]] <- freq
}


# Calculate the mean allele frequency
mean_allele_freqSLOV <- mean(unlist(allele_frequencies))
mean_allele_freqnGE <- mean(unlist(allele_frequencies))
mean_allele_freqGE <- mean(unlist(allele_frequencies))



# Convert the list to a data frame for easier viewing
allele_freq_df_Slov <- do.call(rbind, lapply(allele_frequencies, as.data.frame))
allele_freq_df_nGE <- do.call(rbind, lapply(allele_frequencies, as.data.frame))
allele_freq_df_GE <- do.call(rbind, lapply(allele_frequencies, as.data.frame))






# Function to plot histogram for allele frequencies with non-missing values
plot_allele_histogram <- function(allele_freq_df, title) {
  # Remove rows with missing genotypes (coded as 0)
  allele_freq_df <- allele_freq_df[allele_freq_df$alleles != 0, ]
  
  # Melt the data frame for easy plotting
  melted_df <- reshape2::melt(allele_freq_df)
  
  # Plot histogram with increased number of bins
  ggplot(melted_df, aes(x = value, fill = variable)) +
    geom_histogram(binwidth = 0.01, position = "identity", alpha = 0.5) +  # Increase the number of bins
    facet_wrap(~variable, scales = "free", ncol = 2) +
    labs(title = title, x = "Allele Frequency", y = "Frequency") +
    theme_minimal()
}

# Plot histograms for both data frames with non-missing genotype frequencies and increased number of bins
plot_1 <- plot_allele_histogram(allele_freq_df_Slov, "Allele Frequencies - Slov")
plot_2 <- plot_allele_histogram(allele_freq_df_nGE, "Allele Frequencies - nGE")
plot_3 <- plot_allele_histogram(allele_freq_df_GE, "Allele Frequencies - GE")

# Display the plots
plot_1
plot_2
plot_3







# Function to plot density plot for allele frequencies with non-missing values
plot_allele_density <- function(allele_freq_df, title) {

  # Melt the data frame for easy plotting
  melted_df <- reshape2::melt(allele_freq_df)
  
  # Plot density plot
  ggplot(melted_df, aes(x = value, fill = variable)) +
    geom_density(alpha = 0.5) +  # Use density plot
    labs(title = title, x = "Allele Frequency", y = "Density") +
    theme_minimal()
}

# Plot density plots for both data frames with non-missing genotype frequencies
plot_1 <- plot_allele_density(allele_freq_df_Slov, "Allele Frequencies - Slov")
plot_2 <- plot_allele_density(allele_freq_df_Sim, "Allele Frequencies - Sim")

# Display the plots
plot_1
plot_2





# Extract the minor allele frequencies
maf_nGE <- maf(nGE_vcf)
maf_Slov <- maf(Slov_vcf)
maf_GE <- maf(GE_vcf)

# Create data frames for ggplot
df_nGE <- data.frame(MAF = maf_nGE)
df_Slov <- data.frame(MAF = maf_Slov)
df_GE <- data.frame(MAF = maf_GE)

# Plot for Big_vcf
plot_nGE <- ggplot(df_nGE, aes(x = MAF.Frequency)) +
  geom_density() +
  labs(title = "nGE ", x = "Minor Allele Frequency", y = "Density")

# Plot for Slov_vcf
plot_Slov <- ggplot(df_Slov, aes(x = MAF.Frequency)) +
  geom_density() +
  labs(title = "Slov ", x = "Minor Allele Frequency", y = "Density")

# Plot for SIM_vcf
plot_GE <- ggplot(df_GE, aes(x = MAF.Frequency)) +
  geom_density() +
  labs(title = "GE ", x = "Minor Allele Frequency", y = "Density")


# Arrange plots side by side
library(gridExtra)
grid.arrange(plot_Slov, plot_nGE, plot_GE , ncol =3)




#Slov freq density 
Slov_SNPfreq <- read.table("Slov_SNPdensity.frq", header = T)
plot(1:nrow(Slov_SNPfreq), Slov_SNPfreq$MAF, xlab = "SNP Index", ylab = "SNP Density", main = "SNP Density Plot")

