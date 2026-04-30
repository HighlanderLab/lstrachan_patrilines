#******************************************************************************
# Running and analysing KING pedigree reconstruction software - REAL DATA
#******************************************************************************
rm(list = ls())

args = commandArgs(trailingOnly=TRUE)
workingDir = args[1]
softwareDir = args[2]
pathToPlink <- softwareDir
pathToBeagle <- softwareDir


# Set working directory
setwd(paste0(workingDir, "/Real_data/"))


library(ggplot2)

dir.create("Outputs/KING", showWarnings = F)
system(paste0(pathToPlink,"/plink2 --bfile Data/Slov_fM_QC --make-king-table --out Outputs/KING/Real_KINGoutput")) 

#process the output file kin.0
Real_KINGoutput <- read.table("Outputs/KING/Real_KINGoutput.kin0")
Real_KINGoutput$Test <- "Real"
colnames(Real_KINGoutput) <- c("FID1","IID1",	"FID2",	"IID2",	"NSNP",	"HETHET",	"IBS0",	"KINSHIP", "Test")

plot_ibs0_kinship <- function(data) {
  # Check if the necessary columns are present in the dataframe
  if (!all(c("IBS0", "KINSHIP") %in% colnames(data))) {
    stop("Dataframe must contain 'IBS0' and 'KINSHIP' columns.")
  }
  
  # Create the scatter plot using ggplot2
  p <- ggplot(data, aes(x = KINSHIP, y = IBS0)) +
    geom_point() +
    labs(title = "Scatter Plot of IBS0 vs KINSHIP",
         x = "KINSHIP",
         y = "IBS0") +
    theme_minimal() +
    scale_y_continuous(expand = c(0, 0)) # Ensure the y-axis starts at 0
  
  # Print the plot
  print(p)
}

plot_ibs0_kinship(Real_KINGoutput)

# Select those with a Kinship >0.2 to more accuratetly find the fathers and IBS < 0.005
Real_KING_0.2 <- Real_KINGoutput[Real_KINGoutput$KINSHIP >= 0.2 & Real_KINGoutput$IBS0 < 0.005,]
nsires_assigned <- nrow(Real_KING_0.2)
plot_ibs0_kinship(Real_KING_0.2)

#Load real pedigree to get dpc ids 
samples <- read.csv("Data/SNP_samples_2022.csv")
dpc_ids <- samples$snp_id[samples$biotype == "dpc"]

worker_ids <- samples$snp_id[samples$biotype == "worker"]

#make sure the first ids are workers only 
Real_KING_filter <- Real_KING_0.2[Real_KING_0.2$IID1 %in% worker_ids,]

Real_KING_dpcs <- Real_KING_filter[Real_KING_filter$IID2 %in% dpc_ids,]
ndpcs_assigned <- nrow(Real_KING_dpcs)

df <- data.frame(
  Test = "Real",
  nDPQs_assigned = ndpcs_assigned,
  Multi_perWorker = length(unique(Real_KING_dpcs$IID1)) != nrow(Real_KING_dpcs),
  Software = "KING"
)

write.csv(df, file = "Outputs/KING/KING_summary.txt", quote=F, row.names=F)
save.image(file = "Pipeline/KINGoutput.Rdata")
