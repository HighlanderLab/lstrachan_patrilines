


#******************************************************************************
# Prepare files for AlphaAssign pedigree reconstruction software
#******************************************************************************

#                     1. SIMULATED DATA 
#******************************************************************************

##############.  Prepping the input files ##################
rm(list = ls()) #Lets clear workspace just incase something sneaks in
############### Processing AlphaAssign output #################

Process_AlphaAssign_output <- function(GE = NULL, True_pedigree = NULL, nOffspring = NULL) {
  
  if (isTRUE(GE)) {

    test <- "WithGE"
    filetype <- "_WithGE"
  } else if (isFALSE(GE)) {

    test <- "NoGE"
    filetype <- "_NoGE"
  } else {
    stop("GE must be TRUE or FALSE")
  }
  
  results_list <- list()
  
  for (n in 1:5){
    
    # FIX: proper file name construction
    file_name <- paste0("AlphaGeno_SNP", n, filetype, ".sires")
    
    Alpha_output <- read.table(file_name, header = TRUE)
    
    Sires_assigned <- Alpha_output[Alpha_output$chosen == 1, ]
    nSires_assigned <- nrow(Sires_assigned)
    
    Offspring_and_candidateParent <- Sires_assigned[, c(1,2)]
    colnames(Offspring_and_candidateParent) <- c("id", "sire")
    
    compare_sires <- merge(
      Offspring_and_candidateParent,
      True_pedigree,
      by = "id",
      suffixes = c("_assigned", "_known")
    )
    
    nCorrect_sires <- sum(
      compare_sires$sire_assigned == compare_sires$sire_known,
      na.rm = TRUE
    )
    
    df <- data.frame(
      Test = test,
      nOffspring = length(unique(Alpha_output$id)),
      nDpc = length(unique(Alpha_output$candidate)),
      SNP_group = n,
      nSires_assigned = nSires_assigned,
      nCorrect_sires = nCorrect_sires,
      Software = "AlphaAssign"
    )
    
    results_list[[n]] <- df
  }
  
  # bind all iterations into one dataframe
  final_df <- do.call(rbind, results_list)
  
  return(final_df)
}
ped_to_raw <- function(ped_file, map_file, output_file) {
  # Read the .ped file
  ped_data <- read.table(ped_file, header = FALSE, stringsAsFactors = FALSE)
  map_data <- read.table(map_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Extract the first 6 columns (family ID, individual ID, paternal ID, maternal ID, sex, phenotype)
  ped_info <- ped_data[, 1:6] 
  colnames(ped_info) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
  ped_snps <- ped_data[, -c(1:6)] # Extract the SNP genotype data
  



  # Create new dataframe with paired pasted columns
  ped_snps_raw <- data.frame(
    lapply(seq(1, ncol(ped_snps), by = 2), function(i) {
      paste(ped_snps[[i]], ped_snps[[i+1]], sep = "")
    })
  )
  colnames(ped_snps_raw) <- map_data$V2

  map <- c("11" = 0, "22" = 2, "12" = 1, "21" = 1, "00" = 9)

  ped_snps_raw[] <- lapply(ped_snps_raw, function(col) {
    map[col]
  })

  ped_snps_raw_plink <- cbind(ped_info, ped_snps_raw)
  write.table(ped_snps_raw_plink, file = output_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
}


args = commandArgs(trailingOnly=TRUE)
Rep = args[1]
workingDir = args[2]
softwareDir = args[3]

repDir = paste0(workingDir, "/SimRep", Rep, "/")

# Set working directory
setwd(repDir)

pathToPlink <- softwareDir

#Get known DPQs from pedigree 
pedigree_file <- read.csv("Data/worker_pedigree.csv")

#Known Dpcs
Worker_known <- data.frame(id = pedigree_file$id,
                          dam = pedigree_file$mother,
                          sire = pedigree_file$dpc)
Mothers_known <- data.frame(id = unique(pedigree_file$mother),
                            dam = rep(NA, length(unique(pedigree_file$mother))),
                            sire = rep(NA, length(unique(pedigree_file$mother))))

Dpc_known <- data.frame(id = unique(pedigree_file$dpc),
                        dam = rep(NA, length(unique(pedigree_file$dpc))),
                        sire = rep(NA, length(unique(pedigree_file$dpc))))

True_pedigree <- rbind(Mothers_known, Dpc_known, Worker_known)
dir.create(paste0(repDir, "/Data/AlphaAssign"), showWarnings = FALSE)
write.table(True_pedigree, file = paste0("Data/AlphaAssign/SimData_TruePedigree.csv"), sep = ",", quote = F, col.names = T, row.names = F)


no_individuals <- list()
#With genotyping errors: 
for (method in c("NoGE", "WithGE")) {
  no_individuals_methods <- list()
  for (n in 1:5){ 
    setwd(paste0("Data/Sim_", method))

    print(n)
    ped_to_raw(ped_file = paste0("SNP_",n,"_", method, "_QC.ped"), map_file = paste0("SNP_",n,"_", method, "_QC.map"), output_file = paste0("SNP_",n,"_", method, "_QC_RecodeA.raw"))

    AlphaPed <- read.table(paste0("SNP_",n,"_", method, "_QC_RecodeA.raw"), header=TRUE)

    # Do all this within the loop since some queens and DPCs get removed in the QC
    # If a queen is removed, remove her offspring
    pedigree_file <- read.csv(paste0(repDir, "/Data/worker_pedigree.csv"))
    pedigree_file <- pedigree_file[pedigree_file$id %in% AlphaPed$IID,]
    pedigree_file <- pedigree_file[pedigree_file$mother %in% AlphaPed$IID,]
    pedigree_file <- pedigree_file[pedigree_file$dpc %in% AlphaPed$IID,]

    AlphaPed <- AlphaPed[AlphaPed$IID %in% c(pedigree_file$id, pedigree_file$dpc, pedigree_file$mother),]
    AlphaGeno <- AlphaPed[,7:ncol(AlphaPed)]; AlphaGeno[is.na(AlphaGeno)] <- 9
    AlphaGeno_id <- cbind(AlphaPed$IID, AlphaGeno)
    
    write.table(AlphaGeno_id, file=paste0(repDir, "Data/AlphaAssign/AlphaGeno_SNP",n,"_", method, ".txt"), sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE) 
  

    no_individuals_methods[[n]] <- c("NoDams" = length(unique(pedigree_file$mother)),
                                     "NoDpc" = length(unique(pedigree_file$dpc)),
                                     "NoWorkers" = length(unique(pedigree_file$id)))
   
    Alpha_pedigree_sim <- data.frame(id = pedigree_file$id,
                                sire = rep(0, length(pedigree_file$id)),
                                dam = pedigree_file$mother)
    
    setwd(repDir)
    f = data.frame(id = unique(Alpha_pedigree_sim$sire), sire = 0, dam = 0)
    m = data.frame(id = unique(Alpha_pedigree_sim$dam), sire = 0, dam = 0)
    Alpha_pedigree_sim_full <- rbind(f, m, Alpha_pedigree_sim)
    write.table(Alpha_pedigree_sim_full, file = paste0(repDir, "Data/SimData_Pedigree_Full_Maternal_", method,"_SNP", n, ".txt"), quote = F, row.names = F, col.names = F)

    dir.create("Data/AlphaAssign", showWarnings = FALSE)
    write.table(Alpha_pedigree_sim, file = paste0(repDir, "/Data/AlphaAssign/SimData_Pedigree_", method, "_SNP", n, ".txt"), sep = " ", quote = F, col.names = F, row.names = F)

    if (length(unique(pedigree_file$dpc)) != 4) {
      stop("There should be exactly 4 unique DPCs in the pedigree file.")
    }
    dpcs <- unique(pedigree_file$dpc)

    # Create a data frame with id
    Potential_fathers <- data.frame(id = pedigree_file$id)

    # Add one column per DPC dynamically
    for (i in seq_along(dpcs)) {
      Potential_fathers[[paste0("Dpc", i)]] <- rep(dpcs[i], nrow(pedigree_file))
    }

    write.table(Potential_fathers, file = paste0("Data/AlphaAssign/SimData_PotentialFathers_", method, "SNP_", n, ".list"), sep = " ",  quote = F, col.names = F, row.names = F)
  }
  no_individuals[[method]] <- no_individuals_methods
}

tmp = unlist(no_individuals)
no_individuals_df = data.frame(Name = names(tmp), No = tmp)
no_individuals_df$Method <- sapply(strsplit(no_individuals_df$Name, "\\."), "[", 1)
no_individuals_df$Caste <- sapply(strsplit(no_individuals_df$Name, "\\."), "[", 2)
no_individuals_df$SNP_group <- rep(rep(1:5, each = 3), 2)

############### Running AlphAssign in terminal################
  #create a for loop here for all of the SNP array sizes <-------------- 
print("Creating directories for AlphaAssign outputs")
dir.create("Outputs/", showWarnings = FALSE)
dir.create("Outputs/AlphaAssign", showWarnings = FALSE)

print("Running AlphaAssign")

for (method in c("NoGE", "WithGE")) {
  for (n in 1:5){
    print(method)
    print(n)

    system(paste0("AlphaAssign -genotypes ", repDir, "/Data/AlphaAssign/AlphaGeno_SNP", n, "_", method, ".txt ",
    "-potentialsires ", repDir, "/Data/AlphaAssign/SimData_PotentialFathers_", method, "SNP_", n, ".list ", 
    "-pedigree ", repDir, "/Data/AlphaAssign/SimData_Pedigree_", method, "_SNP", n, ".txt ", 
    "-out ", repDir, "/Outputs/AlphaAssign/AlphaGeno_SNP", n, "_", method, " -runtype likelihood"))
  }
}

print("Done running AlphaAssign")

print("Processing AlphaAssign output")
setwd(paste0(repDir, "/Outputs/AlphaAssign/"))
NoGE_Alpha_output <- Process_AlphaAssign_output(GE = FALSE, True_pedigree = True_pedigree, nOffspring = nrow(Worker_known))
WithGE_Alpha_output <- Process_AlphaAssign_output(GE = TRUE, True_pedigree = True_pedigree, nOffspring = nrow(Worker_known))


print("Updating pedigrees with AlphaAssign results")
# Updating the pedigree 
#Original pedigree = Alpha_pedigree
# NoGE updated pedigree ----------------------------------------------------
#Using the 1.7k SNP and the 50k SNP
NoGE_2k_pedigree <- read.table("AlphaGeno_SNP4_NoGE.sires", header=T)
NoGE_2k_pedigree <- NoGE_2k_pedigree[NoGE_2k_pedigree$chosen == 1, ]
NoGE_2k_pedigree <- NoGE_2k_pedigree[, c(1,2)]
colnames(NoGE_2k_pedigree) <- c("id", "sire")

# Match ids to get corresponding sires from NoGE_2k_pedigree
idx <- match(Alpha_pedigree_sim$id, NoGE_2k_pedigree$id)
# Replace only where sire == 0 and a match exists
update_sires <- Alpha_pedigree_sim$sire == 0 & !is.na(idx)
Alpha_pedigree_sim_2K <- Alpha_pedigree_sim
Alpha_pedigree_sim_2K$sire[update_sires] <- NoGE_2k_pedigree$sire[idx[update_sires]]
write.table(Alpha_pedigree_sim_2K, file = paste0(repDir, "Outputs/AlphaAssign/Alpha_pedigree_2k_NoGE.txt"), sep = " ", quote = F, col.names = F, row.names = F)


NoGE_50k_pedigree <- read.table("AlphaGeno_SNP5_NoGE.sires", header=T)
NoGE_50k_pedigree <- NoGE_50k_pedigree[NoGE_50k_pedigree$chosen == 1, ]
NoGE_50k_pedigree <- NoGE_50k_pedigree[, c(1,2)]
colnames(NoGE_50k_pedigree) <- c("id", "sire")

# Match ids to get corresponding sires from NoGE_50k_pedigree
idx <- match(Alpha_pedigree_sim$id, NoGE_50k_pedigree$id)
# Replace only where sire == 0 and a match exists
update_sires <- Alpha_pedigree_sim$sire == 0 & !is.na(idx)
Alpha_pedigree_sim_50K <- Alpha_pedigree_sim
Alpha_pedigree_sim_50K$sire[update_sires] <- NoGE_50k_pedigree$sire[idx[update_sires]]
write.table(Alpha_pedigree_sim_50K, file = paste0(repDir, "Outputs/AlphaAssign/Alpha_pedigree_50k_NoGE.txt"), sep = " ", quote = F, col.names = F, row.names = F)


# WithGE updated pedigree ----------------------------------------------------
#Using the 1.7k SNP and the 50k SNP
WithGE_2k_pedigree <- read.table("AlphaGeno_SNP4_WithGE.sires", header=T)
WithGE_2k_pedigree <- WithGE_2k_pedigree[WithGE_2k_pedigree$chosen == 1, ]
WithGE_2k_pedigree <- WithGE_2k_pedigree[, c(1,2)]
colnames(WithGE_2k_pedigree) <- c("id", "sire")

# Match ids to get corresponding sires from WithGE_2k_pedigree
idx <- match(Alpha_pedigree_sim$id, WithGE_2k_pedigree$id)
# Replace only where sire == 0 and a match exists
Alpha_pedigree_sim_2K_withGE <- Alpha_pedigree_sim
update_sires <- Alpha_pedigree_sim$sire == 0 & !is.na(idx)
Alpha_pedigree_sim_2K_withGE$sire[update_sires] <- WithGE_2k_pedigree$sire[idx[update_sires]]
write.table(Alpha_pedigree_sim_2K_withGE, file = paste0(repDir, "Outputs/AlphaAssign/Alpha_pedigree_2k_WithGE.txt"), sep = " ", quote = F, col.names = F, row.names = F)


WithGE_50k_pedigree <- read.table("AlphaGeno_SNP5_WithGE.sires", header=T)
WithGE_50k_pedigree <- WithGE_50k_pedigree[WithGE_50k_pedigree$chosen == 1, ]
WithGE_50k_pedigree <- WithGE_50k_pedigree[, c(1,2)]
colnames(WithGE_50k_pedigree) <- c("id", "sire")

# Match ids to get corresponding sires from WithGE_50k_pedigree
idx <- match(Alpha_pedigree_sim$id, WithGE_50k_pedigree$id)
# Replace only where sire == 0 and a match exists
update_sires <- Alpha_pedigree_sim$sire == 0 & !is.na(idx)
Alpha_pedigree_sim_50K_withGE <- Alpha_pedigree_sim
Alpha_pedigree_sim_50K_withGE$sire[update_sires] <- WithGE_50k_pedigree$sire[idx[update_sires]]
write.table(Alpha_pedigree_sim_50K_withGE, file = paste0(repDir, "Outputs/AlphaAssign/Alpha_pedigree_50k_WithGE.txt"), sep = " ", quote = F, col.names = F, row.names = F)


AlphaAssign_Summary <- rbind(NoGE_Alpha_output, WithGE_Alpha_output)
write.table(AlphaAssign_Summary, file = paste0(repDir, "Outputs/AlphaAssign/Alpha_summary.txt"), sep = " ", quote = F, col.names = T, row.names = F)

  save.image(file = paste0(repDir, "/Pipeline/2_AlphaAssign.Rdata"))
