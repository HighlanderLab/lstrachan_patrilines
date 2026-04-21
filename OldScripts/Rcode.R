library(AlphaSimR)
library(SIMplyBee)
library(readr)
library(ggplot2)
library(nadiv)
library(viridisLite)

#load 
load("Pop_withNoFathers.Rdata")
load("SP_object.Rdata")

#NestedSNP4nGE
GenoMatrix <- read.table("AlphaGenoSNP4_nGE.txt", row.names = 1)
GenoMatrix <- as.matrix(GenoMatrix)


sex_all <- getCasteSex(PopMerged_noFathers, caste = "all", collapse = TRUE)

AlleleFreq <- calcBeeAlleleFreq(x = GenoMatrix, sex = sex_all)

GRMIbs <- calcBeeGRMIbs(x = GenoMatrix, sex = sex_all, alleleFreq = AlleleFreq)

#get haplotypes for IBD calculation
Haplotypes <- pullSnpHaplo(PopMerged_noFathers, simParam = SP)
HaploMatrix <- as.matrix(Haplotypes)

GRMIbd <- calcBeeGRMIbd(x = HaploMatrix)
GRMIbd <- GRMIbd$indiv


computeRelationship_pedigree <- function(pedigree) {
  pedigree <- as.data.frame(pedigree)
  # nadiv needs missing as NA
  pedigree$mother[pedigree$mother == 0] <- NA
  pedigree$father[pedigree$father == 0] <- NA
  pedigree$ID <- rownames(pedigree)
  
  # Females are 1, males are 0
  colnames(pedigree) <- c("Dam", "Sire", "Sex", "ID")
  # The order for nadiv should be ID, Dam, Sire, Sex
  pedigree <- pedigree[, c("ID", "Dam", "Sire", "Sex")]
  pedigree$Sire[pedigree$Sex == 1] <- NA
  
  tmp <- makeS(pedigree = pedigree, heterogametic = "1", returnS = FALSE)
  IBDe <- tmp #$S
  #  dimnames(IBDe) <- list(rownames(pedigree), rownames(pedigree))
  return(IBDe)
}

IBDe <- computeRelationship_pedigree(SP$pedigree)
Sinv <- IBDe$Sinv

#Look at opposing homozygotes 
Pedigree <- read_csv("worker_pedigree.csv")

OffspringID <- as.character(Pedigree$id)
QueenID <- unique(as.character(Pedigree$mother))
DpcID <- unique(as.character(Pedigree$dpc))

x_worker <- GenoMatrix[OffspringID,]
x_dpc <- GenoMatrix[DpcID,]
x_queen <- GenoMatrix[QueenID,]

# Look at the diagonal and off-diagonal 

getS <- function(Sinv, ids, with = ids, diagOnly = FALSE, vector = FALSE) {
  ids <- as.numeric(ids)
  with <- as.numeric(with)
  x <- sparseMatrix(i = ids, j = 1:length(ids), dims = c(nrow(Sinv), length(ids)))
  M1 <- as(x, "dMatrix")
  Sids <- as.matrix(solve(Sinv, M1)[with,])
  if (dim(Sids)[1] == length(ids)){
    rownames(Sids) <- ids
  } else {
    rownames(Sids) <- with
  }
  
  if (dim(Sids)[2] == length(ids)){
    colnames(Sids) <- ids
  } else {
    colnames(Sids) <- with
  }
  if (diagOnly) {
    Sids <- diag(Sids)
  }
  if (vector) {
    Sids <- c(as.matrix(Sids))
  }
  return(Sids)
}

prepareData <- function(GRMIbs = NULL, GRMIbd = NULL, IBDe = NULL, OffspringID = NULL, DpcID = NULL, SNP_group = NULL, Data_group = NULL, Test = NULL) {
  
    print("GRM_IBS")
      x <- GRMIbs[OffspringID, DpcID]
      ind <- which(lower.tri(x , diag = FALSE) , arr.ind = TRUE )
      ibs_nonDiagonal <- data.frame(DpcID = dimnames(x)[[2]][ind[,2]] ,
                                    OffspringID = dimnames(x)[[1]][ind[,1]] ,
                                  Value = x[ind],
                                  Type = "ibs_nonDiagonal",
                                  SNP_group = SNP_group,
                                  Data_group = Data_group,
                                  Test = Test)
  
      print("GRM_IBDr")

      x <- GRMIbs[OffspringID, DpcID]
      ind <- which(lower.tri(x , diag = FALSE) , arr.ind = TRUE )
      ibdr_nonDiagonal <- data.frame(DpcID = dimnames(x)[[2]][ind[,2]] ,
                                     OffspringID = dimnames(x)[[1]][ind[,1]] ,
                                    Value = x[ind],
                                    Type = "ibdr_nonDiagonal",
                                    SNP_group = SNP_group,
                                    Data_group = Data_group,
                                    Test = Test)

    print("IBDe")
      x <- getS(Sinv, ids = DpcID, with = OffspringID, vector = FALSE, diagOnly = FALSE)
      ind <- which(lower.tri(x, diag = FALSE), arr.ind = TRUE)
      ibde_nonDiagonal <- data.frame(DpcID = dimnames(x)[[2]][ind[,2]] ,
                                     OffspringID = dimnames(x)[[1]][ind[,1]] ,
                                     Value = x[ind],
                                     Type = "ibde_nonDiagonal",
                                     SNP_group = SNP_group,
                                     Data_group = Data_group,
                                     Test = Test)
   
 
  return(
      rbind(ibs_nonDiagonal,ibdr_nonDiagonal,ibde_nonDiagonal))
}

data <- prepareData(GRMIbs = GRMIbs, GRMIbd = GRMIbd, IBDe = IBDe, OffspringID = OffspringID, DpcID = DpcID, SNP_group = 4, Data_group = "Nested", Test = "No_GenoErr")



save(data, file = "GRM_data_Nested_SNP4_nGE.Rdata")

#Plot 

paletteViridis <- plasma(n = 5, begin = 0.3, end = 0.9) 

plot_selected_types <- function(data, selected_types) {
  # Filter data for selected types
  filtered_data <- data %>% filter(Type %in% selected_types)
  
  # Create a combined plot with different colors for each selected type
  combined_plot <- ggplot(filtered_data, aes(x = Value, fill = Type)) +
    geom_histogram(bins = 30, color = "white", position = "identity") + #add alpha = 0.6 if you want the colours slighlty transparent 
    scale_fill_manual(values = paletteViridis[1:length(selected_types)]) +
    facet_grid(Type ~ ., scales = "free_y") + #remove TYPE from here if you want to ge rid of the rows 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "top", legend.title = element_blank(), 
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.margin = margin(20, 20, 20, 20)) +
    scale_x_continuous(name = "GRM relatedness value" , breaks=seq(-2, 10, 0.25)) +
    scale_y_continuous(name = "Frequency" , breaks=seq(-2, 10, 0.25))
  
  
  return(combined_plot)
}

plot_selected_types(data=data, selected_types = c("ibs_nonDiagonal", "ibs_diagonal"))



#Putting all of the nest and unNested GRM data together 
#nested
load("~/Desktop/Slovenia data/Simulated /nested/16SNP/data/GRM_data_Nested16.Rdata")
SNP_16 <- data

load("~/Desktop/Slovenia data/Simulated /nested/160SNP/data/GRM_data_Nested160.Rdata")
SNP_160 <- data

load("~/Desktop/Slovenia data/Simulated /nested/800SNP/data/GRM_data_Nested800.Rdata")
SNP_800<- data

load("~/Desktop/Slovenia data/Simulated /nested/1776SNP/data/GRM_data_Nested1776.Rdata")
SNP_1776<- data

load("~/Desktop/Slovenia data/Simulated /nested/4176SNP/data/GRM_data_Nested4176.Rdata")
SNP_4176<- data

NestedGRM_all <- rbind(SNP_16, SNP_160, SNP_800, SNP_1776, SNP_4176)

load("~/Desktop/Slovenia data/Simulated /un-nested/16SNPs/data/GRM_data_unNested16.Rdata")
unSNP_16<- data
load("~/Desktop/Slovenia data/Simulated /un-nested/160SNP/data/GRM_data_unNested160.Rdata")
unSNP_160<- data
load("~/Desktop/Slovenia data/Simulated /un-nested/800SNP/data/GRM_data_unNested800.Rdata")
unSNP_800<- data
load("~/Desktop/Slovenia data/Simulated /un-nested/1776SNP/data/GRM_data_unNested1776.Rdata")
unSNP_1776<- data
load("~/Desktop/Slovenia data/Simulated /un-nested/4176SNP/data/GRM_data_unNested4176.Rdata")
unSNP_4176<- data

UnNestedGRM_all <- rbind(unSNP_16, unSNP_160, unSNP_800, unSNP_1776, unSNP_4176)

load("~/Desktop/Slovenia data/Slov_data/R_code_GRM/GRM_data_Slov1722.Rdata")
Slov_GRM <- data
Slov_GRM <- Slov_GRM[Slov_GRM$Type == c("ibs_nonDiagonal"), ]
Slov_GRM$Data_group <- rep("Slov", length(Slov_GRM$Value))

Summary_GRM <- bind_rows(Slov_GRM, UnNestedGRM_all, NestedGRM_all)

save(Summary_GRM, file = "GRM_summary.Rdata")





calculate_opposing_homozygotes <- function(Geno1, Geno2) {
  
  # Initialize an empty data frame to store results
  results <- data.frame(Geno1_ID = character(), Geno2_ID = character(), pOppHom = numeric(),
                        nOppHom = numeric(), stringsAsFactors = FALSE)
  
  # Loop over all individuals in Geno1
  for (i in 1:nrow(Geno1)) {
    
    cat(paste0("Row: ", i, "/", nrow(Geno1), "\n"))
    # Extract the SNP genotypes of the current individual in Geno1
    Geno1_SNP <- Geno1[i, ]
    
    # Loop over all individuals in Geno2
    for (j in 1:nrow(Geno2)) {
      # Extract the SNP genotypes of the current individual in Geno2
      Geno2_SNP <- Geno2[j,]
      
      # Exclude NA values and heterozygous genotypes (1 represents heterozygous)
      selGeno1 = !(is.na(Geno1_SNP) | Geno1_SNP == 1)
      selGeno2 = !(is.na(Geno2_SNP) | Geno2_SNP == 1)
      
      # Calculate opposing homozygotes proportion
      diff = Geno1_SNP - Geno2_SNP
      diff = diff[selGeno1 & selGeno2]
      pOppHom <- sum(diff !=0) / length(diff)
      
      # Calculate number of opposing homozygotes
      nOppHom <- sum(diff != 0)
      
      # Store results in the data frame
      results <- rbind(results, data.frame(Geno1_ID = rownames(Geno1)[i],
                                           Geno2_ID = rownames(Geno2)[j],
                                           pOppHom = pOppHom,
                                           nOppHom = nOppHom))
    }
  }
  
  # Sort the data frame by Geno1 ID and pOppHom
  results <- results[order(results$Geno1_ID, results$pOppHom), ]
  
  return(results)
}


#Compare workers with dpcs 
OppHomo <- calculate_opposing_homozygotes(x_worker, x_dpc)

ggplot(data = OppHomo, aes(x = Geno2_ID, y = Geno1_ID, fill = pOppHom)) +
  geom_tile() +
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(title = "Proportion of Opposing Homozygotes",
       x = "Geno2 ID", y = "Geno1 ID",
       fill = "Opposing Homozygotes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

write_csv(OppHomo, file = "OppHomo.csv", col_names = TRUE)


