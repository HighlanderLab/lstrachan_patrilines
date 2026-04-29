###### Phasing prep #########

rm(list = ls())

library(Eagle)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
workingDir = args[1]
softwareDir = args[2]
realDataDir = args[3]
pathToPlink <- softwareDir

#pathToPlink <- "~/Desktop/PLINK/./"
#workingDir = "~/Desktop/lstrachan_patrilines"
setwd(workingDir)
source("ScriptsApril26/AB_to_12.R")


dir.create("Real_data", showWarnings = FALSE)
setwd(paste0(workingDir, "/Real_data/"))

dir.create("Data", showWarnings = FALSE)
dir.create("Outputs", showWarnings = FALSE)
dir.create("Pipeline", showWarnings = FALSE)

Slov_raw_map <- read.table(paste0(realDataDir, "/newPat.map"))
Slov_raw_ped <- read.table(paste0(realDataDir, "/newPat.ped"))


#Fix the pedigree *********************************
Real_SNP_samples <- read.csv(paste0(realDataDir, "/SNP_samples_2022.csv"))
write.csv(Real_SNP_samples, "Data/SNP_samples_2022.csv", row.names = FALSE) #Just to make sure the file is in the right place and can be read in later
queen_ids <- Real_SNP_samples[Real_SNP_samples$biotype == "queen",]
queen_ids <- queen_ids$snp_id
dpc_ids <- Real_SNP_samples[Real_SNP_samples$biotype == "dpc",]
unknown_workers <- Real_SNP_samples[Real_SNP_samples$microlocation == "na",]
unknown_workers_id <- unknown_workers$snp_id

Real_data_pedigree <- read.csv(paste0(realDataDir, "/Real_Data_pedigree.csv"))
write.csv(Real_data_pedigree, "Data/Real_Data_pedigree.csv", row.names = FALSE) #Just to make sure the file is in the right place and can be read in later
Slov_raw_ped$V4 <- 0
# Match rows from pedigree to PED by ID
match_ids <- match(Slov_raw_ped$V2, Real_data_pedigree$ID)
# Fill sire (V3) where a match exists
Slov_raw_ped$V3[!is.na(match_ids)] <- Real_data_pedigree$FID[match_ids[!is.na(match_ids)]]
# Fill dam (V4) where a match exists
Slov_raw_ped$V4[!is.na(match_ids)] <- Real_data_pedigree$MID[match_ids[!is.na(match_ids)]]
#remove the unknown workers 
Slov_raw_ped <- Slov_raw_ped[!Slov_raw_ped$V2 %in% unknown_workers_id, ]

# Format the map file ********************************
Fix_map <- Slov_raw_map[,2]

# Function to extract the desired parts and create a data frame
extract_parts <- function(x) {
  parts <- unlist(strsplit(x, "[._]"))
  chromosome <- paste(parts[2], paste(parts[3], parts[4], sep="."), sep="_")
  position <- as.numeric(parts[5])
  ref_allele <- parts[6]
  data.frame(Chromosome = chromosome, Position = position, RefAllele = ref_allele, stringsAsFactors = FALSE)
}

# Organise the map (we don't have the Genetic distance colummn yet and will add that at the end)
Slov_map_refAll <- do.call(rbind, lapply(Fix_map, extract_parts))
Slov_map_refAll$markerID <- Fix_map
write.table(Slov_map_refAll, file = "Data/Slov_map_refAllele.map", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = " ")

AB_to_12(ped_file = paste0(realDataDir, "/newPat.ped"), map_ref_file = "Data/Slov_map_refAllele.map", output_file_prefix = "Data/newPat_12")

Slov_map <- read.table("Data/newPat_12.map")
colnames(Slov_map) <- c("Chromosome", "markerID", "GeneticDistance", "Position")
Slov_ped <- read.table("Data/newPat_12.ped")


# Next we only want chromosomes beginning with NC (NW are those with unknown chromomsomes)
filter_chromosome_nc <- function(df) {
  # Filter the dataframe to keep only rows where the Chromosome column starts with "NC"
  filtered_df <- df[grep("^NC", df$Chromosome), ]
  return(filtered_df)
}

Slov_map_filtered <- filter_chromosome_nc(Slov_map)
# Change the chromosome number
#Slov_map_filtered$Chromosome <- gsub()
nrow(Slov_map_filtered)


#Slov_ped_filtered <- filter_columns_nc(Slov_raw_ped)
#ncol(Slov_ped_filtered)
write.table(unique(Slov_map_filtered$markerID), "Data/Slov_map_filtered_markers.txt", row.names = FALSE, quote=F, col.names = FALSE)

system(paste0(pathToPlink, "/plink --file Data/newPat_12 --extract Data/Slov_map_filtered_markers.txt --recode --allow-extra-chr --out Data/Slov_fM"))


# ---  Map Chromosomes numerically ---

Slov_map_filtered_chrom <- Slov_map_filtered %>%
  dplyr::mutate(Chromosome = case_when(
    Chromosome == "NC_007070.3" ~ 1,
    Chromosome == "NC_007071.3" ~ 2,
    Chromosome == "NC_007072.3" ~ 3,
    Chromosome == "NC_007073.3" ~ 4,
    Chromosome == "NC_007074.3" ~ 5,
    Chromosome == "NC_007075.3" ~ 6,
    Chromosome == "NC_007076.3" ~ 7,
    Chromosome == "NC_007077.3" ~ 8,
    Chromosome == "NC_007078.3" ~ 9,
    Chromosome == "NC_007079.3" ~ 10,
    Chromosome == "NC_007080.3" ~ 11,
    Chromosome == "NC_007081.3" ~ 12,
    Chromosome == "NC_007082.3" ~ 13,
    Chromosome == "NC_007083.3" ~ 14,
    Chromosome == "NC_007084.3" ~ 15,
    Chromosome == "NC_007085.3" ~ 16,
    TRUE ~ 0
  ))


Slov_map_filtered_chrom <- Slov_map_filtered_chrom[, c("Chromosome", "markerID", "GeneticDistance", "Position")]

# --- Save Final Map File --- (must be the same name as the ped file to work in PLINK)
write.table(Slov_map_filtered_chrom, file = "Data/Slov_fM.map", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")



## ped file is in 0/A/B format which can't be read by Beagle. Needs to be changed to 0/A/C

# Function to convert allele codes
# convert_ped_AB_to_AC <- function(ped_df) {
#   ped_df[, 7:ncol(ped_df)] <- apply(
#     ped_df[, 7:ncol(ped_df)],
#     2,
#     function(col) {
#       col[col == "B"] <- "C"
#       col
#     }
#   )
#   return(ped_df)
# }


# Slov_ped_filtered_AC <- convert_ped_AB_to_AC(Slov_ped_filtered)
# write.table(Slov_ped_filtered_AC, file = "Data/Real_data/Slov_fM_AC.ped", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")


#### ---- PLINK quality control ###################################################################################################
setwd("Data")

# Step 1: Quality control with PLINK
system(paste0(pathToPlink, "/plink --file Slov_fM --make-bed --geno 0.1 --mind 0.1 --maf 0.01 --out Slov_fM_QC"))

system(paste0(pathToPlink, "/plink --bfile Slov_fM_QC --recode --out Slov_fM_QC"))
Slov_fm_QC_ped <- read.table("Slov_fM_QC.ped")
#Identify which queen was removed 
ped_id <- Slov_fm_QC_ped$V2
#which queen_id is not in ped_id?
removed_queen <- setdiff(queen_ids, ped_id)
#Which workers are the offspring of the removed queen? 
real_ped = read.csv(paste0(realDataDir, "/Real_Data_pedigree.csv"), header = TRUE, stringsAsFactors = FALSE)
workers_to_remove <- data.frame(Fam = "AMEL", ID = real_ped$ID[real_ped$MID %in% removed_queen])
write.table(workers_to_remove, "workers_to_remove.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep=" ")
system(paste0(pathToPlink, "/plink --bfile Slov_fM_QC --remove workers_to_remove.txt --recode --out Slov_fM_QC"))   

# Step 2: Convert to VCF
# system(paste0(pathToPlink, "plink --bfile Slov_fM_QC --recode vcf --out Slov_fM_QC"))

# # Step 3: Sort and index the VCF
# system("bcftools sort Slov_fM_AC_QC.vcf -Oz -o Slov_fM_AC_QC_sorted.vcf.gz")
# system("tabix -p vcf Slov_fM_AC_QC_sorted.vcf.gz")

# # Step 3.5: Check for duplicate positions
# system("bcftools query -f '%CHROM\\t%POS\\n' Slov_fM_AC_QC_sorted.vcf.gz | sort | uniq -d")

# # Step 3.6: Normalize and remove duplicate records
# system("bcftools norm -m -any Slov_fM_AC_QC_sorted.vcf.gz -Oz -o Slov_fM_AC_QC_biallelic.vcf.gz")
# system("tabix -p vcf Slov_fM_AC_QC_biallelic.vcf.gz")
# system("bcftools norm -d all Slov_fM_AC_QC_biallelic.vcf.gz -Oz -o Slov_fM_AC_QC_noDupPos.vcf.gz")
# system("tabix -p vcf Slov_fM_AC_QC_noDupPos.vcf.gz")
#TODO: Add the manual transformation from VCF to Ped and Map



#Step 4 Remove the workers who's queen was deleted during QC
# workers_ids <- paste0("AMEL_",workers_to_remove$V2)
# write.table(workers_ids, "workers_to_remove.txt",
#             quote = FALSE, row.names = FALSE, col.names = FALSE)
# system("bcftools view -S ^workers_to_remove.txt Slov_fM_AC_QC_noDupPos.vcf.gz -Oz -o Slov_fM_AC_QC_filtered.vcf.gz")
# system("tabix -p vcf Slov_fM_AC_QC_filtered.vcf.gz")

# #Step 4.1 Check the samples were removed
# before_n <- as.integer(system("bcftools query -l Slov_fM_AC_QC_noDupPos.vcf.gz | wc -l", intern = TRUE))
# after_n  <- as.integer(system("bcftools query -l Slov_fM_AC_QC_filtered.vcf.gz | wc -l", intern = TRUE))
# cat("Removed:", before_n - after_n, "\n")

#FINAL VCF.GZ FILE TO GO FORWARD IS NAMES Slov_fm_AC_QC_filtered.vcf.gz
  



#******* Don't know what this is yet - May make sense during phasing *************
#Haploid genome data not working. Try formatting haplotypes as 0|/ 1| 

# Define the input and output file paths
# vcf_file <- "Genome_SloDrones.vcf"
# output_file <- "Genome_SloDrones_haploid_fixed.vcf"

# # Open the input and output files
# infile <- file(vcf_file, "r")
# outfile <- file(output_file, "w")

# # Process the file line by line
# while(TRUE) {
#   line <- readLines(infile, n = 1)
  
#   # Stop if we reach the end of the file
#   if(length(line) == 0) {
#     break
#   }
  
#   # Write header lines directly to the output file
#   if(grepl("^#", line)) {
#     writeLines(line, outfile)
#   } else {
#     fields <- strsplit(line, "\t")[[1]]
    
#     # Check if the format field is GT (Genotype)
#     if(fields[9] == "GT") {
#       genotypes <- fields[10:length(fields)]
      
#       # Add a separator to each haploid genotype
#       fixed_genotypes <- sapply(genotypes, function(g) {
#         if(g %in% c("0", "1")) {
#           return(paste0(g, "|"))
#         } else {
#           return(g)
#         }
#       })
      
#       # Write the modified line to the output file
#       writeLines(paste(c(fields[1:9], fixed_genotypes), collapse = "\t"), outfile)
#     } else {
#       writeLines(line, outfile)
#     }
#   }
# }

# Close the input and output files
# close(infile)
# close(outfile)

# cat("Fixed VCF file saved to", output_file, "\n")

# phased_Slov_vcf <- read.vcfR("phased_Slov_filtered_Beagle4.vcf")

# phased_Slov_haplo <- extract.gt(phased_Slov_vcf)




save.image(file = paste0(workingDir, "/Real_data/Pipeline/1_RealData_prepared.Rdata"))
