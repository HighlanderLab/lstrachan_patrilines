# Clean workspace
rm(list = ls())

#library
library(AlphaSimR)
library(SIMplyBee)
library(readr)
library(genio)
library(AGHmatrix)
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)
library(future)
library(future.apply)
library(doParallel)

#Modified functions from Audrey's nesting functions
createArray= function(array_name, array_number, nChr, segSites, nSNPPerChr, pop){
  snpArray = vector("list", nChr)
  for (chr in 1:nChr) {
    #get the haplotypes for all sites within the pop
    x = pullSegSiteHaplo(pop, chr = chr)
    #pop allele frequency at all sites
    alleleFreq = apply(X = x, MARGIN = 2, FUN = mean)
    #Will pick SNP based on allele frequency distribution
    # tmp = runif_from_nonunif(x = data.frame(id = 1:segSites[[chr]], value = alleleFreq), n = nSNPPerChr[[chr]]) #run with runif_from_nonunif

    tmp_data <- data.frame(id = 1:length(alleleFreq), value = alleleFreq)
    #Select SNP at random
    tmp <- sample(nrow(tmp_data), size = nSNPPerChr[[chr]], replace = FALSE) # run without runif_from_nonunif
    #save selected SNP
    sel <- tmp_data[tmp, "id"] #without runif_from_nonunif
    #sel = tmp$id #with runif_from_nonunif
    #save SNP in array
    snpArray[[chr]] = sort(sel) # Must be sorted
  }

  snpArray = do.call("c", snpArray) # Collapse list to vector
  snpArray = new(
    Class = "LociMap",
    nLoci = sum(as.integer(nSNPPerChr)),
    lociPerChr = as.integer(nSNPPerChr),
    lociLoc = snpArray,
    name = array_name
  )
  SP$snpChips[[array_number]] = snpArray

}

runif_from_nonunif <- function(x, n, n_bins = 100) {
  samples_min <- min(x$value)
  samples_max <- max(x$value)
  bin_size <- (samples_max - samples_min) / n_bins
  bin_seq <- seq(from = samples_min, to = samples_max, by = bin_size)
  x$bin <- cut(x = x$value, breaks = bin_seq)
  bin_freq <- as.data.frame(table(x$bin))
  colnames(bin_freq) <- c("bin", "freq")
  x <- merge(x = x, y = bin_freq)
  # Sample without replacement and up-weight low frequency values so that
  # once we sample these out, we can then move to more common & high frequency
  # values
  # TODO: maybe this should be done differently/better by sampling bins at
  # random and then randomly within a bin?
  sel <- sample.int(n = nrow(x), size = n, prob = 1 / x$freq, replace = FALSE)
  return(x[sel, c("id", "value")])
}





# Create/download the founder genomes from SIMplyBee
founderGenomes <- simulateHoneyBeeGenomes(nCar = 30,
                                          nChr = 16,
                                          nSegSites = 3200,
                                          Ne = 3000)

#Set up the SP ####
SP <- SimParamBee$new(founderGenomes, csdChr = ifelse(16 >= 3, 3, 1), nCsdAlleles = 128)

SP$nWorkers <- 30
# Track the pedigree
SP$setTrackPed(TRUE)
# Track the recombination
SP$setTrackRec(TRUE)


# Add a SNP chip (Audrey's)
createArray(array_name = 'BigArray', array_number = 1, nChr = 16, segSites = rep(3200, 16), nSNPPerChr = rep(3200, 16), pop = founderGenomes)
#save array
array = colnames(pullSnpGeno(pop = founderGenomes, snpChip = 1, simParam = SP))
write.table(array, "Bigarray.txt", sep = " ", na = "NA", quote = F, row.names = FALSE, col.names = FALSE)

# define csd chromomsome
csdChr <- SP$csdChr
# Add traits - taken from the QuantGen vignettte (find on simplybee.info)
mean <- c(20, 0)
varA <- c(1, 1 / SP$nWorkers)
corA <- matrix(data = c( 1.0, -0.5,
                         -0.5,  1.0), nrow = 2, byrow = TRUE)
SP$addTraitA(nQtlPerChr = 100, mean = mean, var = varA, corA = corA,
             name = c("queenTrait", "workersTrait"))

varE <- c(3, 3 / SP$nWorkers)
corE <- matrix(data = c(1.0, 0.3,
                        0.3, 1.0), nrow = 2, byrow = TRUE)
SP$setVarE(varE = varE, corE = corE)


#Create vector for array size 1728, 800, 160,16
nSNP_array <- rbind(c(3125, 5L), c(108, 4L), c(50, 3L), c(10, 2L), c(1,1L))

#nest arrayssss
for (n in 2:6){ #we create the first array (larger one) outside of the loop and subset it to create the others
  print(n)
  tmp1 = vector()
  for (chr in 1:16){
    tmp2 = sample(colnames(pullSnpGeno(pop = founderGenomes, snpChip = (n-1), chr = chr, simParam = SP)), size = nSNP_array[(n-1),1], replace = FALSE)

    tmp1 = c(tmp1, tmp2)
  }
  SP$addSnpChipByName(tmp1, name=paste0('SNP_', nSNP_array[n-1,2]))
  write.table(tmp1, paste0("Data/SNP_", nSNP_array[(n-1),2], "_array.txt"), sep = " ", na = "NA", quote = F, row.names = FALSE, col.names = FALSE)  
}


##################################################################################
##################################################################################


# Check the relatedness in the real data
realGeno_meta = read.csv("/home/jana/Documents/1Projects/HoneybeeParentage/data/SNP_samples_2022.csv")
realGeno = read.table("/home/jana/Documents/1Projects/HoneybeeParentage/data/newPat_CleanIndsMarkers.raw",
                      header=T)
table(realGeno$IID %in% realGeno_meta$snp_id)
realGeno_meta = realGeno_meta %>% dplyr::filter(snp_id %in% realGeno$IID)
table(realGeno_meta$biotype)

realGeno_IID = realGeno$IID
realGeno_geno = realGeno[, 7:ncol(realGeno)]
realGeno_snpnames = colnames(realGeno)
realGeno_geno[is.na(realGeno_geno)] = -9
rownames(realGeno_geno) = realGeno_IID
realGeno_G = as.data.frame(Gmatrix(SNPmatrix = as.matrix(realGeno_geno),
                                   missing = -9,
                                   ploidy = 2)) %>%
  rownames_to_column("IID") %>%
  pivot_longer(-IID, names_to = "IID2")

realGeno_G %>% filter(IID == IID2) %>%
  ggplot() + geom_histogram(aes(x = value - 1))

realGeno_G %>% filter(IID == IID2) %>% summarise(mean(value - 1))

#Only queens
realGeno_G %>% filter(IID == IID2) %>% filter(IID %in% realGeno_meta$snp_id[realGeno_meta$biotype == "queen"]) %>%
  ggplot() + geom_histogram(aes(x = value - 1))

realGeno_G %>% filter(IID == IID2) %>% filter(IID %in% realGeno_meta$snp_id[realGeno_meta$biotype == "queen"]) %>%
  summarise(mean(value-1))

#Only DPQ
realGeno_G %>% filter(IID == IID2) %>% filter(IID %in% realGeno_meta$snp_id[realGeno_meta$biotype == "dpc"]) %>%
  ggplot() + geom_histogram(aes(x = value - 1))

realGeno_G %>% filter(IID == IID2) %>% filter(IID %in% realGeno_meta$snp_id[realGeno_meta$biotype == "dpc"]) %>%
  summarise(mean(value-1))

realGeno_G$Pop1 = sapply(realGeno_G$IID, FUN = function(x) realGeno_meta$biotype[realGeno_meta$snp_id == x])
realGeno_G$Pop2 = sapply(realGeno_G$IID2, FUN = function(x) realGeno_meta$biotype[realGeno_meta$snp_id == x])
realGeno_G$Data = "Real"
realGeno_G$Gen = "Real"
realGeno_G$Pop = ifelse(realGeno_G$Pop1 == "queen" & realGeno_G$Pop2 == "queen", "queen",
                        ifelse(realGeno_G$Pop1 == "dpc" & realGeno_G$Pop2 == "dpc", "dpc",
                               ifelse(realGeno_G$Pop1 == "worker" & realGeno_G$Pop2 == "worker", "worker", "across")))
realGeno_G$Mating = "Real"



inbreeding = data.frame()
inbreeding = rbind(inbreeding, realGeno_G)
##################################################################################
##################################################################################

# Create X rounds of inbreeding, check relatedness after each - NON BEE VERSION
nQueens = data.frame()
SP$nThreads = detectCores()
plan(sequential)
oopts <- options(future.globals.maxSize = +Inf)  ## 1.0 GB
on.exit(options(oopts))
f <- future({ expr })  ## Launch a future with large objects

for (gen in 1:10) {
  print(paste0("Base generation is ", gen))
  if (gen == 1) {
    #Create base population#realGenoCreate base population
    #basePop = newPop(founderGenomes[1:30])


    basePop = newPop(founderGenomes[1:10])
    print(paste0(nInd(basePop), " base queens."))
    baseAF = calcBeeAlleleFreq(pullSegSiteGeno(basePop), sex = rep("F", basePop@nInd))
  } else {
    # A round of inbreeding
    basePop = randCross(basePop[1:10], nCrosses = basePop[1:10]@nInd)
  }
  # Compute relatedness
  basePop_geno = pullSegSiteGeno(basePop)
  basePop_G= as.data.frame(Gmatrix(SNPmatrix = basePop_geno,
                                   ploidy = 2,
                                   maf = baseAF)) %>%
    rownames_to_column("IID") %>%
    pivot_longer(-IID, names_to = "IID2")

  basePop_G$Pop1 = paste0("SimFounders", gen)
  basePop_G$Pop2 = paste0("SimFounders", gen)
  basePop_G$Data = "Sim"
  basePop_G$Gen = gen
  basePop_G$Pop = "queen"
  basePop_G$Mating = "NonBee_version"
  inbreeding = rbind(inbreeding, basePop_G)
  nQueens = rbind(nQueens, data.frame(Gen = gen, nQueens = basePop@nInd, Mating = "NonBee"))

  if (gen == 10) {
    # Create bee virgin queens from pop
    founderPop = new(Class = "MapPop",
                     nInd = basePop@nInd,
                     nChr = basePop@nChr,
                     ploidy = basePop@ploidy,
                     nLoci = basePop@nLoci,
                     geno = basePop@geno,
                     genMap = founderGenomes@genMap,
                     centromere = founderGenomes@centromere,
                     inbred = founderGenomes@inbred)

    # Create DPQ and inspect their relatedness
    bee_basePop = createVirginQueens(founderPop)
    pFathers <- nFathersPoisson
    nDronesPerQueen = 50
    drones <- createDrones(bee_basePop[1], nInd = nDronesPerQueen)
    fathers <- pullDroneGroupsFromDCA(drones, n = 1, nDrones = pFathers)

    #Create a colony and cross it
    colony1 <- createColony(x = bee_basePop[2])
    colony1 <- cross(colony1, drones = fathers[[1]])



    #create Mating station Drones where the dpcs are
    DPQs <- createVirginQueens(colony1, nInd = 4)

    # Compute relatedness of DPQs
    DPQ_geno = pullSegSiteGeno(DPQs)
    DPQ_G= as.data.frame(Gmatrix(SNPmatrix = DPQ_geno,
                                 ploidy = 2,
                                 maf = baseAF)) %>%
      rownames_to_column("IID") %>%
      pivot_longer(-IID, names_to = "IID2")

    DPQ_G$Pop1 = paste0("SimFounders", gen)
    DPQ_G$Pop2 = paste0("SimFounders", gen)
    DPQ_G$Data = "Sim"
    DPQ_G$Gen = gen
    DPQ_G$Pop <- "dpc"
    DPQ_G$Mating = "NonBee_version"
    inbreeding = rbind(inbreeding, DPQ_G)
  }
}


#### BEE VERSION
founderGenomes = founderGenomes[1:30]
for (gen in 1:10) {
  print(paste0("Base generation is ", gen))
  if (gen == 1) {
    #Create base population#realGenoCreate base population
    #basePop <- createVirginQueens(founderGenomes[1:30])
    basePop <- createVirginQueens(founderGenomes)
    baseAF = calcBeeAlleleFreq(pullSegSiteGeno(basePop), sex = rep("F", basePop@nInd))
  } else {
    # A round of inbreeding
    # This is bee version
    print("Creating drones")
    drones = createDrones(basePop, nInd = 100)
    print("Crossing queens")
    queens = cross(basePop, drones = drones, crossPlan = "create")
    print("Create colonies")
    baseColonies <- createMultiColony(queens)
    print("Merging pops")
    basePop = mergePops(createVirginQueens(baseColonies, nInd = 1))
  }

  # Compute relatedness
  basePop_geno = pullSegSiteGeno(basePop)
  basePop_G= as.data.frame(Gmatrix(SNPmatrix = basePop_geno,
                                   ploidy = 2,
                                   maf = baseAF)) %>%
    rownames_to_column("IID") %>%
    pivot_longer(-IID, names_to = "IID2")

  basePop_G$Pop1 = paste0("queen")
  basePop_G$Pop2 = paste0("queen")
  basePop_G$Data = "Sim"
  basePop_G$Gen = gen
  basePop_G$Pop <- "queen"
  basePop_G$Mating = "Bee_version"
  inbreeding = rbind(inbreeding, basePop_G)
  nQueens = rbind(nQueens, data.frame(Gen = gen, nQueens = basePop@nInd, Mating = "Bee"))

  if (gen == 10) {
    # Create DPQ and inspect their relatedness
    pFathers <- nFathersPoisson
    nDronesPerQueen = 50
    drones <- createDrones(basePop[1], nInd = nDronesPerQueen)
    fathers <- pullDroneGroupsFromDCA(drones, n = 1, nDrones = pFathers)

    #Create a colony and cross it
    colony1 <- createColony(x = basePop[2])
    colony1 <- cross(colony1, drones = fathers[[1]])

    #create Mating station Drones where the dpcs are
    DPQs <- createVirginQueens(colony1, nInd = 10)

    # Compute relatedness of DPQs
    DPQ_geno = pullSegSiteGeno(DPQs)
    DPQ_G= as.data.frame(Gmatrix(SNPmatrix = DPQ_geno,
                                 ploidy = 2,
                                 maf = baseAF)) %>%
      rownames_to_column("IID") %>%
      pivot_longer(-IID, names_to = "IID2")

    DPQ_G$Pop1 = paste0("dpc")
    DPQ_G$Pop2 = paste0("dpc")
    DPQ_G$Data = "Sim"
    DPQ_G$Gen = gen
    DPQ_G$Pop <- "dpc"
    DPQ_G$Mating = "Bee_version"
    inbreeding = rbind(inbreeding, DPQ_G)
  }
}

save.image(file = "Data/2000NE_HBGenome_CsdOn_10gen_PlusArrays.RData")

##################################################################################
# PLOT
##################################################################################
meanRel_queen = inbreeding %>%
  dplyr::mutate(SelfRel = ifelse((IID == IID2), "SelfRel", "Rel")) %>%
  filter(Data == "Sim" | (Data == "Real" & Pop %in% c("queen", "dpc"))) %>%
  group_by(Pop, SelfRel, Data, Gen, Mating) %>%
  summarise(meanRel = mean(value), varRel = var(value))

inbreeding$PopGen <- paste0(inbreeding$Pop, inbreeding$Gen)
meanRel_queen$PopGen <- paste0(meanRel_queen$Pop, meanRel_queen$Gen)
meanRel_queen %>% ggplot() +
  geom_col(aes(x = PopGen, y = meanRel, fill = SelfRel), position = "dodge") +
  facet_grid(. ~ Mating)


plot_step=5
# Plot for self relatedness
# Real data
inbreeding %>%
  dplyr::mutate(SelfRel = (IID == IID2)) %>%
  filter(Data == "Real" & Pop %in% c("queen", "dpc")) %>%
  filter(SelfRel) %>%
  ggplot() + geom_histogram(aes(x = value)) +
  theme_bw(base_size=16) +
  facet_grid(. ~ Pop)


# Simulated
inbreeding %>%
  dplyr::mutate(SelfRel = (IID == IID2)) %>%
  filter(Data == "Sim") %>%
  filter(SelfRel) %>%
  #filter(Gen %in% c(paste0(seq(1, 50, plot_step)))) %>%
  ggplot() + geom_histogram(aes(x = value)) +
  facet_grid(rows = vars(PopGen), cols = vars(Mating)) +
  theme_bw(base_size=16)

# Plot for across relatedness - queens
#Real
#Real
inbreeding %>%
  dplyr::mutate(SelfRel = (IID == IID2)) %>%
  filter((Data == "Real" & Pop %in% c("queen", "dpc"))) %>%
  filter(SelfRel != "SelfRel") %>%
  ggplot() + geom_histogram(aes(x = value)) +
  facet_wrap(. ~ Pop, nrow = 2) +
  theme_bw(base_size=16)

# Simulated
inbreeding %>%
  dplyr::mutate(SelfRel = (IID == IID2)) %>%
  filter(Data == "Sim") %>%
  filter(!SelfRel) %>%
  #filter(Gen %in% seq(1, 50, plot_step)) %>%
  ggplot() + geom_histogram(aes(x = value), bins = 100) +
  facet_grid(rows = vars(PopGen), cols = vars(Mating)) +

  theme_bw(base_size=16)

meanRel_queen %>%
  filter(Data == "Sim") %>%
  filter(Pop == "queen") %>%
  ggplot() +
  geom_line(aes(x = as.numeric(Gen), y = meanRel, colour = Mating)) +
  facet_wrap(. ~ SelfRel, scales="free", nrow=2) +
  xlab("Mean relatedeness") + ylab("Generation")

meanRel_queen %>% filter(Gen %in%  c(5, "Real")) %>%
  dplyr::mutate(PlotPop = paste0(Data, "_", Pop)) %>%
  ggplot(aes(x = PlotPop, y = meanRel))+
  geom_col() +
  facet_grid(. ~SelfRel)


ggplot(data = nQueens) + geom_line(aes(x = as.numeric(Gen), y =nQueens, colour = Mating)) +
  ylim(c(0, max(nQueens$nQueens)))


