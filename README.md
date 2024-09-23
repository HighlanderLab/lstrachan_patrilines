# Optimizing Genotype Phasing and Patriline Determination in Honeybees Through Pedigree Reconstruction

## Introduction
This repository contain R scripts that optimize genotype phasing and patriline determination in honeybees through pedigree reconstruction using [SIMplyBee](https://cran.r-project.org/web/packages/SIMplyBee/index.html). 

If you have no experience with SIMplyBee, we suggest you go to the software's website [simplybee.info](http://www.simplybee.info) and study the tutorials. This repository uses SIMplyBee to simulate honeybee populations as described in the below manuscript and to calculate and summarise relatedness between individual honeybees.

## Respository contents
The repository contains five scripts:
* ```SimulatedSNPArrays.R``` runs a SNP array with 4 SNP chip sizes, mimicking a mating station with 2 generations of honeybees. Two SNP versions are creates, those without genotyping errors and those with added genotyping errors.
*  ```RealData_processing.R``` describes the processing of the real datafiles in preparation for pedigree reconstuction and further analysis
* ```HaplotypePO_Patrilines.R``` determines the haplotype parental origins, calculates gametic Mendelian sampling values and determines patriline numbers from a vcf file.
* ```Plotting.R``` processes the output data and creates plots to visualise of all Pedigree reconstruction software and Patriline numbers.
