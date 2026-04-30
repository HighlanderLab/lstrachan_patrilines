#!/bin/bash
#
#SBATCH --job-name=3_RunBeagle_simData__REP_
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32GB
#SBATCH --time=200:00:00
#SBATCH --output=log/%x.out
#SBATCH --error=log/%x.err


# Initialise the environment modules
. /etc/profile.d/modules.sh
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate AlphaAssign_SIMplyBee


# Run the program
Rscript /home/jo319kis/Projects/Honeybee_patriline_paper/ScriptsApril26/3_RunBeagle_simData.R _REP_ "/home/jo319kis/Projects/Honeybee_patriline_paper" "/home/jo319kis/bin"

Rscript /home/jo319kis/Projects/Honeybee_patriline_paper/ScriptsApril26/4_Converting_PhasedVCF_simData.R _REP_ "/home/jo319kis/Projects/Honeybee_patriline_paper" "/home/jo319kis/bin"

Rscript /home/jo319kis/Projects/Honeybee_patriline_paper/ScriptsApril26/5_Haplotype_ParentAssignments_simData.R _REP_ "/home/jo319kis/Projects/Honeybee_patriline_paper" "/home/jo319kis/bin"

Rscript /home/jo319kis/Projects/Honeybee_patriline_paper/ScriptsApril26/6_Checking_Gametic_MendelianSampling_Values_simData.R _REP_ "/home/jo319kis/Projects/Honeybee_patriline_paper" "/home/jo319kis/bin"

