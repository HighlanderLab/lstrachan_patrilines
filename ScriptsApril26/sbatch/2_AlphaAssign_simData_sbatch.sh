#!/bin/bash
#
#SBATCH --job-name=2_AlphaAssign_simData__REP_
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --output=log/%x.out
#SBATCH --error=log/%x.err


# Initialise the environment modules
. /etc/profile.d/modules.sh
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate AlphaAssign3.10


# Run the program
Rscript /home/jo319kis/Projects/Honeybee_patriline_paper/ScriptsApril26/2_AlphaAssign_simData.R _REP_ "/home/jo319kis/Projects/Honeybee_patriline_paper" "/home/jo319kis/bin"
