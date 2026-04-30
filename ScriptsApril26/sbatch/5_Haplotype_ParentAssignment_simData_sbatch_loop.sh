nRep=10

for rep in $(seq 1 ${nRep}) 
do
    sed "s/_REP_/${rep}/g" 5_Haplotype_ParentAssignment_simData_sbatch.sh > sbatch_reps/5_Haplotype_ParentAssignment_simData_sbatch_${rep}.sh
    sbatch sbatch_reps/5_Haplotype_ParentAssignment_simData_sbatch_${rep}.sh
done