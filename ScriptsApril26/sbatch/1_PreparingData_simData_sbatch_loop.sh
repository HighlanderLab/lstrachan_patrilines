nRep=10

for rep in $(seq 1 ${nRep}) 
do
    sed "s/_REP_/${rep}/g" 1_PreparingData_simData_sbatch.sh > sbatch_reps/1_PreparingData_simData_sbatch_${rep}.sh
    sbatch sbatch_reps/1_PreparingData_simData_sbatch_${rep}.sh
done