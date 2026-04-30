nRep=10

for rep in $(seq 1 ${nRep}) 
do
    sed "s/_REP_/${rep}/g" 3_RunBeagle_simData_sbatch.sh > sbatch_reps/3_RunBeagle_simData_sbatch_${rep}.sh
    sbatch sbatch_reps/3_RunBeagle_simData_sbatch_${rep}.sh
done