#Run AlphaAssign.sh

snp_array=$1
in_dir=$2
out_dir=$3


source ~/bin/anaconda3/etc/profile.d/conda.sh
conda activate AlphaAssign
AlphaAssign -genotypes ${in_dir}/${snp_array}.txt -potentialsires ${in_dir}/SimData_PotentialFathers.list -pedigree ${in_dir}/SimData_Pedigree.txt -out ${out_dir}/${snp_array} -runtype likelihood