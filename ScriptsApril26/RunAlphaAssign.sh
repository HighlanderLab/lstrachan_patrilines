#Run AlphaAssign.sh

snp_array=$1

source ~/Desktop/Slovenia_data/April26/Pedigree_assignment_softwares/AlphaAssign/AlphaAssign/bin/activate
AlphaAssign -genotypes ${snp_array}.txt -potentialsires PotentialFathers.list -pedigree Pedigree.txt -out ${snp_array} -runtype likelihood