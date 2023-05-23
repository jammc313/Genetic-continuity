#This script iterates 'Continuity_Anchor_Test_PART1.py' across all 22 chr and the given anchor individuals. 

#declare -a anchorind=( flo011 Neanderthal Kotias sf12 UstIshim )
declare -a anchorind=( LBK Loschbour Bichon SIII )

for chr in {1..22}; do
    for (( i = 0; i < ${#anchorind[@]}; ++i )); do
        (echo '#!/bin/bash -l'
echo "module load conda"
echo "source conda_init.sh"
echo "conda activate poptrees"
echo "python Continuity_Anchor_Test_PART1.py ${chr} ${anchorind[i]}
exit 0") | sbatch -p core -n 1 -t 24:00:00 -A p2018003 -J ${chr}_${anchorind[i]} -o ${chr}_${anchorind[i]}.output -e DIR_error/${chr}_${anchorind[i]}.error
    done
done
