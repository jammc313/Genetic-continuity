# This script iterates 'Continuity_Anchor_Test_PART1.py' across all 22 chr and the given anchor individuals. It searches through their vcfs, identifying reliable heterozygote ancestral/derived sites.
# The jobs are run for each individual and each chromosome in parallel. 

declare -a anchorind=( LBK Loschbour Bichon SIII sf12 UstIshim Neanderthal )

for chr in {1..22}; do
    for (( i = 0; i < ${#anchorind[@]}; ++i )); do
        (echo '#!/bin/bash -l'
echo "python Continuity_Anchor_Test_PART1.py ${chr} ${anchorind[i]}
exit 0") | sbatch -p core -n 1 -t 24:00:00 -A p2018 -J ${chr}_${anchorind[i]} -o ${chr}_${anchorind[i]}.output -e DIR_error/${chr}_${anchorind[i]}.error
    done
done
