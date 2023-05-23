#This script iterates 'Continuity_Anchor_Test_PART2.py' across all 22 chr, for each of the vcfs in array.
#The heterozygite sites in each anchor individual are found in each "recent individual", and the numbers of ancestral and derived counted at those conditioned sites. 

declare -a anchorind=( flo011 Neanderthal Kotias UstIhim sf12 SIII )
declare -a recentind=( buk002 buk002dr buk003 buk003dr buk004 buk004dr buk010 buk010dr buk012 buk013 buk018 buk019 buk019dr buk022 buk022dr buk023 buk023dr buk029 buk029dr buk031 buk031dr buk033 buk033dr buk040 lbk101 lbk102 lbk104 lbk138 poz120 poz120dr poz121 poz177 poz236 poz236dr poz252 poz264 poz275 poz297 poz297dr poz375 poz503 poz503dr rom002_rom003_rom057_rom058 rom011 rom046dr )

for chr in {1..22}; do
    for (( i = 0; i < ${#anchorind[@]}; ++i )); do
        for (( j = 0; j < ${#recentind[@]}; ++j )); do
            (echo '#!/bin/bash -l'
echo "python Continuity_Anchor_Test_PART2.py ${chr} ${anchorind[i]} ${recentind[j]}
exit 0") | sbatch -p core -t 03:00:00 -A p2018 -J ${chr}_${anchorind[i]}_${recentind[j]} -o ${chr}_${anchorind[i]}_${recentind[j]}.output -e DIR_error/${chr}_${anchorind[i]}_${recentind[j]}.error
        done
    done
done

