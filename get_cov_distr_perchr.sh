# This SLURM submission script implements "get_coverage_distr.py" to output coverage distributions for each chr of given individuals into
# DIR_cov_distr/.
# Once complete, run "merge_cov_distr.py" to merge across chromosomes, and plot results with "plot_cov_distr.R".
# The goal is to check reasonable site cov thresholds to use when filtering sites. Set these individually, rather than with a hard default such as 10<X<500. 

declare -a ind_1=( LBK Loschbour Bichon SIII UstIshim Kotias sf12 Neanderthal flo011 )

for chr in {1..22}; do
    for (( i = 0; i < ${#ind_1[@]}; ++i )); do
#        echo "${chr} ${ind_1[i]} ${ind_2[j]}"
	(echo '#!/bin/bash -l'
echo "python get_coverage_distr.py ${chr} ${ind_1[i]} 
exit 0") | sbatch -p core -n 1 -t 1-00:00:00 -A p2018 -J ${chr}_${ind_1[i]} -e DIR_error/${chr}_${ind_1[i]}.error
    done
done


