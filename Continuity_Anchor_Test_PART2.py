# High Coverage Anchor Continuity Test PART2. Testing for a deviation from expected freq of derived sites in a more recent sample, from those in a more ancient sample, conditioning upon heterozygosity of site in the ancient sample. Idea being that drift alone will not affect the expected frequency of derived sites in this case.
# This version takes an anchor individual and a single vcf (more recent individual) as argument to command line.
# The iteration across chromosomes is implemented within an external bash script. 
# A position is only called heterozygous if it has genotype 0/1 or 1/0, if one of the REF or ALT alleles is present with full support as ancestral allele in outgroups.
# This version is the same as 'Continuity_Anchor_Test_PART2.py' except it differentiates between how the transition and transversion het sites in the anchor are found as ancestral, derived and other. 
# Arguments to this python script are supplied from arrays in 'ScriptLoop_PART2.sh'

import sys
import gzip
import get_file_name
import random

################################ Functions ###################################
def check_if_missingness(full_info):
    a_list=full_info[8:]
    [info,x]=a_list
    d=x.split(':')
    if d[0]=='./.':
        return 1
    else:
        return 0

def get_haploid_genotype(full_info):
    ref_nt=full_info[3]
    alt_nt=full_info[4]
    a_list=full_info[8:]
    b_geno=''
    coverage=0
    [info,x]=a_list
    d=x.split(':')
    if len(d)>5: # a snpAD vcf has "GT:DP:A:C:G:T:PP:GQ" in sample field whether homozygote or heterozygote
        coverage=int(d[1])
        if d[0]=='0/0': # site is homozygote reference
            return [coverage,ref_nt]
        elif d[0]=='1/1':
            return [coverage,alt_nt]
        else: # you have a heterozygote
            nt_pos={'A':2,'C':3,'G':4,'T':5} # dict to index nucleotide coverages in list
            ref_cov = sum(int(x) for x in d[nt_pos[ref_nt]] if x.isdigit()) # sum reads supporting reference allele
            if len(alt_nt)==1: # you have only a single alt allele
                alt_cov = sum(int(x) for x in d[nt_pos[alt_nt]] if x.isdigit())
                return [coverage,random.choice([ref_nt,alt_nt])]
            else:	#You have >1 alt allele; multiallelic sites should have been filtered before
                print('need to be worried?',info,x,full_info[1])
                input()

    else:  # it is a GATK called vcf with "GT:DP" (homozygote or missing), or "GT:AD:DP:GQ:PL" (heterozygote) in sample field
        if not info=='GT:AD:DP:GQ:PL':  #This depends on GATK pipeline used
            if info=='GT:DP' and d[0]=='0/0':
                return [int(d[1]),ref_nt]
            elif info=='GT':
                return [0,'']
            else:
                print('need to be worried?',info,x,full_info[1])
                input()
        else:   # You have >=1 alt allele
            temp=d[1].split(',') # split the AD field to get covs of ref and alt nucleotides
            if len(temp)==2:
                [ref_cov,alt_cov]=temp
                ref_cov=int(ref_cov)
                alt_cov=int(alt_cov)
                coverage=int(d[2])
                if ref_cov==0:
                    return [coverage,alt_nt]
                elif alt_cov==0:
                    return [coverage,ref_nt]
                else:
                    return [coverage,random.choice([ref_nt,alt_nt])]

            else:	#You have >1 alt allele; multiallelic sites should have been filtered before
                print('need to be worried?',info,x,full_info[1])
                input()


def orient_and_get_count_haploid(haploid_gt,anc_nt,der_nt):
    #print('genotypes:',haploid_genotype,anc_nt,der_nt)
    if haploid_gt==anc_nt:
        if anc_nt+der_nt==('AG' or 'GA' or 'CT' or 'TC'):
            return [1,0,0,0,0,0] # Transition
        else:
            return [0,1,0,0,0,0] # Transversion
    elif haploid_gt==der_nt:
        if anc_nt+der_nt==('AG' or 'GA' or 'CT' or 'TC'):
            return [0,0,1,0,0,0] # Transition
        else:
            return [0,0,0,1,0,0] # Transversion
    elif ((haploid_gt!=anc_nt) and (haploid_gt!=der_nt) and (haploid_gt in nt_set)):
        if anc_nt+der_nt==('AG' or 'GA' or 'CT' or 'TC'):
            return [0,0,0,0,1,0] # Transition
        else:
            return [0,0,0,0,0,1] # Transversion
    else:
        return [0,0,0,0]


def check_if_biallelic(anc_nt,der_nt,ref_nt,alt_nt):
    """
    This function can be used to ensure a site is biallelic with either alt or ref matching the ancestral and/or derived allele
    """
    set1=set([anc_nt,der_nt,ref_nt,alt_nt]).difference('.')
    if set1.issubset(nt_set):
        if len(set1)==2: # this condition means we only consider biallelic sites, where ref or alt allele is same as ancestral or derived
            return 1
        else:
            return 0


def make_out_str(a_list):
    b_str=''
    for x in a_list:
        b_str+=str(x)+'\t'
    return b_str[:-1]

#############################################################################
#############################################################################
#############################################################################

# Comparing recent individuals' vcfs against the HETPOS file created above

hetpath='/domus/h1/jamesm/Tiina_collab/anchor_method/DIR_anchorind_hetcounts_perchr'
vcf_path='/crex/proj/snic2020-2-10/uppstore2018169/private/MERGE_FILES/ALLSITE_FILES/vcf_files/PER_CHR_VCFS'

#Throw an error unless you have at least two arguments to command line
if len(sys.argv)<4:
    sys.exit('Input Error: The Anchor continuity test PART2 requires a command line input of format: python script.py the_chr anchor_ind_vcf recent_ind_vcf')

arg_list=sys.argv
the_chr=arg_list[1]
anchorind=arg_list[2]
recentind=arg_list[3]
#print(the_chr)
#print(recentind)
#input()
###########################################################################################
###########################################################################################

#Get vcf file name
file_dict = get_file_name.get_name_file_dict()
recentind_vcf = vcf_path+'/'+file_dict[recentind]
vcf_in = recentind_vcf.split('.vc')
recent_vcf = vcf_in[0] + the_chr + '.vc' + vcf_in[1]

NUCL=['A','C','G','T']
nt_set=set(NUCL)

OKCHR=[]
for i in range(22):
    OKCHR.append(str(i+1))
OK_CHR=set(OKCHR)

#########################################
#########################################
#########################################

count_dict={the_chr: {}} # Initialize dict

#Create some win_steps to see effect of continuity across and within chromosomes, and also to enable wbj later
win_start=0
win_step=5000000    #(corresponds to about 5 cM)
win_end=win_start+win_step

print(the_chr,win_start,win_end)    #prints to slurm outfile
count_dict[the_chr].update({(win_start,win_end): [0,0,0,0,0,0]}) # In each window total number matching ancestral, derived, and other

# Begin reading HET and VCF file, ensure lines in-sink
with open(hetpath+'/chr'+the_chr+'_'+anchorind+'_HetPos.txt','r') as het_file:
    with gzip.open(recent_vcf,'rt',encoding='utf-8') as RecentInd:
        l='##'
        while l[0]=='#':
            l=RecentInd.readline()
        het_l=het_file.readline()
        while l and het_l:
            vcf_data=l.strip().split()
            if not vcf_data[0] in OK_CHR:
                break
            het_d=het_l.strip().split()
            vcf_pos=vcf_data[1]
            het_pos=het_d[0]

            while not vcf_pos==het_pos:  # hashed loop to sync two files, each of which may have missing positions
                if int(vcf_pos) == min(int(vcf_pos), int(het_pos)):
                    l=RecentInd.readline()
                elif int(het_pos) == min(int(vcf_pos), int(het_pos)):
                    het_l = het_file.readline()
                if l and het_l:
                    vcf_data = l.strip().split()
                    vcf_pos = vcf_data[1]
                    het_d = het_l.strip().split()
                    het_pos = het_d[0]
                else:
                    break
            if not het_pos == vcf_pos:
                break

            het_anc=het_d[1]
            het_der=het_d[2]
            vcf_ref=vcf_data[3]
            vcf_alt=vcf_data[4]
            at_pos=int(float(het_pos))

            while at_pos>win_end: # Update dict with new window as key with empty value
                win_start+=win_step
                win_end+=win_step
                print(the_chr,win_start,win_end)
                count_dict[the_chr].update({(win_start,win_end): [0,0,0,0,0,0]})
            qual=vcf_data[5]
            if not qual=='.':
#                if int(float(qual)) >= 30: # removed this threshold because biasing against low coverage sample heterozygote sites
                flag=vcf_data[6]
                flag_list = ['FAIL','FAIL1','FAIL2','FAIL3']
                if flag not in flag_list:
                    if not check_if_missingness(vcf_data):
                        if check_if_biallelic(het_anc,het_der,vcf_ref,vcf_alt):
                            coverage,haploid_genotype=get_haploid_genotype(vcf_data)    #all correct, extract haploid allele call
                            match_count=orient_and_get_count_haploid(haploid_genotype,het_anc,het_der)  #Update dict based on which of ancestral, derived, other
                            count_dict[the_chr][(win_start,win_end)]=[count_dict[the_chr][(win_start,win_end)][i] + match_count[i] for i in range(6)]

            l=RecentInd.readline()
            het_l=het_file.readline()

with open('DIR_counts_per_5cm/chr'+the_chr+'_'+anchorind+'_vs_'+recentind+'.txt','w') as outf:
    for a_tuple in sorted(count_dict[the_chr].keys()):
        out_str=str(a_tuple[0])+','+str(a_tuple[1])+'\t'+make_out_str(count_dict[the_chr][a_tuple])
        outf.write(out_str+'\n')



