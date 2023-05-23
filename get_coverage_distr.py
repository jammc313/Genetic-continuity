import sys
import gzip
import get_file_name
from zipfile import ZipFile

def get_var_form(a_list):
    for x in a_list:
        d=x.split('=')
        if d[0]=='VariantType':
            return d[1]
    return ''

def check_if_missingness(a_list):
    for x in a_list:
        d=x.split(':')
        if d[0]=='./.':
            return 'missingness'
    return ''

def parse_var_genotypes(a_list):
    b_list=[]
    for x in a_list:
        d=x.split(':')
        if d[0]=='./.':
            return 'missingness'
        b_list.append(d[0])
    return b_list


def get_genotype(a_list):
    b_geno=''
    coverage=0
    for x in a_list:
        d=x.split(':')
        if d[0]=='./.':
            return [0,'']
        b_geno=d[0]
        if len(d)>5: # a snpAD vcf has "GT:DP:A:C:G:T:PP:GQ" in sample field
            coverage+=int(d[1])
        else: # it is a GATK called vcf with "GT:DP" (homozygote or missing), or "GT:AD:DP:GQ:PL" (heterozygote) in sample field 
            if len(d)>2:
                if d[2]=='.':
                    return [0,'']
                coverage+=int(d[2])
            else:
                coverage+=int(d[1])
    return [coverage,b_geno]



def check_if_ok_and_get_var_form(anc_nt,ref_nt1,alt_nt1):
    set1=set([anc_nt,ref_nt1,alt_nt1]).difference('.')
    if set1.issubset(nt_set):
        if len(set1)==1:
            return 'OK_NO_VARIATION'
        if len(set1)==2:
            if alt_nt1=='.':
                return 'OK_NO_VARIATION'
            else:
                return 'OK_POLY'
    return ''


def make_out_str(a_list):
    b_str=''
    for x in a_list:
        b_str+=str(x)+','
    return b_str[:-1]
	
ancPath='/proj/snic2020/private/Analyses/per_chr_vcfs'
vcf_path='/proj/snic2020/ALLSITE_FILES/vcf_files/PER_CHR_VCFS'

arg_list=sys.argv
the_chr=arg_list[1]
ind=arg_list[2]

##########################
##########################
##########################
##########################

#Get file name function to use ind1 ind2 etc 
file_dict=get_file_name.get_name_file_dict()
vcf_file_one=vcf_path+'/'+file_dict[ind]
vcf_fileOne=vcf_file_one.split('.vc')
vcf_file1=vcf_fileOne[0] + the_chr + '.vc' + vcf_fileOne[1]

#sys.exit(outgrp_vcf)
##########################
##########################
##########################
##########################


NUCL=['A','C','G','T']
nt_set=set(NUCL)
ANCESTRAL_FILTER=['A','C','G','T']


win_start=0
win_step=5000000 #(corresponds to about 5 cM)
win_end=win_start+win_step

					
cov_distr={}

with ZipFile(ancPath+'/Ancestral_states.zip', 'r') as z:
    with z.open('Ancestral_states/chr'+the_chr+'.txt','r') as anc_file:
        with gzip.open(vcf_file1,'rt',encoding='utf-8') as myf1:
            if 1==1:
                l1='##'   
                while l1[0]=='#':
                    l1=myf1.readline()   	
                anc_l=anc_file.readline().decode('utf-8')
                while l1 and anc_l:
                    vcf_data=l1.strip().split()
                    anc_d=anc_l.strip().split()
                    vcf_pos=vcf_data[1]
                    anc_pos=anc_d[0]

                    while not vcf_pos == anc_pos:
                        if int(vcf_pos) == min(int(vcf_pos), int(anc_pos)):
                            l1 = myf1.readline()
                        elif int(anc_pos) == min(int(vcf_pos), int(anc_pos)):
                            anc_l = anc_file.readline().decode('utf-8')
                        if l1 and anc_l:
                            vcf_data = l1.strip().split()
                            vcf_pos = vcf_data[1]
                            anc_d = anc_l.strip().split()
                            anc_pos = anc_d[0]
                        else:
                            break
                    if not anc_pos == vcf_pos:   # ugly way to catch mismatched positions at end of files and allow to finish without error
                        break

                    at_pos=int(float(anc_pos))
                    while at_pos>win_end:
                        win_start+=win_step
                        win_end+=win_step
                    anc_support=anc_d[2]
                    if anc_support=='3':
                        qual=vcf_data[5]
                        if not qual=='.' and int(float(qual)) >= 30:			
                            flag=vcf_data[6]
	##############################CONSIDER IF YOU WANT OTHER FILTERS#############################
                            if not (flag in ['FAIL','FAIL1','FAIL2','FAIL3']):
	###########################################################
                                anc_nt=anc_d[1]
                                ref_nt=vcf_data[3]
                                alt_nt=vcf_data[4]
                                [coverage,genotype]=get_genotype(vcf_data[9:])
	##############################CONSIDER IF YOU WANT OTHER FILTERS#############################
                                if 1==1:
	###########################################################
                                    if anc_nt in ANCESTRAL_FILTER:
                                        var_form=check_if_ok_and_get_var_form(anc_nt,ref_nt,alt_nt)
                                        if not var_form=='':
                                            if not coverage in cov_distr.keys():
                                                cov_distr.update({coverage:0})
                                            cov_distr[coverage]+=1
                    l1=myf1.readline()   
                    anc_l=anc_file.readline().decode('utf-8')

outf=open('DIR_cov_distr/'+ind+'_chr'+the_chr+'.txt','w')
for i in range(max(cov_distr.keys())+1):
    if i in cov_distr.keys():
        outf.write(str(i)+' '+str(cov_distr[i])+'\n')
    else:
        outf.write(str(i)+' 0\n')
outf.close()
