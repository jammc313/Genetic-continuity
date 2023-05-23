import os

the_Path='DIR_cov_distr'

file_dict={}
for x in os.listdir(the_Path):
	d=x.split('_')
	a_name=d[0]
	if not a_name in file_dict.keys():
		file_dict.update({a_name:[]})
	file_dict[a_name].append(x)


#for x in file_dict.keys():
#	print(x,len(file_dict[x]))


for x in file_dict.keys():
	print(x)
	cov_distr={}
	for a_file in file_dict[x]:
		f=open(the_Path+'/'+a_file,'r')
		for l in f:
			d=l.split()
			cov=int(d[0])
			count=int(d[1])
			if not cov in cov_distr.keys():
				cov_distr.update({cov:0})
			cov_distr[cov]+=count
	outf=open(the_Path+'/'+x,'w')
	for i in range(max(cov_distr.keys())+1):
		outf.write(str(i)+'\t'+str(cov_distr[i])+'\n')
	outf.close()
