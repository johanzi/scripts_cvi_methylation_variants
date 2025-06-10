#!/usr/bin/python


'''
		 FILE: convertKmers4mapping_methylome_SA.py
		USAGE: ---
  DESCRIPTION: map kmers to references
	  OPTIONS: ---
 REQUIREMENTS: ---
		 BUGS: ---
		NOTES: ---
	   AUTHOR: Emmanuel Tergemina
 ORGANIZATION: Department of Plant Developmental Biology, Max Planck Institute for Plant Breeding Research
	  VERSION: 1.0
	  CREATED: 20220825
	 REVISION: ---
'''

import os
import glob
import sys
import re
import pandas as pd

path2output='/srv/netscratch/irg/grp_hancock/johan/kmerGWAS/output_dir/'
path2ref='/srv/netscratch/irg/grp_hancock/johan/kmerGWAS/reference/'


for folder in glob.glob(path2output + '*'):
	basename = os.path.basename(folder)
	print(basename)
	with open(folder + '/kmers/kmers.fa','w') as fp:
		for line in open(folder + '/kmers/pass_threshold_5per','r'):
			if not 'chr' in line:
				tmp=line.strip().split()
				tmp1=tmp[1].split('_')
				fp.write('>{}\n{}\n'.format(tmp1[1],tmp1[0]))
				# print '>{}\n{}'.format(tmp1[1],tmp1[0])

	##### Map to references
	def map2ref(var):
		os.system('bowtie2 -f -x ' + path2ref + '{}/{}'.format(var,var) +' -U ' + folder + '/kmers/kmers.fa -S ' + folder + '/kmers/' + var +'_aln')
	map2ref('TAIR10')
	map2ref('S1-1')
	map2ref('S15-3')
	map2ref('S5-10')
	map2ref('CVI-0')
	map2ref('Cvi') #Assembly from Korbinian

	##### Add ref positions
	def add_positions(var):
		list_var=[]
		for line in open(folder + '/kmers/' + var + '_aln', 'r'):
			if not '@' in line:
				tmp=line.strip().split()
				list_var.append(tmp[:5])
		return(list_var)
	TAIR10=add_positions('TAIR10')
	S1_1=add_positions('S1-1')
	S15_3=add_positions('S15-3')
	S5_10=add_positions('S5-10')
	CVI_0=add_positions('CVI-0')
	Cvi=add_positions('Cvi')

	n = 0
	with open(folder + '/kmers/' + basename + '.kmers.assoc.txt','w') as fp:
		for line in open(folder + '/kmers/pass_threshold_10per','r'):
			if 'chr' in line:
				fp.write('{}\tChr_TAIR10\tPos_TAIR10\tMP_TAIR10\tChr_S1-1\tPos_S1-1\tMQ_S1-1\tChr_S15_3\tPos_S15_3\tMQ_S15_3\tChr_S5_10\tPos_S5_10\tMQ_S5_10\tChr_Cvi-0\tPos_Cvi-0\tMQ_Cvi-0\tChr_Cvi\tPos_Cvi\tMQ_Cvi\n'.format(line.strip()))
			else:
				TAIR10info=''
				S1_1info=''
				S15_3info=''
				S5_10info=''
				CVI_0info=''
				Cvi_info=''
				tmp=line.strip().split()
				kmer=tmp[1].split('_')[1]
				for i in range(len(TAIR10)):
					if kmer == TAIR10[i][0]:
						TAIR10info=TAIR10[i][2:5]
					if kmer == S1_1[i][0]:
						S1_1info=S1_1[i][2:5]
					if kmer == S15_3[i][0]:
						S15_3info=S15_3[i][2:5]
					if kmer == S5_10[i][0]:
						S5_10info=S5_10[i][2:5]
					if kmer == CVI_0[i][0]:
						CVI_0info=CVI_0[i][2:5]
					if kmer == Cvi[i][0]:
						Cvi_info=Cvi[i][2:5]
						break
				fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(line.strip(),
					'\t'.join(TAIR10info),
					'\t'.join(S1_1info),
					'\t'.join(S15_3info),
					'\t'.join(S5_10info),
					'\t'.join(CVI_0info),
					'\t'.join(Cvi_info)))

		del TAIR10
		del S1_1
		del S15_3
		del S5_10
		del CVI_0
		del Cvi_info






