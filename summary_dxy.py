#!/usr/bin/python
from Bio import SeqIO
import os
import sys
import re
import shutil
import seaborn as sns 
from matplotlib import pyplot as plt
import itertools 
import pandas as pd
import numpy as np
from amas import AMAS

fasta_directory = "/ysm-gpfs/special/cgab/tortoise_fasta/science_fastas/all_torts/all_100kb_w100kb_filtered_fastas/"
name = 'distances_all_100kb_w100kb'
txt_name = str(name+'.txt')
plot_name = str(name+'.pdf')
exclude_samples = ['Target_PNT_LG','Target_Chaco']

meta = pd.read_csv(sys.argv[1], delimiter='\t')

def compseq(seq1, seq2):
	n = 0
	count = 0
	for a, b in zip(seq1, seq2):
		if a != 'N' and b != 'N':
			n += 1
			if a != b:
				count += 1
	if n > 100:
		pair_dif = count/n
		return pair_dif

def Reverse(tuples): 
    new_tup = () 
    for k in reversed(tuples): 
        new_tup = new_tup + (k,) 
    return new_tup

dist_dict = {}
avg_dict = {}
f = 0

for filename in os.listdir(fasta_directory):
	print(f)
	if filename.endswith(".fasta"): 
		full_file = fasta_directory+filename
		parts = filename.split('.fasta')
		meta_aln = AMAS.MetaAlignment(in_files=[full_file], data_type="dna",in_format="fasta", cores=1)
		aln_dicts = meta_aln.get_parsed_alignments()
		combs = list(itertools.combinations(aln_dicts[0], 2))
		f += 1
		for item in combs:
			#rev_item = Reverse(item)
			seq1 = aln_dicts[0][item[0]]
			seq2 = aln_dicts[0][item[1]]
			pdist = compseq(seq1, seq2)
			if item in dist_dict:
				dist_dict[item].append(pdist)
				#dist_dict[rev_item].append(pdist)
			else:	
				dist_dict[item] = [pdist]
				#dist_dict[rev_item] = [pdist]
		
for key in dist_dict:
	f = (sum(x is not None for x in dist_dict[key]))
	avg_dict[key] = sum(filter(None, dist_dict[key]))/f

cat_distances = []

with open(txt_name, 'w') as f_out:
	for key in avg_dict:
		if key[0] not in exclude_samples and key[1] not in exclude_samples:
			pair = ','.join(key)
			avg_het = (float(meta[meta.ID==key[0]]['HET'].to_string(index=False))+float(meta[meta.ID==key[1]]['HET'].to_string(index=False)))/2
			normalized_dist = avg_dict[key]/avg_het
			f_out.write(pair + '\t' + str(normalized_dist) + '\n')

			if meta[meta.ID==key[0]]['SHELL'].values != meta[meta.ID==key[1]]['SHELL'].values:
				cat_distances.append(['SHELL', normalized_dist])
			elif meta[meta.ID==key[0]]['ISLA'].values != meta[meta.ID==key[1]]['ISLA'].values:
				cat_distances.append(['ISLAND', normalized_dist])
			elif meta[meta.ID==key[0]]['SP'].values != meta[meta.ID==key[1]]['SP'].values:
				cat_distances.append(['BETWEEN_SPECIES', normalized_dist])
			elif meta[meta.ID==key[0]]['SP'].values == meta[meta.ID==key[1]]['SP'].values:
				cat_distances.append(['WITHIN_SPECIES', normalized_dist])

dist_df = pd.DataFrame(cat_distances, columns=['category','normalized genetic distance'])
plt.figure(figsize=(8,5))
print(dist_df)
sns.catplot(x="category", y="normalized genetic distance", kind="swarm", data=dist_df)
#plt.ylim(0,0.0013)
#plt.ylabel("normalized genetic distance")
#plt.legend()
plt.savefig(plot_name)
plt.close()


