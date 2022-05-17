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

name = 'distances_all_100kb_w100kb'
txt_name = str(name+'.txt')
plot_name = str(name+'.pdf')
exclude_samples = ['Target_PNT_LG','Target_Chaco']
target_samples = ['Target_FERN_07','Target_FERN_8101']

meta = pd.read_csv(sys.argv[1], delimiter='\t')

avg_dict = {}

with open(sys.argv[2], 'r') as distance_f:
	for line in distance_f:
		parts = line.strip().split()
		avg_dict[parts[0]] = float(parts[1])

cat_distances = []

for key in avg_dict:
	key_p = key.split(',')
	if key_p[0] not in exclude_samples and key_p[1] not in exclude_samples:
		if key_p[0] not in target_samples and key_p[1] not in target_samples:
			pair = key
			#if meta[meta.ID==key_p[0]]['SHELL'].values != meta[meta.ID==key_p[1]]['SHELL'].values:
			#	cat_distances.append(['SHELL', normalized_dist, 'red', 'none'])		
			#elif meta[meta.ID==key_p[0]]['ISLA'].values != meta[meta.ID==key_p[1]]['ISLA'].values:			
			#	cat_distances.append(['ISLAND', normalized_dist, 'green', 'none'])		
			if meta[meta.ID==key_p[0]]['SP'].values != meta[meta.ID==key_p[1]]['SP'].values:			
				cat_distances.append(['Between Species', avg_dict[key], 'blue', 'none'])		
			elif meta[meta.ID==key_p[0]]['SP'].values == meta[meta.ID==key_p[1]]['SP'].values:
				cat_distances.append(['Within Species', avg_dict[key], 'pink', 'none'])
		else:
			pair = key
			#if meta[meta.ID==key_p[0]]['SHELL'].values != meta[meta.ID==key_p[1]]['SHELL'].values:
			#	cat_distances.append(['SHELL', normalized_dist, 'red', 'red'])		
			#elif meta[meta.ID==key_p[0]]['ISLA'].values != meta[meta.ID==key_p[1]]['ISLA'].values:			
			#	cat_distances.append(['ISLAND', normalized_dist, 'green','green'])		
			if meta[meta.ID==key_p[0]]['SP'].values != meta[meta.ID==key_p[1]]['SP'].values:			
				cat_distances.append(['Between Species', avg_dict[key], 'blue', 'blue'])		
			elif meta[meta.ID==key_p[0]]['SP'].values == meta[meta.ID==key_p[1]]['SP'].values:
				cat_distances.append(['Within Species', avg_dict[key], 'pink', 'pink'])

		if key_p[0] == 'Target_FERN_07':
			cat_distances.append(['FERN07', avg_dict[key], 'red', 'none'])
		if key_p[1] == 'Target_FERN_07':
			cat_distances.append(['FERN07', avg_dict[key], 'red', 'none'])

dist_df = pd.DataFrame(cat_distances, columns=['category','normalized genetic distance', 'color', 'fc'])
plt.figure(figsize=(8,5))
sns.catplot(x="category", y="normalized genetic distance", palette="ocean", order=['FERN07','Within Species','Between Species'], kind="strip", data=dist_df)
#sns.displot(dist_df, x="normalized genetic distance", palette="ocean", hue="category", kind="kde")
#ax = sns.violinplot(x="category", y="normalized genetic distance", data=dist_df, order=['Within Species','Between Species'], inner=None, color=".8")
#ax = sns.stripplot(x="category", y="normalized genetic distance", data=dist_df, palette="ocean", order=['Within Species','Between Species'])

plt.tight_layout()
#plt.ylim(0,0.0013)
#plt.ylabel("normalized genetic distance")
#plt.legend()
plt.savefig(plot_name)
plt.close()


