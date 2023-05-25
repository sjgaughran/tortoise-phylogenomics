#!/usr/bin/python
from Bio import SeqIO
import os
import sys
import re
import shutil
from amas import AMAS

# This command filters fasta alignments and moves a set that meet certain criteria to a new directory ("filt_directory").
# It requires all fasta alignments to be in the "raw_directory", and to end with the suffix ".fasta"

# Change the "if" statement to adjust criteria for number of parsimonious sites ("pars"), number of missing sites ("miss"),
# and GC content ("gc"). 

raw_directory = "/raw/fasta/dir/"
filt_directory = "/filtered/fasta/dir/"
moved_files = []

i=0

for filename in os.listdir(raw_directory):
	if filename.endswith(".fasta"): 
		full_file = raw_directory+filename
		parts = filename.split('.fasta')
		meta_aln = AMAS.MetaAlignment(in_files=[full_file], data_type="dna",in_format="fasta", cores=1)
		summaries = meta_aln.get_summaries()
		gc = float(summaries[1][0][11])
		print("GC: " + str(gc))
		miss = float(summaries[1][0][5])
		print("Miss: " + str(miss))
		pars = int(summaries[1][0][8])
		print("Pars: " + str(pars))
		var = int(summaries[1][0][6])
		print("Vars: " + str(var))
		i += 1
		# The below loop can be edited change filters. For example, "and pars > 2" could be added to only include loci with more than 
		# 2 parsimony informative sites. 
		if var >1 and miss < 10 and 0.3 < gc < 0.7: 
			print(filename)
			shutil.copy(full_file, filt_directory)
			moved_files.append(filename)

# A list of loci that meet the filter criteria is written to this file.
with open("moved_files.txt", "w") as out_file:
	for item in moved_files:
		out_file.write("%s\n" % item)
