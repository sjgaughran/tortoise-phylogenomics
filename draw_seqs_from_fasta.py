#!/usr/bin/python

# This script allows the user to draw sequences from set of genomes in fasta format. The genomes
# should have identical coordinates produced by mapping to the same reference assembly. 
# For Jensen, Gaughran et al. 2022 we used ANGSD -doFasta 2 to create each consensus fasta.

import sys
import os
import re
from random import randint
import Bio.SeqIO as Seq

# Specify the directory where the genomic fastas are kept:
directory = "/fasta/dir/"

#Specify the directory where the individual aligned loci will be written:
gene_dir = "/individual/loci/dir/"

# fasta_list_name = sys.argv[1]
# fasta_list_f = open (fasta_list_name, "r")
# fasta_list = fasta_list_f.readlines()

seq_size = 1000 # Size of your sequence. 500bp-1kb is good for BPP, while larger sizes are appropriate for other analyses.
step_size = 100000 # Size of interval between sequences. Make sure it's enough to allow independence/unlinked loci.
start_point = randint(10000,200000) # Pick a random place some distance from the begining of the contig.


# Thie script iterates through each genomic fasta in a directory. The fasta files should end in ".fa" (can be changed below).
# The fifth and sixth lines of this loop draw out the sample name. In our case, our files were named things like "AGO_08_masked.fa"
# so we pull out target_parts[0] (i.e. "AGO") and target_parts[1] (i.e. "08") to create the sample name. 
# Output files are written as an "alignment" to a new file in gene_dir. Note that there is not an actual alignment step, and that this script instead
# relies on identical coordinates across genomes. 
for filename in os.listdir(directory):
	if filename.endswith('.fa'):
		parts = filename.split('.')
		target_sample = parts[0]
		target_parts = target_sample.split('_')
		sample = target_parts[0]+target_parts[1]
		fasta_seq = Seq.parse(directory+filename, "fasta")
		for contig in fasta_seq:
			for i in range(start_point,len(contig.seq),step_size):
				start = i
				end = i+seq_size
				limit = len(contig.seq)-start_point
				out_name = gene_dir+contig.id+'_'+str(start)+'_'+str(end)+'.fasta'
		
				if os.path.exists(out_name):
					append_or_write = 'a' # append if already exists
				else:
					append_or_write = 'w' # make a new file if not
		
				if end < limit:
					target_seq = str(contig.seq[start:end])
					with open(out_name, append_or_write) as out_f:
						out_f.write('>'+sample+'\n'+target_seq+'\n')
				else:
					pass
