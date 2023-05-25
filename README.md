# Galapagos Tortoise Phylogenomics

These scripts detail the computational methods used in Jensen, Gaughran et al. 2022.

## Drawing random sequences and filtering loci

We generated consensus fasta files using `angsd -minQ 35 -minMapQ 20 -uniqueOnly 1 -setMaxDepth 150 -C 50 -doFasta 2 -doCounts 1`. With those genomic fastas in a directory, `draw_seqs_from_fasta.py` can then be run by specifying the input and output directories. The size of the drawn sequence and the interval between sequences can be specified in the script. 

The drawn loci can then be filtered using `loci_filtering.py`. Filtering options and thresholds are specified within the script. Loci meeting the criteria are written into a new directory.

## dxy calculation and plotting
`summary_dxy.py` outputs the Dxy estimate between pairs. `plot_dxy.py` plots those estimates according to species, island, and morphology. `tortoise_meta_data.txt` provides metadata for each tortoise specimen. 

Note: scripts will need to be edited for use in other systems, according to the metadata columns used. 
