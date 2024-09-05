#!/usr/bin/python3
# python3 MiSeq_pipeline_V4_SWARM_part3HPC_noTree.py ExampleFiles/List_samples.txt ExampleFiles/RawData/

#### TO DO BEFORE RUNNING THE SCRIPT ###


__author__ = "Jean-David Grattepanche"
__version__ = "3.02, September 27, 2022"
__email__ = "jeandavid.grattepanche@gmail.com"


import string
import re
import sys
import os
from Bio import SeqIO
seqlist = []
OTUlist = []
duplicatelist = []



def removeoutgroup(outputpath, listsample):
#	os.system('cp '+outputpath+'/chimeras/SWARM_postout_reduced.fas '+outputpath+'outgroup_removal/SWARM_postout_nosingleton_nochimeras_in_only.fasta')
#	os.system('cp '+outputpath+'/OTUs/SWARM_postout.txt '+outputpath+'outgroup_removal/SWARM_postout_nosingleton_nochimeras_in_only.txt')
	os.system('python Miseq_scripts/9_randomly_subsample_ingroup_notree_faster.py '+outputpath+'taxonomic_assignment/Seq_reads_inGroup.fasta '+outputpath+'taxonomic_assignment/Seq_map_inGroup.txt '+ listsample +' 5')
	os.system('python Miseq_scripts/11_makeOTUtable_ingroup_v2.py '+outputpath+'taxonomic_assignment/Seq_map_inGroup.txt '+outputpath+'taxonomic_assignment/subsampled.txt '+ listsample)
	os.system('python Miseq_scripts/12_createFinalfiles_diff_notree.py '+outputpath+'taxonomic_assignment/Seq_reads_inGroup.fasta '+outputpath+'VsearchBLAST.tsv '+outputpath+'OTUs_ingroup/Seq_map_inGroup_subsampled.txt '+listsample)		
	
def main():
	b = sys.argv[1]
	folderraw = sys.argv[2]
	listsamp = []
	outputpath = os.getcwd()+"/"+folderraw.split('/')[0]+ '/outputs/'
	temppath = os.getcwd()+"/"+folderraw.split('/')[0]+'/temp/'
	try:
		listsample = b
	except ValueError:
		b = ""	
	if b == "":
		print ('Your input is empty.  Try again. ')
	for samp in open(listsample,'r'):
		if samp.split('\t')[0] not in listsamp:
			listsamp.append(samp.split('\t')[0])
	if not os.path.exists(outputpath): 
		print('ERROR! No output folder. run part 2!')
		main()
	removeoutgroup(outputpath, listsample)
main()
