#!/usr/bin/python3
#python3 MiSeq_pipeline_V4_SWARM_part3HPC.py RWS_0001-0096_List.txt RWS_0001-0096/RawData

#### TO DO BEFORE RUNNING THE SCRIPT ###
# update path L36-37

__author__ = "Jean-David Grattepanche"
__version__ = "1.01, June 29, 2021"
__email__ = "jeandavid.grattepanche@gmail.com"

import string
import re
import sys
import os
from Bio import SeqIO
seqlist = []
OTUlist = []
duplicatelist = []

def removeoutgroup(outputpath, outgrouptree, tree, listsample):
	print("Remove OTUs based on outgroup tree\n")
	if not os.path.exists(outputpath + 'outgroup_removal/'): 
		print('ERROR! No outgroup_removal folder. Run part 2!')	
		main()
 	os.system('python3 Miseq_scripts/8_remove_outgroup_from_tree.py '+outputpath+'outgroup_removal/' +outgrouptree + ' '+outputpath+'outgroup_removal/' +tree + ' '+outputpath+'/OTUs/SWARM_postout.txt '+outputpath+'/chimeras/Seq_reads_nochimera_nosingleton_renamed_nocont.fasta')
	os.system('cp '+outputpath+'/chimeras/Seq_reads_nochimera_nosingleton_renamed_nocont.fasta '+outputpath+'outgroup_removal/SWARM_postout_nosingleton_nochimeras_in_only.fasta')
	os.system('cp '+outputpath+'/OTUs/SWARM_postout.txt '+outputpath+'outgroup_removal/SWARM_postout_nosingleton_nochimeras_in_only.txt')
	os.system('python3 Miseq_scripts/9_randomly_subsample_ingroup_notree.py '+outputpath+'outgroup_removal/SWARM_postout_nosingleton_nochimeras_in_only.txt '+ listsample)
	os.system('python3 Miseq_scripts/11_makeOTUtable_ingroup_v2.py '+outputpath+'outgroup_removal/SWARM_postout_nosingleton_nochimeras_in_only.txt '+outputpath+'outgroup_removal/subsampled.txt '+ listsample)
	os.system('python3 Miseq_scripts/12_createFinalfiles_diff_notree.py '+outputpath+'outgroup_removal/SWARM_postout_nosingleton_nochimeras_in_only.fasta '+outputpath+'VsearchBLAST.tsv '+outputpath+'OTUs_ingroup/SWARM_postout_nosingleton_nochimeras_in_only_subsampled.txt '+listsample)		
	
def main():
	b = sys.argv[1]
	folderraw = sys.argv[2]
	listsamp = []
	outputpath = "/home/tuk61790/" +folderraw.split('/')[0]+ '/outputs/'
	temppath = "/home/tuk61790/" +folderraw.split('/')[0]+ '/temp/'
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
	removeoutgroup(outputpath, outgrouptree, tree, listsample)
main()
