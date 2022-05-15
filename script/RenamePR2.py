#!/usr/bin/python3
# python RenamePR2.py pr2_version_4.14.0_SSU_UTAX.fasta 4
#=> 4 is the number of taxonomic ranks to include in addition to the species name and accession number 
#### TO DO BEFORE RUNNING THE SCRIPT ###


__author__ = "Jean-David Grattepanche"
__version__ = "1, April 5, 2022"
__email__ = "jeandavid.grattepanche@gmail.com"


import string
import re
import sys
import os
from Bio import SeqIO
from sys import argv



def renamingPr2(seqf, rankn):
	output = open(('_').join(seqf.split('_')[:-1])+"_renamed.fasta","w+")
	for seq in SeqIO.parse(open(seqf,"r"),"fasta"):
		accession = seq.description.split('.')[0]
		spname = seq.description.split(',s:')[1].replace('.','').replace(':',"*").replace(' ','').replace('(','').replace(')','')
		rank = ('_').join(x[2:5] for x in seq.description.split(',g:')[0].split(',')[1:int(rankn+1)])
		print(rank,"_",spname,"_",accession)
		output.write('>'+rank.replace(' ','').replace('(','').replace(')','')+'_'+spname.replace(' ','').replace('(','').replace(')','')+'_'+accession.replace(' ','').replace('(','').replace(')','')+'\n'+str(seq.seq)+'\n')
	output.close()

def main():
	script, seqfile, ranknum = argv
	renamingPr2(seqfile, int(ranknum))
main()

