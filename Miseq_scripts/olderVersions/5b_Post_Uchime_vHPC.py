#!/usr/bin/python2

__author__ = "Jean-David Grattepanche"
__version__ = "6, April 21, 2016"
__email__ = "jeandavid.grattepanche@gmail.com"


import re,os, sys
from Bio import SeqIO
from Bio import Phylo
from sys import argv

# outpath = "/home/tuk61790/v4_SWARM/outputs/"

def main():
	script, seqfile = argv
	fastadict = {}; i =0
	outpath = seqfile.split('chimeras/')[0]
	out= open(outpath + 'chimeras/Seq_reads_nochimera_nosingleton_renamed.fas','w+')
	for Seq in SeqIO.parse(seqfile,'fasta'):
		out.write('>'+Seq.description.replace(';size=','_').replace(';','r')+'\n'+str(Seq.seq)+'\n')
		
	out.close()		
		
main()