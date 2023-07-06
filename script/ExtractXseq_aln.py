#!/usr/bin/python3
# python ExtractSeq_aln.py pr2_aln.fasta 100
#=> 4 is the number of sequences to include in new file 
#### TO DO BEFORE RUNNING THE SCRIPT ###


__author__ = "Jean-David Grattepanche"
__version__ = "1, May 9, 2023"
__email__ = "jeandavid.grattepanche@gmail.com"


import string
import re
import sys
import os
from Bio import SeqIO
from sys import argv



def ExtractSeq_aln(seqf, numbseq):
	i =0; seqlist=[]
	output = open(seqf.split('.')[0]+"_reduced.fasta","w+")
	for seq in SeqIO.parse(open(seqf,"r"),"fasta"):
		while i < numbseq and seq.description not in seqlist:
			seqlist.append(seq.description)
			print(seq.description)
			output = open(seqf.split('.')[0]+"_reduced.fasta","a")
			output.write('>'+seq.description+'\n'+str(seq.seq)+'\n')
			output.close()
			i+=1
		if i >= numbseq:
			break
def main():
	script, seqfile, seqX = argv
	ExtractSeq_aln(seqfile, int(seqX))
main()

