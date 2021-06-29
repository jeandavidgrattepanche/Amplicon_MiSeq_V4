#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "1, July 8,2020"
__email__ = "jeandavid.grattepanche@gmail.com"


import string
import re
import sys
import os
from sys import argv
from Bio import SeqIO

seqdict = {}
Directory = os.path.abspath(os.getcwd())

def main():
	script, seqfile = argv
	for seq in SeqIO.parse(open(seqfile,'r'),'fasta'):
		newseq=str(seq.seq)
		if "N" in str(seq.seq):
# 			print(seq.description, str(seq.seq))
			while newseq[-1] == 'N':
				print('end:',seq.id,newseq[-1],len(newseq),end='\r')
				newseq = newseq[:-1]
			while newseq[0] == 'N':
				print('beg:',seq.id,newseq[0],len(newseq),end='\r')
				newseq = newseq[1:]
			else:
# 				print('\r')
				if 'N' in newseq:
					out = open(seqfile.split('.fasta')[0] + "_withNv2.fas",'a')
					out.write('>' + seq.description + '\n' + newseq +'\n')
					out.close()
				else:
					out1 = open(seqfile.split('.fasta')[0] + "_noN.fas",'a')
					out1.write('>' + seq.description + '\n' + newseq +'\n')
					out1.close()
		else:
			out1 = open(seqfile.split('.fasta')[0] + "_noN.fas",'a')
			out1.write('>' + seq.description + '\n' + newseq +'\n')
			out1.close()
		
main()
