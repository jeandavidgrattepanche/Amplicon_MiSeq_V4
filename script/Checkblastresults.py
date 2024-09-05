#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "1, August 7, 2023"
__email__ = "jeandavid.grattepanche@gmail.com"


import re,os, sys
from Bio import SeqIO
from sys import argv


def checkseq(seqfile, BLASTfile):
	out = open("seq_toBLASTonline.fasta",'w+'); OTUlist = []; i = 0; seqdict = {}
	for OTU in open(BLASTfile, 'r'):
		name = OTU.split('\t')[0]
		sim = OTU.split('\t')[2]
		if float(sim) < 90:
			OTUlist.append((name,int(name.split('_')[1].replace('r',''))))
		else:
			print(name, "\t", sim , " BLAST is not OKAY", end="\r")
			i +=1
	print("\n\n", int(i), " OTUs ready")
	print(len(OTUlist)," OTUs to reBlast")
	OTUlist2 = sorted(OTUlist, key= lambda x :x[1], reverse=True)
	for seq in SeqIO.parse(seqfile,'fasta'):
		seqdict[seq.description]= str(seq.seq)
		print(seq.id, seq.description, end = '\r')
		
	for OTU2 in OTUlist2:
		out = open("seq_toBLASTonline.fasta",'a')
		out.write(">" + OTU2[0] + '\n' + seqdict[OTU2[0]]+'\n')
		out.close
		
		

def main():
	script, seqfile, BLASTfile = argv	
	checkseq(seqfile, BLASTfile)	
main()