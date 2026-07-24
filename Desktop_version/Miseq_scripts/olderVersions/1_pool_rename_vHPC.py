#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "2.01, June 29, 2021"
__email__ = "jeandavid.grattepanche@gmail.com"


import string
import re
import sys
import os
from sys import argv
from Bio import SeqIO
seqlist = []
OTUlist = []
duplicatelist = []
sampledict = {}

def main():
	script, seqfile,  Listsample = argv
	for samplecode in open(Listsample,'r'):
# 		if samplecode.startswith("RWS"): #RWS should match the beginning of your sample or update the line
		sampledict[samplecode.split('\t')[0]] = samplecode.split('\t')[1].split('\n')[0]
#			print(samplecode.split('\t')[0], samplecode.split('\t')[1].split('\n')[0])
	for seq in SeqIO.parse(open(seqfile,'r'),'fasta'):
		sampleRSW = seqfile.split('/')[-1].split('_')[0]
#		print(seqfile, sampleRSW)
		samplename = sampledict[sampleRSW]			
		out = open(('/').join(seqfile.split('/')[:-1]) + "/readpooled.fas",'a')
		out.write('>' + samplename + '_' + seq.id + '\n' + str(seq.seq) +'\n')
		out.close()
main()
