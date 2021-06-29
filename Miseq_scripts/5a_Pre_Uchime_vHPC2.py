#!/usr/bin/python2

__author__ = "Jean-David Grattepanche"
__version__ = "7.01, June 29, 2021"
__email__ = "jeandavid.grattepanche@gmail.com"


import re, os, sys
from Bio import SeqIO
from Bio import Phylo
from sys import argv
import time

t0 = time.time() 
def main():
	script, seqfile, OTUfile, readcutoff = argv
	IDlist = [];OTUin= ''; IDolist=[]; IDdict={}; fastadict = {}; abunlist = []; abddict = {}
	outpath = OTUfile.split('/OTUs')[0]
	out= open(outpath + 'chimeras/Seq_reads_test.fas','w+')
	for OTU in open(OTUfile,'r'):
		if "OTU" in OTU.split('\t')[0]:
# 			IDlist.append(OTU.split('\t')[0])
			timea = (time.time()-t0)
			day = timea // (24 * 3600)
			timea = timea % (24 * 3600)
			hour = timea // 3600
			timea %= 3600
			minutes = timea // 60
			timea %= 60
			seconds = timea
			print(int(day),'d', int(hour),'h', int(minutes),'min', int(seconds), 'sec \t', OTUin,end='\r')
			if int(OTU.split('\t')[2]) not in abunlist and int(OTU.split('\t')[2]) >= int(readcutoff):
				OTUin=OTU.split('\t')[0]
				IDdict[OTU.split('\t')[0]]=OTU.split('\t')[2]
				IDolist.append(OTU.split('\t')[0])
				abunlist.append(int(OTU.split('\t')[2]))
				abddict.setdefault(str(OTU.split('\t')[2]),[]) 
				abddict[str(OTU.split('\t')[2])].append(OTU.split('\t')[0])
# 				print(abddict)
			else:
				if int(OTU.split('\t')[2]) >= int(readcutoff):
					OTUin=OTU.split('\t')[0]
					IDolist.append(OTU.split('\t')[0])
					IDdict[OTU.split('\t')[0]]=OTU.split('\t')[2]
					abddict[str(OTU.split('\t')[2])].append(OTU.split('\t')[0])
	print('\n',len(IDolist), len(abddict))
	for Seq in SeqIO.parse(seqfile,'fasta'):
		if Seq.id in IDolist:
			fastadict[Seq.id] = Seq.seq
			timea = (time.time()-t0)
			day = timea // (24 * 3600)
			timea = timea % (24 * 3600)
			hour = timea // 3600
			timea %= 3600
			minutes = timea // 60
			timea %= 60
			seconds = timea
			print("Creating seq dictionary\t:", int(day),'d', int(hour),'h', int(minutes),'min', int(seconds), 'sec \t', Seq.id,end='\r')
	for size in sorted(abunlist, reverse=True):
		print("Sorting abd dictionary\t:", int(day),'d', int(hour),'h', int(minutes),'min', int(seconds), 'sec \t',end='\r')
		for OTU in abddict[str(size)]:
			out.write('>'+OTU +';size='+IDdict[OTU]+';\n'+str(fastadict[OTU])+'\n')
			timea = (time.time()-t0)
			day = timea // (24 * 3600)
			timea = timea % (24 * 3600)
			hour = timea // 3600
			timea %= 3600
			minutes = timea // 60
			timea %= 60
			seconds = timea
			print("Writing seq file\t:",int(day),'d', int(hour),'h', int(minutes),'min', int(seconds), 'sec \t', OTU,end='\r')
	out.close()		
		
main()
