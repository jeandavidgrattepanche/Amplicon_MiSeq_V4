#!/usr/bin/python 

__author__ = "Jean-David Grattepanche"
__version__ = "1, March 15,2017"
__email__ = "jeandavid.grattepanche@gmail.com"



import sys
import os
import re
import time
import string
import os.path
from Bio import SeqIO
from sys import argv


def main():
	script, treefileout, treefileall, otu_read, seqnosingleton = argv
	nosingletonlist = []; seqdict= {}
	for Seq in SeqIO.parse(open(seqnosingleton,'r'),'fasta'):
		nosingletonlist.append(Seq.id.split('_')[0])
		seqdict[Seq.id.split('_')[0]] =  [Seq.description, str(Seq.seq)]
	print(len(nosingletonlist))
# 	folder = treefileall.split('/')[0]+'/'+treefileall.split('/')[1]
	folder = ('/').join(treefileall.split('/')[:-1])
	print(treefileall, folder)
	out = open(folder+'/'+otu_read.split('/')[-1].split('.')[0] + '_nosingleton_nochimeras_in_only.txt','w+')
	outseq = open(folder+'/'+otu_read.split('/')[-1].split('.')[0] + '_nosingleton_nochimeras_in_only.fasta','w+')
	OTUReaddict = {}; readlist = []
	treeout = open(treefileout,'r').readline()
	treeall = open(treefileall,'r').readline()
	for line in open(otu_read,'r'):
		allread = []
		OTUname = 'QUERY___' +line.split('\t')[0] + "_"
# 		print(OTUname)
		if OTUname in treeall:
			if OTUname in treeout:
				print('outgroup OTU: ', line.split('\t')[0])	
			else:
				if line.split('\t')[0] in nosingletonlist:
					outseq.write('>' + seqdict[line.split('\t')[0]][0] + '\n' + seqdict[line.split('\t')[0]][1] + '\n')
					for read in line.split('\t'  )[1:]:
						samplename = ('-').join(read.split('_')[:-1])+'_'+read.split('_')[-1]
# 						samplename = read.replace('_','-').split('.')[0]
						allread.extend([samplename for x in range(int(read.split(';size=')[1]))])
# 						print(samplename, " has" ,allread)
					out.write(line.split('\t')[0]+ '\t'+ str(allread).replace(',','\t').replace('[','').replace(']','') +  '\n')
					print(line.split('\t')[0],' has ',len(allread),' reads')
		else:
			if line.split('\t')[0] in nosingletonlist:
				print(line.split('\t')[0], ' is not in the tree. TAKE A look')
			
	out.close()

main()