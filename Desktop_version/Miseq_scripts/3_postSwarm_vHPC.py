#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "2.01, June 29, 2021"
__email__ = "jeandavid.grattepanche@gmail.com"



import sys
import os
import re
import time
from Bio import SeqIO
from sys import argv

seqlist = []
unidict = {}
readdict = {}
start_time = time.time()
totalread = 0


def countread(swarmout, derepmapfile, derepmapprimer, derepseqprimer):
	outputpath = derepmapfile.split('/OTUs')[0]
	outseq = open(outputpath+"/OTUs/SWARM_postout.fas",'w+')
	outlist = open(outputpath+"/OTUs/SWARM_postout.txt",'w+')
	for line in open(derepmapfile,'r'):
		if line.split('\t')[0] == 'S':
			unidict.setdefault(line.split('\t')[8].split(';')[0], [line.split('\t')[8].replace('\n','')] )
		if line.split('\t')[0] == 'H':
			unidict[line.split('\t')[9].split(';')[0]].append(line.split('\t')[8].replace('\n',''))
	for linep in open(derepmapprimer,'r'):
		if linep.split('\t')[0] == 'H':
			unidict[linep.split('\t')[9].split(';')[0]].append(unidict[linep.split('\t')[8].split(';')[0]])

	for Seq in SeqIO.parse(open(derepseqprimer),'fasta'):
		readdict[Seq.description] = str(Seq.seq)
	i = 0
	totalread = 0
	for SWARM in open(swarmout, 'r'):
		i = i + 1
		seqlist = []
		representative = SWARM.split('\n')[0].split(" ")[0]
		outseq = open(outputpath+"OTUs/SWARM_postout.fas",'a')
		outseq.write(">OTU" + str(i)+"\t"+ representative + "\n" + str(readdict[representative]) + "\n")
		outseq.close()
		numberofread = 0
		for amplicon in SWARM.split('\n')[0].split(" "):
			seqlist.append(unidict[amplicon.split(';')[0]])
		outlist = open(outputpath+"OTUs/SWARM_postout.txt",'a')
		outlist.write("OTU" + str(i) + "\t" + str(seqlist).replace(',','\t').replace('[','').replace(']','').replace("'","") + "\n")
		outlist.close() 
		
def countreadnp(swarmout, derepmapfile, derepseq):
	outputpath = derepmapfile.split('/OTUs')[0]
	outseq = open(outputpath+"/OTUs/SWARM_postout.fas",'w+')
	outlist = open(outputpath+"/OTUs/SWARM_postout.txt",'w+')
	for line in open(derepmapfile,'r'):
		if line.split('\t')[0] == 'S':
			unidict.setdefault(line.split('\t')[8].split(';')[0], [line.split('\t')[8].replace('\n','')] )
		if line.split('\t')[0] == 'H':
			unidict[line.split('\t')[9].split(';')[0]].append(line.split('\t')[8].replace('\n',''))

	for Seq in SeqIO.parse(open(derepseq),'fasta'):
		readdict[Seq.description] = str(Seq.seq)
	i = 0
	totalread = 0
	for SWARM in open(swarmout, 'r'):
		i = i + 1
		seqlist = []
		representative = SWARM.split('\n')[0].split(" ")[0]
		outseq = open(outputpath+"/OTUs/SWARM_postout.fas",'a')
		outseq.write(">OTU" + str(i)+"\t"+ representative + "\n" + str(readdict[representative]) + "\n")
		outseq.close()
		numberofread = 0
		for amplicon in SWARM.split('\n')[0].split(" "):
			seqlist.append(unidict[amplicon.split(';')[0]])
		outlist = open(outputpath+"/OTUs/SWARM_postout.txt",'a')
		outlist.write("OTU" + str(i) + "\t" + str(seqlist).replace(',','\t').replace('[','').replace(']','').replace("'","") + "\n")
		outlist.close() 
	
def main():
	try:
		script, swarmout, derepmapfile, derepmapprimer, derepseqprimer = argv
		countread(swarmout, derepmapfile, derepmapprimer , derepseqprimer)
	except:
	
		script, swarmout, derepmapfile, derepseq = argv
		countreadnp(swarmout, derepmapfile, derepseq)
		
main()