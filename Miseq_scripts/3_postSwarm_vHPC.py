#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "2, October 12, 2016"
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

# outputpath = '/home/tuk61790/v4_SWARM/outputs/'
# 
# outseq = open(outputpath+"OTUs/SWARM_postout.fas",'w+')
# outlist = open(outputpath+"OTUs/SWARM_postout.txt",'w+')

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
#			unidict[linep.split('\t')[9].replace('\n','')].append(linep.split('\t')[8].replace('\n',''))
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
# 			print(numberofread, len(unidict[amplicon]), unidict[amplicon])
# 			numberofread = numberofread + len(unidict[amplicon]) #seqdict[amplicon.split(";")[0]].count(runname)
# 		print('OTU',str(i)," is represented by ", representative,  "and has ", str(numberofread), "total reads")
# 		totalread = totalread + numberofread
		outlist = open(outputpath+"OTUs/SWARM_postout.txt",'a')
		outlist.write("OTU" + str(i) + "\t" + str(seqlist).replace(',','\t').replace('[','').replace(']','').replace("'","") + "\n")
		outlist.close() 
# 	print("Total number of reads: " , str(totalread))

# 	i = 0
# 	for SWARM in open(statSWARM, 'r'):
# 		i = i + 1
# 		seqlist = []
# 		representative = SWARM.split('\t')[2]
# 		outseq = open("outputs/OTUs/SWARM_postout.fas",'a')
# 		outseq.write(">OTU" + str(i)+ '-'+SWARM.split('\t')[1]+ "r\n" + str(readdict[representative]) + "\n")
# 		outseq.close()
# 		
# 		
# 	outunique = open("outputs/OTUs/dereplicated_listunique.txt",'a')			
		
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
# 			print(numberofread, len(unidict[amplicon]), unidict[amplicon])
# 			numberofread = numberofread + len(unidict[amplicon]) #seqdict[amplicon.split(";")[0]].count(runname)
# 		print('OTU',str(i)," is represented by ", representative,  "and has ", str(numberofread), "total reads")
# 		totalread = totalread + numberofread
		outlist = open(outputpath+"/OTUs/SWARM_postout.txt",'a')
		outlist.write("OTU" + str(i) + "\t" + str(seqlist).replace(',','\t').replace('[','').replace(']','').replace("'","") + "\n")
		outlist.close() 
# 	print("Total number of reads: " , str(totalread))

# 	i = 0
# 	for SWARM in open(statSWARM, 'r'):
# 		i = i + 1
# 		seqlist = []
# 		representative = SWARM.split('\t')[2]
# 		outseq = open("outputs/OTUs/SWARM_postout.fas",'a')
# 		outseq.write(">OTU" + str(i)+ '-'+SWARM.split('\t')[1]+ "r\n" + str(readdict[representative]) + "\n")
# 		outseq.close()
# 		
# 		
# 	outunique = open("outputs/OTUs/dereplicated_listunique.txt",'a')			
	
def main():
	try:
		script, swarmout, derepmapfile, derepmapprimer, derepseqprimer = argv
		countread(swarmout, derepmapfile, derepmapprimer , derepseqprimer)
	except:
	
		script, swarmout, derepmapfile, derepseq = argv
		countreadnp(swarmout, derepmapfile, derepseq)
		
main()