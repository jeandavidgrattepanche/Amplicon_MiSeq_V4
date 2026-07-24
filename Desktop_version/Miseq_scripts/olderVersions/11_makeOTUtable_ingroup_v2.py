#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "5.01, June 29, 2021"
__email__ = "jeandavid.grattepanche@gmail.com"


import sys
import os
import re
from Bio import SeqIO
from sys import argv
from random import randrange

#samples = []; readpersamplesdict= {}; listreaddict= {}
readsk = []; g = 0 ; readtokeepdict={}; OTUlist=[]; OTUtot=0
def countread(readmap, readtokeep, samplelist): 	
	folder = ('/').join(readtokeep.split('/')[:-1])
	folderout = folder.split('/outgroup')[0] + '/OTUs_ingroup/'
# 	print('test',folderout)
	
	if not os.path.exists(folderout):
		os.makedirs(folderout)
	output = open(folderout + readmap.split('/')[-1].split('.')[0] + '_subsampled.txt','w+')
	g = 0; tot = 0
	tot = sum(1 for line in open(readtokeep, 'r'))
	OTUlist=[]
	for read in open(readtokeep,'r'):
		if read.split(';')[0] not in OTUlist:
			OTUlist.append(read.split(';')[0])
			OTUtot=len(OTUlist)
	print("there is ",tot," reads and ", OTUtot," OTUs.")
	OTUlist=[]
	for read in open(readtokeep,'r'):
#		readsk.append(read.split('\n')[0])
		g += 1
		print(g,'reads (',round((g/tot)*100,2),"%) and added =>", len(OTUlist) ,'OTUs (',round((len(OTUlist)/OTUtot)*100,2),"%) ", end = '\r', flush=True)
		OTUIDs = read.split(';')[0]
		readname =read.split(';')[1].split('\n')[0]
		readtokeepdict.setdefault(OTUIDs,[])
		readtokeepdict[OTUIDs].append(readname)
		if read.split(';')[0] not in OTUlist:
			OTUlist.append(read.split(';')[0])
	print("\n[Done].\n\n")
	g=0
	for OTUID in OTUlist:
		g +=1
		output= open(folderout + readmap.split('/')[-1].split('.')[0]+'_subsampled.txt','a')
		print(OTUID, "has", len(readtokeepdict[OTUID]), 'reads (', round((g/OTUtot)*100,2),"%)", end="\r", flush=True)
		output.write(OTUID +'\t'+ str(readtokeepdict[OTUID]).replace("'",'').replace('[','').replace(']','').replace(', ','\t') + '\n')
		output.close()

	print("\n \n List of OTU and read done \n")
#	for line in open(readmap,'r'):
#		OTUID = line.split('\t')[0]
#		print(OTUID)
#		readlist = ''
#		for read in line.split('\t'  )[1:]:
#			if read.replace('\n','').replace(' ','') in readsk:
#				readlist = readlist + '\t' + read.replace('\n','').replace(' ','')	
#				print(OTUID, "has", str(readlist.count("-")), "reads")
#		print(OTUID, "is represented by a total of ",str(readlist.count("-")), "reads")
#		if readlist != '': 
#			output.write(OTUID + readlist + '\n')
	
def main():
	script, readmap, subsamplefile, listofsample = argv 
	countread(readmap, subsamplefile, listofsample) 
main()
