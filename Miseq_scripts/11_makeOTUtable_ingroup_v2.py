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
readsk = []; g = 0 ; readtokeepdict={}; OTUlist=[]
def countread(readmap, readtokeep, samplelist): 	
	folder = ('/').join(readtokeep.split('/')[:-1])
	folderout = folder.split('/outgroup')[0] + '/OTUs_ingroup/'
	print('test',folderout)
	
	if not os.path.exists(folderout):
		os.makedirs(folderout)
	output = open(folderout + readmap.split('/')[-1].split('.')[0] + '_subsampled.txt','w+')
	g = 0; tot = 0
	tot = sum(1 for line in open(readtokeep, 'r'))
	print(tot)
	for read in open(readtokeep,'r'):
#		readsk.append(read.split('\n')[0])
		g += 1
		print(g,'(',round((g/tot)*100,2),"%) and added =>", len(OTUlist) ,'(',round((len(OTUlist)/tot)*100,2),"%) ", end = '\r')
		OTUIDs = read.split(';')[0]
		readname =read.split(';')[1].split('\n')[0]
		readtokeepdict.setdefault(OTUIDs,[])
		readtokeepdict[OTUIDs].append(readname)
		if OTUIDs not in OTUlist:
			OTUlist.append(OTUIDs)
	print(len(OTUlist))

	for OTUID in OTUlist:
		output= open(folderout + readmap.split('/')[-1].split('.')[0]+'_subsampled.txt','a')
		print(OTUID, "has", len(readtokeepdict[OTUID]), 'reads')
		output.write(OTUID +'\t'+ str(readtokeepdict[OTUID]).replace("'",'').replace('[','').replace(']','').replace(', ','\t') + '\n')
		output.close()

	print("List of OTU and read done \n")
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
	script, otufile, subsamplefile, listofsample = argv 
	countread(otufile, subsamplefile, listofsample) 
main()
