#!/usr/bin/python3
#if many more OTU than 100,000 this version should be faster
__author__ = "Jean-David Grattepanche"
__version__ = "5.2, September 8,2022"
__email__ = "jeandavid.grattepanche@gmail.com"


import sys
import os
import re
from Bio import SeqIO
from sys import argv
from random import randrange

u=0;samples = []; tokeep = []; readpersamplesdict= {}; listreaddict= {}

def countread(seqfile, readmap,samplelist): 	
	u = 0 
	for seq in SeqIO.parse(open(seqfile,'r'),'fasta'):
		tokeep.append(seq.id.split('_')[0])
		u+=1
		print(int(u), end ='\r')
	print('\n\n\n',len(tokeep),'\n\n\n')
	outlog = open('readpersample_cleaned_test.txt','w+')
	folder = ('/').join(readmap.split('/')[:-1])
	print("saving folder: ", folder)
	maxread = 0
	for sample in open(samplelist,'r'):
		samples.append(sample.replace('_','-').split('\t')[1].split('\n')[0])
		listreaddict.setdefault(sample.replace('_','-').split('\t')[1].split('\n')[0], [])
	print(samples)
	for line in open(readmap,'r'):
		OTUID = line.split('\t')[0]
		print('\n',OTUID)
		if OTUID in tokeep:
			for read in line.split('\t'  )[1:]:
				samplename = ("-").join(read.replace(" ","").replace("'","").split('_')[:-1])
				print(samplename,end='\t')
				if samplename in samples:
					print("added",end='\t')
					listreaddict[samplename].append(OTUID+";"+read)
# 			else:
# 				print(OTUID, ' was removed in previous step')

	readnumber= 0 
	for sampler in open(samplelist,'r'):
		readnumber= len(listreaddict[sampler.replace('_','-').split('\t')[1].split('\n')[0]])

		readpersamplesdict[sampler.replace('_','-').split('\t')[1].split('\n')[0]] = str(readnumber)
		print(sampler.replace('_','-').split('\t')[1].split('\n')[0], " has ", str(readnumber), " reads.") 
		if readnumber > maxread:
			maxread = readnumber
# 			print(maxread)
		outlog = open('readpersample_cleaned_test.txt','a')
		outlog.write( sampler.replace('_','-').split('\t')[1].split('\n')[0]+ "\t"+ str(readnumber)+'\n')
		outlog.close()
		readnumber= 0
	outlog.close()	

	subsamp = input('Do you need to subsample? (yes or no) \n')
	if subsamp[0] == 'y':
		j = input('How many reads for each file? ( hit return for default of 10,000) \n')
		try:
			randomnum = int(j) + 1
		except TypeError:
			print ('Your input must be a number.  Try again. ')
			main()
		except ValueError:
			j = ""	
		if j == "":
			randomnum = 10000
		else:
			randomnum = int(j)
	elif subsamp[0] == 'n':
		randomnum = maxread
		print(randomnum)
	else:
		print ('Please answer yes or no. ')
		main()

	outfile = open(folder+'/subsampled_test.txt','w+')
	for sample2 in open(samplelist,'r'):
		samplelistname = sample2.replace('_','-').split('\t')[1].split('\n')[0]	
		if subsamp[0] == 'n':
			print(samplelistname, "no subsampling")
			for name in listreaddict[samplelistname]:
				outfile = open(folder+'/subsampled.txt','a')
				outfile.write(name.replace('\n','').replace(' ','') + '\n')
		if subsamp[0] == 'y':
			if int(readpersamplesdict[samplelistname]) >= int(randomnum):
				print('Start ',randomnum,' random picking ',samplelistname,' of ', readpersamplesdict[samplelistname])
				read_to_keep = []
				count = 0
				while count < int(randomnum):
					i = randrange(0, int(readpersamplesdict[samplelistname]))
					if i not in read_to_keep:
						read_to_keep.append(i)
						count = count + 1
				for k in read_to_keep:
					name = listreaddict[samplelistname][k]
					outfile = open(folder+'/subsampled_test.txt','a')
					outfile.write(name.replace('\n','').replace(' ','') + '\n')
			else:
				print(samplelistname, " should not be included because it contains few reads (", readpersamplesdict[samplelistname], ")")
				for k in range(0, int(readpersamplesdict[samplelistname])):
					name = listreaddict[samplelistname][k]
					outfile = open(folder+'/subsampled_test.txt','a')
					outfile.write(name.replace('\n','').replace(' ','') + '\n')
		outfile.close()		
	
def main():
	script, seqfile, otufile, listofsample = argv 
	countread(seqfile, otufile,listofsample) 
main()
