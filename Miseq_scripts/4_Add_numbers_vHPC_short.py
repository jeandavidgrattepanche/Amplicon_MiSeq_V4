#!/usr/bin/python3
#not seqfile here
#old	os.system('python Miseq_scripts/4_Add_numbers_vHPC.py ' + outputpath + '/OTUs/SWARM_postout.fas ' + outputpath + '/OTUs/SWARM_postout.txt '+listsample)# +' '+ dataname) #  '+ str(runref))
#correct	os.system('python Miseq_scripts/4_Add_numbers_vHPC_short.py '  + outputpath + '/OTUs/SWARM_postout.txt '+listsample)

__author__ = "Jean-David Grattepanche"
__version__ = "6.01, June 29, 2021"
__email__ = "jeandavid.grattepanche@gmail.com"


import sys
import os
import re
from sys import argv

seqlist = {}

def countread(otufile,samplelist) : 	
	outpath = otufile.split('/OTUs')[0]
	outfile = open(outpath + '/OTUs/OTUtable_temp.txt','w')
	outfile.write('SWARM\toccurrence\treadnumber\t' + str(samplelist).replace("', '",'\t').replace("[",'').replace("]",'').replace("'","") + '\n') #add the heading row with samples name
	outfile.close()

	for line in open(otufile,'r'):
		OTUID = line.split('\t')[0]
		allread = []; occlist = []; abundance = []; readnumber= 0; r34=0
		for read in line.split('\t'  )[1:]:
			samplename = ('_').join(read.replace(' ','').split('_')[:-1])
#			print(samplename)
			if samplename in samplelist:
				toadd= samplename + ' ,'
				allread.extend([samplename for x in range(int(read.split(';size=')[1]))])
				r34 += int(read.split(';size=')[1])
				if samplename not in occlist:
					occlist.append(samplename)
# 			else:
# 				print("ERROR in list")
		occurrence = len(occlist)
		totalread = len(allread)
		print(OTUID, " has occurred in ", occurrence, " samples and is represented by ", totalread, "or", int(r34), " reads.", end='\r',flush=True) 
		if totalread != r34:
			break
		if totalread > 1:
			if occurrence > 0:
		
				for R in samplelist:
					readnumber = allread.count(R)
					abundance.append(readnumber)
		
		
#				print(OTUID, " represents", totalread , " for ", occurrence, "samples")
				outfile = open(outpath + '/OTUs/OTUtable_temp.txt','a')
				outfile.write( OTUID +'\t' + str(occurrence) + '\t' + str(totalread)  + '\t'+ str(abundance).replace(',','\t').replace('[','').replace(']','') +  '\n')
				outfile.close()
	
	
				
		
	
def main():
	script, otufile, listofsample = argv 
	samplefile = open(listofsample,'r')
	samplelist= []
	for sample in samplefile:
		if sample.split('\n')[0] != "":
			samplelist.append(sample.split('\n')[0].split('\t')[1])
	print(samplelist)
	countread(otufile,samplelist) 
main()
