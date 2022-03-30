#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "5.01, June 29, 2021"
__email__ = "jeandavid.grattepanche@gmail.com"


import sys
import os
import re
from Bio import SeqIO
from sys import argv

seqlist = {}; BLAST= {}; TA_tree= {}

def countread(seqfile,BLASTtsv, otufile,samplelist): 	
	outpath = ('/').join(BLASTtsv.split('/')[:-1])
	for Seq in SeqIO.parse(open(seqfile),'fasta'):
		seqlist[Seq.id.split('_')[0].split('-')[0]] = [Seq.id,str(Seq.seq)]

	for record in open(BLASTtsv,'r'):
*		BLASTresults = record.split('\n')[0].split('\t')[-1].replace(' ','-')+';'+record.split('\t')[6]+';'+record.split('\t')[7]+';'+record.split('\t')[8]
*		print(BLASTresults)
*		BLAST[record.split('_')[0].split('\t')[0]] = BLASTresults
		
	outfile = open(outpath+'/OTUs_ingroup/OTUtable_ingroup_noTtree.txt','w')
	outfile.write('OTU\tBtaxo_rank1\tBtaxo_rank2\tBtaxo_rank3\tBtaxo_rank4\tBsp\tBacc_number\tid%\tEvalue\toccurrence\treadnumber\t' + str(samplelist).replace("', '",'\t').replace("[",'').replace("]",'').replace("'","") + '\n') #add the heading row with samples name
	outfile.close()
	outseq = open(outpath+'/OTUs_ingroup/' +seqfile.split('/')[-1].split(".")[0]+ "_norare.fas",'w+')

	for line in open(otufile,'r'):
		OTUID = line.split('\t')[0]
		allread = []
		occlist = []
		abundance = []
		readnumber= 0
		for read in line.split('\n')[0].split('\t'  )[1:]:
			samplename = ('-').join(read.replace(" ","").replace('"',"").replace('"','').split('_')[0:3])
			if samplename in samplelist:
				allread.append(samplename)
				if samplename not in occlist:
					occlist.append(samplename)
			else:
				print("ERROR in list", samplename)
		occurrence = len(occlist)
		totalread = len(allread)
		print(OTUID, " has occurred in ", occurrence, " samples and is represented by ", totalread, " reads.") 
		if totalread > 0:
			if occurrence > 0:
		
				for R in samplelist:
					readnumber = allread.count(R)
					abundance.append(readnumber)
		
		
#				print(OTUID, " represents", totalread , " for ", occurrence, "samples")
				outfile = open(outpath+'/OTUs_ingroup/OTUtable_ingroup.txt','a')
# 				if seqlist[OTUID][0].split('_')[2] != 'No':
# 					print(seqlist[OTUID][0].split('_')[0],' ',seqlist[OTUID][0].split('_')[1])
# 					nameB = seqlist[OTUID][0].split('_')[0] + '\t' +seqlist[OTUID][0].split('_')[2]+ '\t' +seqlist[OTUID][0].split('_')[3]+ '\t' +seqlist[OTUID][0].split('_')[4]+ '\t' +seqlist[OTUID][0].split('_')[5]+ '\t' +seqlist[OTUID][0].split('_')[6]+ '\t' +seqlist[OTUID][0].split('_')[7]+ '_' +seqlist[OTUID][0].split('_')[8] +'\t' +seqlist[OTUID][0].split('_')[-3]+'\t' +seqlist[OTUID][0].split('_')[-2] +'\t' +seqlist[OTUID][0].split('_')[-1]
# 				else:
# 					nameB = seqlist[OTUID][0].split('_')[0] + '\t' +seqlist[OTUID][0].split('_')[1]+ '\t\t\t\t\t\t\t\t'
				nameB = BLAST[OTUID].replace(';','\t').replace('_','\t')
				outfile.write(OTUID+'\t'+ nameB + '\t'+ str(occurrence) + '\t' + str(totalread)  + '\t'+ str(abundance).replace(',','\t').replace('[','').replace(']','') +  '\n')
				outfile.close()
				outseq = open(outpath+'/OTUs_ingroup/' +seqfile.split('/')[-1].split(".")[0]+ "_norare.fas",'a')
				outseq.write(">" + seqlist[OTUID][0]+ '\n' + seqlist[OTUID][1] + '\n')
				outseq.close()
	
	
				
		
	
def main():
	script, seqfile, BLASTtsv, otufile, listofsample = argv #, dataname = argv
	samplelist= []
	for sample in open(listofsample,'r'):
		if sample.split('\n')[0] != "":
			samplelist.append(sample.replace("_","-").split('\n')[0].split('\t')[1])
	print(samplelist)
	countread(seqfile,BLASTtsv, otufile,samplelist) #,dataname)
main()
