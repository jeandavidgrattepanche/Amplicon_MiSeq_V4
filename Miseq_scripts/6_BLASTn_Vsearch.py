#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "3.11, October 3, 2022"
__email__ = "jeandavid.grattepanche@gmail.com"



import sys, os, re, time, string, os.path
from distutils import spawn
from Bio import SeqIO
from sys import argv

vsearch_path = spawn.find_executable("vsearch")
SSU_db = "db_v4/pr2_version_4.14.0_SSU_UTAX.fasta" ##change to the correct database
blastdict = {} 

def getBLAST( NGSfile, otumap, idmin, qcov, Taxa): #, readcutoff):
	print("start BLAST SSU_Euk_pr2_version_4.14.0")
	outputpath = NGSfile.split('chimeras')[0]
	outblast = open(outputpath+'/VsearchBLAST.tsv','w+')
	ublast_self = vsearch_path + ' --threads 64 --usearch_global '+NGSfile.split(".fas")[0]+'_reduced.fas --db '+SSU_db+ ' --strand both --id '+str(idmin/100)+' --query_cov '+ str(qcov/100)+' --blast6out '+outputpath+'/VsearchBLAST.tsv ' ## No -evalue 1e-15 as usearch
	print(ublast_self)
	os.system(ublast_self)
# 	os.system("head -1000 "+outputpath+"/VsearchBLAST.tsv > "+outputpath+"/VsearchBLAST_test.tsv")
# 	### Replace VsearchBLAST.tsv by VsearchBLAST_test.tsv for testing
	for blast_record in open(outputpath+'/VsearchBLAST.tsv','r'):
		try:
			print(float(blast_record.split('\t')[2]), " <> ",float(blastdict[blast_record.split('\t')[0]].split('\t')[2]))
			if float(blast_record.split('\t')[2]) >  float(blastdict[blast_record.split('\t')[0]].split('\t')[2]):
				blastdict[blast_record.split('\t')[0]] = blast_record.split('\n')[0]
			else:
				print(blast_record, "duplicated")
		except:
			blastdict[blast_record.split('\t')[0]] = blast_record.split('\n')[0]
			print(blast_record," added")
	outblastc = open(outputpath+'/VsearchBLAST_clean.tsv','w+')
	outblastsp = open(outputpath+'/VsearchBLAST_inGroup.tsv','w+')
	tokeep = []
#	print(blastdict.items())
	for key, value in blastdict.items():
		outblastc.write(blastdict[key]+'\n')
		if Taxa in blastdict[key]:
			outblastsp.write(blastdict[key]+'\n')
			tokeep.append(key)
	outblastc.close()

	outseqIN = open(outputpath+'taxonomic_assignment/Seq_reads_inGroup.fasta','w+')
	for seqi in SeqIO.parse(NGSfile.split(".fas")[0]+'_reduced.fas','fasta'):
		try:
			blastdict[seqi.id]
			if Taxa in blastdict[seqi.id]: #ID.split('p:')[1].split(',')[0] == str(Taxa):# or ID.split('_')[1] == Taxa (need to check PR@ format):
				print(seqi.id, 'blasted with', blastdict[seqi.id] )
				outseqIN = open(outputpath+'taxonomic_assignment/Seq_reads_inGroup.fasta','a')
				outseqIN.write('>'+seqi.description+'\n'+str(seqi.seq) + '\n')
				outseqIN.close()
		except:
			print(seqi, " not in ??")
	s=0; t=0
	outmapIN = open(outputpath+'taxonomic_assignment/Seq_map_inGroup.txt','w+')
	for line in open(otumap, 'r'):
		t+=1
		if line.split('\t')[0] in tokeep:
			s+=1
			outmapIN.write(line)
		print(line.split('\t')[0], round((s/t)*100,1),"%",end='\r')
def main():
	script, NGSfile, otumap, idminy, qcovz, Taxa, readcutoff = argv
	outseq = open(NGSfile.split(".fas")[0]+'_reduced.fas','w+')
	for seq in SeqIO.parse(NGSfile,'fasta'):
		if int(seq.description.split(";size=")[1]) > int(readcutoff):
			outseq.write('>'+seq.description+ '\n'+str(seq.seq) + '\n')
	outseq.close()
	getBLAST(NGSfile, otumap, float(idminy),float(qcovz), Taxa)
main()
