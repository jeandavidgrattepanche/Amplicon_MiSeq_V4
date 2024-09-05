#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "3.12, June 21, 2023"
__email__ = "jeandavid.grattepanche@gmail.com"



import sys, os, re, time, string, os.path
from distutils import spawn
from Bio import SeqIO
from sys import argv

vsearch_path = spawn.find_executable("vsearch")
SSU_db = "db_v4/SILVA_138.2_SSURef_NR99_tax_silva.fasta" ##change to the correct database
blastdict = {} ; taxdict = {}

def getBLAST( NGSfile, otumap, idmin, qcov, Taxa): #, readcutoff):
	for ref in SeqIO.parse(SSU_db,"fasta"):
		taxdict[ref.id] = ('_').join(ref.description.split('\t')[:1])
		
	print("start BLAST SILVA SSU NR")
	outputpath = NGSfile.split('chimeras')[0]
	if not os.path.exists(outputpath + 'taxonomic_assignment/'): 
		os.makedirs(outputpath + 'taxonomic_assignment/') 

	outblast = open(outputpath+'/VsearchBLAST_SILVA.tsv','w+')
	ublast_self = vsearch_path + ' --threads 64 --usearch_global '+NGSfile.split(".fas")[0]+'_SILVA_reduced.fas --db '+SSU_db+ ' --strand both --id '+str(idmin/100)+' --query_cov '+ str(qcov/100)+' --blast6out '+outputpath+'/VsearchBLAST_SILVA.tsv ' ## No -evalue 1e-15 as usearch
	print(ublast_self)
	os.system(ublast_self)
# 	os.system("head -1000 "+outputpath+"/VsearchBLAST.tsv > "+outputpath+"/VsearchBLAST_test.tsv")
# 	### Replace VsearchBLAST.tsv by VsearchBLAST_test.tsv for testing
	for blast_record in open(outputpath+'/VsearchBLAST_SILVA.tsv','r'):
		try:
			print(float(blast_record.split('\t')[2]), " <> ",float(blastdict[blast_record.split('\t')[0]].split('\t')[2]), end='\r', flush=True)
			if float(blast_record.split('\t')[2]) >  float(blastdict[blast_record.split('\t')[0]].split('\t')[2]):
				blastdict[blast_record.split('\t')[0]] = taxdict[blast_record.split('\t')[1]]+'\t'+blast_record.split('\n')[0]
			else:
				print(blast_record, "duplicated")
		except:
			blastdict[blast_record.split('\t')[0]] = taxdict[blast_record.split('\t')[1]]+'\t'+blast_record.split('\n')[0]
			print(blast_record," added")
	outblastc = open(outputpath+'/VsearchBLAST_SILVA_clean.tsv','w+')
	outblastsp = open(outputpath+'/VsearchBLAST_SILVA_inGroup.tsv','w+')
	tokeep = []
#	print(blastdict.items())
	for key, value in blastdict.items():
		outblastc.write(blastdict[key]+'\n')
		if Taxa in blastdict[key]:
			outblastsp.write(blastdict[key]+'\n')
			tokeep.append(key.split('_')[0])
	outblastc.close()

	outseqIN = open(outputpath+'taxonomic_assignment/Seq_reads_inGroup_SILVA.fasta','w+')
	for seqi in SeqIO.parse(NGSfile.split(".fas")[0]+'_SILVA_reduced.fas','fasta'):
		try:
			blastdict[seqi.id]
			if Taxa in blastdict[seqi.id]: #ID.split('p:')[1].split(',')[0] == str(Taxa):# or ID.split('_')[1] == Taxa (need to check PR@ format):
				print(seqi.id, 'blasted with', blastdict[seqi.id], end='\r', flush=True )
				outseqIN = open(outputpath+'taxonomic_assignment/Seq_reads_inGroup_SILVA.fasta','a')
				outseqIN.write('>'+seqi.description+'\n'+str(seqi.seq) + '\n')
				outseqIN.close()
		except:
			print(seqi, " not in ??")
	s=0; t=0
	outmapIN = open(outputpath+'taxonomic_assignment/Seq_map_inGroup_SILVA.txt','w+')
	for line in open(otumap, 'r'):
		t+=1
		if line.split('\t')[0] in tokeep:
			s+=1
			outmapIN.write(line)
		print(line.split('\t')[0], round((s/t)*100,1),"%",end='\r')
def main():
	script, NGSfile, otumap, idminy, qcovz, Taxa, readcutoff = argv
	outseq = open(NGSfile.split(".fas")[0]+'_SILVA_reduced.fas','w+')
	for seq in SeqIO.parse(NGSfile,'fasta'):
		if int(seq.description.split("_")[1].replace('r','')) > int(readcutoff):
			outseq.write('>'+seq.description+ '\n'+str(seq.seq) + '\n')
	outseq.close()
	getBLAST(NGSfile, otumap, float(idminy),float(qcovz), Taxa)
main()
