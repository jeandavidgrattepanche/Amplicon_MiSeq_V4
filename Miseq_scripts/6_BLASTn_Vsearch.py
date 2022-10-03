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

def getBLAST( NGSfile, idmin, qcov, Taxa): #, readcutoff):
	print("start BLAST SSU_Euk_pr2_version_4.14.0")
	outputpath = NGSfile.split('chimeras')[0]
# 	outblast = open(outputpath+'/VsearchBLAST.tsv','w+')
	ublast_self = vsearch_path + ' --threads 64 --usearch_global '+NGSfile.split(".fas")[0]+'_reduced.fas --db '+SSU_db+ ' --strand both --id '+str(idmin/100)+' --query_cov '+ str(qcov/100)+' --blast6out '+outputpath+'/VsearchBLAST.tsv ' ## No -evalue 1e-15 as usearch
	print(ublast_self)
# 	os.system(ublast_self)
# 	os.system("head -1000 "+outputpath+"/VsearchBLAST.tsv > "+outputpath+"/VsearchBLAST_test.tsv")
# 	###
# 	### Replace VsearchBLAST.tsv by VsearchBLAST_test.tsv for testing
# 	###
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
	outblastsp = open(outputpath+'/VsearchBLAST_sp.tsv','w+')
#	print(blastdict.items())
	for key, value in blastdict.items():
#		print(key)
#		print(blastdict[key])
		outblastc.write(blastdict[key]+'\n')
		if Taxa in blastdict[key]:
			outblastsp.write(blastdict[key]+'\n')
	outblastc.close()

	outseqIN = open(outputpath+'taxonomic_assignment/Seq_reads_inGroup_vsearch_pr2_4.14.0.fasta','w+')
	for seq in SeqIO.parse(NGSfile.split(".fas")[0]+'_reduced.fas','fasta'):
		try:
			blastdict[seq.id]

			if Taxa in blastdict[seq.id]: #ID.split('p:')[1].split(',')[0] == str(Taxa):# or ID.split('_')[1] == Taxa (need to check PR@ format):
				print(seq.id, 'blasted with', ID.split(';size=')[0] , " at ", ident , "% and coverage:", cov )
				outseqIN = open(outputpath+'taxonomic_assignment/Seq_reads_inGroup_vsearch_pr2_4.14.0.fasta','a')
				outseqIN.write('>'+seq.description+ '_'+ ID.split('_rid_')[0] + '_' +str(cov)+'_'+ str(Sim) + '%\n'+str(seq.seq) + '\n')
				outseqIN.close()

def main():
	script,  NGSfile, idminy, qcovz, Taxa, readcutoff = argv
	outseq = open(NGSfile.split(".fas")[0]+'_reduced.fas','w+')
	for seq in SeqIO.parse(NGSfile,'fasta'):
		if int(seq.description.split(";size=")[1]) > int(readcutoff):
			outseq.write('>'+seq.description+ '\n'+str(seq.seq) + '\n')
	outseq.close()		
	getBLAST(NGSfile, float(idminy),float(qcovz), Taxa)
main()
