#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "3, March 2, 2017"
__email__ = "jeandavid.grattepanche@gmail.com"



import sys, os, re, time, string, os.path
from distutils import spawn
from Bio import SeqIO
from sys import argv

vsearch_path = spawn.find_executable("vsearch")
SSU_db = "db_v4/SSU_Euk_PR2plus_V16_d.fasta" ##change to the correct database
blastdict = {}

def getBLAST( NGSfile, idmin, qcov, Taxa, readcutoff):
	print("start BLAST with custom dabase + PR2")
	outputpath = NGSfile.split('outgroup')[0]
	outblast = open(outputpath+'/VsearchBLAST_4.tsv','w+')
	ublast_self = vsearch_path + ' --thread 32 --usearch_global '+NGSfile+' --db '+SSU_db+ ' --strand both --id '+str(idmin/100)+' --query_cov '+ str(qcov/100)+' --blast6out '+outputpath+'/VsearchBLAST_4.tsv ' ## No -evalue 1e-15 as usearch
	print(ublast_self)
	os.system(ublast_self)
	print('Vsearch with ', SSU_db, ' and ', NGSfile, " done. \n Creating a dictionnary ... it takes a while!!")
	outblast.close()
	qb =0
	tot = sum(0.5 for line in open(NGSfile,'r'))
	os.system("awk '!seen[$1]++' "+outputpath+'/VsearchBLAST_4.tsv > '+outputpath+'/VsearchBLAST_4sorted.tsv')
	for blast_record in open(outputpath+'/VsearchBLAST_4sorted.tsv','r'):
		print(int(qb), '(', round((int(qb)/tot)*100,2), '%)', end = '\r')
#		if blast_record.split('\t')[0] not in blastdict.values():
		qb += 1
		blastdict[blast_record.split('\t')[0]] = blast_record.split('\n')[0]
#		else:
#			print(blast_record, "duplicated")
	print("Dictionary done! Merging sequence and taxonomic data... it takes a while too !")
	outseq = open(outputpath+'taxonomic_assignment/Seq_reads_nochimera_nosingleton_vsearch_cdb.fasta','w+')
	outseqSAR = open(outputpath+'taxonomic_assignment/Seq_reads_nochimera_nosingleton_specTaxa_vsearch_cdb.fasta','w+')
	for seq in SeqIO.parse(NGSfile,'fasta'):
		try:
			blastdict[seq.id]
			ident =  blastdict[seq.id].split('\t')[2]
			seqmatch = float(blastdict[seq.id].split('\t')[3]) -float(blastdict[seq.id].split('\t')[5])
			cnt = blastdict[seq.id].split('\t')[3] 
			Sim = round(float(seqmatch-int(blastdict[seq.id].split('\t')[4])) / float(cnt) * 100)
			ID = blastdict[seq.id].split('\t')[1]
			seqused = 1 + int(blastdict[seq.id].split('\t')[7]) - int(blastdict[seq.id].split('\t')[6])
			cov = round(float(seqused) / float(len(seq.seq)) * 100)
			outseq = open(outputpath+'taxonomic_assignment/Seq_reads_nochimera_nosingleton_vsearch_cdb.fasta','a')
			if int(seq.description.split('_')[1].replace('r','')) > (int(readcutoff)-1):
				outseq.write('>'+seq.description+ '_'+ ID.split('_rid_')[0] + '_' +str(cov)+'_'+ str(Sim) + '%\n'+str(seq.seq) + '\n')
				outseq.close()

			if ID.split('p:')[1].split(',')[0] == str(Taxa):# or ID.split('_')[1] == Taxa (need to check PR@ format):
				if int(seq.description.split('=')[1].replace('r','')) > (int(readcutoff)-1):
					print(seq.id, 'blasted with', ID.split(';size=')[0] , " at ", ident , "% and coverage:", cov, end = '\r' )
					outseqSAR = open(outputpath+'taxonomic_assignment/Seq_reads_nochimera_nosingleton_specTaxa_vsearch_cdb.fasta','a')
					outseqSAR.write('>'+seq.description+ '_'+ ID.split('_rid_')[0] + '_' +str(cov)+'_'+ str(Sim) + '%\n'+str(seq.seq) + '\n')
					outseqSAR.close()
		except:
			print("NO BLAST for ",seq.id, end = '\r')
			outseq = open(outputpath+'taxonomic_assignment/Seq_reads_nochimera_nosingleton_vsearch_cdb.fasta','a')
			outblast = open(outputpath+'/VsearchBLAST_4sorted.tsv','a')
			if int(seq.description.split('=')[1].replace('r','')) > (int(readcutoff)-1):
				outseq.write('>'+seq.description+ '_No_BLASTrecord\n'+str(seq.seq) + '\n')
				outseq.close()
			outblast.write(seq.description+'\tNO_BLAST\t\t\t\t\t\t\t\t\t\t\n')
			outblast.close()
	print("CONGRATULATION the taxonomic assignment is done!")
def main():
	script,  NGSfile, idminy, qcovz, Taxa, readcutoff = argv
	getBLAST(NGSfile, float(idminy),float(qcovz), Taxa, readcutoff)
main()
