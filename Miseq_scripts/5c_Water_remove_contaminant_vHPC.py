#!/usr/bin/python2

__author__ = "Jean-David Grattepanche"
__version__ = "4.01, June 29, 2021"
__email__ = "jeandavid.grattepanche@gmail.com"
from Bio import AlignIO
import re,os, sys
from Bio import SeqIO
from Bio import Phylo
from sys import argv




def test(seqfile, seqdict):
	outpath = seqfile.split('chimeras/')[0]
	flag = 1
	midline = ""
	gap = 0.0
	ident = 3.0
	sim = 3.0
	diff = 0.0
	outfinal = open( seqfile.split('.')[0]+'_nocont.fasta','a')
	outcont = open(seqfile.split('.')[0]+'_cont.fasta','a')
	outscore = open(seqfile.split('.')[0]+'_pairwise_out_scores.csv','a')
	
	for line in open(outpath + 'chimeras/water.txt','r'):
#		outfinal.write(line)
		if re.search('# 1: ', line):
			seq1 = str(line.split(':')[1].strip())
		if re.search('# 2: ', line):
			seq2 = str(line.split(':')[1].strip())		
		if re.search('# Length:',line):
			tlengthlist = re.findall('\d+', line)
			tlength = float(tlengthlist[0])
			
		#print line[0]
		if line.strip() == '#=======================================':
			flag = 0
		
		if line.strip() == '#---------------------------------------':
			flag = 1
			
		if flag == 0:
			if line[0] == ' ':
				midline=midline + line.strip()
#	outfinal.write('\n')	
	midline = midline + '*'
#	print(midline)
	midline2 = ""
	id = 0 
	gp = 0
	for char in midline:
		if id < 4:	
			if char == '|':
				id = id + 1
			else: 
				id = 0
		else:
			midline2 = midline2 + char
	
	for char in midline2:
	
		if gp < 4:	
			if char == ' ':
#				gp = gp + 1
				gap = gap + 1
			else:
				gp = 0
			if char == '|':
				ident = ident + 1
				sim = sim + 1
			if char == ':':
				sim = sim + 1
			if char == '.':
				diff = diff + 1
			if char == '*':
				gp = 10
	len = gap + sim + diff
	
	outscore.write(seq1 + ',' + seq2 + ', len: ' + str(len) + ', tlength: ' + str(tlength) + ', identity: ' + str(ident) + ', similarity: ' + str(sim) + ', len/tlength: ' + str(len/tlength) + ', sim/len: ' + str(sim/len) + '\n')
	print(seq2, len, round(((float(sim)/float(len))*100),2))
	if float(sim)/float(len) * 100 < float(50):# and sim/ident > 1.02: #len/tlength > 0.67 and sim/len > 0.75 and 
		print("too different", seq1, seq2)
		outcont.write('>'+seq2+'\n'+seqdict[seq2]+'\n')
	else:
		outfinal.write('>'+seq2+'\n'+seqdict[seq2]+'\n')
		
	outfinal.close()		
def main():
	script, seqfile = argv
	outpath = seqfile.split('chimeras/')[0]
	seqdict= {}; seqlist= []; readmax = 0
	outfinal = open( seqfile.split('.')[0]+'_nocont.fasta','w+')
	outcont = open(seqfile.split('.')[0]+'_cont.fasta','w+')
	outscore = open(seqfile.split('.')[0]+'_pairwise_out_scores.csv','w+')
	for Seq in SeqIO.parse(seqfile,'fasta'):
		seqdict[Seq.id]=str(Seq.seq)
		readnumber = Seq.id.split('_')[1].split('r')[0]
		if int(readnumber) > int(readmax):
			readmax = readnumber
	for Seq in SeqIO.parse(seqfile,'fasta'):
		readnumber = Seq.id.split('r')[0].split('_')[1]
		if int(readnumber) == int(readmax):
			refseq=open(outpath+'ref_for_water.fasta','w+')
			refseq.write('>'+Seq.id+'\n'+str(Seq.seq)+'\n')
			refseq.close()
	
	for Seq in SeqIO.parse(seqfile,'fasta'):
		if Seq.id not in seqlist:
			seqlist.append(Seq.id)
			outseq=open(outpath + 'chimeras/seqtemp.fas','w')
			outseq.write('>'+Seq.id+'\n'+str(Seq.seq)+'\n')
			outseq.close()
			cline = os.getcwd()+'/software/EMBOSS-6.6.0/emboss/water -asequence='+ outpath + 'ref_for_water.fasta -bsequence='+ outpath + 'chimeras/seqtemp.fas -gapopen=10 -gapextend=0.5 -outfile='+outpath + 'chimeras/water.txt'		
#			print(cline)
			os.system(cline)
			test(seqfile, seqdict)
main()
