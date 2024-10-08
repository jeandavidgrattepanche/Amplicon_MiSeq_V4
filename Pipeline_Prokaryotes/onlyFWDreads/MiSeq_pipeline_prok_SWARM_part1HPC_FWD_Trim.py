#!/usr/bin/python3
# this script will assemble your MiSeq PE read using bbmap and vsearch toolkits
# python3 MiSeq_pipeline_V4_SWARM_part1HPC.py ExampleFiles/List_samples.txt ExampleFiles/RawData/
# then reply to prompts (read guide before)

#### TO DO BEFORE RUNNING THE SCRIPT ###



# update L94 and L120 for length cut-off, update the values of bbmerge (L86) to match the stringency you want. 
 
__author__ = "Jean-David Grattepanche"
__version__ = "2, June 29, 2021"
__email__ = "jeandavid.grattepanche@gmail.com"

import string
import re
import sys
import os
from Bio import SeqIO
from sys import argv
	
def main():
	samplefile = sys.argv[1]
	listsamp = []
	try:
		listsample = samplefile
	except ValueError:
		samplefile = ""	
	if samplefile == "":
		print ('Your input samplefile is empty.  Try again. ')
	for samp in open(listsample,'r'):
		if samp.split('\t')[0] not in listsamp:
			listsamp.append(samp.split('\t')[0])


	folderraw = sys.argv[2]
	pathA = os.getcwd()+"/"+ folderraw #where your folder are located
	bbmappath = os.path.abspath(os.path.join(os.getcwd(), '..')) +"/software/bbmap/" #where bbmap is installed
	os.path.abspath(os.path.join(os.getcwd(), '..'))
	outputpath = os.getcwd()+"/"+folderraw.split('/')[0]+ '/outputsFWDonly/'
	if not os.path.exists(outputpath): 
		os.makedirs(outputpath) 	
	temppath = os.getcwd() +"/"+folderraw.split('/')[0]+  '/temp/'
	if not os.path.exists(temppath): 
		os.makedirs(temppath) 	
	mergepath = outputpath + 'merge/'
	if not os.path.exists(mergepath): 
		os.makedirs(mergepath) 	
	unmergepath = outputpath + 'unmerged/'
	if not os.path.exists(unmergepath): 
		os.makedirs(unmergepath) 	
	histpath = outputpath + 'histo/'
	if not os.path.exists(histpath): 
		os.makedirs(histpath) 	
	derepApath = outputpath + 'derepA/'
	if not os.path.exists(derepApath): 
		os.makedirs(derepApath) 			
	derepBpath = outputpath + 'derepB/'
	if not os.path.exists(derepBpath): 
		os.makedirs(derepBpath) 		
	Qlenpath = outputpath + 'Qlen/'
	if not os.path.exists(Qlenpath): 
		os.makedirs(Qlenpath) 		
	SWARMpath = outputpath + 'SWARM/'
	if not os.path.exists(SWARMpath): 
		os.makedirs(SWARMpath) 		
	statSWARMpath = outputpath + 'statSWARM/'
	if not os.path.exists(statSWARMpath): 
		os.makedirs(statSWARMpath)
	resultfile = open("SWARM_sample_FWDP.txt","w+") 
	resultfile.write("Sample\treads\tcleanreads\tuniquereads\tSWARM\tSWARM10\tSWARM100\n")
	resultfile.close()	
	os.system("module load java")		
	for sample in listsamp:
		print(sample)#+"_S")
		sample = sample+"_"
		for rawfile in os.listdir(pathA):
# 				print(rawfile,pathA)
			num_SWARM=0;num_reads=0;numUreads=0;numrawreads=0;SWARMnr=0;SWARMnr2=0
			if rawfile.startswith(sample) and rawfile.endswith(".fastq.gz"):
				print(rawfile)
				if "R1_" in rawfile:
					FWD_reads=pathA+rawfile
					REV_reads=pathA+rawfile.replace("R1","R2")
#					print(FWD_reads,REV_reads)
					os.system("vsearch --fastx_uniques "+FWD_reads +" --sizeout --fasta_width 0 --fastaout "+derepApath+sample+"_derepAp.fasta")
					for linec in open(derepApath+sample+"_derepAp.fasta",'r'):
						if linec.startswith('>'):
							numrawreads += int(linec.split(";size=")[1])
					os.system("python script/ext_remove_N_in_seqfile_v2.py "+derepApath+sample+"_derepAp.fasta")
					os.system("sh "+bbmappath+"bbduk.sh in="+derepApath+sample+"_derepAp_noN.fas out="+Qlenpath+sample+"_Qlenp.fasta qtrim=r trimq=20 ftr=140 minlen=140")
					##ftr=140 for trimming 140bp length
					os.system("vsearch --derep_fulllength "+Qlenpath+sample+"_Qlenp.fasta --sizein --sizeout --fasta_width 0 --output "+derepBpath+sample+"_derepBp.fasta")
					os.system("swarm -s "+statSWARMpath+sample +".Stat -d 1 -z "+derepBpath+sample+"_derepBp.fasta > "+SWARMpath+sample+".swarm")
					num_SWARM = sum(1 for linev in open(statSWARMpath+sample +".Stat"))
					for line2 in open(statSWARMpath+sample +".Stat",'r'):
						num_reads += int(line2.split('\t')[1])
						numUreads += int(line2.split('\t')[0])
						if int(line2.split('\t')[1]) > 9:
							SWARMnr += 1
						if int(line2.split('\t')[1]) > 99:
							SWARMnr2 += 1
					resultfile = open("SWARM_sample_FWDP.txt","a") 
					resultfile.write(sample+'\t'+ str(numrawreads)+'\t'+str(num_reads)+'\t'+str(numUreads)+'\t'+str(num_SWARM)+'\t'+str(SWARMnr)+'\t'+str(SWARMnr2)+'\n')
					resultfile.close()
				elif "_1.fastq" in rawfile:
					FWD_reads=pathA+rawfile
					REV_reads=pathA+rawfile.replace("_1.fastq","_2.fastq")
# 					print(FWD_reads,REV_reads)
# 					print(bbmappath+"bbmerge-auto.sh in1="+FWD_reads+" in2="+REV_reads+" out="+mergepath+sample+"_merge.fastq outu="+unmergepath+sample+"_unmerge.fastq  ihist="+histpath+sample+"_ihist.txt ecct extend2=150 loose iterations=5")
#					os.system(bbmappath+"bbmerge-auto.sh ecct extend2=150 loose iterations=5 in1="+FWD_reads+" in2="+REV_reads+" out="+mergepath+sample+"_merge.fastq")# outu="+unmergepath+sample+"_unmerge.fastq ihist="+histpath+sample+"_ihist.txt")
					# ecct = error correction by Kmer, extend2 = length to add after failed merging, iterations = number of failed allowed. loose = strictness ( from strict to ultraloose and fast)
					os.system("vsearch --fastx_uniques "+FWD_reads+" --sizeout --fasta_width 0 --fastaout "+derepApath+sample+"_derepAp.fasta")
					for linec in open(derepApath+sample+"_derepAp.fasta",'r'):
						if linec.startswith('>'):
							numrawreads += int(linec.split(";size=")[1])
					os.system("python script/ext_remove_N_in_seqfile_v2.py "+derepApath+sample+"_derepAp.fasta")
					os.system("sh "+bbmappath+"bbduk.sh in="+derepApath+sample+"_derepAp_noN.fas out="+Qlenpath+sample+"_Qlenp.fasta qtrim=r trimq=20 ftr=140 minlen=140")
					##ftr=140 for trimming 140bp length
					os.system("vsearch --derep_fulllength "+Qlenpath+sample+"_Qlenp.fasta --sizein --sizeout --fasta_width 0 --output "+derepBpath+sample+"_derepBp.fasta")
					os.system("swarm -s "+statSWARMpath+sample +".Stat -d 1 -z "+derepBpath+sample+"_derepBp.fasta > "+SWARMpath+sample+".swarm")
					num_SWARM = sum(1 for linev in open(statSWARMpath+sample +".Stat"))
					for line2 in open(statSWARMpath+sample +".Stat",'r'):
						num_reads += int(line2.split('\t')[1])
						numUreads += int(line2.split('\t')[0])
						if int(line2.split('\t')[1]) > 9:
							SWARMnr += 1
						if int(line2.split('\t')[1]) > 99:
							SWARMnr2 += 1
					resultfile = open("SWARM_sample_FWDP.txt","a") 
					resultfile.write(sample+'\t'+ str(numrawreads)+'\t'+str(num_reads)+'\t'+str(numUreads)+'\t'+str(num_SWARM)+'\t'+str(SWARMnr)+'\t'+str(SWARMnr2)+'\n')
					resultfile.close()
				else:
					print(rawfile, "Not a FWD file")
main()

