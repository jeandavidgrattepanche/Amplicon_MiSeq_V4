#!/usr/bin/python3
# this script will assemble your MiSeq PE read using bbmap and vsearch toolkits
# python3 MiSeq_pipeline_V4_SWARM_part1HPC.py ExampleFiles/RawData/
# then reply to prompts (read guide before)

#### TO DO BEFORE RUNNING THE SCRIPT ###

# update first few letter of file (RWS00 in my case) L64 and L96
# update the number of samples to L61 (range of sample from 5 to 14 [in my case RWS0005 to RWS0014] should be range(5, 15))
# update L74 and L106, update the values of bbmerge to match the stringency you want. 
 
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
	folderraw = sys.argv[1]
	pathA = os.getcwd()+"/"+ folderraw #where your folder are located
	bbmappath = os.path.abspath(os.path.join(os.getcwd(), '..')) +"/software/bbmap/" #where bbmap is installed
	os.path.abspath(os.path.join(os.getcwd(), '..'))
	outputpath = os.getcwd()+"/"+folderraw.split('/')[0]+ '/outputs/'
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
	resultfile = open("SWARM_sample.txt","w+") 
	resultfile.write("Sample\treads\tcleanreads\tuniquereads\tSWARM\tSWARM10\tSWARM100\n")
	resultfile.close()	
	os.system("module load java")		
	for i in range(1,97):
		print(i)
		if i < 10:
			sample="RWS000"+str(i) #RWS should match the beginning of your sample or update the line
			print(sample)
			for rawfile in os.listdir(pathA):
# 				print(rawfile,pathA)
				num_SWARM=0;num_reads=0;numUreads=0;numrawreads=0;SWARMnr=0;SWARMnr2=0
				if rawfile.startswith(sample+'_S') and rawfile.endswith(".fastq.gz"):
					if "R1" in rawfile:
						FWD_reads=pathA+rawfile
						REV_reads=pathA+rawfile.replace("R1","R2")
#						print(FWD_reads,REV_reads)
						print(bbmappath+"bbmerge-auto.sh in1="+FWD_reads+" in2="+REV_reads+" out="+mergepath+sample+"_merge.fastq outu="+unmergepath+sample+"_unmerge.fastq  ihist="+histpath+sample+"_ihist.txt ecct extend2=150 loose iterations=5")
						os.system(bbmappath+"bbmerge-auto.sh ecct extend2=150 loose iterations=5 in1="+FWD_reads+" in2="+REV_reads+" out="+mergepath+sample+"_merge.fastq")# outu="+unmergepath+sample+"_unmerge.fastq ihist="+histpath+sample+"_ihist.txt")
						# ecct = error correction by Kmer, extend2 = length to add after failed merging, iterations = number of failed allowed. loose = strictness ( from strict to ultraloose and fast)
						for linec in open(mergepath+sample+"_merge.fastq",'r'):
							if linec.startswith('@'):
								numrawreads += 1
						os.system("vsearch --derep_fulllength "+mergepath+sample+"_merge.fastq --sizeout --fasta_width 0 --output "+derepApath+sample+"_derepA.fasta")
						os.system("python script/ext_remove_N_in_seqfile_v2.py "+derepApath+sample+"_derepA.fasta")
						os.system("vsearch --derep_fulllength "+derepApath+sample+"_derepA_noN.fas --sizein --sizeout --fasta_width 0 --output "+derepBpath+sample+"_derepB.fasta")
						os.system("sh "+bbmappath+"bbduk.sh in="+derepBpath+sample+"_derepB.fasta out="+Qlenpath+sample+"_Qlen.fasta minlen=400")
						os.system("swarm -s "+statSWARMpath+sample +".Stat -d 1 -z "+Qlenpath+sample+"_Qlen.fasta > "+SWARMpath+sample+".swarm")
						num_SWARM = sum(1 for linev in open(statSWARMpath+sample +".Stat"))
						for line2 in open(statSWARMpath+sample +".Stat",'r'):
							num_reads += int(line2.split('\t')[1])
							numUreads += int(line2.split('\t')[0])
							if int(line2.split('\t')[1]) > 9:
								SWARMnr += 1
							if int(line2.split('\t')[1]) > 99:
								SWARMnr2 += 1
						resultfile = open("SWARM_sample.txt","a") 
						resultfile.write(sample+'\t'+ str(numrawreads)+'\t'+str(num_reads)+'\t'+str(numUreads)+'\t'+str(num_SWARM)+'\t'+str(SWARMnr)+'\t'+str(SWARMnr2)+'\n')
						resultfile.close()
		else:
			sample="RWS00"+str(i) #RWS should match the beginning of your sample or update the line
			print(sample)
			for rawfile in os.listdir(pathA):
				print(rawfile,pathA)
				num_SWARM=0;num_reads=0;numUreads=0;numrawreads=0;SWARMnr=0;SWARMnr2=0
				if rawfile.startswith(sample+'_S') and rawfile.endswith(".fastq.gz"):
					if "R1" in rawfile:
						FWD_reads=pathA+rawfile
						REV_reads=pathA+rawfile.replace("R1","R2")
						print(FWD_reads,REV_reads)
						os.system(bbmappath+"bbmerge-auto.sh in1="+FWD_reads+" in2="+REV_reads+" out="+mergepath+sample+"_merge.fastq outu="+unmergepath+sample+"_unmerge.fastq  ihist="+histpath+sample+"_ihist.txt ecct extend2=150 loose iterations=5")
						for linec in open(mergepath+sample+"_merge.fastq",'r'):
							if linec.startswith('@'):
								numrawreads += 1
						os.system("vsearch --derep_fulllength "+mergepath+sample+"_merge.fastq --sizeout --fasta_width 0 --output "+derepApath+sample+"_derepA.fasta")
						os.system("python script/ext_remove_N_in_seqfile_v2.py "+derepApath+sample+"_derepA.fasta")
						os.system("vsearch --derep_fulllength "+derepApath+sample+"_derepA_noN.fas --sizein --sizeout --fasta_width 0 --output "+derepBpath+sample+"_derepB.fasta")
						os.system("sh "+bbmappath+"bbduk.sh in="+derepBpath+sample+"_derepB.fasta out="+Qlenpath+sample+"_Qlen.fasta minlen=400")
						os.system("swarm -s "+statSWARMpath+sample +".Stat -d 1 -z "+Qlenpath+sample+"_Qlen.fasta > "+SWARMpath+sample+".swarm")
						num_SWARM = sum(1 for linev in open(statSWARMpath+sample +".Stat"))
						for line2 in open(statSWARMpath+sample +".Stat",'r'):
							num_reads += int(line2.split('\t')[1])
							numUreads += int(line2.split('\t')[0])
							if int(line2.split('\t')[1]) > 9:
								SWARMnr += 1
							if int(line2.split('\t')[1]) > 99:
								SWARMnr2 += 1
						resultfile = open("SWARM_sample.txt","a") 
						resultfile.write(sample+'\t'+ str(numrawreads)+'\t'+str(num_reads)+'\t'+str(numUreads)+'\t'+str(num_SWARM)+'\t'+str(SWARMnr)+'\t'+str(SWARMnr2)+'\n')
						resultfile.close()
main()
