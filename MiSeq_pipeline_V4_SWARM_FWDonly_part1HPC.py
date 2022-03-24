#!/usr/bin/python3
# this script will not assemble your MiSeq PE read using bbmap and used only FWD read and vsearch toolkits

# python3 MiSeq_pipeline_v4_SWARM_part1HPC.py folder
# then reply to prompts (read guide before)

#### TO DO BEFORE RUNNING THE SCRIPT ###
# L62 and 86 + L59 (range if sample from 5 to 14 should be (5, 11))
  
__author__ = "Jean-David Grattepanche"
__version__ = "1.01, June 29, 2021"
__email__ = "jeandavid.grattepanche@gmail.com"

import string
import re
import sys
import os
from Bio import SeqIO
from sys import argv
	
def main():
	folderraw = sys.argv[1]
	pathA = os.getcwd() + folderraw #RWS_0001-0096/RawData
	bbmappath = os.getcwd() +"/software/bbmap/"
	outputpath = os.getcwd() +folderraw.split('/')[0]+ '/outputsFWDonly'
	if not os.path.exists(outputpath): 
		os.makedirs(outputpath) 	
	temppath = os.getcwd() +folderraw.split('/')[0]+  '/temp/'
	if not os.path.exists(temppath): 
		os.makedirs(temppath) 	
	mergepath = outputpath + '/merge/'
	if not os.path.exists(mergepath): 
		os.makedirs(mergepath) 	
	unmergepath = outputpath + '/unmerged/'
	if not os.path.exists(unmergepath): 
		os.makedirs(unmergepath) 	
	histpath = outputpath + '/histo/'
	if not os.path.exists(histpath): 
		os.makedirs(histpath) 	
	derepApath = outputpath + '/derepA/'
	if not os.path.exists(derepApath): 
		os.makedirs(derepApath) 			
	derepBpath = outputpath + '/derepB/'
	if not os.path.exists(derepBpath): 
		os.makedirs(derepBpath) 		
	Qlenpath = outputpath + '/Qlen/'
	if not os.path.exists(Qlenpath): 
		os.makedirs(Qlenpath) 		
	SWARMpath = outputpath + '/SWARM/'
	if not os.path.exists(SWARMpath): 
		os.makedirs(SWARMpath) 		
	statSWARMpath = outputpath + '/statSWARM/'
	if not os.path.exists(statSWARMpath): 
		os.makedirs(statSWARMpath)
	resultfile = open("SWARM_sample.txt","w+") 
	resultfile.write("Sample\treads\tcleanreads\tuniquereads\tSWARM\tSWARM10\tSWARM100\n")
	resultfile.close()	
	os.system("module load java")		
	for i in range(1,97):
		print(i)
		if i < 10:
			sample="RWS000"+str(i)
			print(sample)
			for rawfile in os.listdir(pathA):
				num_SWARM=0;num_reads=0;numUreads=0;numrawreads=0;SWARMnr=0;SWARMnr2=0
				if rawfile.startswith(sample+'_S') and rawfile.endswith(".fastq.gz"):
					if "R1" in rawfile:
						FWD_reads=pathA+rawfile
						os.system("vsearch --derep_fulllength "+FWD_reads+" --sizeout --fasta_width 0 --output "+derepApath+sample+"_derepA.fasta")
						os.system("python3 script/ext_remove_N_in_seqfile_v2.py "+derepApath+sample+"_derepA.fasta")
						os.system("vsearch --derep_fulllength "+derepApath+sample+"_derepA_noN.fas --sizein --sizeout --fasta_width 0 --output "+derepBpath+sample+"_derepB.fasta")
						os.system("sh "+bbmappath+"bbduk.sh in="+derepBpath+sample+"_derepB.fasta out="+Qlenpath+sample+"_Qlen.fasta qtrim=r trimq=20 minlen=200")
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
			sample="RWS00"+str(i)
			print(sample)
			for rawfile in os.listdir(pathA):
				num_SWARM=0;num_reads=0;numUreads=0;numrawreads=0;SWARMnr=0;SWARMnr2=0
				if rawfile.startswith(sample+'_S') and rawfile.endswith(".fastq.gz"):
					if "R1" in rawfile:
						FWD_reads=pathA+rawfile
						os.system("vsearch --derep_fulllength "+FWD_reads+" --sizeout --fasta_width 0 --output "+derepApath+sample+"_derepA.fasta")
						os.system("python3 script/ext_remove_N_in_seqfile_v2.py "+derepApath+sample+"_derepA.fasta")
						os.system("vsearch --derep_fulllength "+derepApath+sample+"_derepA_noN.fas --sizein --sizeout --fasta_width 0 --output "+derepBpath+sample+"_derepB.fasta")
						os.system("sh "+bbmappath+"bbduk.sh in="+derepBpath+sample+"_derepB.fasta out="+Qlenpath+sample+"_Qlen.fasta qtrim=r trimq=20 minlen=200")
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
