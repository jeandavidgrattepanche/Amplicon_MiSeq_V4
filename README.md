# Amplicon_MiSeq_V4
updated version of my pipeline using the v4 primers for eukaryotic SSU and running on a HPC linux server.


# Amplicon_MiSeq_pipeline
This repository contains scripts and use PR2 database to analyse MiSeq data for Amplicon generate with eukaryotic V4 SSU primers.

# Pipeline Guide for eukV4SSU

A whole pipeline (scripts and folders structure) for MiSeq analysis.

Prepare your data and folders:

0- keep the same folders and files structure from the repository or the scripts will crash!

1- Create a folder named Rawdata with all your MiSeq sequence files (e.g. HTS1_test1_S1_L001_R1.fastq.gz, HTS1_test1_S1_L001_R2.fastq.gz)

	You can use the script movefile.py to create this folder

2- Create a file with your sample code and sample name (there should be a file named List_samples.txt containing: RWS## (tab) samplename ) \n '

	You can use excel to create this file and save as a Tab Delimited Text (.txt) file
	
	example of list: 

		HTS1	test1
	
		HTS2	test2
	
		HTS3	test3

Note: the first element is a part of name of the file generated from HTS and the second is the name you want the sample to be labelled
	
3- Copy the script folder, the 3 scripts named MiSeq_pipeline_V4_SWARM_part(1,2 and 3)_HPC.py from this repository and add all in the folder where you save the List_samples.txt and the rawdata folder (suggestion MiSeq_folder).

4- primer file should be in the db folder (with the PR2 UTAX.fasta). Add a fasta file with your primers such as: 

	'>'V4_f
	
	CCAGCASCYGCGGTAATTCC
	
	'>'V4.1_f
	
	CCAGCAGCCGCGGTAATTCC
	
	'>'V4.2_f
	
	CCAGCAGCTGCGGTAATTCC
	
	'>'V4.3_f
	
	CCAGCACCCGCGGTAATTCC
	
	'>'V4.4_f
	
	CCAGCACCTGCGGTAATTCC
	
	'>'V4_r
	
	ACTTTCGTTCTTGATYRA
	
	'>'V4.1_r
	
	ACTTTCGTTCTTGATCAA
	
	'>'V4.2_r
	
	ACTTTCGTTCTTGATTAA

	'>'V4.3_r
	
	ACTTTCGTTCTTGATCGA
	
	'>'V4.4_r
	
	ACTTTCGTTCTTGATTGA
	
	

Note: if you have degenerated bases, you should add the various combination as the script do not do it for you. (remove the ' to create your file)
**IMPORTANT: the name of the primer must be as following: text_r or text_f. r and f for reverse and forward and text should not contain any underscores.**


More descriptions are available in Guide_MiSeqPipeline_2018.txt (not up to date yet).

An undergraduate proofread guide can be shared on request.

# Amplicon_MiSeq_pipeline for other lineages
This set of scripts can be edited to use another set of primers/taxa.

part 1: Do not forget to update PEAR parameters in MiSeq_pipeline_V4_SWARM_part1HPC.py (see here for parameters https://sco.h-its.org/exelixis/web/software/pear/doc.html)

part2: build your own database folder (see https://github.com/jeandavidgrattepanche/SSU_DataBase_builiding for scripts and more information) and replacing the database in the script MiSeq_pipeline_V4_SWARM_part2HPC.py and script 6 (look for SSU_db and replace by the corresponding files). 

Where to update your reference database:
- MiSeq_pipeline_V4_SWARM_part2HPC.py  

	*71 and 78: replace file after '--mapout' by your reference sequence alignment
	
	*72 and 79: replace the last argument by your list of column with missing data
	
	*76 and 83/84: replace the "-t" argument by your reference tree
	
- script 6: line 15 

	line 15: replace the SSU_db value by your db e.g. the PR2database (do not forget to add the database folder)


Create also your primer file. **IMPORTANT: the name of the primer must be as following: text_r or text_f. r and f for reverse and forward and text should not contain any underscores.** see above for more details


