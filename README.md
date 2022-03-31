# Amplicon_MiSeq_V4
updated version of my pipeline using the v4 primers for eukaryotic SSU and running on a HPC linux server.


# Amplicon_MiSeq_pipeline
This repository contains scripts and use PR2 database to analyse MiSeq data for Amplicon generate with eukaryotic V4 SSU primers.

# Pipeline Guide for eukV4SSU

A whole pipeline (scripts and folders structure) for MiSeq analysis.

Prepare your data and folders (see **ExampleFile**):

0- keep the same folders and files structure from the repository or the scripts will crash!

1- Create a folder named Rawdata with all your MiSeq sequence files (e.g. HTS1_test1_S1_L001_R1.fastq.gz, HTS1_test1_S1_L001_R2.fastq.gz)


You can use the script movefile.py to create this folder
	
=> check ExampleFile for an example

2- Create a file with your sample code and sample name (there should be a file named List_samples.txt containing: RWS## (tab) samplename ) \n '

You can use excel to create this file and save as a Tab Delimited Text (.txt) file

example of list: 

		HTS1	test1	
		HTS2	test2	
		HTS3	test3

Note: the first element is a part of name of the file generated from HTS and the second is the name you want the sample to be labelled => check ExampleFile for an example
	
3- Copy the script folder, the 3 scripts named MiSeq_pipeline_V4_SWARM_part(1,2 and 3)_HPC.py from this repository and add all in the folder where you save the List_samples.txt and the rawdata folder (suggestion MiSeq_folder).

(facultative) 4- primer file should be in the db folder (with the PR2 UTAX.fasta). Add a fasta file with your primers such as: 

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

# Run on Temple server
Check the ToInstall.txt which describes the softwares required by the pipeline and the best way to install them on a linux server.

```
$ module load python/3.x.x (python you used to install biopython and other softwares)
$ module load java
$ python3 MiSeq_pipeline_V4_SWARM_part1.py ExampleFile/RawData/
```

Results for the 4 samples attached to this pipeline:
SWARM_sample.txt 



|  Sample |  reads | cleanreads | uniquereads | SWARM | SWARM10 | SWARM100 |
|---------:|--------:|------------:|-------------:|-------:|---------:|----------:|
| RWS0001 | 105398 |   86663    |       65786 | 42301 |     191 |       44 | 
| RWS0002 | 134835 |   49025    |       43050 | 32866 |      93 |       13 | 
| RWS0003 |  99197 |   77539    |       67826 | 51393 |     184 |       31 | 
| RWS0004 |  87184 |   72912    |       59482 | 41748 |     178 |       40 | 


```
$ python3 MiSeq_pipeline_V4_SWARM_part2.py ExampleFile/List_samples.txt ExampleFile/RawData/
```
**working in the last part of the pipeline to reomve outgroup based on tree or to skip this step and produce the final table.**


