import string
import sys
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sys import argv
from Bio import AlignIO
from Bio import Align
from Bio.Align import MultipleSeqAlignment

def removecolumnfrommask(seqfile,filetype,maskstart, maskend):
	outFile = open(seqfile.split('.')[0] + '_trimmed.fas','w+')
	alignment = AlignIO.read(seqfile, filetype)
	trimAlign = MultipleSeqAlignment([])
	numCol = alignment.get_alignment_length()
	colToKeep=[]; coltoremove=[]
# 
# 	for k in open(mask,'r'):
# 		coltoremove.append(int(k.split('\n')[0]))	
# 	print(len(coltoremove))

	for i in range(numCol):
		if i >= int(maskstart) and i<= int(maskend):
			colToKeep.append(i)
	print(len(colToKeep))
	for record in alignment:
		newseq = ""
		for j in colToKeep:
			newseq= newseq + (record[j])
			
		newRecord = SeqRecord(Seq(newseq), id=record.id)
		trimAlign.append(newRecord)
		if 'SWARM' in record.id:
			outFile.write('>' + record.id.split('_')[0] + '\n' + newseq + '\n')
		else:
			outFile.write('>' + record.description + '\n' + newseq + '\n')		
	outFile.close()
	print("Total number of columns remaining: %i" % trimAlign.get_alignment_length())


def main():
	print("*************************************************************************************************")
	print("This script will take an alignment and a mask for removing gapped columns." )
	print("Usage is 'python mask_gaps.py <inputAlignment> <typeofalignment>fasta,phylip,nexus <maskst>start and <maskend> end of fragment to keep'")
	print("*************************************************************************************************\n\n")
	
	script, alignedaddfrag, type, maskstart, maskend = argv 
	removecolumnfrommask(alignedaddfrag,type,maskstart, maskend) 
main()
