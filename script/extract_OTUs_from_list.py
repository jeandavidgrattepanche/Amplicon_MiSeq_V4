import sys
import os
import re
import time
import string
import os.path
from sys import argv
import numpy

#occurence = 13 (only number to sum after)

def Merged_taxa(OTUfile, listtokeep):
	i = 0; listunique = []; listtowrite= []
	with open(listtokeep) as lk:
		for line in lk:
			if line.startswith("OTU"):
				listunique.append(line.split(' ')[0])
	print("number of OTUs to keep: ", len(listunique))
	
	with open(OTUfile) as f:
		header = f.readline()
# 	print(header)
	with open(OTUfile) as f:
		next(f)
		for elt1 in f:
			if elt1.split('\t')[0] in listunique:
				listtowrite.append(elt1)
	print("number of OTUs kept: ", len(listtowrite))
	if len(listunique) == len(listtowrite):
		out = open("indVal_table.txt", "w+")
		out.write(header)
		for elt2 in listtowrite:
			out.write(elt2)
		out.close()
	else:
		print(listunique[1], "\n\n\n", elt1.split('\t')[0])
def main():
	script, OTUfile, listtokeep = argv
	Merged_taxa(OTUfile, listtokeep)
main()