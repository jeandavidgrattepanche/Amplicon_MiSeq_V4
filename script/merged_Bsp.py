import sys
import os
import re
import time
import string
import os.path
from sys import argv
import numpy

#occurence = 13 (only number to sum after)

def Merged_taxa(OTUfile, taxarank):
	with open(OTUfile) as f:
		header = f.readline()
	print(header)
	i = 0; listunique = []; dicttaxa = {}; dictname = {}
	for elt in header.split('\t'):
		i += 1
		if elt == taxarank:	
			taxarp = i - 1
			print(taxarp, taxarank, header.split('\t')[taxarp])
	with open(OTUfile) as f:
		next(f)
		for elt1 in f:
# 			print(elt1.split('\t')[taxarp])
# 			print(len(elt1.split('\n')[0].split('\t')[13:]))
			dictname[elt1.split('\t')[taxarp]] = ('\t').join(elt1.split('\n')[0].split('\t')[1:13])
			if not elt1.split('\t')[taxarp] in dicttaxa:
				dicttaxa.setdefault(elt1.split('\t')[taxarp], [])
				dicttaxa[elt1.split('\t')[taxarp]].append([int(t) for t in elt1.split('\n')[0].split('\t')[14:]])
				listunique.append(elt1.split('\t')[taxarp])
			else:
				dicttaxa[elt1.split('\t')[taxarp]].append([int(t) for t in elt1.split('\n')[0].split('\t')[14:]])
	print(len(listunique), listunique)
	out = open("SPtable.txt", "w+")
	out.write(('\t').join(header.split('\n')[0].split('\t')[1:])+'\n')
	list3 = []
	for elt2 in listunique:
#		print(dicttaxa[listunique[0]])
		my_list = dicttaxa[elt2]
#		print(my_list)
		my_list = sum(map(numpy.array, my_list))
		list3 = my_list.tolist()
# 		for elt3 in my_list:
# 			list3.append(elt3)
		out = open("SPtable.txt","a")
		out.write(elt2+'\t'+dictname[elt2]+'\t'+ str(list3).replace(',','\t').replace('[','').replace(']','')+'\n')
		out.close()
def main():
	script, OTUfile, taxarank = argv
	Merged_taxa(OTUfile, taxarank)
main()