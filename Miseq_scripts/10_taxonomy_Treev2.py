#!/usr/bin/python2

__author__ = "Jean-David Grattepanche"
__version__ = "6.01, June 29, 2021"
__email__ = "jeandavid.grattepanche@gmail.com"



import re,os, sys

sys.setrecursionlimit(10**7)

from Bio import SeqIO
from Bio import Phylo
from sys import argv



#take a leaf, compare it to its sister
def get_parent(tree, child_clade, i):
		#############################################################
		#############################################################
	node_path = tree.get_path(child_clade)
#	print(node_path)
	try:
		return node_path[-i]				
	except:
		sys.exit("Error in path")	
	
			
def checkClade(tree, treefile):
	print("start walking in tree")
	seen = []
	outpath = treefile.split('/outgroup_removal/')[0]
	out = open(outpath +'/taxonomic_assignment/taxonomy_by_Tree_r.txt','w+')
	print('ouput clear') ; j = 0
	for clade in tree.get_terminals():
		i = 0 ; j+=1
		print(str(j), "clade")
		if "OTU" in clade.name: 
#			print(clade.name, "checkclade")
			while clade.name not in seen:
				print('start' , i)
				i = i +1
				parent = get_parent(tree,clade, i)
				for leaf in parent.get_terminals():
#					print(leaf.name)
					if 'OTU' not in leaf.name and clade.name not in seen:
						print(clade.name, 'closely related at', str(i),'to ',leaf.name, 'distance=', tree.distance(clade,leaf))
						out = open(outpath +'/taxonomic_assignment/taxonomy_by_Tree_r.txt','a')
						out.write(clade.name+'\t'+ str(i) +'\t'+str(tree.distance(clade,leaf))+ '\t'+leaf.name+'\n')
						seen.append(clade.name)
						out.close()
#							else:
#								print(clade.name, 'nothing close at 2 node aways')
						#parent = get_parent(tree,leaf)
						#print(parent.get_terminals())
def main():
	script, treefile = argv	
# 	if not os.path.exists('outputs/taxonomic_assignment/'): 
# 		os.makedirs('outputs/taxonomic_assignment/') 	


	tree = Phylo.read(treefile,'newick')
	print('tree read')
	checkClade(tree, treefile)	
main()
