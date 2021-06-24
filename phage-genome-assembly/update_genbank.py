from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import argparse
import sys 
import os 
import itertools
from collections import defaultdict

'''
This code takes the available genbank file and updates the file with blast hits. 
The blast output should be in a tabular format, output format 6. The script updates the
genbank files and adds the annotation only if the gene is annotated as hypothetical 
protein

Running code
python update_genbank.py -g genbank -b blast_output.tsv -o updated_genbank
'''

def prepend(list, str):
	# Using format()
	str += '{0}'
	list = [str.format(i) for i in list]
	return(list)
  
def update_genbank(gbk, blast, gbk_out):
	'''
	This bit  is reading the genbank file and counting the number of hypothetical proteins
	'''
	count=0
	list=[]
	for record in SeqIO.parse(gbk, "genbank"):
		for f in record.features:
			#print (f.qualifiers)
			#print (f.type)
			if f.type =="CDS" and "product" in f.qualifiers:
				pdt=f.qualifiers["product"][0]
				if pdt=="hypothetical protein":
					count=count+1
				pdt=f.qualifiers["product"][0]
				list.append(f.qualifiers["db_xref"])
					
	#print (list)

	'''
	Taking a look at the tabular format from blast output, matching the hypoethical proteins to blast output
	'''
	tsv=open(blast, 'r')
	star="RAST2:"
	rast=[]
	name=[]
	for line in tsv:
		l=line.split(',')
		li=l[1].strip()
		rast.append(li)
		protein=l[13].rstrip("\n")
		name.append(protein)
	new_rast=prepend(rast, star)
	#print (new_rast)

	'''
	Comparing both the genbank and the blastp hit list to grab the annotations that are not novel
	'''
	fixed={}
	for i in list:
		index=0
		for j in new_rast:
			if (i[0]==j):
				fixed[i[0]]=name[index]
				#fixed.append(name[index])
			index=index+1		
	#print (fixed)

	'''
	Adding the new product name to genbank
	'''
	output_file=open(gbk_out, 'w')
	for record in SeqIO.parse(gbk, "genbank"):
		for f in record.features:
			if f.type =="CDS" and "product" in f.qualifiers:
				for k in fixed:
					if (f.qualifiers["db_xref"][0]==k):
						f.qualifiers["note"]="BLASTP"
						f.qualifiers["product"]=fixed[k]
			#print (f)
		SeqIO.write(record, output_file, "genbank")
							
if __name__=='__main__':
	parser=argparse.ArgumentParser(description="Update genbank file")
	parser.add_argument ('-g', dest="genbank", help="Enter the genbank filename")
	parser.add_argument ('-b', dest="blast", help="Enter the blast output filename")
	parser.add_argument ('-o', dest="output", help="Enter the output filename")
	results=parser.parse_args()
	update_genbank(results.genbank, results.blast, results.output)
