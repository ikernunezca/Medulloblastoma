# This converts the second column of the
# file multi.layer.net (provided by CURIE) to Entrez IDs.
# When the conversion is not possible, the original label is kept.
# A vocabulary between Entrez IDs and original labels is also generated.

import pandas as pd
import mygene
import shutil
import urllib.parse
import urllib.request
from urllib.request import urlopen
from contextlib import closing
import csv
from tqdm import tqdm
import re, os
import gzip
import os.path
import requests
import io
import subprocess

#load the file
az_df = pd.read_csv('multi.layer.net',sep=' ')
original = az_df.copy()

#create a vocabulay of the original entities before alterating them
az_dict = dict(enumerate(az_df['to']))

#remove redundant labels in column 'to'

#it semms that labels of the Protein layer are UniProt IDs devoid of the '_HUMAN' sufix
az_df['to'] = az_df['to'].str.replace(".Protein", "_HUMAN.Protein").tolist()

def clean(text):
	for ch in ['.Methylation','.Phospho','.Protein','.RNAseq']:
		if ch in text:
			text = text.replace(ch,'')
	return(text)
az_df['to'] = [clean(i) for i in az_df['to'].tolist()]

def extractUniprot(text):
	m = re.search('(^.*?_HUMAN)', text)
	if m:
		found = m.group(1)
	else:
		found = text
	return(found)
az_df['to'] = [extractUniprot(i) for i in az_df['to'].tolist()]

uniprot = az_df[az_df['to'].str.contains("_HUMAN")]['to']
genename = az_df[~az_df['to'].str.contains("_HUMAN")]['to']

#convert UniProt ids to Entrez ids

def myuniprot(lst): # a list of uniprot ids
	url = 'https://www.uniprot.org/uploadlists/'

	params = {
		'from': 'ACC+ID',
		'to': 'P_ENTREZGENEID',
		'format': 'tab',
		'query': lst
	}

	data = urllib.parse.urlencode(params)
	data = data.encode('utf-8')
	req = urllib.request.Request(url, data)
	with urllib.request.urlopen(req) as f:
		response = f.read()
	out = response
	out = [ i.split('\\t') for i in str(out).split('\\n') ]
	out = list(filter (lambda x: len(x) == 2 , out))
	up_dict = {d[0]: d[1] for d in out}
	return(up_dict)

print('converting from uniprot to entrez (1/2)')
sect = uniprot
lst = ' '.join(list(set(sect)))
d = myuniprot(lst)
sect = sect.map(d).fillna(sect)
sect = dict(zip(sect.index,sect))
az_dict.update(sect)
uniprot = uniprot[~uniprot.isin(d.keys())]

#search the remaining ones directly
def fromUniprot(j):
	url = 'https://www.uniprot.org/uniprot/'+j+'.txt'
	urlData = urllib.request.urlopen(url)
	lst = str(urlData.read()).split('\\n')
	entrez = re.compile('DR\s+GeneID;\s+(\d+);')
	for i in lst:
		if(entrez.findall(i)):
			return(entrez.findall(i)[0])

print('converting from uniprot to entrez (2/2)')
sect = uniprot
lst = list(set(sect))
d = {}
for j in tqdm(lst):
	d[j] = fromUniprot(j)
sect = sect.map(d).fillna(sect)
sect = dict(zip(sect.index,sect))
az_dict.update(sect)

#convert gene symbols to Entrez ids

#download gene2accession conversion table from NCBI and extract Homo sapiens (if not there already)
pipe = subprocess.Popen(["perl", "extract_Hsap.pl"], stdin=subprocess.PIPE)
pipe.stdin.close()
file = 'gene2accession_Hsap'
df = pd.read_csv(file, sep='\t', dtype=str)

print('converting from symbols to entrez (1/3)')
sect = genename
d = {}
for i in ['RNA_nucleotide_accession.version','protein_accession.version','genomic_nucleotide_accession.version']:
	b = dict(zip(df[i],df['GeneID']))
	d.update(b)
sect = sect.map(d).fillna(sect)
sect = dict(zip(sect.index,sect))
az_dict.update(sect)
genename = genename[~genename.isin(d.keys())]

#download Entrez conversion table and load it to a dataframe
file = 'Homo_sapiens.gene_info.gz'
fileout = os.path.splitext(file)[0]+'.txt'
url = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/'+file
with closing(urllib.request.urlopen(url)) as r:
	with open(file, 'wb') as f:
		shutil.copyfileobj(r, f)
with gzip.open(file, 'rb') as f_in:
	with open(fileout, 'wb') as f_out:
		shutil.copyfileobj(f_in, f_out)
os.remove(file)
df = pd.read_csv(fileout,sep='\t', dtype=str)
df = df[df['#tax_id']=='9606']
df = df.iloc[:,[1,2,4]] # GeneID, Symbol, Synonyms
os.remove(fileout)

print('converting from symbols to entrez (2/3)')
sect = genename
lst = set(sect)
d = {}
b = dict(zip(df.Symbol,df.GeneID))
for k,v in b.items():
	if k in lst:
		d[k] = v
b = dict(zip(df.Synonyms,df.GeneID))
for k,v in b.items():
	if v != '\-':
		for i in k.split('|'):
			if i in lst:
				d[i] = v
sect = sect.map(d).fillna(sect)
sect = dict(zip(sect.index,sect))
az_dict.update(sect)
genename = genename[~genename.isin(d.keys())]

#check the remaining ones in mygene
mg = mygene.MyGeneInfo()
def mygene(lst): # a list of gene symbols
	mg_dict = {}
	out = mg.querymany(lst, scopes='symbol', species=9606)
	for i in out:
		try:
			mg_dict[i['query']] = i['entrezgene']
		except:
			pass
	return(mg_dict)

print('converting from symbols to entrez (3/3)')
sect = genename
lst = set(sect)
d = mygene(lst)
sect = sect.map(d).fillna(sect)
sect = dict(zip(sect.index,sect))
az_dict.update(sect)

#save the updated edgelist and the vocabulary
print('creating output edgelist')
lst = []
for k in tqdm(list(az_dict.keys())):
	lst.append([az_df.loc[k]['from'],az_dict[k],az_df.loc[k]['layer'],original.loc[k]['to']])

df = pd.DataFrame.from_records(lst)
df.columns = ['from','to','layer','to_original']
#df.to_csv('multi.layer.net.check',sep='\t',index=False,header=False)

df[['from','to','layer']].drop_duplicates().to_csv('multi.layer.net.gr',sep='\t',index=False,header=False)
df[['to','to_original']].drop_duplicates().to_csv('multi.layer.net.vocab',sep='\t',index=False,header=False)
