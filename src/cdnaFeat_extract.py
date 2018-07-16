#!/usr/bin/env python

# Take Biomart File for getiing cDna Features Table

# Written by Arnold Franz AKE
#aerod7710@gmail.com in April 2018
#@UMR7216 Paris Diderot

# -*- coding: utf-8 *-*-

#Loading packages
from Bio import SeqIO
import pandas as pd
import numpy as np
import sys
import argparse

'''
#Options/Arguments parser
parser = argparse.ArgumentParser(description = "Getting Features from a cDNAtable")
#Positionnal Arguments
parser.add_argument('infile', metavar = "input_file",
 type = argparse.FileType('r'), help = 'Path (name) of the cDna_table.')
parser.add_argument('outfile', nargs = "?", default = sys.stdout, metavar = "output_file",
 type = argparse.FileType('w'), help = 'Path (name) of the cdnaFeaturesTable. Default: writing to the Stdout')


args = parser.parse_args()
'''

list_param = []
list_seq = []
with open("../data/BertrandandRandom_Mart_Export.txt", "rU") as handle:
	for record in SeqIO.parse(handle, "fasta"):
		if(len(record.id.split("|"))==13):
			list_param.append(str(record.id))
			list_seq.append(str(record.seq))
		

cdna_seq = pd.Series(list_seq)
ensembl_gene_id = []
ensembl_transcript_id = []
external_gene_name = []
transcript_start = []
transcript_end = []
utr5_start = []
utr5_end = []
utr3_start = []
utr3_end = []
transcription_start_site = []
transcript_biotype = []
cdna_coding_start = []
cdna_coding_end = []


for elem in list_param:
	elem=elem.split("|")
	ensembl_gene_id.append(elem[0])
	ensembl_transcript_id.append(elem[1])
	external_gene_name.append(elem[2])
	transcript_start.append(elem[3])
	transcript_end.append(elem[4])
	utr5_start.append(elem[5])
	utr5_end.append(elem[6])
	utr3_start.append(elem[7])
	utr3_end.append(elem[8])
	transcription_start_site.append(elem[9])
	transcript_biotype.append(elem[10])
	cdna_coding_start.append(elem[11])
	cdna_coding_end.append(elem[12])


cdna_feat = pd.DataFrame(cdna_seq)
cdna_feat[1] = pd.Series(ensembl_gene_id)
cdna_feat[2] = pd.Series(ensembl_transcript_id)
cdna_feat[3] = pd.Series(external_gene_name)
cdna_feat[4] = pd.Series(transcript_start)
cdna_feat[5] = pd.Series(transcript_end)
cdna_feat[6] = pd.Series(utr5_start)
cdna_feat[7] = pd.Series(utr5_end)
cdna_feat[8] = pd.Series(utr3_start)
cdna_feat[9] = pd.Series(utr3_end)
cdna_feat[10] = pd.Series(transcription_start_site)
cdna_feat[11] = pd.Series(transcript_biotype)
cdna_feat[12] = pd.Series(cdna_coding_start)
cdna_feat[13] = pd.Series(cdna_coding_end)
cdna_feat.columns=["cdna_seq", "ensembl_gene_id", "ensembl_transcript_id", "external_gene_name", "transcript_start",
 "transcript_end", "utr5_start", "utr5_end", "utr3_start", "utr3_end", "transcription_start_site", "transcript_biotype",
  "cdna_coding_start", "cdna_coding_end"]

print cdna_feat.shape

def get_utr5MAX(cdna_feat_row):
	print("getting UtrMax...")
	#take the cdna_feat-row with multiples utrs

	utr5s_start = cdna_feat_row["utr5_start"].values[0].split(";")
	utr5s_end = cdna_feat_row["utr5_end"].values[0].split(";")

	size_liste = []
	for i in range(len(utr5s_start)):
		size = int(utr5s_end[i])-int(utr5s_start[i])
		size_liste.append(size)
	indice_max = size_liste.index(max(size_liste))
	max_utr5_start = int(utr5s_start[indice_max])
	max_utr5_end = int(utr5s_end[indice_max])

	return([max_utr5_start,max_utr5_end])


def get_utr3MAX(cdna_feat_row):
	print("getting UtrMax...")
	#take the cdna_feat-row with multiples utrs

	utr3s_start = cdna_feat_row["utr3_start"].values[0].split(";")
	utr3s_end = cdna_feat_row["utr3_end"].values[0].split(";")

	size_liste = []
	for i in range(len(utr3s_start)):
		size = int(utr3s_end[i])-int(utr3s_start[i])
		size_liste.append(size)
	indice_max = size_liste.index(max(size_liste))
	max_utr3_start = int(utr3s_start[indice_max])
	max_utr3_end = int(utr3s_end[indice_max])

	return([max_utr3_start,max_utr3_end])

def get_cdnaMAX(cdna_feat_row):
	print("getting cDNA_MAX...")
	#take the cdna_feat_row with multiples cdnas

	cdnas_start = cdna_feat_row["cdna_coding_start"].values[0].split(";")
	cdnas_end = cdna_feat_row["cdna_coding_end"].values[0].split(";")

	size_liste = []
	for i in range(len(cdnas_start)):
		size = int(cdnas_end[i])-int(cdnas_start[i])
		size_liste.append(size)
	indice_max = size_liste.index(max(size_liste))
	max_cdna_start = int(cdnas_start[indice_max])
	max_cdna_end = int(cdnas_end[indice_max])

	return([max_cdna_start,max_cdna_end])


#getting the Longest UTR & CDNAs
for i in range(cdna_feat.shape[0]):
	indices=0
	if(str(cdna_feat.iloc[i:i+1,:]["utr5_start"].values[0]) != "nan" and 
		str(cdna_feat.iloc[i:i+1,:]["utr5_start"].values[0]) != "NaN"):
		if(len(cdna_feat.iloc[i:i+1,:]["utr5_start"].values[0].split(";")) > 1):
			indices = get_utr5MAX(cdna_feat.iloc[i:i+1,:])
			cdna_feat.iloc[i:i+1,6:7] = indices[0]
			cdna_feat.iloc[i:i+1,7:8] = indices[1]

	if(str(cdna_feat.iloc[i:i+1,:]["utr3_start"].values[0]) != "nan" and 
		str(cdna_feat.iloc[i:i+1,:]["utr3_start"].values[0]) != "NaN"):
		if(len(cdna_feat.iloc[i:i+1,:]["utr3_start"].values[0].split(";")) > 1):
			indices = get_utr3MAX(cdna_feat.iloc[i:i+1,:])
			cdna_feat.iloc[i:i+1,8:9] = indices[0]
			cdna_feat.iloc[i:i+1,9:10] = indices[1]

	if(str(cdna_feat.iloc[i:i+1,:]["cdna_coding_start"].values[0]) != "nan" and 
		str(cdna_feat.iloc[i:i+1,:]["cdna_coding_start"].values[0]) != "NaN"):
		if( len(cdna_feat.iloc[i:i+1,:]["cdna_coding_start"].values[0].split(";")) > 1):
			indices = get_cdnaMAX(cdna_feat.iloc[i:i+1,:])
			cdna_feat.iloc[i:i+1,12:13] = indices[0]
			cdna_feat.iloc[i:i+1,13:14] = indices[1]


	
#Remplacer les Champs vides par Nan et Delete it
cdna_feat = cdna_feat.replace('',np.NaN)
cdna_feat = cdna_feat.dropna()


cdna_feat = cdna_feat.reset_index()



tmp =cdna_feat.loc[:,['transcript_start', 'transcript_end']].astype(int)

trlength = pd.DataFrame(tmp.transcript_end-tmp.transcript_start,
 columns=["trlength"])

tmp2 = pd.DataFrame()
tmp2[0] = cdna_feat['ensembl_gene_id']
tmp2[1] = trlength
tmp2.columns= ['ensembl_gene_id', 'trlength']


new_tab = pd.DataFrame()
for geneid in np.unique(tmp2.ensembl_gene_id):
	ok1 = tmp2[tmp2.ensembl_gene_id==geneid]
	if(ok1.shape[0]==1):
		new_tab=pd.concat([new_tab,ok1])
		tmp2 = tmp2.drop(ok1.index)
	elif(ok1.shape[0]>1):
		ok2 = ok1[ok1.trlength==max(ok1.trlength)]
		if(ok2.shape[0]==1):
			new_tab=pd.concat([new_tab,ok2])
		else:
			new_tab=pd.concat([new_tab,ok2.iloc[0:1,:]])
		tmp2 = tmp2.drop(ok1.index)

	

cdna_feat = cdna_feat.iloc[new_tab.index,:].reset_index().iloc[:,2:]

'''
with open("bioMartGRCh38_92_RBPx_Random.fasta","w") as handle2:
	for i in range(cdna_feat.shape[0]):	
		handle2.write(">%s|%s\n%s\n" % (cdna_feat['ensembl_gene_id'][i],
			cdna_feat['ensembl_transcript_id'][i],cdna_feat['cdna_seq'][i]))
'''
cdna_feat.to_csv("BertrandandRandom_Mart_Export_EXTRACT.txt",
 sep= "\t", index=False)
