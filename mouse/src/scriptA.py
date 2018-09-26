#!/usr/bin/env python
# coding: utf-8

# From a list of Gene_ID ..getting Feature's Table

# Author: AKE Franz-Arnold - Antoine LU
# mail: aerod7710@gmail.com
# 17 Juillet 2018
# UMR7216 Paris Diderot


# Loading Packages
import subprocess
import argparse
import sys
import pandas as pd
import library
from biomart import BiomartServer


# Options/Arguments parser
parser = argparse.ArgumentParser(description='Fetch features for an ensembl gene ID list and output a multi FASTA file',  epilog="Author: Franz-Arnold Ake and Antoine Lu, 2018")
# Positional Arguments
parser.add_argument('infile', metavar="input_file", type=argparse.FileType('r'), help='Path for Ensembl Gene_ID file')
parser.add_argument('outfile', nargs="?", default=sys.stdout, metavar="output_file", type=argparse.FileType('w'), help='Path for saving Output Fasta file (Default: output in Stdout)')

# Optionnal Arguments
parser.add_argument('-d, --dataset', nargs="?", default='mmusculus_gene_ensembl', metavar="Dataset for collecting informations", type=str, help="Choice of Ensembl Dataset -- hsapiens_gene_ensembl Or mmusculus_gene_ensembl -- Default : Mouse_dataset")


args = parser.parse_args()


# Variables
listID = args.infile.read().splitlines()
listAttrib = ['ensembl_gene_id', 'ensembl_transcript_id',
              'external_gene_name', 'transcript_start', 'transcript_end',
              '5_utr_start', '5_utr_end', '3_utr_start', '3_utr_end',
              'transcription_start_site', 'transcript_biotype',
              'cdna_coding_start', 'cdna_coding_end', 'start_position', 'cdna']

# Step_1
print "\nChargement des Données sur Biomart ..."
# Loading Data from Biomart
server = BiomartServer("http://www.ensembl.org/biomart/")
hsapiens_dataset = server.datasets[args.dataset]
#FIXME is it human or mouse...???

response = hsapiens_dataset.search({
    'filters': {listAttrib[0]: listID},
    'attributes': listAttrib}, header=1)
#FIXME convert directly the database responce object to a pandas data table.
#FIXME no need of a TMP file then.
# A améliorer !
################################################################
'''l'idéal serait de ne pas stocker dans un fichier
puis relire pour avoir les datas
mais de le réaliser directement but how ?
'''
filout = open("temp.tsv", "w")
filout.write(response.text)
filout.close()
print "Done"
#FIXME what is the need for a tmp file?

# Step_2
print "\nRecupération des Données..."
# Reading temp_1 for getting preloaded data
with open("temp.tsv", "r") as preloaded_data:
    cdna = pd.read_csv(preloaded_data, sep="\t")
    subprocess.call("rm -f", stdin=preloaded_data, shell=True)
    print "Done"
subprocess.call("rm -f temp.tsv", shell=True)
##################################################################

print("Dimension_Data :", cdna.shape)

print "Elimination des NaN values..."
cdna = cdna.dropna()
cdna = cdna.reset_index()
print "Done"

print("Dimension_Data :", cdna.shape)

# Step3
print "\nElimination des 5utr start et End multiples..."
print "Elimination des Cdna coding start et end multiples..."
for i in range(cdna.shape[0]):
    ligne = pd.DataFrame(cdna.loc[i, :]).transpose()

    # For 5_UTRs
    if len(cdna.iloc[i:i+1, :]["5' UTR start"].values[0].split(";")) > 1:
        indices = library.get_utr5MAX(ligne)
        cdna.iloc[i:i+1, 7:8] = indices[1]
        cdna.iloc[i:i+1, 8:9] = indices[0]

    # For 3_UTRs
    if len(str(cdna.iloc[i:i+1, :]["3' UTR start"].values[0]).split(";")) > 1:
        indices = library.get_utr3MAX(ligne)
        cdna.iloc[i:i+1, 9:10] = indices[1]
        cdna.iloc[i:i+1, 10:11] = indices[0]

    # For cDNA
    if len(cdna.iloc[i:i+1, :]["cDNA coding start"].values[0].split(";")) > 1:
        indices = library.get_cDNAMIN(ligne)
        cdna.iloc[i:i+1, 13:14] = indices[0]
        cdna.iloc[i:i+1, 14:15] = indices[1]
        print "Done"

# Step4
print "\nSauvegarde fasta_file..."
library.txt2fasta(cdna, args.outfile)
print "Done\n"
