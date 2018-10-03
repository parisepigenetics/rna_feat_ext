#!/usr/bin/python3
# -*- coding: UTF-8 -*-

"""Calculate feature's Table from a fasta file with the appropriate header.

Authors: AKE Franz-Arnold - LU Antoine
mail: aerod7710@gmail.com lu.zhao.antoine@gmail.com
July 2018
UMR7216 Paris Diderot"""

__version__ = "0.3a01"

import argparse
import pandas as pd
from Bio import SeqIO

import rnaFeaturesLib as rnalib


parser = argparse.ArgumentParser(prog='fasta2table', description='Calculate features from a fasta file (ENSEML header) and output a featurs table for clustering', epilog="Authors: Arnold-Franz AKE, Antoine LU, 2018")
parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=__version__))
parser.add_argument("infile", help="input FASTA file", type=str)
parser.add_argument("outfile", help="output CSV filename; add .csv", type=str)
parser.add_argument("motifs_file", help="MEME motifs file", default="", type=str)

args = parser.parse_args()

# Populate the Seq.Record generator.
seqRecs = SeqIO.parse(args.infile, "fasta")
# Instantiate the ensebl class.
ensRecs = rnalib.ENSEMBLSeqs(seqRecs).bioSeqRecs
# Extract the features from ENSEMBL.
ensFeat = rnalib.FeaturesExtract(ensRecs)
de = ensFeat.collectFeatures()
# Calculate features by using external programs.
dc = ensFeat.calculateFeatures()
# Concatenate the results
dd = pd.concat([de, dc], axis=1, sort=False)
# Re-arrange the columns.
dd = dd[['ensembl_gene_id', 'gene_name', 'coding_len', '5pUTR_len', '5pUTR_GC', '5pUTR_MFE', '5pUTR_MfeBP', '3pUTR_len', '3pUTR_GC', '3pUTR_MFE', '3pUTR_MfeBP', 'Kozak_Sequence', 'Kozak_Context']]
# Sort and print to csv file.
dd.sort_values(by=['ensembl_gene_id', 'coding_len'])
dd.to_csv(args.outfile, sep=";")