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


parser = argparse.ArgumentParser(prog='fasta2table', description='Calculate features from a fasta file (ENSEMBL header) and output a featurs table for clustering', epilog="Authors: Arnold-Franz AKE, Antoine LU & Costas Bouyioukos, UMR7216, Jul-Oct 2018")

parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=__version__))
parser.add_argument("infile", nargs='?', default='-', type=argparse.FileType('r'), metavar="input_file", help="Path to input FASTA file. (or STDIN).")
parser.add_argument("outfile", nargs='?', default='-', type=argparse.FileType('w'), metavar="output_file", help="Path to output CSV filename. (or STDOUT).")
parser.add_argument('-e', '--expressed-transcripts', help="An expressed transcripts file. It can contain an arbitrary number of columns but the first MUST be the transcript IDs. (Default=None).", type=argparse.FileType('r'), default=None, dest="exprTrans")
parser.add_argument('-l', '--length-3pUTR', help="The maximum allowed length of a 3'UTR. (Default=8000)", type=int, default=10000, dest="utr3len")
parser.add_argument('-u', '--utr-files', nargs=2, help="Return two files containing the 5' and 3' UTRs. (Default=None)", type=str, dest="utrFiles")
#TODO add FIMO MEME motifs. parser.add_argument("motifs_file", help="MEME motifs file", default="", type=str)

args = parser.parse_args()

# Populate the Seq.Record generator.
seqRecs = SeqIO.parse(args.infile, "fasta")

# Instantiate the ensebl class.
ensRecs = rnalib.ENSEMBLSeqs(seqRecs, args.exprTrans).bioSeqRecs

# Extract the features from ENSEMBL.
ensFeat = rnalib.FeaturesExtract(ensRecs, args.utr3len)
de = ensFeat.collectFeatures()

# Calculate features by using external programs.
dc = ensFeat.calculateFeatures()

#TODO include all these last steps to a ninalise function in the FeaturesExtract class.
# Concatenate the results
dd = pd.concat([de, dc], axis=1)

# Re-arrange data frame columns.
dd = dd[['ensembl_gene_id', 'gene_name', 'coding_len', '5pUTR_len', '5pUTR_GC', '5pUTR_MFE', '5pUTR_MfeBP', '3pUTR_len', '3pUTR_GC', '3pUTR_MFE', '3pUTR_MfeBP', 'Kozak_Sequence', 'Kozak_Context']]

# Sort and print to csv file.
dd.sort_values(by=['ensembl_gene_id', 'coding_len'])
dd.to_csv(args.outfile, sep=";")
