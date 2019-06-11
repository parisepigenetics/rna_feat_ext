#!/usr/bin/python3
# -*- coding: UTF-8 -*-

"""Calculate features' table from a fasta file with the appropriate header.

Authors: Costas Bouyioukos, AKE Franz-Arnold and LU Antoine
mail: costas.bouyioukso@univ-paris-diderot.fr aerod7710@gmail.com lu.zhao.antoine@gmail.com
2018-19
@UMR7216 Paris Diderot
"""

__version__ = "0.4a01"

import sys
import argparse
import datetime
import pandas as pd
from Bio import SeqIO

import rnaFeaturesLib as rnalib

parser = argparse.ArgumentParser(prog='fasta2table', description="Calculate features from a transcripts fasta file (ENSEMBL header) and return a transcript features' table", epilog="Authors: Costas Bouyioukos, Franz-Arnold Ake and Antoine Lu, 2018-19, Paris UMR7216.")

parser.add_argument('-v', '--version', action='version', version='%(prog)s  v. {version}'.format(version=__version__))
parser.add_argument("infile", nargs='?', default='-', type=argparse.FileType('r'), metavar="input_file", help="Path to input FASTA file. (or STDIN).")
parser.add_argument("outfile", nargs='?', default='-', type=argparse.FileType('w'), metavar="output_file", help="Path to output CSV file. (or STDOUT).")
parser.add_argument('-l', '--length-3pUTR', help="The maximum allowed length of a 3'UTR. (Default=5000).", type=int, default=20000, dest="utr3len", metavar="3'UTRLength")
parser.add_argument('-u', '--utr-files', nargs=2, help="Return two seperate fasta files containing the 5' and 3' UTRs. (Default=None).", type=str, dest="utrFiles", metavar=("5'UTRFile", "3'UTRFile"))
parser.add_argument('-c', '--clip', help="The 5'UTR segment size to calulate the TOP mRNA local score. (Default=20).", type=int, default=20, dest="clip", metavar="5'UTRclip")
# TODO add FIMO MEME motifs. parser.add_argument("motifs_file", help="MEME motifs file", default="", type=str)

# Parse the command line arguments.
optArgs = parser.parse_args()

# Populate the Seq.Record generator.
seqRecs = SeqIO.parse(optArgs.infile, "fasta")

# Instantiate the ensebl class.
ensRecs = rnalib.ENSEMBLSeqs(seqRecs).bioSeqRecs

# Extract the features from ENSEMBL.
ensFeat = rnalib.FeaturesExtract(ensRecs, optArgs)
de = ensFeat.collect_features()

# Calculate features by using external programs.
dc = ensFeat.calculate_features()
# TODO include all these last steps to an ininalise function in the FeaturesExtract class.

# Concatenate the results
dd = pd.concat([de, dc], axis=1, sort=False)

# Re-arrange data frame columns.
dd = dd[['ensembl_gene_id', 'gene_name', 'coding_len', 'GC', '5pUTR_len', '5pUTR_GC', '5pUTR_MFE', '5pUTR_MfeBP', '3pUTR_len', '3pUTR_GC', '3pUTR_MFE', '3pUTR_MfeBP', 'TOP_localScore', 'CAI', 'Kozak_Sequence', 'Kozak_Context']]

# Sort and print to csv file.
dd.sort_values(by=['ensembl_gene_id', 'coding_len'])
dd.to_csv(optArgs.outfile, sep=";")
# Print the command line arguments in the csv file.
optArgs.outfile.write('# {}\n# {}\n'.format(str(sys.argv), datetime.datetime.now().strftime("%d/%m/%Y at %H:%M:%S")))
