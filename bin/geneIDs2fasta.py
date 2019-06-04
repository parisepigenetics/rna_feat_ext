#!/usr/bin/python3
# -*- coding: UTF-8 -*-

"""From a list of ENSEMBL Gene_IDs get a fasta file with metadata as header.

Authors: Costas Bouyioukos, AKE Franz-Arnold and LU Antoine
mail: costas.bouyioukso@univ-paris-diderot.fr aerod7710@gmail.com lu.zhao.antoine@gmail.com
2018-19
@UMR7216 Paris Diderot
"""

__version__ = "0.3a08"

import argparse
import rnaFeaturesLib

parser = argparse.ArgumentParser(prog='geneIDs2Fasta', description='Fetch transcript features (sequence and data) from ENSEMBL, for each gene ID in a list and output a multi FASTA file.', epilog="Authors: Costas Bouyioukos, Franz-Arnold Ake and Antoine Lu, 2018-19, Paris UMR7216.")
parser.add_argument('infile', nargs='?', default='-', type=argparse.FileType('r'), metavar="input_file", help='Path to the Ensembl Gene_ID list file. (or STDIN).')
parser.add_argument("outfile", nargs='?', default='-', type=argparse.FileType('w'), metavar='output_file', help="Path to output FASTA file. (or STDOUT).")
parser.add_argument('-d', '--dataset', nargs="?", default='hsapiens_gene_ensembl', metavar="ENSEMBL Dataset Name", type=str, help="Choise of the Ensembl Dataset, taken from the web API of ENSEMBL. (Default='hsapiens_gene_ensembl').")
parser.add_argument('-e', '--expressed-transcripts', help="An expressed transcripts file. It can contain an arbitrary number of columns but the first MUST be the gene name, the second the transcript ID and the third the transcription level estimate of the transcript. (Default=None).", type=argparse.FileType('r'), default=None, dest="exprTrans", metavar="EXPTRANSFile")
parser.add_argument('-v', '--version', action='version', version='%(prog)s  v. {version}'.format(version=__version__))

# Parse the command line arguments.
optArgs = parser.parse_args()

# Quickly take the genes of interest from the file.
listID = optArgs.infile.read().splitlines()

#print(rnaFeaturesLib.parse_transcripts_expression(optArgs.exprTrans))

# Connect to ENSEBL and select sequences and data.
transcripts = rnaFeaturesLib.get_ENSEMBL_data(listID, optArgs.dataset, optArgs.exprTrans)

# Print transcript sequences and data to a FASTA file.
rnaFeaturesLib.txt2fasta(transcripts, optArgs.outfile)
